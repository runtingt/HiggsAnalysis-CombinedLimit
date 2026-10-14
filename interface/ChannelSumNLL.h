#ifndef HiggsAnalysis_CombinedLimit_ChannelSumNLL_h
#define HiggsAnalysis_CombinedLimit_ChannelSumNLL_h

#include <ROOT/RConfig.hxx>  // for ROOT_VERSION

// The RooFit evaluation backends that this class builds on exist since ROOT 6.30.
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)

#include "SimNLLBase.h"

#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooCategory.h>
#include <RooDataSet.h>
#include <RooListProxy.h>
#include <RooRealVar.h>
#include <RooSetProxy.h>
#include <RooSimultaneous.h>

class RooMultiPdf;

#include <memory>
#include <string>
#include <vector>

namespace cacheutils {

  /// Likelihood of a RooSimultaneous assembled from one RooFit likelihood per
  /// channel, plus one for the constraint terms.
  ///
  /// Each term is created with RooAbsPdf::createNLL() and the requested
  /// RooFit evaluation backend ("cpu", "codegen", ...), so every term owns its
  /// own RooFit::Evaluator. The sum over the terms is done here, which is what
  /// makes channel masking possible on these backends: a masked channel is
  /// left out of the sum and of getParameters(), and its evaluator never runs.
  /// Inside a single evaluator this is impossible, because the evaluator always
  /// recomputes every dirty node.
  ///
  /// The constraint terms are factorized out of the channel pdfs and evaluated
  /// by a separate likelihood that is never masked. It is created with the same
  /// backend, through createNLL() with the constraint pdfs as external
  /// constraints. The pdf of that likelihood is a uniform distribution of a
  /// dummy observable with a one-entry dataset, which contributes exactly zero,
  /// so only the constraint sum remains, normalized like in the monolithic
  /// RooFit likelihood.
  ///
  /// Like CachingSimNLL, this class adds the penalty term of every RooMultiPdf
  /// (RooMultiPdf::getCorrection(), which depends on the active function) to
  /// the channel it belongs to. RooFit's createNLL() only knows the correction
  /// of the top-level pdf, which is zero when the RooMultiPdf is a component.
  ///
  /// With the "codegen" backend every term provides an analytic gradient, and
  /// this class sums them: hasGradient() is true only if all terms have one,
  /// there is no numerical fallback for individual terms.
  ///
  /// The RooFit likelihoods keep their data by pointer, so this class owns the
  /// per-channel datasets created by splitting the input data. setData()
  /// splits the new dataset and hands the pieces to the channel likelihoods.
  ///
  /// A RooFit::Evaluator stores a data token in every node of its graph, so a
  /// node can't be part of two evaluators, except for RooRealVar parameters.
  /// createNLL() compiles a private copy of the pdf graph, but the fundamental
  /// nodes stay shared, and the RooCategory of a RooMultiPdf is one of them.
  /// Channels that depend on such a category therefore get a cloned pdf in
  /// which the categories are replaced by private copies. The copies are kept
  /// in sync with the originals before every evaluation, and getParameters()
  /// reports the originals, so nothing else in Combine needs to know.
  class ChannelSumNLL : public SimNLLBase {
  public:
    /// \param simPdf     The simultaneous pdf, usually a RooSimultaneousOpt
    ///                   whose channel pdfs still contain the constraint terms.
    /// \param data       The combined dataset, containing the index category.
    /// \param constrain  Parameters to constrain (the nuisance parameters), or
    ///                   nullptr for RooFit's automatic selection.
    /// \param evalBackend  Name of the RooFit evaluation backend.
    ChannelSumNLL(const char *name,
                  const char *title,
                  RooSimultaneous &simPdf,
                  RooAbsData &data,
                  const RooArgSet *constrain,
                  std::string const &evalBackend);
    ChannelSumNLL(const ChannelSumNLL &other, const char *name = nullptr);
    ChannelSumNLL *clone(const char *name = nullptr) const override { return new ChannelSumNLL(*this, name); }

    bool isDerived() const override { return true; }
    double defaultErrorLevel() const override { return 0.5; }
    /// Always evaluate: the terms track their own dirty state, and the private
    /// category copies are synchronized in evaluate().
    double getValV(const RooArgSet *) const override { return evaluate(); }

    /// Replace the data of all channels. The new dataset must have the same
    /// structure as the original one (same columns and index category).
    bool setData(RooAbsData &data, bool cloneData = true) override;
    bool getParameters(const RooArgSet *observables,
                       RooArgSet &outputSet,
                       bool stripDisconnected = true) const override;
    void applyWeightSquared(bool flag) override;

    bool hasGradient() const override { return hasGradient_; }
    /// Gradient with respect to the parameters returned by
    /// getParameters(nullptr, set, /*stripDisconnected=*/false), in their
    /// order, like RooMinimizer expects it. Masked channels don't contribute.
    void gradient(double *out) const override;

    void setZeroPoint() override;
    void clearZeroPoint() override;
    /// There is no constant zero point in this class, nothing to do.
    void clearConstantZeroPoint() override {}
    void setMaskNonDiscreteChannels(bool mask) override;

    std::size_t numChannels() const { return channelNLLs_.size(); }
    std::string const &channelName(std::size_t i) const { return channelNames_[i]; }
    RooAbsReal const &channelNLL(std::size_t i) const { return *channelNLLs_[i]; }
    bool isChannelMasked(std::size_t i) const { return !internalMasks_.empty() && !internalMasks_[i]; }
    /// The likelihood of the constraint terms, or nullptr if there are none.
    RooAbsReal const *constraintNLL() const { return constraintNLL_.get(); }
    std::string const &evalBackend() const { return evalBackend_; }

  protected:
    double evaluate() const override;

  private:
    struct CategorySync {
      RooCategory *original = nullptr;
      RooCategory *clone = nullptr;
    };

    void setup_(RooAbsData &data);
    std::vector<std::unique_ptr<RooAbsData>> splitData_(RooAbsData &data) const;
    /// The pdf to build the channel likelihood from: the channel pdf itself, or
    /// a clone with private copies of the categories it depends on.
    RooAbsPdf &channelPdfForNLL_(RooAbsPdf &pdf, std::vector<CategorySync> &syncs);
    void syncCategories_() const;
    void replaceClonesByOriginals_(RooArgSet &params) const;
    void updateGradientMaps_() const;

    RooSimultaneous *simPdf_ = nullptr;  // not owned
    RooAbsData *data_ = nullptr;         // not owned, kept for cloning
    bool hasConstrain_ = false;
    RooArgSet constrain_;
    std::string evalBackend_;

    std::unique_ptr<RooAbsPdf> factorizedPdfOwner_;  // set if the factorization created a new pdf
    RooSimultaneous *factorizedPdf_ = nullptr;       // channel pdfs without constraint terms
    RooArgList constraints_;

    std::vector<std::string> channelNames_;
    std::vector<std::unique_ptr<RooAbsPdf>> privateChannelPdfs_;
    std::vector<std::unique_ptr<RooCategory>> categoryClones_;
    std::vector<std::vector<CategorySync>> categorySyncs_;
    RooSetProxy originalCategories_;  // the original categories as servers, for graph traversals
    std::vector<std::unique_ptr<RooAbsData>> datasets_;
    std::vector<std::unique_ptr<RooAbsReal>> channelNLLs_;
    std::vector<RooArgSet> channelParams_;                            // all parameters of each channel likelihood
    std::vector<std::vector<RooMultiPdf const *>> channelMultiPdfs_;  // for the penalty terms
    std::vector<double> zeroPoints_;

    std::unique_ptr<RooRealVar> constraintDummyObs_;
    std::unique_ptr<RooAbsPdf> constraintPdf_;
    std::unique_ptr<RooDataSet> constraintData_;
    std::unique_ptr<RooAbsReal> constraintNLL_;
    double constraintZeroPoint_ = 0;

    // All terms as servers, for dirty-state propagation and graph traversal.
    RooListProxy components_;

    bool hasGradient_ = false;

    // Channel masking, see setMaskNonDiscreteChannels(). Empty if no mask is
    // set, otherwise true for the channels that are kept.
    std::vector<bool> internalMasks_;
    RooArgSet activeParameters_;
    double maskingOffset_ = 0;      // keeps the value continuous when channels are masked
    double maskingOffsetZero_ = 0;  // and the corresponding zero point

    // Maps from the parameter index of each term to the index in our own
    // parameter list, -1 if the parameter is not in our list (masked).
    mutable bool gradientMapsDirty_ = true;
    mutable std::size_t nGradientParams_ = 0;
    mutable std::vector<std::vector<int>> gradientMaps_;
    mutable std::vector<double> gradientBuffer_;
  };

}  // namespace cacheutils

#endif  // ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)

#endif
