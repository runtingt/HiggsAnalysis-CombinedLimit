#ifndef HiggsAnalysis_CombinedLimit_SimNLLBase_h
#define HiggsAnalysis_CombinedLimit_SimNLLBase_h

#include <RooAbsReal.h>

namespace cacheutils {

  /// Interface shared by Combine's likelihood classes for simultaneous pdfs,
  /// CachingSimNLL and ChannelSumNLL.
  ///
  /// The CascadeMinimizer talks to the likelihood through this interface for
  /// the operations that don't depend on how the channel terms are evaluated:
  /// the zero-point subtraction that keeps the function values seen by Minuit
  /// close to zero, and the channel masking of the discrete profiling scans.
  class SimNLLBase : public RooAbsReal {
  public:
    SimNLLBase() = default;
    SimNLLBase(const char *name, const char *title) : RooAbsReal(name, title) {}
    SimNLLBase(const SimNLLBase &other, const char *name = nullptr) : RooAbsReal(other, name) {}

    /// Subtract the current value of every term from the likelihood, so that
    /// the minimizer works with numbers close to zero.
    virtual void setZeroPoint() = 0;
    /// Undo setZeroPoint().
    virtual void clearZeroPoint() = 0;
    void updateZeroPoint() {
      clearZeroPoint();
      setZeroPoint();
    }
    /// Remove any constant offset applied at construction time, so that the
    /// value is the full negative log-likelihood (needed when likelihoods of
    /// different models are compared, like in the goodness-of-fit tests).
    virtual void clearConstantZeroPoint() = 0;

    /// Keep only the channels whose parameters include a floating RooCategory,
    /// i.e. the ones taking part in a discrete profiling scan, in the
    /// likelihood sum and in getParameters(). The value of the likelihood is
    /// kept continuous by an offset, so that masked and unmasked values can be
    /// compared. Calling it with `mask = false` restores all channels.
    virtual void setMaskNonDiscreteChannels(bool mask) = 0;
  };

}  // namespace cacheutils

#endif
