#include "../interface/ChannelSumNLL.h"

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)

#include "../interface/CombineLogger.h"
#include "../interface/RooMultiPdfCombine.h"
#include "../interface/utils.h"

#include <RooAbsCategoryLValue.h>
#include <RooCmdArg.h>
#include <RooGlobalFunc.h>
#include <RooUniform.h>

#include <Math/Util.h>
#include <TString.h>

#include <algorithm>
#include <stdexcept>
#include <unordered_map>

namespace {

  std::string subName(RooAbsArg const &owner, const char *suffix) {
    return std::string(owner.GetName()) + "_" + suffix;
  }

}  // namespace

cacheutils::ChannelSumNLL::ChannelSumNLL(const char *name,
                                         const char *title,
                                         RooSimultaneous &simPdf,
                                         RooAbsData &data,
                                         const RooArgSet *constrain,
                                         std::string const &evalBackend)
    : SimNLLBase(name, title),
      simPdf_(&simPdf),
      data_(&data),
      hasConstrain_(constrain != nullptr),
      evalBackend_(evalBackend),
      originalCategories_("originalCategories", "Categories replaced by private copies in the channel pdfs", this),
      components_("components", "Likelihood terms", this) {
  if (constrain)
    constrain_.add(*constrain);
  setup_(data);
}

cacheutils::ChannelSumNLL::ChannelSumNLL(const ChannelSumNLL &other, const char *name)
    : SimNLLBase(other, name),
      simPdf_(other.simPdf_),
      data_(other.data_),
      hasConstrain_(other.hasConstrain_),
      constrain_(other.constrain_),
      evalBackend_(other.evalBackend_),
      originalCategories_("originalCategories", "Categories replaced by private copies in the channel pdfs", this),
      components_("components", "Likelihood terms", this) {
  // The RooFit likelihoods can't be shared between two sums, so a copy
  // rebuilds them from the same pdf and data. Masks and zero points are not
  // copied, the copy starts in the unmasked state.
  setup_(*data_);
}

void cacheutils::ChannelSumNLL::setup_(RooAbsData &data) {
  // Separate the constraint terms from the channel pdfs. The constraints are
  // deduplicated across channels, and the extra constraints of a
  // RooSimultaneousOpt are included.
  RooArgList constraints;
  RooAbsPdf *factorized = utils::factorizePdf(*data.get(), *simPdf_, constraints);
  factorizedPdf_ = dynamic_cast<RooSimultaneous *>(factorized);
  if (!factorizedPdf_) {
    throw std::runtime_error(std::string("ChannelSumNLL: the pdf ") + simPdf_->GetName() +
                             " did not factorize into a RooSimultaneous");
  }
  if (factorized != simPdf_)
    factorizedPdfOwner_.reset(factorized);
  constraints_.add(constraints);

  RooAbsCategoryLValue const &indexCat = factorizedPdf_->indexCat();
  // A RooSimultaneousOpt without channels is used for likelihoods made of
  // constraint terms only (the prefit nuisance fit of FitDiagnostics). There
  // is nothing to split then, and the category may not even have states.
  if (hasChannels_(indexCat))
    datasets_ = splitData_(data);

  RooFit::EvalBackend backendArg(evalBackend_);

  // One RooFit likelihood per channel. The channel pdfs contain no
  // constraint terms anymore, so createNLL() finds none to add.
  for (auto const &state : indexCat) {
    std::string const &label = state.first;
    RooAbsPdf *pdf = factorizedPdf_->getPdf(label.c_str());
    if (!pdf)
      continue;
    auto found =
        std::find_if(datasets_.begin(), datasets_.end(), [&](auto const &ds) { return label == ds->GetName(); });
    if (found == datasets_.end()) {
      throw std::logic_error("ChannelSumNLL: no data for channel " + label);
    }
    std::vector<CategorySync> syncs;
    RooAbsPdf &nllPdf = channelPdfForNLL_(*pdf, syncs);
    std::unique_ptr<RooAbsReal> nll{nllPdf.createNLL(**found, backendArg)};
    if (!nll) {
      throw std::runtime_error("ChannelSumNLL: failed to create the likelihood of channel " + label);
    }
    nll->SetName(subName(*this, label.c_str()).c_str());
    channelNames_.push_back(label);
    categorySyncs_.push_back(syncs);
    // The RooMultiPdfs of the original channel pdf: their penalty follows the
    // index of the original category.
    std::vector<RooMultiPdf const *> multiPdfs;
    std::unique_ptr<RooArgSet> components{pdf->getComponents()};
    for (RooAbsArg const *component : *components) {
      if (auto const *multiPdf = dynamic_cast<RooMultiPdf const *>(component))
        multiPdfs.push_back(multiPdf);
    }
    channelMultiPdfs_.push_back(multiPdfs);
    RooArgSet params;
    nll->getParameters(nullptr, params, /*stripDisconnected=*/false);
    replaceClonesByOriginals_(params);
    channelParams_.push_back(params);
    components_.add(*nll);
    channelNLLs_.push_back(std::move(nll));
  }
  zeroPoints_.assign(channelNLLs_.size(), 0.);

  // The constraint terms are evaluated by a likelihood of their own. Its pdf
  // is a uniform distribution of a dummy observable with a one-entry dataset,
  // whose likelihood is exactly zero, and the constraint pdfs are passed as
  // external constraints: they enter through the constraint sum like in the
  // monolithic RooFit likelihood. Passing them explicitly makes sure that all
  // of them are used, also the ones on parameters outside of the "constrain"
  // set (e.g. the ones of the "constr" directive).
  if (!constraints_.empty()) {
    constraintDummyObs_ = std::make_unique<RooRealVar>(subName(*this, "constraintDummyObs").c_str(), "", 0.5, 0., 1.);
    constraintPdf_ = std::make_unique<RooUniform>(subName(*this, "constraintPdf").c_str(), "", *constraintDummyObs_);
    constraintData_ =
        std::make_unique<RooDataSet>(subName(*this, "constraintData").c_str(), "", RooArgSet(*constraintDummyObs_));
    constraintData_->add(RooArgSet(*constraintDummyObs_));
    RooArgSet constraintSet(constraints_);
    RooCmdArg constrainArg = hasConstrain_ ? RooFit::Constrain(constrain_) : RooCmdArg::none();
    constraintNLL_ = std::unique_ptr<RooAbsReal>{constraintPdf_->createNLL(
        *constraintData_, backendArg, RooFit::ExternalConstraints(constraintSet), constrainArg)};
    if (!constraintNLL_) {
      throw std::runtime_error("ChannelSumNLL: failed to create the likelihood of the constraint terms");
    }
    constraintNLL_->SetName(subName(*this, "constraints").c_str());
    components_.add(*constraintNLL_);
  }

  hasGradient_ = !channelNLLs_.empty();
  for (auto const &nll : channelNLLs_) {
    hasGradient_ = hasGradient_ && nll->hasGradient();
  }
  if (constraintNLL_) {
    hasGradient_ = hasGradient_ && constraintNLL_->hasGradient();
  }

  CombineLogger::instance().log(
      "ChannelSumNLL.cc",
      __LINE__,
      std::string(
          Form("Created %s with %d channels and %d constraint terms, RooFit backend \"%s\", analytic gradient: %s",
               GetName(),
               int(channelNLLs_.size()),
               int(constraints_.size()),
               evalBackend_.c_str(),
               hasGradient_ ? "yes" : "no")),
      __func__);

  gradientMapsDirty_ = true;
  setValueDirty();
}

RooAbsPdf &cacheutils::ChannelSumNLL::channelPdfForNLL_(RooAbsPdf &pdf, std::vector<CategorySync> &syncs) {
  // All parameters, also those of the inactive pdfs of a RooMultiPdf.
  std::unique_ptr<RooArgSet> allParams{utils::getAllParameters(pdf)};
  RooArgSet categories;
  for (RooAbsArg *arg : *allParams) {
    if (dynamic_cast<RooCategory *>(arg))
      categories.add(*arg);
  }
  if (categories.empty())
    return pdf;

  RooArgSet privateCategories;
  for (RooAbsArg *arg : categories) {
    auto *original = static_cast<RooCategory *>(arg);
    // The copy keeps the name, so that the servers can be redirected by name
    // and the parameters of the term can be matched to the originals.
    auto clone = std::make_unique<RooCategory>(*original);
    syncs.push_back({original, clone.get()});
    privateCategories.add(*clone);
    categoryClones_.push_back(std::move(clone));
    if (!originalCategories_.containsInstance(*original))
      originalCategories_.add(*original);
  }

  // cloneTree() copies all nodes. Attach the real-valued parameters of the
  // copy back to the originals, which are shared with the other channels,
  // and only keep the categories private.
  RooArgSet sharedParams(*allParams);
  sharedParams.remove(categories);
  auto *clonedPdf = dynamic_cast<RooAbsPdf *>(pdf.cloneTree());
  if (!clonedPdf) {
    throw std::runtime_error(std::string("ChannelSumNLL: failed to clone the pdf ") + pdf.GetName());
  }
  privateChannelPdfs_.emplace_back(clonedPdf);
  clonedPdf->recursiveRedirectServers(sharedParams);
  clonedPdf->recursiveRedirectServers(privateCategories);
  return *clonedPdf;
}

void cacheutils::ChannelSumNLL::syncCategories_() const {
  for (auto const &syncs : categorySyncs_) {
    for (CategorySync const &sync : syncs) {
      if (sync.clone->getCurrentIndex() != sync.original->getCurrentIndex()) {
        sync.clone->setIndex(sync.original->getCurrentIndex());
      }
    }
  }
}

void cacheutils::ChannelSumNLL::replaceClonesByOriginals_(RooArgSet &params) const {
  for (auto const &syncs : categorySyncs_) {
    for (CategorySync const &sync : syncs) {
      if (params.containsInstance(*sync.clone)) {
        params.remove(*sync.clone);
        params.add(*sync.original, /*silent=*/true);
      }
    }
  }
}

std::vector<std::unique_ptr<RooAbsData>> cacheutils::ChannelSumNLL::splitData_(RooAbsData &data) const {
  // Split by the index category of the simultaneous pdf. Each piece contains
  // only the observables of its channel, keeps the weights, including the
  // zero-weight entries that RooFit needs for binned likelihoods, and pieces
  // are created also for channels without entries.
  return data.split(*factorizedPdf_, /*createEmptyDataSets=*/true);
}

bool cacheutils::ChannelSumNLL::hasChannels_(RooAbsCategoryLValue const &indexCat) const {
  for (auto const &state : indexCat) {
    if (factorizedPdf_->getPdf(state.first.c_str()))
      return true;
  }
  return false;
}

bool cacheutils::ChannelSumNLL::setData(RooAbsData &data, bool /*cloneData*/) {
  if (channelNLLs_.empty()) {
    data_ = &data;
    return true;
  }
  auto newDatasets = splitData_(data);
  for (std::size_t i = 0; i < channelNLLs_.size(); ++i) {
    std::string const &label = channelNames_[i];
    auto found =
        std::find_if(newDatasets.begin(), newDatasets.end(), [&](auto const &ds) { return label == ds->GetName(); });
    if (found == newDatasets.end()) {
      throw std::logic_error("ChannelSumNLL::setData: no data for channel " + label);
    }
    // The channel likelihoods keep a pointer to the data, we own it.
    channelNLLs_[i]->setData(**found, /*cloneData=*/false);
  }
  // Only release the old datasets after all terms point to the new ones.
  datasets_ = std::move(newDatasets);
  data_ = &data;
  setValueDirty();
  return true;
}

double cacheutils::ChannelSumNLL::evaluate() const {
  syncCategories_();
  ROOT::Math::KahanSum<double> sum;
  for (std::size_t i = 0; i < channelNLLs_.size(); ++i) {
    if (isChannelMasked(i))
      continue;
    sum += channelNLLs_[i]->getVal() - zeroPoints_[i];
    for (RooMultiPdf const *multiPdf : channelMultiPdfs_[i]) {
      sum += multiPdf->getCorrection();
    }
  }
  if (constraintNLL_) {
    sum += constraintNLL_->getVal() - constraintZeroPoint_;
  }
  sum += maskingOffset_ - maskingOffsetZero_;
  return sum.Sum();
}

bool cacheutils::ChannelSumNLL::getParameters(const RooArgSet *observables,
                                              RooArgSet &outputSet,
                                              bool stripDisconnected) const {
  if (internalMasks_.empty()) {
    // Collect per term and map the category copies to the originals before
    // merging: the copies of different channels share their names.
    auto addParamsOf = [&](RooAbsReal const &term) {
      RooArgSet termParams;
      term.getParameters(observables, termParams, stripDisconnected);
      replaceClonesByOriginals_(termParams);
      outputSet.add(termParams, /*silent=*/true);
    };
    for (auto const &nll : channelNLLs_) {
      addParamsOf(*nll);
    }
    if (constraintNLL_) {
      addParamsOf(*constraintNLL_);
    }
  } else {
    // Only the parameters of the active channels. Parameters that appear
    // only in masked channels or in the constraint terms stay out, which
    // reduces the dimension of the fits in the discrete profiling scans.
    outputSet.add(activeParameters_, /*silent=*/true);
    if (observables) {
      outputSet.remove(*observables, /*silent=*/true, /*matchByNameOnly=*/true);
    }
  }
  replaceClonesByOriginals_(outputSet);
  // Same ordering convention as RooAbsArg::getParameters().
  outputSet.sort();
  return false;
}

void cacheutils::ChannelSumNLL::applyWeightSquared(bool flag) {
  for (auto const &nll : channelNLLs_) {
    nll->applyWeightSquared(flag);
  }
}

void cacheutils::ChannelSumNLL::updateGradientMaps_() const {
  // The gradient of each term is laid out like the term's own unstripped
  // parameter list, see RooMinimizerFcn. Map these positions into the layout
  // of our own list, which is what the minimizer asks us for.
  RooArgSet params;
  getParameters(nullptr, params, /*stripDisconnected=*/false);
  nGradientParams_ = params.size();
  std::unordered_map<std::string, int> indexByName;
  int idx = 0;
  for (RooAbsArg const *arg : params) {
    indexByName[arg->GetName()] = idx++;
  }

  auto makeMap = [&](RooAbsReal const &term) {
    RooArgSet termParams;
    term.getParameters(nullptr, termParams, /*stripDisconnected=*/false);
    std::vector<int> map;
    map.reserve(termParams.size());
    for (RooAbsArg const *arg : termParams) {
      auto found = indexByName.find(arg->GetName());
      map.push_back(found == indexByName.end() ? -1 : found->second);
    }
    return map;
  };

  gradientMaps_.clear();
  for (auto const &nll : channelNLLs_) {
    gradientMaps_.push_back(makeMap(*nll));
  }
  if (constraintNLL_) {
    gradientMaps_.push_back(makeMap(*constraintNLL_));
  }
  gradientMapsDirty_ = false;
}

void cacheutils::ChannelSumNLL::gradient(double *out) const {
  if (!hasGradient_) {
    throw std::logic_error("ChannelSumNLL::gradient() called, but not all likelihood terms provide a gradient");
  }
  if (gradientMapsDirty_)
    updateGradientMaps_();
  syncCategories_();
  std::fill(out, out + nGradientParams_, 0.0);

  auto accumulate = [&](RooAbsReal const &term, std::vector<int> const &map) {
    gradientBuffer_.assign(map.size(), 0.0);
    term.gradient(gradientBuffer_.data());
    for (std::size_t i = 0; i < map.size(); ++i) {
      if (map[i] >= 0)
        out[map[i]] += gradientBuffer_[i];
    }
  };

  for (std::size_t i = 0; i < channelNLLs_.size(); ++i) {
    if (isChannelMasked(i))
      continue;
    accumulate(*channelNLLs_[i], gradientMaps_[i]);
  }
  if (constraintNLL_) {
    accumulate(*constraintNLL_, gradientMaps_.back());
  }
}

void cacheutils::ChannelSumNLL::setZeroPoint() {
  for (std::size_t i = 0; i < channelNLLs_.size(); ++i) {
    zeroPoints_[i] = channelNLLs_[i]->getVal();
  }
  if (constraintNLL_) {
    constraintZeroPoint_ = constraintNLL_->getVal();
  }
  maskingOffsetZero_ = maskingOffset_;
  setValueDirty();
}

void cacheutils::ChannelSumNLL::clearZeroPoint() {
  std::fill(zeroPoints_.begin(), zeroPoints_.end(), 0.0);
  constraintZeroPoint_ = 0;
  maskingOffsetZero_ = 0;
  setValueDirty();
}

void cacheutils::ChannelSumNLL::setMaskNonDiscreteChannels(bool mask) {
  double nllBefore = evaluate();
  internalMasks_.clear();
  activeParameters_.removeAll();
  if (mask) {
    internalMasks_.resize(channelNLLs_.size(), false);
    for (std::size_t i = 0; i < channelNLLs_.size(); ++i) {
      for (RooAbsArg *param : channelParams_[i]) {
        auto *cat = dynamic_cast<RooCategory *>(param);
        if (cat && !cat->isConstant()) {
          internalMasks_[i] = true;
          activeParameters_.add(channelParams_[i], /*silent=*/true);
          CombineLogger::instance().log("ChannelSumNLL.cc",
                                        __LINE__,
                                        std::string(Form("Enabling channel %s that depends on non-constant category %s",
                                                         channelNames_[i].c_str(),
                                                         cat->GetName())),
                                        __func__);
          break;
        }
      }
    }
  }
  gradientMapsDirty_ = true;
  setValueDirty();
  double nllAfter = evaluate();
  maskingOffset_ += (nllBefore - nllAfter);
  setValueDirty();
}

#endif  // ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)
