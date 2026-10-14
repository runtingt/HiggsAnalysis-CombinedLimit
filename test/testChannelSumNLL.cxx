// Tests of cacheutils::ChannelSumNLL, the likelihood that Combine builds for
// a RooSimultaneousOpt on the RooFit evaluation backends: one RooFit
// likelihood per channel plus one for the constraint terms, summed by
// Combine so that channels can be masked in the discrete profiling scans.
//
// The first test uses the RooMultiPdf workspace of the CI
// (ws_RooMultiPdf.root, produced by the "multipdf" text2workspace test). The
// others use a small toy model built here, because the CI workspace contains
// a signal pdf (RooDoubleCBFast) that the "codegen" backend doesn't support.

#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "TError.h"
#include "TFile.h"

#include "Math/MinimizerOptions.h"
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooMsgService.h"
#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"

#include "../interface/CachingNLL.h"
#include "../interface/CascadeMinimizer.h"
#include "../interface/ChannelSumNLL.h"
#include "../interface/Combine.h"
#include "../interface/ProfilingTools.h"
#include "../interface/RooMultiPdfCombine.h"
#include "../interface/RooSimultaneousOpt.h"

#include <gtest/gtest.h>

namespace {

  // The runtime defines that bin/combine.cpp sets by default and that change
  // how CachingSimNLL evaluates the likelihood.
  void setCombineRuntimeDefaults() {
    for (const char *flag : {"OPTIMIZE_BOUNDS",
                             "ADDNLL_RECURSIVE",
                             "ADDNLL_GAUSSNLL",
                             "ADDNLL_HISTNLL",
                             "ADDNLL_CBNLL",
                             "TMCSO_AdaptivePseudoAsimov",
                             "MINIMIZER_rooFitOffset",
                             "ADDNLL_ROOREALSUM_FACTOR",
                             "ADDNLL_ROOREALSUM_NONORM",
                             "ADDNLL_ROOREALSUM_BASICINT",
                             "ADDNLL_ROOREALSUM_KEEPZEROS",
                             "ADDNLL_PRODNLL",
                             "ADDNLL_HFNLL",
                             "ADDNLL_HISTFUNCNLL",
                             "ADDNLL_ROOREALSUM_CHEAPPROD"}) {
      runtimedef::set(flag, 1);
    }
    runtimedef::set("MINIMIZER_optimizeConst", 2);
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  }

  std::unique_ptr<RooAbsReal> makeNLL(RooAbsPdf &pdf,
                                      RooAbsData &data,
                                      const RooArgSet *nuisances,
                                      std::string const &backend) {
    Combine::setNllBackend(backend);
    return combineCreateNLL(pdf, data, nuisances, /*offset=*/true);
  }

  // Floating real-valued parameters of a likelihood, in the order of its
  // gradient (the unstripped, name-sorted parameter list).
  struct GradientLayout {
    RooArgSet allParams;
    std::vector<std::pair<std::size_t, RooRealVar *>> floating;  // (index in the gradient, parameter)
  };

  GradientLayout gradientLayout(RooAbsReal const &nll) {
    GradientLayout layout;
    nll.getParameters(nullptr, layout.allParams, /*stripDisconnected=*/false);
    std::size_t i = 0;
    for (RooAbsArg *arg : layout.allParams) {
      if (auto *var = dynamic_cast<RooRealVar *>(arg); var && !var->isConstant())
        layout.floating.emplace_back(i, var);
      ++i;
    }
    return layout;
  }

  // Compares the analytic gradient of the likelihood with central finite
  // differences for every floating parameter.
  void checkGradientAgainstFiniteDifferences(RooAbsReal &nll, const char *label) {
    ASSERT_TRUE(nll.hasGradient()) << label;
    GradientLayout layout = gradientLayout(nll);
    ASSERT_FALSE(layout.floating.empty()) << label;
    std::vector<double> grad(layout.allParams.size(), 0.);
    nll.gradient(grad.data());
    for (auto const &[index, var] : layout.floating) {
      const double v0 = var->getVal();
      const double h = 1e-4 * (1. + std::abs(v0));
      var->setVal(v0 + h);
      const double up = nll.getVal();
      var->setVal(v0 - h);
      const double down = nll.getVal();
      var->setVal(v0);
      const double fd = (up - down) / (2 * h);
      EXPECT_NEAR(grad[index], fd, 1e-3 * (1. + std::abs(fd))) << label << ": parameter " << var->GetName();
    }
  }

}  // namespace

// On the CI workspace, the sum of the per-channel likelihoods must reproduce
// RooFit's monolithic likelihood of the same backend. They differ by a known
// constant: the monolithic one adds sum(weights) * log(numberOfChannels) in
// every channel (RooNLLVarNew with simCount > 1), the per-channel ones don't.
TEST(ChannelSumNLL, ReproducesMonolithicRooFitLikelihood) {
  setCombineRuntimeDefaults();
  std::unique_ptr<TFile> file{TFile::Open("ws_RooMultiPdf.root", "READ")};
  ASSERT_TRUE(file && !file->IsZombie());
  auto *w = dynamic_cast<RooWorkspace *>(file->Get("w"));
  ASSERT_TRUE(w);
  auto *mc = dynamic_cast<RooStats::ModelConfig *>(w->genobj("ModelConfig"));
  RooAbsData *data = w->data("data_obs");
  auto *pdf = dynamic_cast<RooSimultaneousOpt *>(mc->GetPdf());
  RooCategory *index = w->cat("pdf_index_ggh");
  ASSERT_TRUE(mc && data && pdf && index);
  w->var("MH")->setVal(125.38);
  w->var("MH")->setConstant(true);

  runtimedef::set("SIMNLL_ROOFIT_MONOLITHIC", 1);
  auto nllMonolithic = makeNLL(*pdf, *data, mc->GetNuisanceParameters(), "cpu");
  runtimedef::set("SIMNLL_ROOFIT_MONOLITHIC", 0);
  auto nllSplit = makeNLL(*pdf, *data, mc->GetNuisanceParameters(), "cpu");
  ASSERT_FALSE(dynamic_cast<cacheutils::ChannelSumNLL *>(nllMonolithic.get()));
  auto *sum = dynamic_cast<cacheutils::ChannelSumNLL *>(nllSplit.get());
  ASSERT_TRUE(sum);
  EXPECT_EQ(sum->numChannels(), 5u);
  EXPECT_TRUE(sum->constraintNLL());

  const double expectedOffset = data->sumEntries() * std::log(double(sum->numChannels()));
  for (int i = 0; i < index->numTypes(); ++i) {
    index->setIndex(i);
    const double split = nllSplit->getVal();
    const double monolithic = nllMonolithic->getVal();
    std::cout << "index " << i << ": split " << split << ", monolithic " << monolithic
              << ", difference minus expected offset " << (monolithic - split - expectedOffset) << std::endl;
    EXPECT_NEAR(monolithic - split, expectedOffset, 1e-9 * std::abs(monolithic)) << "index " << i;
  }
  // Both parameter lists must agree (the split one reports the original category).
  RooArgSet paramsSplit;
  RooArgSet paramsMonolithic;
  nllSplit->getParameters(nullptr, paramsSplit, false);
  nllMonolithic->getParameters(nullptr, paramsMonolithic, false);
  EXPECT_TRUE(paramsSplit.equals(paramsMonolithic));
  EXPECT_TRUE(paramsSplit.containsInstance(*index));
}

// The toy model of test/makeToyMultiPdfWorkspace.C (toy_multipdf.root): two
// channels sharing the signal strength "r" and the nuisance "lumi". Channel A
// has a RooMultiPdf background and takes part in the discrete profiling,
// channel B has a fixed background shape and a nuisance "bnorm" of its own.
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 37, 0)  // RooFit's RooMultiPdf, which has codegen support

namespace {

  struct ToyModel {
    std::unique_ptr<TFile> file;
    RooWorkspace *w = nullptr;
    RooSimultaneousOpt *pdf = nullptr;
    RooDataSet *data = nullptr;
    RooCategory *index = nullptr;
    RooMultiPdf *multiPdf = nullptr;
    RooRealVar *r = nullptr;
    RooArgSet nuisances;
    RooArgSet observables;
  };

  ToyModel buildToyModel() {
    ToyModel m;
    m.file.reset(TFile::Open("toy_multipdf.root", "READ"));
    if (!m.file || m.file->IsZombie())
      throw std::runtime_error("could not open toy_multipdf.root");
    m.w = dynamic_cast<RooWorkspace *>(m.file->Get("w"));
    if (!m.w)
      throw std::runtime_error("no workspace 'w' in toy_multipdf.root");
    RooWorkspace &w = *m.w;
    auto *mc = dynamic_cast<RooStats::ModelConfig *>(w.genobj("ModelConfig"));
    m.pdf = dynamic_cast<RooSimultaneousOpt *>(w.pdf("model_s"));
    m.data = dynamic_cast<RooDataSet *>(w.data("data_obs"));
    m.index = w.cat("pdf_index_A");
    m.multiPdf = dynamic_cast<RooMultiPdf *>(w.pdf("bkgA"));
    m.r = w.var("r");
    if (!mc || !m.pdf || !m.data || !m.index || !m.multiPdf || !m.r)
      throw std::runtime_error("toy workspace content not as expected");
    m.nuisances.add(*mc->GetNuisanceParameters());
    m.observables.add(*mc->GetObservables());
    m.index->setIndex(1);
    w.saveSnapshot("start", w.allVars());
    return m;
  }

  void restoreStart(ToyModel &m) {
    m.w->loadSnapshot("start");
    m.index->setIndex(1);
  }

  // A second parameter point, away from the start values.
  void moveToSecondPoint(ToyModel &m) {
    m.w->var("r")->setVal(0.4);
    m.w->var("lumi")->setVal(0.8);
    m.w->var("bnorm")->setVal(-0.5);
    m.w->var("slopeA")->setVal(-0.05);
    m.w->var("cA1")->setVal(-0.2);
    m.w->var("pA1")->setVal(-0.002);
    m.w->var("nbkgA")->setVal(2100.);
    m.w->var("slopeB")->setVal(-0.015);
    m.w->var("nbkgB")->setVal(900.);
  }

  std::vector<std::string> const kBackends{"combine", "cpu", "codegen"};

}  // namespace

// All backends must give the same likelihood up to a constant, at every pdf
// index and parameter point. The RooMultiPdf penalty is included in all of
// them, so the differences between the indices agree too.
TEST(ChannelSumNLLToy, ValuesAgreeBetweenBackends) {
  setCombineRuntimeDefaults();
  ToyModel m = buildToyModel();
  std::map<std::string, std::unique_ptr<RooAbsReal>> nlls;
  for (auto const &backend : kBackends) {
    nlls[backend] = makeNLL(*m.pdf, *m.data, &m.nuisances, backend);
    ASSERT_TRUE(nlls[backend]) << backend;
  }
  ASSERT_TRUE(dynamic_cast<cacheutils::CachingSimNLL *>(nlls["combine"].get()));
  ASSERT_TRUE(dynamic_cast<cacheutils::ChannelSumNLL *>(nlls["cpu"].get()));
  ASSERT_TRUE(dynamic_cast<cacheutils::ChannelSumNLL *>(nlls["codegen"].get()));
  EXPECT_FALSE(nlls["cpu"]->hasGradient());
  EXPECT_TRUE(nlls["codegen"]->hasGradient());

  std::map<std::string, double> offset;
  bool first = true;
  for (int point = 0; point < 2; ++point) {
    restoreStart(m);
    if (point == 1)
      moveToSecondPoint(m);
    for (int i = 0; i < m.index->numTypes(); ++i) {
      m.index->setIndex(i);
      const double reference = nlls["combine"]->getVal();
      const double penalty = m.multiPdf->getCorrection();
      for (auto const &backend : kBackends) {
        const double value = nlls[backend]->getVal();
        if (first)
          offset[backend] = value - reference;
        EXPECT_NEAR(value - reference, offset[backend], 1e-4) << backend << ", point " << point << ", index " << i;
      }
      std::cout << "point " << point << " index " << i << " (penalty " << penalty << "): combine " << reference
                << ", cpu - combine " << nlls["cpu"]->getVal() - reference << ", codegen - combine "
                << nlls["codegen"]->getVal() - reference << std::endl;
      first = false;
    }
  }
  EXPECT_GT(m.multiPdf->getCorrection(), 0.) << "the penalty term should not vanish in this model";
}

// With the codegen backend the summed gradient must match finite differences
// of the summed likelihood, unmasked and with the non-discrete channel masked.
TEST(ChannelSumNLLToy, CodegenGradient) {
  setCombineRuntimeDefaults();
  ToyModel m = buildToyModel();
  auto nll = makeNLL(*m.pdf, *m.data, &m.nuisances, "codegen");
  auto *sum = dynamic_cast<cacheutils::ChannelSumNLL *>(nll.get());
  ASSERT_TRUE(sum);
  for (int point = 0; point < 2; ++point) {
    restoreStart(m);
    if (point == 1)
      moveToSecondPoint(m);
    for (int i = 0; i < m.index->numTypes(); ++i) {
      m.index->setIndex(i);
      checkGradientAgainstFiniteDifferences(*nll, (std::string("unmasked, index ") + std::to_string(i)).c_str());
      sum->setMaskNonDiscreteChannels(true);
      checkGradientAgainstFiniteDifferences(*nll, (std::string("masked, index ") + std::to_string(i)).c_str());
      sum->setMaskNonDiscreteChannels(false);
    }
  }
}

// Masking keeps the value continuous, drops the channel without a floating
// RooCategory from the sum and from the parameter list, and is reversible.
TEST(ChannelSumNLLToy, MaskNonDiscreteChannels) {
  setCombineRuntimeDefaults();
  ToyModel m = buildToyModel();
  for (auto const &backend : {"cpu", "codegen"}) {
    restoreStart(m);
    auto nll = makeNLL(*m.pdf, *m.data, &m.nuisances, backend);
    auto *sum = dynamic_cast<cacheutils::ChannelSumNLL *>(nll.get());
    ASSERT_TRUE(sum) << backend;
    ASSERT_EQ(sum->numChannels(), 2u);
    RooArgSet fullParams;
    nll->getParameters(nullptr, fullParams, false);
    const double before = nll->getVal();

    sum->setMaskNonDiscreteChannels(true);
    EXPECT_FALSE(sum->isChannelMasked(0)) << backend;  // A has the RooMultiPdf
    EXPECT_TRUE(sum->isChannelMasked(1)) << backend;
    EXPECT_NEAR(nll->getVal(), before, 1e-9 * std::abs(before)) << backend;
    RooArgSet maskedParams;
    nll->getParameters(nullptr, maskedParams, false);
    for (const char *name : {"r", "lumi", "slopeA", "cA1", "cA2", "pA1", "nbkgA", "pdf_index_A"})
      EXPECT_TRUE(maskedParams.find(name)) << backend << ": " << name << " should be active";
    for (const char *name : {"bnorm", "slopeB", "nbkgB"})
      EXPECT_FALSE(maskedParams.find(name)) << backend << ": " << name << " should be hidden";
    EXPECT_TRUE(maskedParams.containsInstance(*m.index)) << backend << ": the original category must be reported";

    // The masked channel does not contribute anymore: moving one of its
    // private parameters leaves the value unchanged.
    m.w->var("nbkgB")->setVal(1200.);
    EXPECT_NEAR(nll->getVal(), before, 1e-9 * std::abs(before)) << backend;
    restoreStart(m);

    // With the category constant no channel takes part in the discrete
    // profiling, so all of them are masked (the CascadeMinimizer uses this
    // in mode 2, where the categories of the other channels are constant).
    sum->setMaskNonDiscreteChannels(false);
    m.index->setConstant(true);
    sum->setMaskNonDiscreteChannels(true);
    EXPECT_TRUE(sum->isChannelMasked(0)) << backend;
    EXPECT_TRUE(sum->isChannelMasked(1)) << backend;
    EXPECT_NEAR(nll->getVal(), before, 1e-9 * std::abs(before)) << backend;
    m.index->setConstant(false);

    sum->setMaskNonDiscreteChannels(false);
    EXPECT_NEAR(nll->getVal(), before, 1e-9 * std::abs(before)) << backend;
    RooArgSet restoredParams;
    nll->getParameters(nullptr, restoredParams, false);
    EXPECT_TRUE(restoredParams.equals(fullParams)) << backend;
  }
}

// setData() must give the same likelihood as building a new one on the new
// dataset (used when fitting toys).
TEST(ChannelSumNLLToy, SetData) {
  setCombineRuntimeDefaults();
  ToyModel m = buildToyModel();
  RooRandom::randomGenerator()->SetSeed(12345);
  std::unique_ptr<RooDataSet> toy{m.pdf->generate(m.observables, RooFit::Extended(), RooFit::Name("toy"))};
  for (auto const &backend : {"cpu", "codegen"}) {
    restoreStart(m);
    auto nll = makeNLL(*m.pdf, *m.data, &m.nuisances, backend);
    auto fresh = makeNLL(*m.pdf, *toy, &m.nuisances, backend);
    ASSERT_TRUE(dynamic_cast<cacheutils::ChannelSumNLL *>(nll.get())) << backend;
    EXPECT_TRUE(nll->setData(*toy, false)) << backend;
    for (int i = 0; i < m.index->numTypes(); ++i) {
      m.index->setIndex(i);
      EXPECT_NEAR(nll->getVal(), fresh->getVal(), 1e-9 * std::abs(fresh->getVal())) << backend << ", index " << i;
    }
    moveToSecondPoint(m);
    EXPECT_NEAR(nll->getVal(), fresh->getVal(), 1e-9 * std::abs(fresh->getVal())) << backend;
    if (std::string(backend) == "codegen")
      checkGradientAgainstFiniteDifferences(*nll, "after setData");
  }
}

// The discrete profiling of the CascadeMinimizer with channel masking must
// find the same minimum on every backend.
TEST(ChannelSumNLLToy, DiscreteProfilingFit) {
  setCombineRuntimeDefaults();
  // These are read into static variables on first use, so they have to be
  // set before the first CascadeMinimizer call in this process.
  runtimedef::set("MINIMIZER_freezeDisassociatedParams", 1);
  runtimedef::set("MINIMIZER_multiMin_maskChannels", 2);
  CascadeMinimizer::initOptions();
  gErrorIgnoreLevel = kWarning;  // silence the Minuit2 info messages
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.001);

  ToyModel m = buildToyModel();
  // Note: assigning a RooArgList copies the values of same-named elements,
  // not the membership, hence removeAll() and add().
  auto &cfg = CascadeMinimizerGlobalConfigs::O();
  cfg.pdfCategories.removeAll();
  cfg.pdfCategories.add(*m.index);
  cfg.allRooMultiPdfs.removeAll();
  cfg.allRooMultiPdfs.add(*m.multiPdf);
  cfg.parametersOfInterest.removeAll();
  cfg.parametersOfInterest.add(*m.r);
  cfg.nuisanceParameters.removeAll();
  cfg.nuisanceParameters.add(m.nuisances);
  cfg.allFloatingParameters.removeAll();
  for (RooAbsArg *arg : m.w->allVars()) {
    if (!arg->isConstant())
      cfg.allFloatingParameters.add(*arg);
  }

  struct Result {
    double r = 0;
    double rError = 0;
    int index = -1;
    double nllImprovement = 0;
  };
  std::map<std::string, Result> results;
  for (auto const &backend : kBackends) {
    restoreStart(m);
    auto nll = makeNLL(*m.pdf, *m.data, &m.nuisances, backend);
    const double nllStart = nll->getVal();
    CascadeMinimizer minim(*nll, CascadeMinimizer::Constrained, m.r);
    ASSERT_TRUE(minim.minimize(/*verbose=*/0)) << backend;
    minim.hesse(/*verbose=*/0);
    Result &res = results[backend];
    res.r = m.r->getVal();
    res.rError = m.r->getError();
    res.index = m.index->getIndex();
    res.nllImprovement = nllStart - nll->getVal();
    std::cout << backend << ": r = " << res.r << " +/- " << res.rError << ", index " << res.index
              << ", NLL improvement " << res.nllImprovement << std::endl;
  }
  for (auto const &backend : kBackends) {
    Result const &ref = results["combine"];
    Result const &res = results[backend];
    EXPECT_EQ(res.index, ref.index) << backend;
    EXPECT_NEAR(res.r, ref.r, 2e-3) << backend;
    EXPECT_NEAR(res.rError, ref.rError, 0.02 * ref.rError) << backend;
    EXPECT_NEAR(res.nllImprovement, ref.nllImprovement, 1e-2) << backend;
  }
}

#endif  // ROOT_VERSION_CODE >= ROOT_VERSION(6, 37, 0)
