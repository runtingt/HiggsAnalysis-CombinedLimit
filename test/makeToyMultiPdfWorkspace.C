// Builds a small workspace for testing the discrete profiling of a
// RooMultiPdf on all likelihood backends of Combine, including channel
// masking and the analytic gradient of the "codegen" backend.
//
// Two channels share the signal strength "r" and the nuisance "lumi":
//   - channel A: Gaussian signal plus a RooMultiPdf background with three
//     candidate functions (exponential, Chebychev, polynomial);
//   - channel B: Gaussian signal plus an exponential background, and a
//     nuisance "bnorm" of its own.
// The constraint terms are multiplied into every channel like text2workspace
// does. Only pdf classes with codegen support are used, which is why the
// RooMultiPdf datacard of the CI (with its RooDoubleCBFast signal) is not
// enough for these tests.
//
// Usage: root -l -b -q 'makeToyMultiPdfWorkspace.C("toy_multipdf.root")'

void makeToyMultiPdfWorkspace(const char *outputFile = "toy_multipdf.root") {
  gSystem->Load("libHiggsAnalysisCombinedLimit");

  RooWorkspace w("w");
  w.factory("x[100,180]");
  w.factory("r[1,-5,5]");
  w.factory("lumi[0,-5,5]");
  w.factory("lumi_In[0]");
  w.factory("Gaussian::lumi_Pdf(lumi, lumi_In, 1)");
  w.factory("bnorm[0,-5,5]");
  w.factory("bnorm_In[0]");
  w.factory("Gaussian::bnorm_Pdf(bnorm, bnorm_In, 1)");
  w.var("lumi_In")->setConstant(true);
  w.var("bnorm_In")->setConstant(true);

  // Channel A
  w.factory("Gaussian::sigA(x, 125, 2.0)");
  w.factory("Exponential::bkgA_exp(x, slopeA[-0.03,-0.2,0.0])");
  w.factory("Chebychev::bkgA_cheb(x, {cA1[-0.3,-1,1], cA2[0.1,-1,1]})");
  w.factory("Polynomial::bkgA_pol(x, {pA1[-0.003,-0.005,0.005]})");
  RooCategory index("pdf_index_A", "");
  RooMultiPdf multiPdf("bkgA", "", index, RooArgList(*w.pdf("bkgA_exp"), *w.pdf("bkgA_cheb"), *w.pdf("bkgA_pol")));
  ProcessNormalization nsigA("nsigA", "", 60.);
  nsigA.addLogNormal(1.05, *w.var("lumi"));
  nsigA.addOtherFactor(*w.var("r"));
  w.factory("nbkgA[2000,0,10000]");
  RooAddPdf pdfA("pdfA", "", RooArgList(*w.pdf("sigA"), multiPdf), RooArgList(nsigA, *w.var("nbkgA")));
  RooProdPdf pdfA_nuis("pdfA_nuis", "", RooArgList(pdfA, *w.pdf("lumi_Pdf"), *w.pdf("bnorm_Pdf")));

  // Channel B
  w.factory("Gaussian::sigB(x, 125, 3.0)");
  w.factory("Exponential::bkgB(x, slopeB[-0.02,-0.2,0.0])");
  ProcessNormalization nsigB("nsigB", "", 30.);
  nsigB.addLogNormal(1.05, *w.var("lumi"));
  nsigB.addLogNormal(1.10, *w.var("bnorm"));
  nsigB.addOtherFactor(*w.var("r"));
  w.factory("nbkgB[1000,0,10000]");
  RooAddPdf pdfB("pdfB", "", RooArgList(*w.pdf("sigB"), *w.pdf("bkgB")), RooArgList(nsigB, *w.var("nbkgB")));
  RooProdPdf pdfB_nuis("pdfB_nuis", "", RooArgList(pdfB, *w.pdf("lumi_Pdf"), *w.pdf("bnorm_Pdf")));

  RooCategory channel("CMS_channel", "");
  channel.defineType("A");
  channel.defineType("B");
  RooSimultaneousOpt sim("model_s", "", channel);
  sim.addPdf(pdfA_nuis, "A");
  sim.addPdf(pdfB_nuis, "B");
  w.import(sim, RooFit::RecycleConflictNodes(), RooFit::Silence());

  // A dataset generated with the exponential background, then the start
  // index is set to the Chebychev function so that the discrete profiling has
  // something to do.
  RooArgSet observables(*w.var("x"), *w.cat("CMS_channel"));
  RooRandom::randomGenerator()->SetSeed(4357);
  w.cat("pdf_index_A")->setIndex(0);
  std::unique_ptr<RooDataSet> data{w.pdf("model_s")->generate(observables, RooFit::Extended(), RooFit::Name("data_obs"))};
  w.import(*data);
  w.cat("pdf_index_A")->setIndex(1);

  RooStats::ModelConfig mc("ModelConfig", &w);
  mc.SetPdf(*w.pdf("model_s"));
  mc.SetObservables(observables);
  mc.SetParametersOfInterest(*w.var("r"));
  mc.SetNuisanceParameters(RooArgSet(*w.var("lumi"), *w.var("bnorm")));
  mc.SetGlobalObservables(RooArgSet(*w.var("lumi_In"), *w.var("bnorm_In")));
  w.import(mc);

  // Combine takes the discrete parameters from this set.
  RooArgSet discreteParams("discreteParams");
  discreteParams.add(*w.cat("pdf_index_A"));
  w.import(discreteParams, "discreteParams");

  w.writeToFile(outputFile);
  std::cout << "Wrote " << outputFile << " with " << data->numEntries() << " events" << std::endl;
}
