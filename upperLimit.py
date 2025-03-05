from ROOT import RDataFrame, TChain, TH1D, TCanvas, TFile, RooStats
import ROOT
import numpy as np
import yaml

# columns: { "DecR", "E", "E3H", "EPi", "GenPt", "L", "M", "Pt", "Px", "Px3H", "PxPi", "Py", "Py3H", "PyPi", "Pz", "Pz3H", "PzPi", "ct", "fCentralityFT0A", "fCentralityFT0C", "fCentralityFT0M", "fDca3H", "fDcaPi", "fDcaV0Daug", "fEta3H", "fEtaPi", "fFlags", "fGenEta", "fGenPhi", "fGenPt", "fGenPt3H", "fGenXDecVtx", "fGenYDecVtx", "fGenZDecVtx", "fITSclusterSizes3H", "fITSclusterSizesPi", "fIsMatter", "fIsReco", "fIsSignal", "fMassTrTOF", "fNSigma3H", "fNTPCclus3H", "fNTPCclusPi", "P3H", "fPhi3H", "fPhiPi", "fPt3H", "fPtPi", "fSurvivedEventSelection", "fTPCchi3H", "fTPCmom3H", "fTPCmomPi", "fTPCsignal3H", "fTPCsignalPi", "fXDecVtx", "fXPrimVtx", "fYDecVtx", "fYPrimVtx", "fZDecVtx", "fZPrimVtx", "genMatter", "mTOF3H_centred", "signedP3H" }

with open("configurations.yaml", "r") as stream:
  configurations = yaml.safe_load(stream)

ROOT.EnableImplicitMT()
eventCuts = f"fCentralityFT0C < {configurations['centmax']} && fCentralityFT0C > {configurations['centmin']}"
ptcuts = f"Pt > {configurations['ptmin']} && Pt < {configurations['ptmax']}"

anResultsData = ROOT.TFile("data/AnalysisResults.root")
counterHist = anResultsData.Get("lnn-reco-task/hCentFT0C")
file = ROOT.TFile.Open("bw/H3L_BWFit.root")
bw_0_10 = file.Get("BlastWave_H3L_0_10")
bw_10_30 = file.Get("BlastWave_H3L_10_30")
bw_30_50 = file.Get("BlastWave_H3L_30_50")
bw_shapes = [bw_0_10, bw_10_30, bw_30_50]
br = 0.25
delta_pt = configurations['ptmax'] - configurations['ptmin']
rapidity_window = 2
mat_antimat_fac = 2
n_ev = counterHist.Integral(counterHist.GetXaxis().FindBin(configurations['centmin'] + 0.001), counterHist.GetXaxis().FindBin(configurations['centmax'] - 0.001))

fileData = TFile('data/AO2D.root')
fileMC = TFile('mc/AO2D.root')
chainData = TChain('O2lnncands')
chainMC = TChain('O2mclnncands')
for key in fileData.GetListOfKeys() :
    keyName = key.GetName()
    if 'DF_' in keyName :
      chainData.Add(f'data/AO2D.root/{keyName}/{chainData.GetName()}')
for key in fileMC.GetListOfKeys() :
    keyName = key.GetName()
    if 'DF_' in keyName :
      chainMC.Add(f'mc/AO2D.root/{keyName}/{chainMC.GetName()}')

dfData = ROOT.RDataFrame(chainData).Define("Px3H", "fPt3H * std::cos(fPhi3H)").Define("Py3H", "fPt3H * std::sin(fPhi3H)").Define("Pz3H", "fPt3H * std::sinh(fEta3H)").Define("P3H", "std::hypot(fPt3H, Pz3H)").Define("mTOF3H_centred", "fMassTrTOF - 2.80892113298 * 2.80892113298").Define("signedP3H", "(fIsMatter * 2 - 1) * fTPCmom3H").Define("E3H", "std::hypot(P3H, 2.80892113298)").Define("EP", "std::hypot(P3H, 0.93827208816)").Define("EPiK0s", "std::hypot(P3H, 0.139570)").Define("PxPi", "fPtPi * std::cos(fPhiPi)").Define("PyPi", "fPtPi * std::sin(fPhiPi)").Define("PzPi", "fPtPi * std::sinh(fEtaPi)").Define("PPi", "std::hypot(fPtPi, PzPi)").Define("EPi", "std::hypot(PPi, 0.139570)").Define("Px", "Px3H + PxPi").Define("Py", "Py3H + PyPi").Define("Pz", "Pz3H + PzPi").Define("E", "E3H + EPi").Define("Pt", "std::hypot(Px, Py)").Define("P", "std::hypot(Pt, Pz)").Define("M", "std::sqrt(E * E - P * P)").Define("L", "std::hypot(fXDecVtx, fYDecVtx, fZDecVtx)").Define("ct", "L * 2.991 / std::hypot(Pt, Pz)").Define("DecR", "std::hypot(fXDecVtx + fXPrimVtx, fYDecVtx + fYPrimVtx)").Define("cosPA", "(Px * fXDecVtx + Py * fYDecVtx + Pz * fZDecVtx) / (P * L)").Define("MLambda0", "std::sqrt((EPi + EP) * (EPi + EP) - P * P)").Define("MK0s", "std::sqrt((EPi + EPiK0s) * (EPi + EPiK0s) - P * P)").Define("PxPos", "fIsMatter ? Px3H : PxPi").Define("PyPos", "fIsMatter ? Py3H : PyPi").Define("PzPos", "fIsMatter ? Pz3H : PzPi").Define("PxNeg", "!fIsMatter ? Px3H : PxPi").Define("PyNeg", "!fIsMatter ? Py3H : PyPi").Define("PzNeg", "!fIsMatter ? Pz3H : PzPi").Define("alpha", "(PxPos * Px + PyPos * Py + PzPos * Pz - PxNeg * Px - PyNeg * Py - PzNeg * Pz) / (PxPos * Px + PyPos * Py + PzPos * Pz + PxNeg * Px + PyNeg * Py + PzNeg * Pz)").Define("qt", "std::sqrt(std::pow(fIsMatter ? P3H : PPi ,2) - std::pow((Px * PxPos + Py * PyPos + Pz * PzPos) / P,2))") \
    .Filter(eventCuts)

dfMC = ROOT.RDataFrame(chainMC).Define("Px3H", "fPt3H * std::cos(fPhi3H)").Define("Py3H", "fPt3H * std::sin(fPhi3H)").Define("Pz3H", "fPt3H * std::sinh(fEta3H)").Define("P3H", "std::hypot(fPt3H, Pz3H)").Define("mTOF3H_centred", "fMassTrTOF - 2.80892113298 * 2.80892113298").Define("signedP3H", "(fIsMatter * 2 - 1) * fTPCmom3H").Define("E3H", "std::hypot(P3H, 2.80892113298)").Define("EP", "std::hypot(P3H, 0.93827208816)").Define("EPiK0s", "std::hypot(P3H, 0.139570)").Define("PxPi", "fPtPi * std::cos(fPhiPi)").Define("PyPi", "fPtPi * std::sin(fPhiPi)").Define("PzPi", "fPtPi * std::sinh(fEtaPi)").Define("PPi", "std::hypot(fPtPi, PzPi)").Define("EPi", "std::hypot(PPi, 0.139570)").Define("Px", "Px3H + PxPi").Define("Py", "Py3H + PyPi").Define("Pz", "Pz3H + PzPi").Define("E", "E3H + EPi").Define("Pt", "std::hypot(Px, Py)").Define("P", "std::hypot(Pt, Pz)").Define("M", "std::sqrt(E * E - P * P)").Define("L", "std::hypot(fXDecVtx, fYDecVtx, fZDecVtx)").Define("ct", "L * 2.991 / std::hypot(Pt, Pz)").Define("DecR", "std::hypot(fXDecVtx + fXPrimVtx, fYDecVtx + fYPrimVtx)").Define("cosPA", "(Px * fXDecVtx + Py * fYDecVtx  + Pz * fZDecVtx) / (P * L)").Define("MLambda0", "std::sqrt((EPi + EP) * (EPi + EP) - P * P)").Define("MK0s", "std::sqrt((EPi + EPiK0s) * (EPi + EPiK0s) - P * P)").Define("PxPos", "fIsMatter ? Px3H : PxPi").Define("PyPos", "fIsMatter ? Py3H : PyPi").Define("PzPos", "fIsMatter ? Pz3H : PzPi").Define("PxNeg", "!fIsMatter ? Px3H : PxPi").Define("PyNeg", "!fIsMatter ? Py3H : PyPi").Define("PzNeg", "!fIsMatter ? Pz3H : PzPi").Define("alpha", "(PxPos * Px + PyPos * Py + PzPos * Pz - PxNeg * Px - PyNeg * Py - PzNeg * Pz) / (PxPos * Px + PyPos * Py + PzPos * Pz + PxNeg * Px + PyNeg * Py + PzNeg * Pz)").Define("qt", "std::sqrt(std::pow(fIsMatter ? P3H : PPi ,2) - std::pow((Px * PxPos + Py * PyPos + Pz * PzPos) / P,2))") \
  .Define("genMatter", "fGenPt > 0").Define("GenPt", "std::abs(fGenPt)")

dfRecMC = dfMC.Filter("fIsReco").Filter(eventCuts)
hGenMC = dfMC.Filter(eventCuts).Histo2D(("hMC", "MC", 100, 0, 10.,3, np.array([0, 10, 30, 50], dtype=float)), "GenPt", "fCentralityFT0C")

outFile = TFile("output.root", "RECREATE")

hRecMC = []
hPID3H = []
hRpt = []
hCosPAptMC = []
hCosPApt = []
hMassPtMC = []
hMassLambdaPtMC = []
hMassPt = []
hDCAv0DaughPtMC = []
hMassLambdaPt = []
hDCAv0DaughPt = []
hMassK0sPt = []
hDCAptPi = []
hDCAptPiMC = []
hDCApt3H = []
hDCApt3HMC = []
hArmenterosPodolanski = []
hArmenterosPodolanskiMC = []
hMassFit = []
hMassFitMC = []

upperlimits = []

for dcapi in configurations['dcapi_variations']:
    for dcav0daugters in configurations['dcav0daughters_variations']:
       for cospa in configurations['cospa_variations']:
          cuts = f"fNTPCclusPi > 100 && fNTPCclus3H > 100 && (MLambda0 < 1.1 || MLambda0 > 1.13) && (MK0s < 0.47 || MK0s > 0.53) && std::abs(fDcaPi) > {dcapi} && cosPA > {cospa} && fDcaV0Daug < {dcav0daugters} && std::abs(alpha) < 0.94 && qt > 0.005"
          dfDataFiltered = dfData.Filter(cuts)
          dfRecMCFiltered = dfRecMC.Filter(cuts)
          hRecMC.append(dfRecMCFiltered.Histo2D(("hData", "Data", 100, 0, 10.,3, np.array([0, 10, 30, 50], dtype=float)), "Pt", "fCentralityFT0C"))
          hPID3H.append(dfDataFiltered.Histo3D(("hPID3H", ";#it{p} (GeV/#it{c}) / q;TPC n#sigma; TOF m^{2} - m^{2}_{PDG}", 50, -10, 10., 50, -5, 5, 60, -3, 3), "signedP3H", "fNSigma3H", "mTOF3H_centred"))
          hRpt.append(dfDataFiltered.Histo2D(("hRpt", ";#it{p}_{T} (GeV/#it{c});R", 50, 0, 10., 100, 0, 35), "Pt", "DecR"))
          hCosPAptMC.append(dfRecMCFiltered.Histo2D(("hCosPAptMC", ";#it{p}_{T} (GeV/#it{c});cos#it{#theta}_{P}", 20, 0, 10., 1000, 0.9, 1), "Pt", "cosPA"))
          hCosPApt.append(dfDataFiltered.Histo2D(("hCosPApt", ";#it{p}_{T} (GeV/#it{c});cos#it{#theta}_{P}", 20, 0, 10., 1000, 0.9, 1), "Pt", "cosPA"))
          hMassPtMC.append(dfRecMCFiltered.Histo2D(("hMassPtMC", ";#it{p}_{T} (GeV/#it{c});m_{^{3}H + #pi^{-} + c.c.} (GeV/#it{c}^{2})", 20, 0, 10., 50, 2.97, 3.02), "Pt", "M"))
          hMassLambdaPtMC.append(dfRecMCFiltered.Histo3D(("hMassLambdaPtMC", ";#it{p}_{T} (GeV/#it{c});m_{p + #pi^{-} + c.c.} (GeV/#it{c}^{2});m_{^{3}H + #pi^{-} + c.c.} (GeV/#it{c}^{2})", 20, 0, 10., 425, 1.075, 1.5, 80, 2.95, 3.03), "Pt", "MLambda0", "M"))
          hMassPt.append(dfDataFiltered.Histo2D(("hMassPt", ";#it{p}_{T} (GeV/#it{c});m_{^{3}H + #pi^{-} + c.c.} (GeV/#it{c})^{2}", 20, 0, 10., 50, 2.97, 3.02), "Pt", "M"))
          hDCAv0DaughPtMC.append(dfRecMCFiltered.Histo2D(("hDCAv0DaughPtMC", ";#it{p}_{T} (GeV/#it{c});V0 daughters DCA (cm)", 20, 0, 10., 2000, 0, 1), "Pt", "fDcaV0Daug"))
          hMassLambdaPt.append(dfDataFiltered.Histo3D(("hMassLambdaPt", ";#it{p}_{T} (GeV/#it{c});m_{p + #pi^{-} + c.c.} (GeV/#it{c}^{2});m_{^{3}H + #pi^{-} + c.c.} (GeV/#it{c}^{2})", 20, 0, 10., 425, 1.075, 1.5, 80, 2.95, 3.03), "Pt", "MLambda0", "M"))
          hDCAv0DaughPt.append(dfDataFiltered.Histo2D(("hDCAv0DaughPt", ";#it{p}_{T} (GeV/#it{c});V0 daughters DCA (cm)", 20, 0, 10., 2000, 0, 1), "Pt", "fDcaV0Daug"))
          hMassK0sPt.append(dfDataFiltered.Histo2D(("hMassK0sPt", ";#it{p}_{T} (GeV/#it{c});m_{#pi^{+} + #pi^{-}} (GeV/#it{c})^{2}", 20, 0, 10., 80, 0.455, 0.535), "Pt", "MK0s"))
          hDCAptPi.append(dfDataFiltered.Histo2D(("hDCAptPi", ";#it{p}_{T} (GeV/#it{c});#pi DCA (cm)", 20, 0, 10., 2000, -1, 1), "Pt", "fDcaPi"))
          hDCAptPiMC.append(dfRecMCFiltered.Histo2D(("hDCAptPiMC", ";#it{p}_{T} (GeV/#it{c});#pi DCA (cm)", 20, 0, 10., 2000, -1, 1), "Pt", "fDcaPi"))
          hDCApt3H.append(dfDataFiltered.Histo2D(("hDCApt3H", ";#it{p}_{T} (GeV/#it{c});^{3}H DCA (cm)", 20, 0, 10., 2000, -1, 1), "Pt", "fDca3H"))
          hDCApt3HMC.append(dfRecMCFiltered.Histo2D(("hDCApt3HMC", ";#it{p}_{T} (GeV/#it{c});^{3}H DCA (cm)", 20, 0, 10., 2000, -1, 1), "Pt", "fDca3H"))

          hArmenterosPodolanski.append(dfDataFiltered.Histo2D(("hArmenterosPodolanski", ";#alpha;#it{q}_{T} (GeV/#it{c})", 200, -1, 1, 500, 0, 0.5), "alpha", "qt"))
          hArmenterosPodolanskiMC.append(dfRecMCFiltered.Histo2D(("hArmenterosPodolanskiMC", ";#alpha;#it{q}_{T} (GeV/#it{c})", 200, -1, 1, 500, 0, 0.5), "alpha", "qt"))

          hMassFit.append(dfDataFiltered.Filter(ptcuts).Histo1D(("hMassFit", ";m_{^{3}H + #pi^{-} + c.c.} (GeV/#it{c}^{2})", 45, 2.975, 3.02), "M"))
          hMassFitMC.append(dfRecMCFiltered.Filter(ptcuts).Histo1D(("hMassFitMC", ";m_{^{3}H + #pi^{-} + c.c.} (GeV/#it{c}^{2})", 45, 2.975, 3.02), "M"))


upperLimitsMassHists = []
for iTrial in range(len(hMassFit)):
  outFile.mkdir(f"trial{iTrial}")
  outFile.cd(f"trial{iTrial}")

  effMC = hRecMC[iTrial].Clone(f"effMC{iTrial}")
  effMC.Divide(hRecMC[iTrial].GetPtr(), hGenMC.GetPtr(), 1, 1, "B")

  effBinStart = effMC.GetXaxis().FindBin(configurations['ptmin'] + 0.0001)
  effBinEnd = effMC.GetXaxis().FindBin(configurations['ptmax'] - 0.0001)

  efficiencies = [0, 0, 0]
  for k in range(3):
    bw_cent = bw_shapes[k]
    eff_cent = effMC.ProjectionX(f"eff_{k}", k + 1, k + 1)

    sum_eff = 0
    sum_weights = 0
    for j in range(effBinStart, effBinEnd + 1):
      bw_weight = bw_cent.Integral(effMC.GetXaxis().GetBinLowEdge(j), effMC.GetXaxis().GetBinUpEdge(j))
      sum_eff += eff_cent.GetBinContent(j) * bw_weight
      sum_weights += bw_weight

    efficiencies[k] = sum_eff / sum_weights

  centrality_weights = [counterHist.Integral(counterHist.GetXaxis().FindBin(0.001), counterHist.GetXaxis().FindBin(9.999)), counterHist.Integral(counterHist.GetXaxis().FindBin(10.001), counterHist.GetXaxis().FindBin(29.999)), counterHist.Integral(counterHist.GetXaxis().FindBin(30.001), counterHist.GetXaxis().FindBin(49.999))]
  yield_weights = [bw_0_10.Integral(configurations['ptmin'], configurations['ptmax']), bw_10_30.Integral(configurations['ptmin'], configurations['ptmax']), bw_30_50.Integral(configurations['ptmin'], configurations['ptmax'])]
  print('efficiencies: ', efficiencies)
  print('centrality weights: ', centrality_weights)
  print('yield weights: ', yield_weights)
  eff_reweighted = sum([a*b*c for a,b,c in zip(efficiencies, centrality_weights, yield_weights)]) / sum([a*b for a,b in zip(centrality_weights, yield_weights)])

  mass = ROOT.RooRealVar("mass", "mass", 2.991, hMassFit[iTrial].GetXaxis().GetXmin(), hMassFit[iTrial].GetXaxis().GetXmax())
  mu = ROOT.RooRealVar("mu", "#mu", 2.991, hMassFit[iTrial].GetXaxis().GetXmin(), hMassFit[iTrial].GetXaxis().GetXmax())
  sigma = ROOT.RooRealVar("sigma", "#sigma", 0.002, 0.001, 0.01)
  a1 = ROOT.RooRealVar("a1", "a1", 0., -1., 1.)
  bkg = ROOT.RooPolynomial("bkg", "bkg", mass, ROOT.RooArgList())
  alpha0 = ROOT.RooRealVar("alpha0", "#alpha_{0}", 1., 0., 5.)
  n0 = ROOT.RooRealVar("n0", "n_{0}", 1., 0., 10.)
  alpha1 = ROOT.RooRealVar("alpha1", "#alpha_{1}", 1., 0., 5.)
  n1 = ROOT.RooRealVar("n1", "n_1", 1., 0., 10.)
  cb = ROOT.RooCrystalBall("cb", "cb", mass, mu, sigma, alpha0, n0, alpha1, n1)
  signal = ROOT.RooGaussian("signal", "signal", mass, mu, sigma)
  n_sig = ROOT.RooRealVar("n_sig", "n_sig", 1000, -4000, 4000)
  n_bkg = ROOT.RooRealVar("n_bkg", "n_bkg", 1000, 0, 1e7)
  fulltpl = ROOT.RooAddPdf("fulltpl", "fulltpl", ROOT.RooArgList(cb, bkg), ROOT.RooArgList(n_sig, n_bkg))

  rooHistMC = ROOT.RooDataHist(f"rooHistMC{iTrial}", "rooHistMC", ROOT.RooArgList(mass), hMassFitMC[iTrial].GetPtr())
  cb.fitTo(rooHistMC)
  frameMC = mass.frame()
  rooHistMC.plotOn(frameMC)
  cb.plotOn(frameMC)
  cb.paramOn(frameMC)
  alpha0.setConstant(True)
  n0.setConstant(True)
  alpha1.setConstant(True)
  n1.setConstant(True)
  sigma.setConstant(True)
  mu.setConstant(True)

  rooHistData = ROOT.RooDataHist(f"rooHistData{iTrial}", "rooHistData", ROOT.RooArgList(mass), hMassFit[iTrial].GetPtr())
  bkg.fitTo(rooHistData)
  frame = mass.frame()
  rooHistData.plotOn(frame)
  bkg.plotOn(frame)
  chi2Val = frame.chiSquare()
  label = ROOT.TPaveText(0.6, 0.8, 0.9, 0.85, "NDC")
  label.AddText(f"#chi^{{2}} = {chi2Val:.2f}")
  label.SetFillColor(0)
  label.SetBorderSize(0)
  frame.addObject(label)

  rebinUpperLimit = 3
  upperLimitsMassHists.append(ROOT.TH1D(f"upperLimitsMass{iTrial}", "; m (GeV/#it{c}); Upper limit", int(hMassFit[iTrial].GetNbinsX() / rebinUpperLimit), hMassFit[iTrial].GetXaxis().GetXmin(), hMassFit[iTrial].GetXaxis().GetXmax()))
  for iMass in range(upperLimitsMassHists[iTrial].GetNbinsX()):
    # Create a workspace
    workspace = ROOT.RooWorkspace(f"workspace{iTrial}_{iMass}")
    workspace.Import(rooHistData)
    workspace.Import(fulltpl)

    workspace.Print()
    # Set the confidence level
    confidence_level = 0.95

    # Create a ModelConfig
    model_config = RooStats.ModelConfig(f"model_config{iTrial}_{iMass}", workspace)
    model_config.SetPdf(fulltpl)
    model_config.SetParametersOfInterest(n_sig)
    model_config.SetObservables(ROOT.RooArgSet(mass))
    model_config.SetNuisanceParameters(ROOT.RooArgSet())
    model_config.SetNuisanceParameters(ROOT.RooArgSet(n_bkg))
    oldValue = n_sig.getVal()

    # Create a null (background-only) model
    null_model = RooStats.ModelConfig(f"null_model{iTrial}_{iMass}", workspace)
    null_model.SetPdf(fulltpl)
    n_sig.setVal(0)
    workspace.var("n_sig").setVal(0)
    null_model.SetParametersOfInterest(n_sig)
    null_model.SetObservables(ROOT.RooArgSet(mass))
    null_model.SetNuisanceParameters(ROOT.RooArgSet(mu))
    null_model.SetSnapshot(ROOT.RooArgSet(n_sig))
    null_model.GetParametersOfInterest().first().setVal(oldValue)
    print("Null model snapshot:")
    null_model.GetSnapshot().Print("v")

    # Perform the upper limit calculation
    calculator = RooStats.AsymptoticCalculator(rooHistData, null_model, model_config)
    calculator.SetOneSided(True)
    mu.setVal(upperLimitsMassHists[iTrial].GetXaxis().GetBinCenter(iMass + 1))
    print (f"Calculating upper limit for mass {mu.getVal()}, trial {iTrial}")
    hypo_test_result = calculator.GetHypoTest()
    if not hypo_test_result:
      print("Erro: HypoTest result is None!")

    # Extract the p-value
    invereter = RooStats.HypoTestInverter(calculator)
    invereter.SetConfidenceLevel(confidence_level)
    invereter.UseCLs(True)
    invereter.SetVerbose(True)
    result = invereter.GetInterval()

    # Print the result
    upper_limit = result.UpperLimit()
    error_up = result.UpperLimitEstimatedError()
    error_lower = result.LowerLimitEstimatedError()

    ### upper limit calculation
    upper_limit_corrected = upper_limit / (eff_reweighted * br * rapidity_window * mat_antimat_fac * n_ev * delta_pt)
    upperLimitsMassHists[iTrial].SetBinContent(iMass + 1, upper_limit_corrected)


  frame.Write(f"fitData{iTrial}")
  frameMC.Write(f"fitMC{iTrial}")
  effMC.Write()
  hRpt[iTrial].Write()
  hPID3H[iTrial].Write()
  hCosPAptMC[iTrial].Write()
  hCosPApt[iTrial].Write()
  hMassPtMC[iTrial].Write()
  hMassPt[iTrial].Write()
  hMassLambdaPtMC[iTrial].Write()
  hMassLambdaPt[iTrial].Write()
  hMassK0sPt[iTrial].Write()
  hDCAptPi[iTrial].Write()
  hDCAptPiMC[iTrial].Write()
  hDCApt3H[iTrial].Write()
  hDCApt3HMC[iTrial].Write()
  hDCAv0DaughPt[iTrial].Write()
  hDCAv0DaughPtMC[iTrial].Write()
  hArmenterosPodolanski[iTrial].Write()
  hArmenterosPodolanskiMC[iTrial].Write()
  hMassFit[iTrial].Write()
  upperLimitsMassHists[-1].Write()

outFile.cd()
topLimit = upperLimitsMassHists[0].Clone("topLimit")
for upperLimit in upperLimitsMassHists:
   for i in range(upperLimit.GetNbinsX()):
     if upperLimit.GetBinContent(i + 1) > topLimit.GetBinContent(i + 1):
       topLimit.SetBinContent(i + 1, upperLimit.GetBinContent(i + 1))
topLimit.Write()
shmPred = ROOT.TF1("shmPred", "469.23 * TMath::Exp(-x/0.155)", 2.975, 3.02)
shmPred.Write()
variationCanvas = ROOT.TCanvas("variationCanvas", "variationCanvas", 800, 600)
variationCanvas.DrawFrame(2.975, 1.e-7, 3.02, 2.5e-6, "; m (GeV/#it{c}^{2}); Upper limit")
for upperLimit in upperLimitsMassHists:
  upperLimit.Draw("plc l same")
topLimit.SetLineColor(ROOT.kBlack)
topLimit.SetLineWidth(2)
topLimit.SetMarkerStyle(20)
topLimit.SetMarkerSize(0.5)
topLimit.Draw("lp same")
variationCanvas.Write()

outFile.Close()
