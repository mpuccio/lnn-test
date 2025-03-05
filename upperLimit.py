from ROOT import RDataFrame, TChain, TH1D, TCanvas, TFile, RooStats
import ROOT
import numpy as np

# columns: { "DecR", "E", "E3H", "EPi", "GenPt", "L", "M", "Pt", "Px", "Px3H", "PxPi", "Py", "Py3H", "PyPi", "Pz", "Pz3H", "PzPi", "ct", "fCentralityFT0A", "fCentralityFT0C", "fCentralityFT0M", "fDca3H", "fDcaPi", "fDcaV0Daug", "fEta3H", "fEtaPi", "fFlags", "fGenEta", "fGenPhi", "fGenPt", "fGenPt3H", "fGenXDecVtx", "fGenYDecVtx", "fGenZDecVtx", "fITSclusterSizes3H", "fITSclusterSizesPi", "fIsMatter", "fIsReco", "fIsSignal", "fMassTrTOF", "fNSigma3H", "fNTPCclus3H", "fNTPCclusPi", "P3H", "fPhi3H", "fPhiPi", "fPt3H", "fPtPi", "fSurvivedEventSelection", "fTPCchi3H", "fTPCmom3H", "fTPCmomPi", "fTPCsignal3H", "fTPCsignalPi", "fXDecVtx", "fXPrimVtx", "fYDecVtx", "fYPrimVtx", "fZDecVtx", "fZPrimVtx", "genMatter", "mTOF3H_centred", "signedP3H" }

ROOT.EnableImplicitMT()
eventCuts = "fCentralityFT0C < 50 && fCentralityFT0C > 0"
cuts = "fNTPCclusPi > 60 && fNTPCclus3H > 100 && (MLambda0 < 1.1 || MLambda0 > 1.13) && (MK0s < 0.47 || MK0s > 0.53) && std::abs(fDcaPi) > 0.2 && cosPA > 0.9995 && fDcaV0Daug < 0.2 && std::abs(alpha) < 0.94 && qt > 0.005"
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
    .Filter(cuts).Filter(eventCuts)

dfMC = ROOT.RDataFrame(chainMC).Define("Px3H", "fPt3H * std::cos(fPhi3H)").Define("Py3H", "fPt3H * std::sin(fPhi3H)").Define("Pz3H", "fPt3H * std::sinh(fEta3H)").Define("P3H", "std::hypot(fPt3H, Pz3H)").Define("mTOF3H_centred", "fMassTrTOF - 2.80892113298 * 2.80892113298").Define("signedP3H", "(fIsMatter * 2 - 1) * fTPCmom3H").Define("E3H", "std::hypot(P3H, 2.80892113298)").Define("EP", "std::hypot(P3H, 0.93827208816)").Define("EPiK0s", "std::hypot(P3H, 0.139570)").Define("PxPi", "fPtPi * std::cos(fPhiPi)").Define("PyPi", "fPtPi * std::sin(fPhiPi)").Define("PzPi", "fPtPi * std::sinh(fEtaPi)").Define("PPi", "std::hypot(fPtPi, PzPi)").Define("EPi", "std::hypot(PPi, 0.139570)").Define("Px", "Px3H + PxPi").Define("Py", "Py3H + PyPi").Define("Pz", "Pz3H + PzPi").Define("E", "E3H + EPi").Define("Pt", "std::hypot(Px, Py)").Define("P", "std::hypot(Pt, Pz)").Define("M", "std::sqrt(E * E - P * P)").Define("L", "std::hypot(fXDecVtx, fYDecVtx, fZDecVtx)").Define("ct", "L * 2.991 / std::hypot(Pt, Pz)").Define("DecR", "std::hypot(fXDecVtx + fXPrimVtx, fYDecVtx + fYPrimVtx)").Define("cosPA", "(Px * fXDecVtx + Py * fYDecVtx  + Pz * fZDecVtx) / (P * L)").Define("MLambda0", "std::sqrt((EPi + EP) * (EPi + EP) - P * P)").Define("MK0s", "std::sqrt((EPi + EPiK0s) * (EPi + EPiK0s) - P * P)").Define("PxPos", "fIsMatter ? Px3H : PxPi").Define("PyPos", "fIsMatter ? Py3H : PyPi").Define("PzPos", "fIsMatter ? Pz3H : PzPi").Define("PxNeg", "!fIsMatter ? Px3H : PxPi").Define("PyNeg", "!fIsMatter ? Py3H : PyPi").Define("PzNeg", "!fIsMatter ? Pz3H : PzPi").Define("alpha", "(PxPos * Px + PyPos * Py + PzPos * Pz - PxNeg * Px - PyNeg * Py - PzNeg * Pz) / (PxPos * Px + PyPos * Py + PzPos * Pz + PxNeg * Px + PyNeg * Py + PzNeg * Pz)").Define("qt", "std::sqrt(std::pow(fIsMatter ? P3H : PPi ,2) - std::pow((Px * PxPos + Py * PyPos + Pz * PzPos) / P,2))") \
  .Define("genMatter", "fGenPt > 0").Define("GenPt", "std::abs(fGenPt)")

dfRecMC = dfMC.Filter("fIsReco").Filter(cuts).Filter(eventCuts)

print(dfMC.GetColumnNames())
outFile = TFile("output.root", "RECREATE")
hRecMC = dfRecMC.Histo1D(("hData", "Data", 100, 0, 10.), "Pt")
hGenMC = dfMC.Filter(eventCuts).Histo1D(("hMC", "MC", 100, 0, 10.), "GenPt")
hPID3H = dfData.Histo3D(("hPID3H", ";#it{p} (GeV/#it{c}) / q;TPC n#sigma; TOF m^{2} - m^{2}_{PDG}", 50, -10, 10., 50, -5, 5, 60, -3, 3), "signedP3H", "fNSigma3H", "mTOF3H_centred")
hRpt = dfData.Histo2D(("hRpt", ";#it{p}_{T} (GeV/#it{c});R", 50, 0, 10., 100, 0, 35), "Pt", "DecR")
hCosPAptMC = dfRecMC.Histo2D(("hCosPAptMC", ";#it{p}_{T} (GeV/#it{c});cos#it{#theta}_{P}", 20, 0, 10., 1000, 0.9, 1), "Pt", "cosPA")
hCosPApt = dfData.Histo2D(("hCosPApt", ";#it{p}_{T} (GeV/#it{c});cos#it{#theta}_{P}", 20, 0, 10., 1000, 0.9, 1), "Pt", "cosPA")
hMassPtMC = dfRecMC.Histo2D(("hMassPtMC", ";#it{p}_{T} (GeV/#it{c});m_{^{3}H + #pi^{-} + c.c.} (GeV/#it{c}^{2})", 20, 0, 10., 50, 2.97, 3.02), "Pt", "M")
hMassLambdaPtMC = dfRecMC.Histo3D(("hMassLambdaPtMC", ";#it{p}_{T} (GeV/#it{c});m_{p + #pi^{-} + c.c.} (GeV/#it{c}^{2});m_{^{3}H + #pi^{-} + c.c.} (GeV/#it{c}^{2})", 20, 0, 10., 425, 1.075, 1.5, 80, 2.95, 3.03), "Pt", "MLambda0", "M")
hMassPt = dfData.Histo2D(("hMassPt", ";#it{p}_{T} (GeV/#it{c});m_{^{3}H + #pi^{-} + c.c.} (GeV/#it{c})^{2}", 20, 0, 10., 50, 2.97, 3.02), "Pt", "M")
hDCAv0DaughPtMC = dfRecMC.Histo2D(("hDCAv0DaughPtMC", ";#it{p}_{T} (GeV/#it{c});V0 daughters DCA (cm)", 20, 0, 10., 2000, 0, 1), "Pt", "fDcaV0Daug")
hMassLambdaPt = dfData.Histo3D(("hMassLambdaPt", ";#it{p}_{T} (GeV/#it{c});m_{p + #pi^{-} + c.c.} (GeV/#it{c}^{2});m_{^{3}H + #pi^{-} + c.c.} (GeV/#it{c}^{2})", 20, 0, 10., 425, 1.075, 1.5, 80, 2.95, 3.03), "Pt", "MLambda0", "M")
hDCAv0DaughPt = dfData.Histo2D(("hDCAv0DaughPt", ";#it{p}_{T} (GeV/#it{c});V0 daughters DCA (cm)", 20, 0, 10., 2000, 0, 1), "Pt", "fDcaV0Daug")
hMassK0sPt = dfData.Histo2D(("hMassK0sPt", ";#it{p}_{T} (GeV/#it{c});m_{#pi^{+} + #pi^{-}} (GeV/#it{c})^{2}", 20, 0, 10., 80, 0.455, 0.535), "Pt", "MK0s")
hDCAptPi = dfData.Histo2D(("hDCAptPi", ";#it{p}_{T} (GeV/#it{c});#pi DCA (cm)", 20, 0, 10., 2000, -1, 1), "Pt", "fDcaPi")
hDCAptPiMC = dfRecMC.Histo2D(("hDCAptPiMC", ";#it{p}_{T} (GeV/#it{c});#pi DCA (cm)", 20, 0, 10., 2000, -1, 1), "Pt", "fDcaPi")
hDCApt3H = dfData.Histo2D(("hDCApt3H", ";#it{p}_{T} (GeV/#it{c});^{3}H DCA (cm)", 20, 0, 10., 2000, -1, 1), "Pt", "fDca3H")
hDCApt3HMC = dfRecMC.Histo2D(("hDCApt3HMC", ";#it{p}_{T} (GeV/#it{c});^{3}H DCA (cm)", 20, 0, 10., 2000, -1, 1), "Pt", "fDca3H")

hArmenterosPodolanski = dfData.Histo2D(("hArmenterosPodolanski", ";#alpha;#it{q}_{T} (GeV/#it{c})", 200, -1, 1, 500, 0, 0.5), "alpha", "qt")
hArmenterosPodolanskiMC = dfRecMC.Histo2D(("hArmenterosPodolanskiMC", ";#alpha;#it{q}_{T} (GeV/#it{c})", 200, -1, 1, 500, 0, 0.5), "alpha", "qt")

hMassFit = dfData.Filter("Pt > 4 && Pt < 8").Histo1D(("hMassFit", ";m_{^{3}H + #pi^{-} + c.c.} (GeV/#it{c}^{2})", 45, 2.975, 3.02), "M")
hMassFitMC = dfRecMC.Filter("Pt > 4 && Pt < 8").Histo1D(("hMassFitMC", ";m_{^{3}H + #pi^{-} + c.c.} (GeV/#it{c}^{2})", 45, 2.975, 3.02), "M")

mass = ROOT.RooRealVar("mass", "mass", 2.991, hMassFit.GetXaxis().GetXmin(), hMassFit.GetXaxis().GetXmax())
mu = ROOT.RooRealVar("mu", "#mu", 2.991, 2.989, 2.993)
sigma = ROOT.RooRealVar("sigma", "#sigma", 0.002, 0.001, 0.01)
a1 = ROOT.RooRealVar("a1", "a1", 0., -1., 1.)
bkg = ROOT.RooPolynomial("bkg", "bkg", mass, ROOT.RooArgList(a1))
alpha = ROOT.RooRealVar("alpha", "#alpha", 1., 0., 5.)
n = ROOT.RooRealVar("n", "n", 1., 0., 10.)
cb = ROOT.RooCrystalBall("cb", "cb", mass, mu, sigma, alpha, n, True)
signal = ROOT.RooGaussian("signal", "signal", mass, mu, sigma)
n_sig = ROOT.RooRealVar("n_sig", "n_sig", 1000, 0, 1e7)
n_bkg = ROOT.RooRealVar("n_bkg", "n_bkg", 1000, 0, 1e7)
fulltpl = ROOT.RooAddPdf("fulltpl", "fulltpl", ROOT.RooArgList(cb, bkg), ROOT.RooArgList(n_sig, n_bkg))

rooHistMC = ROOT.RooDataHist("rooHistMC", "rooHistMC", ROOT.RooArgList(mass), hMassFitMC.GetPtr())
cb.fitTo(rooHistMC)
frameMC = mass.frame()
rooHistMC.plotOn(frameMC)
cb.plotOn(frameMC)
alpha.setConstant(True)
n.setConstant(True)
sigma.setConstant(True)
print(f'chi2 MC: {frameMC.chiSquare()}')

rooHistData = ROOT.RooDataHist("rooHistData", "rooHistData", ROOT.RooArgList(mass), hMassFit.GetPtr())
bkg.fitTo(rooHistData)
frame = mass.frame()
rooHistData.plotOn(frame)
bkg.plotOn(frame)
print(f'chi2: {frame.chiSquare()}')


# Create a workspace
workspace = ROOT.RooWorkspace("workspace")
workspace.Import(rooHistData)
workspace.Import(fulltpl)

workspace.Print()
# Set the confidence level
confidence_level = 0.95

# Create a ModelConfig
model_config = RooStats.ModelConfig("model_config", workspace)
model_config.SetPdf(fulltpl)
model_config.SetParametersOfInterest(n_sig)
model_config.SetObservables(ROOT.RooArgSet(mass))
model_config.SetNuisanceParameters(ROOT.RooArgSet(mu))
model_config.SetNuisanceParameters(ROOT.RooArgSet(n_bkg))
oldValue = n_sig.getVal()

# Create a null (background-only) model
null_model = RooStats.ModelConfig("null_model", workspace)
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
low_limit = result.LowerLimit()
error_up = result.UpperLimitEstimatedError()
error_lower = result.LowerLimitEstimatedError()

print(f"Upper limit at {confidence_level*100:.0f}% CL: {upper_limit:.2f} +- {error_up}")



eff = 0.03  ## needs to be rewighted
br = 0.25
delta_pt = 4
rapidity_window = 2
mat_antimat_fac = 2
n_ev = 2.98e9







upper_limit_corrected = upper_limit / (eff * br * rapidity_window * mat_antimat_fac * n_ev * delta_pt)
print(f"Upper limit corrected: {upper_limit_corrected:.2e}")




effMC = hRecMC.Clone("effMC")
effMC.Divide(hRecMC.GetPtr(), hGenMC.GetPtr(), 1., 1., "B")
frame.Write("fitData")
frameMC.Write("fitMC")
hRpt.Write()
effMC.Write()
hPID3H.Write()
hCosPAptMC.Write()
hCosPApt.Write()
hMassPtMC.Write()
hMassPt.Write()
hMassLambdaPtMC.Write()
hMassLambdaPt.Write()
hMassK0sPt.Write()
hDCAptPi.Write()
hDCAptPiMC.Write()
hDCApt3H.Write()
hDCApt3HMC.Write()
hDCAv0DaughPt.Write()
hDCAv0DaughPtMC.Write()
hArmenterosPodolanski.Write()
hArmenterosPodolanskiMC.Write()
hMassFit.Write()
outFile.Close()
