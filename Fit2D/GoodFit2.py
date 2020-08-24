import ROOT
from ROOT import *
import numpy
import sys
import plotFrame as pf

gROOT.SetBatch(ROOT.kTRUE)
resbins = [0.0, 30.0, 50.0, 70.0, 90.0, 110.0, 130.0, 150.0, 170.0, 190.0, 210.0, 230.0, 250.0, 270.0, 290.0, 310.0, 330.0, # uniform bins\ 
	350.0, 370.0, 390.0, 412.0, 433.0, 455.0, 478.0, 502.0, 526.0, 551.0, 577.0, 603.0, 631.0, 659.0, 688.0, 717.0, 748.0, 779.0, 812.0, 845.0, 879.0, 914.0, 950.0, 988.0, 1026.0, 1066.0, 1106.0, 1148.0, 1191.0, 1235.0, 1281.0, 1328.0, 1376.0, 1426.0, 1477.0, 1529.0, 1583.0, # resolution-based binning\
	1640, 1700, 1760, 1820, 1880, 1940, 2000, 2060, 2120, 2180, 2240, 2300, 2360, 2420, 2480, 2540, 2600, 2660, 2720, 2780, 2840, 2900, 2960, 3020, 3060, 3120, 3180, 3240, 3300, 3360, 3420, 3480, 3540, 3600, 3660, 3720, 3780, 3840, 3900, 3960, 4020, 4080, 4140, 4200, 4260, 4320, 4380, 4440, 4500, 4560, 4620, 4680, 4740, 4800, 4860, 4920, 4980, 5040] # uniform bins

##############
## controls ##
##############
fname = "P2"
polN = 3
startCut = 400
makeToy = True
sigNames = ["Diquark_chi500suu2000",\
	"Diquark_chi1000suu4000",\
	"Diquark_chi1800suu8000"]
##############

f = pf.funcForms[fname]
numpar = ROOT.TF1("t",f,1,2).GetNpar()
lumi = pf.lumifb[2000+int(sys.argv[1])]

pf.UseWS("QCD_20"+sys.argv[1]+"_"+fname)

alpha_binning = pf.MakeNBinsFromMinToMax(100, 0.1, 0.4)
Mbins = pf.MakeNBinsFromMinToMax(400, 0.0, 4000.0)
m4j_binning = pf.MakeNBinsFromMinToMax(400, 0., 10000.)

QCD = pf.importTH2("/home/th544/CMSSW_10_6_2/src/tests/AvJJ_alpha_test/QCD_20"+sys.argv[1]+"_rebin.root", "AvJJ_alpha_rebin", logz=True, png=True)
QCD.RebinX(int(sys.argv[3]))
QCD.RebinY(int(sys.argv[4]))
if makeToy: QCD = pf.MakeAToy(QCD)

signals = []
for n in range(len(sigNames)):
	sig = pf.importTH1(pf.sigDir[sigNames[n]], "AvJJ_alpha_rebin", name=sigNames[n], png=True)
	sig.RebinX(int(sys.argv[3]))
	sig.RebinY(int(sys.argv[4]))
	sig.Scale(1/500.)
	sig.SetTitle(sigNames[n])
	signals.append(sig)

norms = []
starts = []
chi2Ndof = []
yalph = []
yalphe = []
CIs = []

for j in range(1,QCD.GetNbinsY()):
	print "|===> Slice "+ str(j)
	Slice = QCD.ProjectionX("q_"+str(j),j,j)
	minAlpha = QCD.GetYaxis().GetBinLowEdge(j)
	maxAlpha = QCD.GetYaxis().GetBinUpEdge(j)
	sigProjections = []
	for n in range(len(signals)): sigProjections.append(signals[n].ProjectionX(sigNames[n],j,j))
	if not (Slice.Integral() > 0): continue
	if QCD.GetYaxis().GetBinCenter(j) > 0.35 or QCD.GetYaxis().GetBinCenter(j) < 0.15: continue
	m = Slice.GetMaximumBin() + int(sys.argv[2])
	start = max(Slice.GetBinLowEdge(m), startCut)
	end = Slice.GetBinCenter(QCD.GetNbinsX())
	fit = ROOT.TF1("tF_"+str(j), f, start, end)
	fFrame = pf.Fit1D(Slice, fit)
	fFrame.makePull("Slice_"+str(j)+"_"+fname, title="", signals=sigProjections, xTitle="#bar{M}_{jj}", yTitle="Number of Events", displayCI=True, legTitle=str(minAlpha) + " < #alpha < " + str(maxAlpha), lumi=lumi, CMSExtra="Preliminary")
	ndof = fit.GetNDF()
	chi2 = fit.GetChisquare()
	if (ndof < 5) or (chi2/ndof > 1.5): continue
	starts.append(start)
	norms.append(Slice.Integral())
	chi2Ndof.append(chi2/ndof)
	yalph.append(QCD.GetYaxis().GetBinCenter(j))
	yalphe.append(QCD.GetYaxis().GetBinWidth(j)/2.)

pf.MakeTG1(yalph, starts, "start", "Starts", xTitle="#alpha", yTitle="#bar{M}_{jj} [GeV]")
pf.MakeTG1(yalph, chi2Ndof, "chi2Ndof", "#chi^{2}/Ndof", xTitle="#alpha")
pf.MakeTG1(yalph, norms, "norm", "Slice Integrals", xTitle="#alpha")

f2Frame = pf.Fit2D(QCD, SF=signals)
f2Frame.SetCMSText(lumi=lumi, CMSExtra="Preliminary")
f2Frame.SliceYFit2D(fname=fname, polN=polN, options="ER0", startbin=sys.argv[2], startcut=startCut)

f2Frame.SetTitleX("#bar{M}_{jj} [GeV]")
f2Frame.SetTitleY("#alpha")
f2Frame.Print2D("G", zrangeP=[-6., 6.])

f2Frame.SetTitleX("#bar{M}_{jj} [GeV]")
f2Frame.SetTitleY("M_{jjjj} [GeV]")
f2Frame.MakeMJJM4J("GT", zrangeP=[-6., 6.])

f2Frame.SetTitleX("#bar{M}_{jj} [GeV]")
f2Frame.SetTitleY("Number of events")
f2Frame.ProjectionX("MJJ_all", "", legTitle="Fit Proj_#bar{M}_{JJ}")
f2Frame.ProjectionX("MJJ_L", "", legTitle="Fit Proj_#bar{M}_{JJ} L", start=0, end=0.33, pullRange=[-3.,3.])
f2Frame.ProjectionX("MJJ_M", "", legTitle="Fit Proj_#bar{M}_{JJ} M", start=0.33, end=0.67, pullRange=[-3.,3.])
f2Frame.ProjectionX("MJJ_H", "", legTitle="Fit Proj_#bar{M}_{JJ} H", start=0.67, end=1, pullRange=[-3.,3.])

f2Frame.SetTitleX("#alpha")
f2Frame.ProjectionY("Alpha_all", "", legTitle="Fit Proj_#alpha")
f2Frame.ProjectionY("Alpha_L", "", legTitle="Fit Proj_#alpha L", start=0, end=0.33)
f2Frame.ProjectionY("Alpha_M", "", legTitle="Fit Proj_#alpha M", start=0.33, end=0.67)
f2Frame.ProjectionY("Alpha_H", "", legTitle="Fit Proj_#alpha H", start=0.67, end=1)

f2Frame.SetTitleX("M_{jjjj} [GeV]")
f2Frame.ProjectionM4J("M4J_all", "", legTitle="Global Fit Projection onto M_{JJJJ}")
f2Frame.ProjectionM4J("M4J_L", "", legTitle="Fit Proj_M_{JJJJ} L", start=0, end=0.33)
f2Frame.ProjectionM4J("M4J_M", "", legTitle="Fit Proj_M_{JJJJ} M", start=0.33, end=0.67)
f2Frame.ProjectionM4J("M4J_H", "", legTitle="Fit Proj_M_{JJJJ} H", start=0.67, end=1)

f2Frame.PrintChi2Ndof()