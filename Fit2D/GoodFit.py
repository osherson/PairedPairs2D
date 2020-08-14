#
import ROOT
from ROOT import *
import numpy
import sys
import GoodFormat
import os
import math
from plotFrame import *

gROOT.SetBatch(ROOT.kTRUE)
resbins = [0.0, 30.0, 50.0, 70.0, 90.0, 110.0, 130.0, 150.0, 170.0, 190.0, 210.0, 230.0, 250.0, 270.0, 290.0, 310.0, 330.0, # uniform bins\ 
	350.0, 370.0, 390.0, 412.0, 433.0, 455.0, 478.0, 502.0, 526.0, 551.0, 577.0, 603.0, 631.0, 659.0, 688.0, 717.0, 748.0, 779.0, 812.0, 845.0, 879.0, 914.0, 950.0, 988.0, 1026.0, 1066.0, 1106.0, 1148.0, 1191.0, 1235.0, 1281.0, 1328.0, 1376.0, 1426.0, 1477.0, 1529.0, 1583.0, # resolution-based binning\
	1640, 1700, 1760, 1820, 1880, 1940, 2000, 2060, 2120, 2180, 2240, 2300, 2360, 2420, 2480, 2540, 2600, 2660, 2720, 2780, 2840, 2900, 2960, 3020, 3060, 3120, 3180, 3240, 3300, 3360, 3420, 3480, 3540, 3600, 3660, 3720, 3780, 3840, 3900, 3960, 4020, 4080, 4140, 4200, 4260, 4320, 4380, 4440, 4500, 4560, 4620, 4680, 4740, 4800, 4860, 4920, 4980, 5040] # uniform bins

def MakeNBinsFromMinToMax(N,Min,Max):
	BINS = []
	for i in range(N+1):
		BINS.append(Min+(i*(Max-Min)/N))
	return BINS

def makeParStrings(npar, polN):
	Ps = []
	for n in range(numpar):
		P = "([" + str(n*(polN + 1)) + "]"
		for m in range(polN):
			P = P + " + [" + str(n*(polN + 1) + m + 1) + "]"
			for o in range(m + 1): P = P + "*y"
		P = P + ")"
		Ps.append(P)
	return Ps

def MakeAToy(Hi):
	R = numpy.random
	Ho = Hi.Clone("H"+str(R.randint(1)))
	for j in range(1,Ho.GetNbinsY()):
		for i in range(1,Ho.GetNbinsX()):
			v = Ho.GetBinContent(i,j)
			e = Ho.GetBinError(i,j)
			n = max(0,numpy.rint(v + (e*R.normal())))
			Ho.SetBinContent(i,j,int(n))
	Ho.Sumw2()
	return Ho

def takeRMS(arr): 
    square, mean, root = 0, 0.0, 0.0
    for i in range(len(arr)): square += (arr[i]**2) 
    mean = square /len(arr)
    root = math.sqrt(mean) 
    return root 

def takeMean(arr):
	avg = 0
	for i in range(len(arr)): avg += arr[i]
	avg = avg/len(arr)
	return avg

numMbins = 400
Mmin, Mmax = 0.0, 4000.0
Mbinwidth = (Mmax-Mmin)/float(numMbins)
Mbins = MakeNBinsFromMinToMax(numMbins, Mmin, Mmax)

##############
## controls ##
##############
constrainedfit = False
fname = "P2"
polN = 3
##############

numpar = 0
tfname = ""
if fname is "P2": tfname = "TMath::Power(1-(x/13000.0),[0])/TMath::Power(x/13000.0,[1])"
if fname is "P3": tfname = "TMath::Power(1.-(x/13000.0),[0])/(TMath::Power(x/13000.0,[1]+([2]*TMath::Log(x/13000.0))))"
if fname is "P4": tfname = "TMath::Power(1.-(x/13000.0),[0])/(TMath::Power(x/13000.0,[1]+([2]*TMath::Log(x/13000.0))+([3]*TMath::Power(TMath::Log(x/13000.0),2))))"
if fname is "CDF,1": tfname = "(1/TMath::Power(x, [0]))*TMath::Power(1-(x/13000.0),[1])"
if fname is "CDF,2": tfname = "(1/TMath::Power(x, [0]))*TMath::Power(1-(x/13000.0)+TMath::Power(x/13000.0,2),[1])"
if fname is "CDF,3": tfname = "(1/TMath::Power(x, [0]))*TMath::Power(1-(x/13000.0)+TMath::Power(x/13000.0,2)-TMath::Power(x/13000.0,3),[1])"
if fname is "expoPoli,1": tfname = "TMath::Exp([0]+([1]*TMath::Log(x)))"
if fname is "expoPoli,2": tfname = "TMath::Exp([0]+([1]*TMath::Log(x))+([2]*TMath::Power(TMath::Log(x),2)))"
if fname is "expoPoli,3": tfname = "TMath::Exp([0]+([1]*TMath::Log(x))+([2]*TMath::Power(TMath::Log(x),2))+([3]*TMath::Power(TMath::Log(x),3)))"
numpar = ROOT.TF1("t",tfname,1,2).GetNpar()
print(numpar)

# workspace
if not os.path.exists("QCD_20"+sys.argv[1]+"_"+fname):
    os.makedirs("QCD_20"+sys.argv[1]+"_"+fname)
os.chdir("QCD_20"+sys.argv[1]+"_"+fname)

numabins = 100
amin, amax = 0.1, 0.4
abinwidth = (amax-amin)/float(numabins)
alpha_binning = MakeNBinsFromMinToMax(numabins, amin, amax)
m4j_binning = MakeNBinsFromMinToMax(60, 0., 10000.)

InFile = TFile("/home/th544/CMSSW_10_6_2/src/test/AvJJ_alpha_test/QCD_20"+sys.argv[1]+"_rebin.root")
QCD = InFile.Get("AvJJ_alpha_rebin")
try: QCD.RebinY(int(sys.argv[3]))
except: print("no rebin")
# QCD = MakeAToy(QCD)
# QCD.RebinX(4)
C_QCD_uncut = TCanvas()
C_QCD_uncut.cd()
C_QCD_uncut.SetLogz()
QCD.Draw("colz")
C_QCD_uncut.Print("QCD_raw.png")

yalph = []
yalphe = []
norm = []
starts = []
chiN = []

ps = [[[],[]] for i in range(numpar)]
# p0 = []
# p0e = []
# p1 = []
# p1e = []
# p2 = []
# p3 = []

Cs = TCanvas()
Cs.cd()
Cs.SetLogy()

y = []
x = []
z = []
ex = []
ey = []
ez = []

pulls = [[] for i in range(QCD.GetNbinsX())]
for j in range(1,QCD.GetNbinsY()):
	print "====---> "+ str(j)
	Q = QCD.ProjectionX("q_"+str(j),j,j)
	# Q = MakeAToy(Q)
	if not (Q.Integral() > 0): continue
	if QCD.GetYaxis().GetBinCenter(j) > 0.35: continue
	m = Q.GetMaximumBin() + int(sys.argv[2])
	for i in range(m,QCD.GetNbinsX()):
		y.append(QCD.GetYaxis().GetBinCenter(j))
		x.append(QCD.GetXaxis().GetBinCenter(i))
		ey.append(QCD.GetYaxis().GetBinWidth(j)/2.)
		ex.append(QCD.GetXaxis().GetBinWidth(i)/2.)
		z.append(Q.GetBinContent(i))
		ez.append(Q.GetBinError(i))
	start = Q.GetBinLowEdge(m)
	# end = Q.GetBinCenter(Mbins[-1])
	end = Q.GetBinCenter(QCD.GetNbinsX())
	fit = TF1("tF", tfname, start, end)
	Q.Fit(fit,"EM0R")
	P = TH1F("P%d" % j, "P%d" % j, len(Mbins)-1, numpy.array(Mbins))
	F = TH1F("F%d" % j, "F%d" % j, len(Mbins)-1, numpy.array(Mbins))

	for n in range(1,len(Mbins)-1):
		p = 0
		f = 0
		if n >= Q.FindBin(start):
			if Q.GetBinContent(n) == 0: Q.SetBinError(n, 1.4)
			p = (Q.GetBinContent(n)-fit.Eval(Q.GetBinCenter(n)))/Q.GetBinError(n)
			f = fit.Eval(Q.GetBinCenter(n))/Q.GetBinWidth(n)
		P.SetBinContent(n, p)
		F.SetBinContent(n, f)
		pulls[n].append(p)
	stylePull(Q, fit, "Slice_"+str(j)+"_"+fname)

	ndof = fit.GetNDF()
	chi2 = fit.GetChisquare()
	if (ndof < 5) or (chi2/ndof > 1.5): continue
	starts.append(start)
	norm.append(Q.Integral())
	chiN.append(chi2/ndof)
	for n in range(numpar):
		ps[n][0].append(fit.GetParameter(n))
		ps[n][1].append(fit.GetParError(n))
	# p0.append(fit.GetParameter(0))
	# p1.append(fit.GetParameter(1))
	# p0e.append(fit.GetParError(0))
	# p1e.append(fit.GetParError(1))
	#p2.append(fit.GetParameter(2))
	#p3.append(fit.GetParameter(3))
	yalph.append(QCD.GetYaxis().GetBinCenter(j))
	yalphe.append(QCD.GetYaxis().GetBinWidth(j)/2.)
	# Q.Draw()
	# fit.Draw("same")
	# Cs.Print("Slice_"+str(j)+"_"+fname+".png")

# print(ps)
# print(yalph)

def MakePlot(V, N):
	g = TGraph(len(V), numpy.array(yalph), numpy.array(V))
	g.SetTitle(N)
	C = TCanvas()
	C.cd()
	g.Draw("AL")
	C.Print(N+".png")
def MakeEPlot(V, N, E, O):
	g = TGraphErrors(len(V), numpy.array(yalph), numpy.array(V), numpy.array(yalphe), numpy.array(E))
	f = TF1("f", "pol"+str(O), min(yalph), max(yalph))
	g.Fit(f)
	g.SetTitle(N)
	C = TCanvas()
	C.cd()
	g.Draw("AL")
	f.Draw("same")
	C.Print(N+".png")
	return f

MakePlot(starts, "start")
MakePlot(chiN, "chi2N")
MakePlot(norm, "norm")
oPs = []
for n in range(numpar): oPs.append(MakeEPlot(ps[n][0], "p"+str(n), ps[n][1], polN))
# oP0 = MakeEPlot(p0, "p0", p0e, polN)
# oP1 = MakeEPlot(p1, "p1", p1e, polN)
# oP2 = MakeEPlot(p2, "p1", p2e, polN)

G = TGraph2DErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(z), numpy.array(ex), numpy.array(ey), numpy.array(ez))
C = TCanvas()
C.cd()
G.Draw("AP")

X = "x/13000.0"
Ps = makeParStrings(numpar, polN)
# P0 = "([0] + [1]*y + [2]*y*y + [3]*y*y*y)"
# P1 = "([4] + [5]*y + [6]*y*y + [7]*y*y*y)"

FitFunc = ""
# fname: "name;N"
if fname is "P2": FitFunc = "TMath::Power((1.-"+X+"),"+Ps[0]+")/TMath::Power("+X+","+Ps[1]+")"
if fname is "P3": FitFunc = "TMath::Power((1.-"+X+"),"+Ps[0]+")/(TMath::Power("+X+","+Ps[1]+"+("+Ps[2]+"*TMath::Log("+X+"))))"
if fname is "P4": FitFunc = "TMath::Power((1.-"+X+"),"+Ps[0]+")/(TMath::Power("+X+","+Ps[1]+"+("+Ps[2]+"*TMath::Log("+X+"))+("+Ps[3]+"*TMath::Power(TMath::Log("+X+"),2))))"
if fname is "CDF,1": FitFunc = "TMath::Power(1-"+X+","+Ps[1]+")/TMath::Power(x, "+Ps[0]+")"
if fname is "CDF,2": FitFunc = "TMath::Power(1-"+X+"+TMath::Power("+X+",2),"+Ps[1]+")/TMath::Power(x, "+Ps[0]+")"
if fname is "CDF,3": FitFunc = "TMath::Power(1-"+X+"+TMath::Power("+X+",2)-TMath::Power("+X+",3),"+Ps[1]+")/TMath::Power(x, "+Ps[0]+")"
if fname is "expoPoli,1": FitFunc = "TMath::Exp("+Ps[0]+"+("+Ps[1]+"*TMath::Log(x)))"
if fname is "expoPoli,2": FitFunc = "TMath::Exp("+Ps[0]+"+("+Ps[1]+"*TMath::Log(x))+("+Ps[2]+"*TMath::Power(TMath::Log(x),2)))"
if fname is "expoPoli,3": FitFunc = "TMath::Exp("+Ps[0]+"+("+Ps[1]+"*TMath::Log(x))+("+Ps[2]+"*TMath::Power(TMath::Log(x),2))+("+Ps[3]+"*TMath::Power(TMath::Log(x),3)))"

F = TF2("hardFit", FitFunc, 0., 3000., 0.1, 0.4)

for n in range(numpar):
	for m in range(polN + 1):
		F.SetParameter(n*(polN + 1) + m, oPs[n].GetParameter(m))
# F.SetParameter(0, oP0.GetParameter(0))
# F.SetParameter(1, oP0.GetParameter(1))
# F.SetParameter(2, oP0.GetParameter(2))
# F.SetParameter(3, oP0.GetParameter(3))
# F.SetParameter(4, oP1.GetParameter(0))
# F.SetParameter(5, oP1.GetParameter(1))
# F.SetParameter(6, oP1.GetParameter(2))
# F.SetParameter(7, oP1.GetParameter(3))

G.Fit(F, "EMR0")
ndof = F.GetNDF()
chi2 = F.GetChisquare()
print("Chi2Ndof = "  + str(chi2)+"/"+str(ndof))

F.Draw("same")

Qout = QCD.Clone("out")
Qout.Add(Qout,-1.)
Qout.SetStats(0)

PullHist = TH1F("PullHist", ";pull;bins", 1200, -150, 150)
globalbkg = ROOT.TH2F("gbkg","Global bkg;#bar M_{jj}[GeV];#alpha", len(Mbins)-1, numpy.array(Mbins), len(alpha_binning)-1, numpy.array(alpha_binning))
resbins = [0.0, 30.0, 50.0, 70.0, 90.0, 110.0, 130.0, 150.0, 170.0, 190.0, 210.0, 230.0, 250.0, 270.0, 290.0, 310.0, 330.0, # uniform bins\ 
	350.0, 370.0, 390.0, 412.0, 433.0, 455.0, 478.0, 502.0, 526.0, 551.0, 577.0, 603.0, 631.0, 659.0, 688.0, 717.0, 748.0, 779.0, 812.0, 845.0, 879.0, 914.0, 950.0, 988.0, 1026.0, 1066.0, 1106.0, 1148.0, 1191.0, 1235.0, 1281.0, 1328.0, 1376.0, 1426.0, 1477.0, 1529.0, 1583.0, # resolution-based binning\
	1640, 1700, 1760, 1820, 1880, 1940, 2000, 2060, 2120, 2180, 2240, 2300, 2360, 2420, 2480, 2540, 2600, 2660, 2720, 2780, 2840, 2900, 2960, 3020, 3060, 3120, 3180, 3240, 3300, 3360, 3420, 3480, 3540, 3600, 3660, 3720, 3780, 3840, 3900, 3960, 4020, 4080, 4140, 4200, 4260, 4320, 4380, 4440, 4500, 4560, 4620, 4680, 4740, 4800, 4860, 4920, 4980, 5040] # uniform bins
globalfit = ROOT.TH2F("gfit","Global fit of " + ("un" if constrainedfit is False else "") + "constrained fit;#bar M_{jj}[GeV];#alpha", len(Mbins)-1, numpy.array(Mbins), len(alpha_binning)-1, numpy.array(alpha_binning))
# initialize these th2's with QCD.Clone()

# paramfile = open('parameters.txt', 'w')
# paramfile.write('[')
Qout.GetZaxis().SetRangeUser(-5.,5.)
globalmax = 0 # used to calibrate the pulls
for j in range(1,QCD.GetNbinsY()):
	Q = QCD.ProjectionX("q_"+str(j),j,j)
	if (Q.Integral() > 0) and (QCD.GetYaxis().GetBinCenter(j)) < 0.35:
		m = Q.GetMaximumBin() + int(sys.argv[2])
		for i in range(m,QCD.GetNbinsX()):
				q = QCD.GetBinContent(i,j)
				eq = Q.GetBinError(i)
				if q == 0: eq = 1.4
				est = F.Eval(Q.GetBinCenter(i), QCD.GetYaxis().GetBinCenter(j))
				P = (q - est)/eq
				Qout.SetBinContent(i,j,P)
				PullHist.Fill(P)
				globalfit.SetBinContent(i,j,est)
				globalbkg.SetBinContent(i,j,q)
				globalmax = max(abs(P), globalmax)
		for i in range(1,m-1):
			# continue
			Qout.SetBinContent(i,j,0)
	else:
		for i in range(1,QCD.GetNbinsX()):
			# continue
			Qout.SetBinContent(i,j,0)
# paramfile.write(']')

CPH = TCanvas()
CPH.cd()
PullHist.Draw("hist")		
		
# globalfit
C_GF = TCanvas()
C_GF.SetLogz()
C_GF.cd()
# globalfit.GetZaxis().SetLimits(0.2, 1000)
globalfit.GetZaxis().SetRangeUser(0.001, 1000)
globalfit.SetStats(0)
globalfit.Draw("colz")
C_GF.Print("GFit.png")

# globalbkg
C_GB = TCanvas()
C_GB.SetLogz()
C_GB.cd()
# globalbkg.GetZaxis().SetLimits(0.2, 1000)
globalbkg.GetZaxis().SetRangeUser(0.00001, 1000)
globalbkg.SetStats(0)
globalbkg.Draw("colz")
C_GB.Print("GBkg.png")

# transformed fits
tpull = ROOT.TH2F("tpull","Global pull of " + ("un" if constrainedfit is False else "") + "constrained fit in #bar M_{jj} #times M_{4j} space;#bar M_{jj} [GeV];M_{4j} [GeV]", len(Mbins)-1, numpy.array(Mbins), len(m4j_binning)-1, numpy.array(m4j_binning))
tbkg = ROOT.TH2F("tbkg","Global bkg in #bar M_{jj} #times M_{4j} space;#bar M_{jj} [GeV];M_{4j} [GeV]", len(Mbins)-1, numpy.array(Mbins), len(m4j_binning)-1, numpy.array(m4j_binning))
tfit = ROOT.TH2F("tfit","Global fit of " + ("un" if constrainedfit is False else "") + "constrained fit in #bar M_{jj} #times M_{4j} space;#bar M_{jj} [GeV];M_{4j} [GeV]", len(Mbins)-1, numpy.array(Mbins), len(m4j_binning)-1, numpy.array(m4j_binning))

for i in range(1,len(m4j_binning)):
	for j in range(1,len(Mbins)):
		# calculate alpha
		M2J = tfit.GetXaxis().GetBinCenter(j)
		M4J = tfit.GetYaxis().GetBinCenter(i)
		a = M2J/M4J
		ybin = globalbkg.GetYaxis().FindBin(a)

		# find corresponding bin
		tpullval = Qout.GetBinContent(j, ybin)
		tbkgval = globalbkg.GetBinContent(j, ybin)
		tfitval = globalfit.GetBinContent(j, ybin)

		tpull.SetBinContent(j, i, tpullval)
		tfit.SetBinContent(j, i, tfitval)
		tbkg.SetBinContent(j, i, tbkgval)

ROOT.gStyle.SetPalette(57)
tc1 = ROOT.TCanvas()
tc1.cd()
ROOT.gPad.SetPhi(150)
# ROOT.gPad.SetTheta(0)
tc1.SetLogz()
tfit.UseCurrentStyle()
tfit.SetStats(0)
tfit.GetZaxis().SetRangeUser(0.00001, 1000)
tfit.Draw("SURF3")
tc1.Print("GTFitSurf3.png")

tc2 = ROOT.TCanvas()
tc2.cd()
ROOT.gPad.SetPhi(150)
tc2.SetLogz()
tbkg.UseCurrentStyle()
tbkg.SetStats(0)
tbkg.GetZaxis().SetRangeUser(0.00001, 1000)
tbkg.Draw("SURF3")
tc2.Print("GTBkgSurf3.png")

mjj_bkg = globalbkg.ProjectionX("mjj_bkg;#bar{M}_{jj};Number of events", 0, -1, "e")
mjj_bkg.SetTitle("#bar{M}_{jj} QCD MC with Estimate by " + fname)
mjjjj_bkg = tbkg.ProjectionY("mjjjj_bkg;#bar{M}_{jjjj};Number of events", 0, -1, "e")
alpha_bkg = globalbkg.ProjectionY("alpha_bkg;#alpha;Number of events", 0, -1, "e")

mjj_est = globalfit.ProjectionX("mjj_est;#bar{M}_{jj};Number of events", 0, -1, "e")
mjjjj_est = tfit.ProjectionY("mjjjj_est;#bar{M}_{jjjj};Number of events", 0, -1, "e")
alpha_est = globalfit.ProjectionY("alpha_est;#alpha;Number of events", 0, -1, "e")

C_mjj_bkg = ROOT.TCanvas()
C_mjj_bkg.cd()
mjj_bkg.Draw()
C_mjj_bkg.Print("mjj_bkg.png")

C_mjjjj_bkg = ROOT.TCanvas()
C_mjjjj_bkg.cd()
mjjjj_bkg.Draw()
C_mjjjj_bkg.Print("mjjjj_bkg.png")

C_alpha_bkg = ROOT.TCanvas()
C_alpha_bkg.cd()
alpha_bkg.Draw()
C_alpha_bkg.Print("alpha_bkg.png")

C_mjj_est = ROOT.TCanvas()
C_mjj_est.cd()
mjj_est.Draw()
C_mjj_est.Print("mjj_est.png")

C_mjjjj_est = ROOT.TCanvas()
C_mjjjj_est.cd()
mjjjj_est.Draw()
C_mjjjj_est.Print("mjjjj_est.png")

C_alpha_est = ROOT.TCanvas()
C_alpha_est.cd()
alpha_est.Draw()
C_alpha_est.Print("alpha_est.png")

# C_mjj_all = ROOT.TCanvas('c_mjj_all',800,800)
# C_mjj_all.SetLogy()
# C_mjj_all.cd()
# C_mjj_p = ROOT.TPad()
# mjj_bkg.SetLineColor(ROOT.kBlack)
# mjj_bkg.SetMarkerStyle(20)
# mjj_bkg.Draw()
# mjj_est.SetLineColor(ROOT.kRed)
# mjj_est.Draw("hist same")
# C_mjj_all.BuildLegend()
# C_mjj_all.GetPad(0).SetTitle("mjj_all;#bar{M}_{jj} [GeV];Number of events")
mjj_pull = ROOT.TH1F("mjj_pull", "", len(Mbins)-1, numpy.array(Mbins))
for i in range(1, len(Mbins)-1):
	# p = sum(pulls[i])/len(pulls[i])
	p = takeMean(pulls[i])
	# if mjj_bkg.GetBinError(i) != 0: p = (mjj_bkg.GetBinContent(i) - mjj_est.GetBinContent(i))/mjj_bkg.GetBinError(i)
	# else: p = (mjj_bkg.GetBinContent(i) - mjj_est.GetBinContent(i))/1.4
	# p = (mjj_bkg.GetBinContent(i) - mjj_est.GetBinContent(i))/max(mjj_bkg.GetBinError(i), 1.4)
	mjj_pull.SetBinContent(i, p)
# stylePull(mjj_bkg, mjj_est, mjj_pull, "mjj_all_fancy_" + fname)
styleRatio(mjj_bkg, mjj_est, "mjj_all_fancy_" + fname)
# C_mjj_all = stylePull(mjj_bkg, mjj_est, mjj_pull, "M_JJ QCD MC 20" + sys.argv[1])
# C_mjj_all.Print("mjj_all.png")

C_mjjjj_all = ROOT.TCanvas()
C_mjjjj_all.SetLogy()
C_mjjjj_all.cd()
mjjjj_bkg.SetLineColor(ROOT.kBlack)
mjjjj_bkg.SetMarkerStyle(20)
mjjjj_bkg.Draw()
mjjjj_est.SetLineColor(ROOT.kRed)
mjjjj_est.Draw("hist same")
# C_mjjjj_all.BuildLegend()
C_mjjjj_all.Print("mjjjj_all.png")

C_alpha_all = ROOT.TCanvas()
# C_alpha_all.SetLogy()
C_alpha_all.cd()
alpha_bkg.SetLineColor(ROOT.kBlue)
alpha_bkg.SetMarkerStyle(20)
alpha_bkg.Draw()
alpha_est.SetLineColor(ROOT.kRed)
alpha_est.Draw("hist same")
# C_alpha_all.BuildLegend()
C_alpha_all.Print("alpha_all.png")


# Qout.GetZaxis().SetRangeUser(-globalmax - 0.2, globalmax+0.2)
GoodFormat.pullFormat(Qout)

Cpull = TCanvas()
Cpull.cd()
Qout.SetStats(0)
Qout.Draw("COLZ")
Cpull.Print("PULL.png")

GoodFormat.pullFormat(tpull)

# tpull.GetZaxis().SetRangeUser(-globalmax - 0.19, globalmax+0.2)
tpull.GetZaxis().SetRangeUser(-6, 6)
tpull.SetStats(0)
tc3 = ROOT.TCanvas()
tc3.cd()
tpull.Draw("COLZ")
tc3.Print("GTPull.png")