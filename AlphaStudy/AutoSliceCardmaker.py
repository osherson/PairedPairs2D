import ROOT
from ROOT import *
import numpy
import os
gStyle.SetOptFit(1111)

def Template_Replace(F, O, R):
	with open(F, 'r') as file :
		filedata = file.read()
	filedata = filedata.replace(O, R)
	with open(F, 'w') as file:
		file.write(filedata)

def MakeAToy(Hi):
	R = numpy.random
	Ho = Hi.Clone("H"+str(R.randint(1)))
	for i in range(1,Ho.GetNbinsX()):
		v = Ho.GetBinContent(i)
		e = Ho.GetBinError(i)
		n = max(0,numpy.rint(v + (e*R.normal())))
		Ho.SetBinContent(i,int(n))
	Ho.Sumw2()
	Ho.SetStats(1)
	return Ho

def AutoSliceCardmaker(iBKG, SIG, FUNC, NAME, START, MC):
	if MC: BKG = MakeAToy(iBKG)
	else: BKG = iBKG
	
	mStart = BKG.GetMaximumBin() + int(START)
	start = BKG.GetBinLowEdge(mStart)
	end = BKG.GetBinCenter(BKG.GetNbinsX())
	ATLASFUNC = "TMath::Power(1.0-TMath::Power(x/13000.0, 1/3.), [0])/TMath::Power(x/13000.0, [1] + [2]*TMath::Log(x))"
	dumbfit = TF1("dumbfit", ATLASFUNC, start,end)
	BKG.Fit(dumbfit,"EM0R")
	sigfit = TF1("sigfit", "gaus", start, end)
	SIG.Fit(sigfit, "EM0R")
	
	x = RooRealVar("x_"+NAME, "Mjj", BKG.GetXaxis().GetBinLowEdge(mStart), BKG.GetXaxis().GetBinUpEdge(BKG.GetNbinsX()))
	background = RooDataHist("data_obs", "data_obs", RooArgList(x), BKG)
	
	p0 = RooRealVar("p0_"+NAME, "p0_"+NAME, dumbfit.GetParameter(0))
	p1 = RooRealVar("p1_"+NAME, "p1_"+NAME, dumbfit.GetParameter(1))
	p2 = RooRealVar("p2_"+NAME, "p2_"+NAME, dumbfit.GetParameter(2))
	mu = RooRealVar("mu_"+NAME, "mu_"+NAME, sigfit.GetParameter(1))
	si = RooRealVar("si_"+NAME, "si_"+NAME, sigfit.GetParameter(2))
	p0.setConstant(kFALSE)
	p1.setConstant(kFALSE)
	p2.setConstant(kFALSE)
	mu.setConstant(kTRUE)
	si.setConstant(kTRUE)
	fit = RooGenericPdf("bkg_"+NAME,"pow(1.0-pow(@0/13000.0, 1/3.), @1)/pow(@0/13000.0, @2 + @3*log(@0))", RooArgList(x,p0,p1,p2))
	signal = RooGaussian("sig", "signal peak", x, mu, si)
	fit.fitTo(background)
	frame = x.frame()
	background.plotOn(frame)
	signal.plotOn(frame)
	fit.plotOn(frame)
	fit.paramOn(frame)
	
	CANVAS = TCanvas(NAME, "", 1150, 450)
	CANVAS.Divide(2,1)
	CANVAS.cd(1)
	gPad.SetLogy()
	BKG.Draw("E")
	SIG.Draw("histsame")
	dumbfit.Draw("same")
	CANVAS.cd(2)
	gPad.SetLogy()
	frame.Draw()
	CANVAS.Print(NAME+".png")
	CANVAS.Print(NAME+".root")
	
	os.system("cp CombineCard.tpl CARD_"+NAME+".txt") 
	Template_Replace("CARD_"+NAME+".txt", "#NAME#", NAME) 
	Template_Replace("CARD_"+NAME+".txt", "#P0#", str(dumbfit.GetParameter(0))) 
	Template_Replace("CARD_"+NAME+".txt", "#P1#", str(dumbfit.GetParameter(1))) 
	Template_Replace("CARD_"+NAME+".txt", "#P2#", str(dumbfit.GetParameter(2))) 
	Template_Replace("CARD_"+NAME+".txt", "#P0E#", str(dumbfit.GetParError(0))) 
	Template_Replace("CARD_"+NAME+".txt", "#P1E#", str(dumbfit.GetParError(1))) 
	Template_Replace("CARD_"+NAME+".txt", "#P2E#", str(dumbfit.GetParError(2))) 
	Template_Replace("CARD_"+NAME+".txt", "#DATAINT#", str(BKG.Integral())) 
	Template_Replace("CARD_"+NAME+".txt", "#SIGINT#", str(SIG.Integral()))
	
	W = RooWorkspace(NAME)
	getattr(W, 'import')(background)
	getattr(W, 'import')(signal)
	getattr(W, 'import')(fit)
	W.writeToFile("workspace_"+NAME+".root")
	
if __name__ == "__main__":
	from optparse import OptionParser
	parser = OptionParser()
	parser.add_option("-b", "--bkgfile", dest="bkgF", help="file containing background distribution", metavar="BKGFILE")
	parser.add_option("-s", "--sigfile", dest="sigF", help="file containing signal disstribution", metavar="SIGFILE")
	parser.add_option("--bkg", dest="bkgN", help="name of background distribution", metavar="BKGNAME")
	parser.add_option("--name", dest="name", help="outputname", metavar="NAME")
	parser.add_option("--sig", dest="sigN", help="name of signal disstribution", metavar="SIGNAME")
	parser.add_option("--mc", action="store_true", dest="mc", help="is this a weighted MC?")
	parser.add_option("--start", dest="start", default = "3", help="fit starting point")
	(o, args) = parser.parse_args()
	
	bF = TFile(o.bkgF)
	B = bF.Get(o.bkgN)
	sF = TFile(o.sigF)
	S = sF.Get(o.sigN)
	AutoSliceCardmaker(B,S,"", o.name, o.start, o.mc)
