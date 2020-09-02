#
import ROOT
from ROOT import *
import numpy
import os
import sys
import Functions as FUNC
RDF = ROOT.RDataFrame
gROOT.SetBatch(kTRUE)

lumi = 59740 # in pb
cuts = "evt_Masym < 0.1 && evt_Deta < 1.1 && dj1_dR < 2.0 && dj2_dR < 2.0 && evt_4JetM > 1530." # phase space of interest

def GetM2jVM4j(F, L, N):
	C = ROOT.TChain("tree_nominal")
	C.Add(F)
	rdf = RDF(C)
	rdf = rdf.Define("total_weight", "weight_xsN*weight_PU*"+str(L))
	Rlazy = rdf.Filter(cuts).Histo2D((N+"temp", ";Average Dijet Mass (GeV);Four Jet Mass (GeV)", 390, FUNC.MakeNBinsFromMinToMax(390, 100.0, 4000.0), 425, FUNC.MakeNBinsFromMinToMax(425, 1500.0, 10000.0)),"evt_2JetM","evt_4JetM","total_weight")
	R = Rlazy.GetValue()
	R.SetStats(0)
	Cv = TCanvas()
	Cv.cd()
	R.Draw("colz")
	Cv.Print("plots/m2m4_"+N+".png")
	Cv.SetLogz()
	Cv.Print("plots/m2m4_logz_"+N+".png")
	Cv.SetLogy()
	Cv.SetLogx()
	Cv.Print("plots/m2m4_logxyz_"+N+".png")
	return R.Clone(N+"m2m4")
	
QCD = GetM2jVM4j("/cms/xaastorage-2/PicoTrees/4JETS/2018/v6/QCD_HT_all_2018.root", lumi, "QCD")
SH = GetM2jVM4j("/cms/xaastorage-2/PicoTrees/4JETS/2016/v6/Diquark_S8400_chi2000_2016.root", lumi, "SH")
SM = GetM2jVM4j("/cms/xaastorage-2/PicoTrees/4JETS/2016/v6/Diquark_S4000_chi1000_2016.root", lumi, "SM")
SL = GetM2jVM4j("/cms/xaastorage-2/PicoTrees/4JETS/2016/v6/Diquark_S2000_chi500_2016.root", lumi, "SL")

def GetAlphaBinning(F, L, N): # here N is not the name, it's the number of bins we want to split up into
	C = ROOT.TChain("tree_nominal")
	C.Add(F)
	rdf = RDF(C)
	rdf = rdf.Define("total_weight", "weight_xsN*weight_PU*"+str(L))
	rdf = rdf.Define("evt_alpha", "evt_2JetM/evt_4JetM")
	Rlazy = rdf.Filter(cuts).Histo1D(("temp", ";#alpha", 100000, FUNC.MakeNBinsFromMinToMax(100000, 0.15, 0.35)),"evt_alpha","total_weight")
	R = Rlazy.GetValue()
	total = R.Integral()
	binedges = [0.15]
	binsum = 0
	for n in range(1, R.GetNbinsX()+1):
		binsum+=R.GetBinContent(n)
		if binsum >= total/N:
			binsum = 0
			binedges.append(R.GetXaxis().GetBinUpEdge(n))
	binedges.append(0.35)
	Rlazy2 = rdf.Filter(cuts).Histo1D(("temp2", ";#alpha", N, numpy.array(binedges)),"evt_alpha","total_weight")
	R2 = Rlazy2.GetValue()
	Cv = TCanvas()
	Cv.cd()
	R2.Draw("hist")
	Cv.Print("plots/AlphaBins.png")
	Cv.SetLogy()
	Cv.Print("plots/AlphaBins_logy.png")
	return numpy.array(binedges)

Bins = GetAlphaBinning("/cms/xaastorage-2/PicoTrees/4JETS/2018/v6/QCD_HT_all_2018.root", lumi, 25)

def GetM2jVvA(F, L, N, B):
	C = ROOT.TChain("tree_nominal")
	C.Add(F)
	rdf = RDF(C)
	rdf = rdf.Define("total_weight", "weight_xsN*weight_PU*"+str(L))
	rdf = rdf.Define("evt_alpha", "evt_2JetM/evt_4JetM")
	Rlazy = rdf.Filter(cuts).Histo2D((N+"temp", ";Average Dijet Mass (GeV); #alpha", 117, FUNC.MakeNBinsFromMinToMax(117, 100.0, 4000.0), len(B)-1, B),"evt_2JetM","evt_alpha","total_weight")
	R = Rlazy.GetValue()
	R.SetStats(0)
	Cv = TCanvas()
	Cv.cd()
	R.Draw("colz")
	Cv.Print("plots/m2va_"+N+".png")
	Cv.SetLogz()
	Cv.Print("plots/m2va_logz_"+N+".png")
	Cv.SetLogy()
	Cv.SetLogx()
	Cv.Print("plots/m2va_logxyz_"+N+".png")
	return R.Clone(N+"m2a_v")
	
def GetM2jVA(F, L, N):
	C = ROOT.TChain("tree_nominal")
	C.Add(F)
	rdf = RDF(C)
	rdf = rdf.Define("total_weight", "weight_xsN*weight_PU*"+str(L))
	rdf = rdf.Define("evt_alpha", "evt_2JetM/evt_4JetM")
	#Rlazy = rdf.Filter(cuts).Histo2D((N+"temp", ";Average Dijet Mass (GeV); #alpha", 117, FUNC.MakeNBinsFromMinToMax(117, 100.0, 4000.0), 25, FUNC.MakeNBinsFromMinToMax(25, 0.15, 0.35)),"evt_2JetM","evt_alpha","total_weight")
	Rlazy = rdf.Filter(cuts).Histo2D((N+"temp", ";Average Dijet Mass (GeV); #alpha", 100, FUNC.MakeNBinsFromMinToMax(100, 100.0, 4000.0), 14, FUNC.MakeNBinsFromMinToMax(14, 0.15, 0.35)),"evt_2JetM","evt_alpha","total_weight")
	R = Rlazy.GetValue()
	R.SetStats(0)
	Cv = TCanvas()
	Cv.cd()
	R.Draw("colz")
	Cv.Print("plots/m2a_"+N+".png")
	Cv.SetLogz()
	Cv.Print("plots/m2a_logz_"+N+".png")
	Cv.SetLogy()
	Cv.SetLogx()
	Cv.Print("plots/m2a_logxyz_"+N+".png")
	return R.Clone(N+"m2a")
	
vQCD = GetM2jVvA("/cms/xaastorage-2/PicoTrees/4JETS/2018/v6/QCD_HT_all_2018.root", lumi, "QCD", Bins)
vSH = GetM2jVvA("/cms/xaastorage-2/PicoTrees/4JETS/2016/v6/Diquark_S8400_chi2000_2016.root", lumi, "SH", Bins)
vSM = GetM2jVvA("/cms/xaastorage-2/PicoTrees/4JETS/2016/v6/Diquark_S4000_chi1000_2016.root", lumi, "SM", Bins)
vSL = GetM2jVvA("/cms/xaastorage-2/PicoTrees/4JETS/2016/v6/Diquark_S2000_chi500_2016.root", lumi, "SL", Bins)
QCD = GetM2jVA("/cms/xaastorage-2/PicoTrees/4JETS/2018/v6/QCD_HT_all_2018.root", lumi, "QCD")
SH = GetM2jVA("/cms/xaastorage-2/PicoTrees/4JETS/2016/v6/Diquark_S8400_chi2000_2016.root", lumi, "SH")
SM = GetM2jVA("/cms/xaastorage-2/PicoTrees/4JETS/2016/v6/Diquark_S4000_chi1000_2016.root", lumi, "SM")
SL = GetM2jVA("/cms/xaastorage-2/PicoTrees/4JETS/2016/v6/Diquark_S2000_chi500_2016.root", lumi, "SL")

oF = TFile("dists/FitHists2DWB.root","recreate")
oF.cd()
vQCD.Write()
vSH.Write()
vSM.Write()
vSL.Write()
QCD.Write()
SH.Write()
SM.Write()
SL.Write()
oF.Write()
oF.Save()
oF.Close()
