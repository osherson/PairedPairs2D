#
import ROOT
from ROOT import *
import numpy
import os
import sys
import Functions as FUNC
RDF = ROOT.RDataFrame

cutsnew = "evt_Masym < 0.1 && evt_Deta < 1.1 && dj1_dR < 2.0 && dj2_dR < 2.0 && evt_HTAK4 > 1200."
cutsold = "evt_Masym < 0.1 && evt_Deta < 1. && dj1_delta > 200. && dj2_delta > 200. && evt_HTAK4 > 1200."

def GetDists(F, L, N):
	C = ROOT.TChain("tree_nominal")
	C.Add(F)
	rdf = RDF(C)
	rdf = rdf.Define("total_weight", "weight_xsN*weight_PU*"+str(L))
	rdf = rdf.Define("evt_alpha", "evt_2JetM/evt_4JetM")
	RlazyNew = rdf.Filter(cutsnew).Histo2D((N+"tNew", ";Average Dijet Mass (GeV); #alpha", 100, FUNC.MakeNBinsFromMinToMax(100, 100.0, 4000.0), 30, FUNC.MakeNBinsFromMinToMax(30, 0.0, 0.5)),"evt_2JetM","evt_alpha","total_weight")
	RlazyOld = rdf.Filter(cutsold).Histo1D((N+"tOld", ";Average Dijet Mass (GeV)", 100, FUNC.MakeNBinsFromMinToMax(100, 100.0, 4000.0)),"evt_2JetM","total_weight")
	ROld = RlazyOld.GetValue()
	RNew = RlazyNew.GetValue()
	print "new/old: " + str(RNew.Integral()/ROld.Integral())
	return (ROld.Clone(N+"old"), RNew.Clone(N+"new"))

(Qo, Qn) = GetDists("/cms/xaastorage-2/PicoTrees/4JETS/2018/v6/QCD_HT_all_2018.root", 59740., "Q")
(S500o, S500n) = GetDists("/cms/xaastorage-2/PicoTrees/4JETS/2017/v6/RPV_M500_2017.root", 59740., "S500")
(S750o, S750n) = GetDists("/cms/xaastorage-2/PicoTrees/4JETS/2017/v6/RPV_M750_2017.root", 59740., "S750")
(S1000o, S1000n) = GetDists("/cms/xaastorage-2/PicoTrees/4JETS/2017/v6/RPV_M1000_2017.root", 59740., "S1000")
(S2000o, S2000n) = GetDists("/cms/xaastorage-2/PicoTrees/4JETS/2017/v6/RPV_M2000_2017.root", 59740., "S2000")

C = TCanvas()
C.Divide(2,1)
C.cd(1)
S500n.Draw("colz")
C.cd(2)
Qn.Draw("colz")
F = TFile("dists/Significance.root", "recreate")
F.cd()
Qo.Write()
Qn.Write()
S500o.Write()
S500n.Write()
S750o.Write()
S750n.Write()
S1000o.Write()
S1000n.Write()
S2000o.Write()
S2000n.Write()
F.Write()
F.Save()



for (o, n, N) in zip([S500o, S750o, S1000o, S2000o],[S500n, S750n, S1000n, S2000n], ["500", "750", "1000", "2000"]):
	srboS = 0.
	srbnS = 0.
	srboB = 0.
	srbnB = 0.
	for i in range(1,Qo.GetNbinsX()+1):
		if numpy.abs(i - o.GetMaximumBin()) < 4.:
			B = float(max(1.,Qo.GetBinContent(i)))
			S = float(o.GetBinContent(i))
			srboS+=S
			srboB+=B
		else: o.SetBinContent(i,0.)
	for i in range(1,Qn.GetNbinsX()+1):
		for j in range(1,Qn.GetNbinsY()+1):
			if numpy.abs(i - o.GetMaximumBin()) < 4.:
				B = float(max(1.,Qn.GetBinContent(i,j)))
				S = float(n.GetBinContent(i,j))
				srbnS+=S
				srbnB+=B
			else: n.SetBinContent(i,j,0.)
	SRBn = srbnS/numpy.sqrt(srbnB)
	SRBo = srboS/numpy.sqrt(srboB)
	print "RPV"+N+" ratio: " + str(SRBn/SRBo)


Qp = Qn.ProjectionX("Qp", 1, Qn.GetNbinsX())
S500p = S500n.ProjectionX("S500p", 1, Qn.GetNbinsX())
S750p = S750n.ProjectionX("S750p", 1, Qn.GetNbinsX())
S1000p = S1000n.ProjectionX("S1000p", 1, Qn.GetNbinsX())
S2000p = S2000n.ProjectionX("S2000p", 1, Qn.GetNbinsX())
	
for h in [Qo,Qp]:
	h.SetLineWidth(2)
	h.SetLineColor(kBlack)
	h.SetStats(0)
for h in [S500p,S500o]:
	h.SetLineWidth(2)
	h.SetLineColor(kViolet)
	h.SetStats(0)
for h in [S750p,S750o]:
	h.SetLineWidth(2)
	h.SetLineColor(kRed)
	h.SetStats(0)
for h in [S1000p,S1000o]:
	h.SetLineWidth(2)
	h.SetLineColor(kBlue)
	h.SetStats(0)
for h in [S2000p,S2000o]:
	h.SetLineWidth(2)
	h.SetLineColor(kMagenta)
	h.SetStats(0)
for h in [Qp, S500p, S750p, S1000p, S2000p]:
	h.SetLineStyle(2)

L = TLegend(0.5,0.5,0.89,0.89)
L.SetLineColor(0)
L.SetFillColor(0)
L.AddEntry(Qo, "QCD (old)", "L")
L.AddEntry(Qp, "QCD (new, projection of #alpha)", "L")
L.AddEntry(S500o, "RPV stop (500GeV) (old)", "L")
L.AddEntry(S500p, "RPV stop (500GeV) (new, projection of #alpha)", "L")
L.AddEntry(S750o, "RPV stop (750GeV) (old)", "L")
L.AddEntry(S750p, "RPV stop (750GeV) (new, projection of #alpha)", "L")
L.AddEntry(S1000o, "RPV stop (1000GeV) (old)", "L")
L.AddEntry(S1000p, "RPV stop (1000GeV) (new, projection of #alpha)", "L")
L.AddEntry(S2000o, "RPV stop (2000GeV) (old)", "L")
L.AddEntry(S2000p, "RPV stop (2000GeV) (new, projection of #alpha)", "L")

C = TCanvas()
C.cd()
C.SetLogy()
Qo.Draw("hist")
S500o.Draw("histsame")
S750o.Draw("histsame")
S1000o.Draw("histsame")
S2000o.Draw("histsame")
Qp.Draw("samehist")
S500p.Draw("histsame")
S750p.Draw("histsame")
S1000p.Draw("histsame")
S2000p.Draw("histsame")
L.Draw("same")



F.Close()
