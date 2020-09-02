import ROOT
from ROOT import *
import numpy
import os
import sys
import Functions as FUNC

iF = TFile(sys.argv[1])
B = iF.Get("BKG")
F = iF.Get("FIT")
S = iF.Get("SLICE")
PULL = B.Clone("PULL")
DIFF = F.Clone("DIFF")
PULL.Add(F, -1.)
DIFF.Add(S, -1.)


for h,n in zip([B,F,S],["background", "fit", "slices"]):
	C = TCanvas()
	C.cd()
	h.Draw("colz")
	C.Print("plots/"+n+".png")
	C.SetLogz()
	C.Print("plots/"+n+"_logz.png")

for h in [PULL, DIFF]:
	for j in range(1,B.GetNbinsY()):
		for i in range(1,B.GetNbinsX()):
			if B.GetBinContent(i,j) > 0:
				e = max(1.,B.GetBinError(i,j))
				h.SetBinContent(i,j,h.GetBinContent(i,j)/e)
			else:
				h.SetBinContent(i,j,0)
	h.GetZaxis().SetRangeUser(-5.,5.)
				
PULL.GetZaxis().SetRangeUser(-5.,5.)

for h,n in zip([PULL, DIFF],["diff", "fitdiff"]):
	C = TCanvas()
	C.cd()
	h.Draw("colz")
	C.Print("plots/"+n+".png")
	
Bx = B.ProjectionX()
By = B.ProjectionY()
Fx = F.ProjectionX()
Fy = F.ProjectionY()

for (b,f,n) in zip([Bx,By],[Fx,Fy],["x","y"]):
	b.SetMarkerStyle(20)
	b.SetMarkerColor(1)
	b.SetLineColor(1)
	f.SetLineColor(kRed)
	C = TCanvas()
	C.cd()
	C.SetLogy()
	b.Draw("e")
	f.Draw("histsame")
	C.Print("plots/Proj"+n+".png")
	
aBl = int(sys.argv[2])
aBh = int(sys.argv[3])


	
Bx = B.ProjectionX("bPx", aBl, aBh)
Fx = F.ProjectionX("fPx", aBl, aBh)

InFile = TFile("dists/FitHists2D.root")
SH = InFile.Get("SHm2a")
SM = InFile.Get("SMm2a")
SL = InFile.Get("SLm2a")

for h in [SH,SM,SL]: h.Scale(100./h.Integral())

SHx = SH.ProjectionX("bSHx", aBl, aBh)
SMx = SM.ProjectionX("bSMx", aBl, aBh)
SLx = SL.ProjectionX("bSLx", aBl, aBh)


Bx.SetMarkerStyle(20)
Bx.SetMarkerColor(1)
Bx.SetLineColor(1)
Fx.SetLineColor(kRed)
SHx.SetLineColor(kBlue)
SMx.SetLineColor(kViolet)
SLx.SetLineColor(kGreen)
C = TCanvas()
C.cd()
C.SetLogy()
Bx.Draw("e")
Fx.Draw("histsame")
SHx.Draw("histsame")
SMx.Draw("histsame")
SLx.Draw("histsame")
C.Print("plots/ProjZoom.png")

BH4 = TH1F("BH4", ";M_{jjjj} (GeV);events", 50, 0, 10000)
FH4 = TH1F("FH4", ";M_{jjjj} (GeV);events", 50, 0, 10000)
BM4 = TH1F("BM4", ";M_{jjjj} (GeV);events", 50, 0, 10000)
FM4 = TH1F("FM4", ";M_{jjjj} (GeV);events", 50, 0, 10000)
BL4 = TH1F("BL4", ";M_{jjjj} (GeV);events", 50, 0, 10000)
FL4 = TH1F("FL4", ";M_{jjjj} (GeV);events", 50, 0, 10000)

SH4 = TH1F("SH4", ";M_{jjjj} (GeV);events", 50, 0, 10000)
SM4 = TH1F("SM4", ";M_{jjjj} (GeV);events", 50, 0, 10000)
SL4 = TH1F("SL4", ";M_{jjjj} (GeV);events", 50, 0, 10000)
SH4.SetLineColor(kBlue)
SM4.SetLineColor(kViolet)
SL4.SetLineColor(kGreen)

for j in range(1,B.GetNbinsY()):
	for i in range(1,B.GetNbinsX()):
		m4 = B.GetXaxis().GetBinCenter(i)/B.GetYaxis().GetBinCenter(j)
		if B.GetXaxis().GetBinLowEdge(i) > 400. and B.GetXaxis().GetBinUpEdge(i) < 600.:
			BL4.Fill(m4, B.GetBinContent(i,j))
			FL4.Fill(m4, F.GetBinContent(i,j))
			SL4.Fill(m4, SL.GetBinContent(i,j))
		if B.GetXaxis().GetBinLowEdge(i) > 700. and B.GetXaxis().GetBinUpEdge(i) < 1100.:
			BM4.Fill(m4, B.GetBinContent(i,j))
			FM4.Fill(m4, F.GetBinContent(i,j))
			SM4.Fill(m4, SM.GetBinContent(i,j))
		if B.GetXaxis().GetBinLowEdge(i) > 1400. and B.GetXaxis().GetBinUpEdge(i) < 2200.:
			BH4.Fill(m4, B.GetBinContent(i,j))
			FH4.Fill(m4, F.GetBinContent(i,j))
			SH4.Fill(m4, SH.GetBinContent(i,j))

FUNC.FindAndSetMax(BH4, FH4, SH4)
FUNC.FindAndSetMax(BM4, FM4, SM4)
FUNC.FindAndSetMax(BL4, FL4, SL4)

for b in [BH4,BM4,BL4]:
	b.SetMarkerStyle(20)
	b.SetMarkerColor(1)
	b.SetLineColor(1)
	for i in range(1,b.GetNbinsX()):
		if b.GetBinContent(i) > 0:
			b.SetBinError(i, numpy.sqrt(b.GetBinContent(i)))
for f in [FH4, FM4, FL4]:
	f.SetLineColor(kRed)

CH = TCanvas()
CH.cd()
CH.SetLogy()
BH4.Draw("e")
FH4.Draw("histsame")
SH4.Draw("histsame")
CH.Print("plots/M4H.png")
CM = TCanvas()
CM.cd()
CM.SetLogy()
BM4.Draw("e")
FM4.Draw("histsame")
SM4.Draw("histsame")
CM.Print("plots/M4M.png")
CL = TCanvas()
CL.cd()
CL.SetLogy()
BL4.Draw("e")
FL4.Draw("histsame")
SL4.Draw("histsame")
CL.Print("plots/M4L.png")

	
