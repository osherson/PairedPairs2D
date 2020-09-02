# 
import ROOT
from ROOT import *
import numpy
import sys

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

InFile = TFile(sys.argv[1])
QCD = InFile.Get("AvJJ_alpha_rebin")
QCD.RebinX(4)
QCD.RebinY(4)
QCD = MakeAToy(QCD)
QCD.SetTitle("")
QCD.SetStats(0)
C_QCD_uncut = TCanvas()
C_QCD_uncut.cd()
C_QCD_uncut.SetLogz()
QCD.Draw("colz")
C_QCD_uncut.Print("QCD_raw.png")

mV = []
m = []
yV = []
y = []

for j in range(1,QCD.GetNbinsY()):
	Q = QCD.ProjectionX("q_"+str(j),j,j)
	if not (Q.Integral() > 0): continue
	m.append(float(Q.GetMaximumBin()))
	mV.append(Q.GetBinLowEdge(Q.GetMaximumBin()))
	yV.append(QCD.GetYaxis().GetBinCenter(j))
	y.append(float(j))
	
Gm = TGraph(len(m), numpy.array(y), numpy.array(m))
GmV = TGraph(len(m), numpy.array(yV), numpy.array(mV))
GmV.SetTitle(";#alpha bin center;m_{jj} in max bin")
FmV = TF1("FmV", "pol2", min(yV), max(yV))
Fm = TF1("Fm", "pol2", min(y), max(y))
GmV.Fit(FmV)
Gm.Fit(Fm)

C_MSV = TCanvas()
C_MSV.cd()
GmV.Draw("APL")
FmV.Draw("same")
C_MSV.Print("maxbin_vs_alpha.png")

C_MS = TCanvas()
C_MS.cd()
Gm.Draw("APL")
Fm.Draw("same")

x = []
y = []
z = []
ex = []
ey = []
ez = []
a = []
ae = []

p0 = []
p0e = []
p1 = []
p1e = []
p2 = []
p2e = []
c2 = []

BKG = QCD.Clone("BKG")
BKG.Scale(0)
FIT = QCD.Clone("FIT")
FIT.Scale(0)
SLICE = QCD.Clone("SLICE")
SLICE.Scale(0)

tC = TCanvas()
tC.cd()
tC.SetLogy()
for j in range(1,QCD.GetNbinsY()):
	Q = QCD.ProjectionX("Q_"+str(j),j,j)
	if not (Q.Integral() > 0): continue
	Q.Draw("e")
	m = max(0, int(sys.argv[2]) + int(numpy.rint(Fm.Eval(j))))
	for i in range(m,QCD.GetNbinsX()):
		y.append(QCD.GetYaxis().GetBinCenter(j))
		x.append(QCD.GetXaxis().GetBinCenter(i))
		#ey.append(QCD.GetYaxis().GetBinWidth(j)/2.)
		#ex.append(QCD.GetXaxis().GetBinWidth(i)/2.)
		ey.append(0.)
		ex.append(0.)
		z.append(Q.GetBinContent(i))
		ez.append(Q.GetBinError(i))
		BKG.SetBinContent(i,j,Q.GetBinContent(i))
		BKG.SetBinError(i,j,Q.GetBinError(i))
	start = Q.GetBinLowEdge(m)
	end = Q.GetBinCenter(QCD.GetNbinsX())
	fit = TF1("tF", "TMath::Power(1-(x/13000.0),[0])/(TMath::Power(x/13000.0,[1]+([2]*TMath::Log(x/13000.0))))", start,end)
	Q.Fit(fit,"EM0R")
	fit.Draw("same")
	tC.Print("S"+str(j)+".png")
	for i in range(m,QCD.GetNbinsX()):
		SLICE.SetBinContent(i,j,fit.Eval(Q.GetBinCenter(i)))
	ndof = fit.GetNDF()
	chi2 = fit.GetChisquare()
	
	if (ndof < 5) or (chi2/ndof > 3.): continue
	c2.append(chi2/ndof)
	p0.append(fit.GetParameter(0))
	p1.append(fit.GetParameter(1))
	p2.append(fit.GetParameter(2))
	p0e.append(fit.GetParError(0))
	p1e.append(fit.GetParError(1))
	p2e.append(fit.GetParError(2))
	a.append(QCD.GetYaxis().GetBinCenter(j))
	ae.append(QCD.GetYaxis().GetBinWidth(j)/2.)

Bx = BKG.ProjectionX()
By = BKG.ProjectionY()
Sx = SLICE.ProjectionX()
Sy = SLICE.ProjectionY()


gC2 = TGraph(len(a), numpy.array(a), numpy.array(c2))
gC2.SetTitle(";#alpha bin center;#chi^{2}/Ndof")
C_chi2 = TCanvas()
C_chi2.cd()
gC2.Draw("APL")
C_chi2.Print("SliceChi2Ndof.png")
def MakeEPlot(V, N, E, O):
	g = TGraphErrors(len(a), numpy.array(a), numpy.array(V), numpy.array(ae), numpy.array(E))
	f = TF1("f", "pol"+str(O), min(a), max(a))
	g.Fit(f)
	g.SetTitle(N)
	C = TCanvas()
	C.cd()
	g.Draw("AL")
	f.Draw("same")
	C.Print(N+".png")
	return f
oP0 = MakeEPlot(p0, "p0", p0e, 3)
oP1 = MakeEPlot(p1, "p1", p1e, 3)
oP2 = MakeEPlot(p2, "p2", p2e, 3)

X = "x/13000.0"
P0 = "([0] + [1]*y + [2]*y*y + [3]*y*y*y)"
P1 = "([4] + [5]*y + [6]*y*y + [7]*y*y*y)"
P2 = "([8] + [9]*y + [10]*y*y + [11]*y*y*y)"
FitFunc = "TMath::Power((1.-"+X+"),"+P0+")/TMath::Power("+X+","+P1+"+("+P2+"*TMath::Log("+X+")))"
F = TF2("hardFit", FitFunc, 0., 3000., 0.1, 0.4)
F.SetParameter(0, oP0.GetParameter(0))
F.SetParameter(1, oP0.GetParameter(1))
F.SetParameter(2, oP0.GetParameter(2))
F.SetParameter(3, oP0.GetParameter(3))
F.SetParameter(4, oP1.GetParameter(0))
F.SetParameter(5, oP1.GetParameter(1))
F.SetParameter(6, oP1.GetParameter(2))
F.SetParameter(7, oP1.GetParameter(3))
F.SetParameter(8, oP2.GetParameter(0))
F.SetParameter(9, oP2.GetParameter(1))
F.SetParameter(10, oP2.GetParameter(2))
F.SetParameter(11, oP2.GetParameter(3))

G = TGraph2DErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(z), numpy.array(ex), numpy.array(ey), numpy.array(ez))
G.Fit(F, "EM0R")
ndof = F.GetNDF()
chi2 = F.GetChisquare()
print "Chi2Ndof = "  + str(chi2)+"/"+str(ndof) + " = " + str(chi2/ndof)

for n in range(len(z)):
		b = BKG.FindBin(x[n],y[n])
		est = F.Eval(x[n],y[n])
		FIT.SetBinContent(b, est)

Fx = FIT.ProjectionX()
Fy = FIT.ProjectionY()



PULL = BKG.Clone("pull")
PULL.Add(FIT,-1.)
for j in range(1,QCD.GetNbinsY()):
	for i in range(1,QCD.GetNbinsX()):
		if FIT.GetBinContent(i,j) > 0:
			e = max(1.,BKG.GetBinError(i,j))
			PULL.SetBinContent(i,j,PULL.GetBinContent(i,j)/e)
PULL.GetZaxis().SetRangeUser(-5.,5.)
			

C = TCanvas()
C.cd()
PULL.Draw("colz")
C.Print("pull.png")

bL = BKG.ProjectionX("bL", 1,int(BKG.GetNbinsY()/3))
bM = BKG.ProjectionX("bM", int(BKG.GetNbinsY()/3),int(2*BKG.GetNbinsY()/3))
bH = BKG.ProjectionX("bH", int(2*BKG.GetNbinsY()/3),BKG.GetNbinsY())
fL = FIT.ProjectionX("fL", 1,int(BKG.GetNbinsY()/3))
fM = FIT.ProjectionX("fM", int(BKG.GetNbinsY()/3),int(2*BKG.GetNbinsY()/3))
fH = FIT.ProjectionX("fH", int(2*BKG.GetNbinsY()/3),BKG.GetNbinsY())
sL = SLICE.ProjectionX("sL", 1,int(BKG.GetNbinsY()/3))
sM = SLICE.ProjectionX("sM", int(BKG.GetNbinsY()/3),int(2*BKG.GetNbinsY()/3))
sH = SLICE.ProjectionX("sH", int(2*BKG.GetNbinsY()/3),BKG.GetNbinsY())

pL = bL.Clone("pL")
pM = bM.Clone("pM")
pH = bH.Clone("pH")
pL.Add(fL, -1.)
pM.Add(fM, -1.)
pH.Add(fH, -1.)
for i in range(1,pL.GetNbinsX()):
	for p,f,b in zip([pL,pM,pH],[fL,fM,fH],[bL,bM,bH]):
		if f.GetBinContent(i) > 0:
			e = max(1.,b.GetBinError(i))
			p.SetBinContent(i,p.GetBinContent(i)/e)
pLs = bL.Clone("pL")
pMs = bM.Clone("pM")
pHs = bH.Clone("pH")
pLs.Add(sL, -1.)
pMs.Add(sM, -1.)
pHs.Add(sH, -1.)
for i in range(1,pL.GetNbinsX()):
	for p,s,b in zip([pLs,pMs,pHs],[sL,sM,sH],[bL,bM,bH]):
		if s.GetBinContent(i) > 0:
			e = max(1.,b.GetBinError(i))
			p.SetBinContent(i,p.GetBinContent(i)/e)

for f in [Fx,Fy, fL, fM, fH, pL, pM, pH]:
	f.SetLineColor(kBlue)
	f.SetLineWidth(2)
for b in [Bx, By, bL, bM, bH]:
	b.SetMarkerStyle(20)
	b.SetMarkerColor(1)
	b.SetLineColor(1)
	b.SetStats(0)
for s in [Sx,Sy, sL, sM, sH, pLs, pMs, pHs]:
	s.SetLineColor(kRed)
	s.SetLineWidth(2)
	
C_slices = TCanvas("C_slices", "", 1200, 550)
C_slices.Divide(2,1)
C_slices.cd(1)
gPad.SetLogy()
Bx.Draw("e")
Sx.Draw("histsame")
Fx.Draw("histsame")
C_slices.cd(2)
By.Draw("e")
Sy.Draw("histsame")
Fy.Draw("histsame")
C_slices.Print("SliceAgree.png")

C3 = TCanvas("C3", "", 1000, 300)
C3.Divide(3,1)
C3.cd(1)
gPad.SetLogy()
bL.Draw("e")
fL.Draw("histsame")
sL.Draw("histsame")
C3.cd(2)
gPad.SetLogy()
bM.Draw("e")
fM.Draw("histsame")
sM.Draw("histsame")
C3.cd(3)
gPad.SetLogy()
bH.Draw("e")
fH.Draw("histsame")
sH.Draw("histsame")
C3.Print("Breakdowns.png")

for p in [pL,pM,pH]:
	p.GetYaxis().SetRangeUser(-5.,5.)
C3p = TCanvas("C3p", "", 1000, 300)
C3p.Divide(3,1)
C3p.cd(1)
pL.Draw("hist")
pLs.Draw("histsame")
C3p.cd(2)
pM.Draw("hist")
pMs.Draw("histsame")
C3p.cd(3)
pH.Draw("hist")
pHs.Draw("histsame")
C3p.Print("BreakdownsPulls.png")
