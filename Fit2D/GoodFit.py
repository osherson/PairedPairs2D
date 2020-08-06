#
import ROOT
from ROOT import *
import numpy
import sys

InFile = TFile(sys.argv[1])
QCD = InFile.Get(sys.argv[2])
QCD.RebinY(int(sys.argv[3]))
QCD.RebinX(4)
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

p0 = []
p0e = []
p1 = []
p1e = []
p2 = []
p3 = []

Cs = TCanvas()
Cs.cd()
Cs.SetLogy()

y = []
x = []
z = []
ex = []
ey = []
ez = []

for j in range(1,QCD.GetNbinsY()):
	print "====---> "+ str(j)
	Q = QCD.ProjectionX("q_"+str(j),j,j)
	if not (Q.Integral() > 0): continue
	if QCD.GetYaxis().GetBinCenter(j) > 0.4: continue
	m = Q.GetMaximumBin() + int(sys.argv[4])
	for i in range(m,QCD.GetNbinsX()):
			y.append(QCD.GetYaxis().GetBinCenter(j))
			x.append(QCD.GetXaxis().GetBinCenter(i))
			ey.append(QCD.GetYaxis().GetBinWidth(j)/2.)
			ex.append(QCD.GetXaxis().GetBinWidth(i)/2.)
			z.append(Q.GetBinContent(i))
			ez.append(Q.GetBinError(i))
	start = Q.GetBinLowEdge(m)
	end = Q.GetBinCenter(QCD.GetNbinsX())
	fit = TF1("tF", "TMath::Power(1-(x/13000.0),[0])/(TMath::Power(x/13000.0,[1]))", start,end)
	Q.Fit(fit,"EM0R")
	ndof = fit.GetNDF()
	chi2 = fit.GetChisquare()
	if (ndof < 5) or (chi2/ndof > 1.5): continue
	starts.append(start)
	norm.append(Q.Integral())
	chiN.append(chi2/ndof)
	p0.append(fit.GetParameter(0))
	p1.append(fit.GetParameter(1))
	p0e.append(fit.GetParError(0))
	p1e.append(fit.GetParError(1))
	#p2.append(fit.GetParameter(2))
	#p3.append(fit.GetParameter(3))
	yalph.append(QCD.GetYaxis().GetBinCenter(j))
	yalphe.append(QCD.GetYaxis().GetBinWidth(j)/2.)
	Q.Draw()
	fit.Draw("same")
	Cs.Print("Slice_"+str(j)+".png")

print yalph

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
	g.Draw("AP")
	f.Draw("same")
	C.Print(N+".png")
	return f
MakePlot(starts, "start")
MakePlot(chiN, "chi2N")
MakePlot(norm, "norm")
oP0 = MakeEPlot(p0, "p0", p0e, 3)
oP1 = MakeEPlot(p1, "p1", p1e, 3)

G = TGraph2DErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(z), numpy.array(ex), numpy.array(ey), numpy.array(ez))


X = "x/13000.0"
P0 = "([0] + [1]*y + [2]*y*y + [3]*y*y*y)"
P1 = "([4] + [5]*y + [6]*y*y + [7]*y*y*y)"

FitFunc = "TMath::Power((1.-"+X+"),"+P0+")/TMath::Power("+X+","+P1+")"

F = TF2("hardFit", FitFunc, 0., 3000., 0.1, 0.4)

F.SetParameter(0, oP0.GetParameter(0))
F.SetParameter(1, oP0.GetParameter(1))
F.SetParameter(2, oP0.GetParameter(2))
F.SetParameter(3, oP0.GetParameter(3))
F.SetParameter(4, oP1.GetParameter(0))
F.SetParameter(5, oP1.GetParameter(1))
F.SetParameter(6, oP1.GetParameter(2))
F.SetParameter(7, oP1.GetParameter(3))

G.Fit(F)
ndof = F.GetNDF()
chi2 = F.GetChisquare()
print "Chi2Ndof = "  + str(chi2)+"/"+str(ndof)

C = TCanvas()
C.cd()
G.Draw("AE")
F.Draw("same")

Qout = QCD.Clone("out")
Qout.Add(Qout,-1.)
Qout.SetStats(0)

PullHist = TH1F("PullHist", ";pull;bins", 1200, -150, 150)

Qout.GetZaxis().SetRangeUser(-6.,6.)
for j in range(1,QCD.GetNbinsY()):
	Q = QCD.ProjectionX("q_"+str(j),j,j)
	if (Q.Integral() > 0) and (QCD.GetYaxis().GetBinCenter(j)) < 0.4:
		m = Q.GetMaximumBin() + int(sys.argv[4])
		for i in range(m,QCD.GetNbinsX()):
				q = Q.GetBinContent(i)
				eq = Q.GetBinError(i)
				if q == 0: eq = 1.4
				est = F.Eval(Q.GetBinCenter(i), QCD.GetYaxis().GetBinCenter(j))
				P = (q - est)/eq
				Qout.SetBinContent(i,j,P)
				PullHist.Fill(P)
		for i in range(1,m-1):
			Qout.SetBinContent(i,j,0)
	else:
		for i in range(1,QCD.GetNbinsX()):
			Qout.SetBinContent(i,j,0)

CPH = TCanvas()
CPH.cd()
PullHist.Draw("hist")		

CFO = TCanvas()
CFO.cd()
F.Draw()

Cpull = TCanvas()
Cpull.cd()
Qout.Draw("colz")
Cpull.Print("PULL.root")

Cpull2 = TCanvas()
Cpull2.Divide(2,1)
Cpull2.cd(1)
Qout.ProjectionX("xp").Draw("hist")
Cpull2.cd(2)
Qout.ProjectionY("yp").Draw("hist")

			



