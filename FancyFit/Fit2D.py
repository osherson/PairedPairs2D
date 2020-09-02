#
import ROOT
from ROOT import *
import numpy
import os
import sys
import Functions as FUNC
from optparse import OptionParser
gROOT.SetBatch(kTRUE)

parser = OptionParser()
parser.add_option("-v", action="store_true", dest="varbins", help="use variable bins")
parser.add_option("-t", action="store_true", dest="toy", help="use toys instead of scaled MC")
parser.add_option("-f", "--fit", dest="Pfit", default = "P2", help="Mjj fit function(s)", metavar="FILE")
parser.add_option("-r", "--reject", dest="reject", default = "10", help="chi2 to reject slice", metavar="FILE")
parser.add_option("-p", "--polfit", dest="Ofit", default = "3", help="alpha dimension poly order", metavar="FILE")
parser.add_option("-b", "--buffer", dest="buffer", default = "3", help="how many bins past the max do we start in a slice", metavar="FILE")
(o, args) = parser.parse_args()

InFile = TFile("dists/FitHists2DWB.root")
if o.toy:
	QCD = FUNC.MakeAToy(InFile.Get("QCDm2a"+("_v" if o.varbins else "")))
	SH = FUNC.MakeAToy(InFile.Get("SHm2a"+("_v" if o.varbins else "")))
	SM = FUNC.MakeAToy(InFile.Get("SMm2a"+("_v" if o.varbins else "")))
	SL = FUNC.MakeAToy(InFile.Get("SLm2a"+("_v" if o.varbins else "")))
else:
	QCD = InFile.Get("QCDm2a"+("_v" if o.varbins else ""))
	SH = InFile.Get("SHm2a"+("_v" if o.varbins else ""))
	SM = InFile.Get("SMm2a"+("_v" if o.varbins else ""))
	SL = InFile.Get("SLm2a"+("_v" if o.varbins else ""))

max2d = QCD.GetBinContent(QCD.GetMaximumBin())

tC = TCanvas()
tC.cd()
tC.SetLogy()

BKG = QCD.Clone("BKG")
BKG.Scale(0)
FIT = QCD.Clone("FIT")
FIT.Scale(0)
SLICE = QCD.Clone("SLICE")
SLICE.Scale(0)

a = [[],[]]
x = [[],[]]
y = [[],[]]
z = [[],[]]
c2 = []
nPar = FUNC.WhichFit(o.Pfit)[1]
Ps = []
for p in range(nPar): Ps.append([[],[]])
# Get the individual slices and fit them with the chosen functions:
for j in range(1,QCD.GetNbinsY()+1):
	Q = QCD.ProjectionX("Q_"+str(j),j,j)
	Q.GetYaxis().SetRangeUser(0.05,max2d*2.5)
	Q.SetStats(0)
	m = Q.GetMaximumBin() + int(o.buffer)
	for i in range(m,QCD.GetNbinsX()+1):
		y[0].append(QCD.GetYaxis().GetBinCenter(j))
		x[0].append(QCD.GetXaxis().GetBinCenter(i))
		y[1].append(0.)
		x[1].append(0.)
		z[0].append(Q.GetBinContent(i))
		z[1].append(Q.GetBinError(i))
		BKG.SetBinContent(i,j,Q.GetBinContent(i))
		BKG.SetBinError(i,j,Q.GetBinError(i))
	start = Q.GetBinLowEdge(m)
	end = Q.GetBinCenter(QCD.GetNbinsX())
	fit = TF1("tF", FUNC.WhichFit(o.Pfit)[0], start,end)
	Q.Fit(fit,"EM0R")
	for i in range(m,QCD.GetNbinsX()+1):
		SLICE.SetBinContent(i,j,fit.Eval(Q.GetBinCenter(i)))
	Q.SetMarkerStyle(20)
	Q.SetMarkerColor(1)
	Q.SetLineColor(1)
	fit.SetLineColor(kRed)
	Q.Draw("e")
	fit.Draw("same")
	L = TLegend(0.5,0.5,0.89,0.89)
	L.SetLineColor(0)
	L.SetFillColor(0)
	L.SetHeader("#alpha#in("+str(QCD.GetYaxis().GetBinLowEdge(j))+","+str(QCD.GetYaxis().GetBinUpEdge(j))+")")
	L.AddEntry(Q, "QCD MC"+(" (toy)" if o.toy else ""), "P")
	L.AddEntry(fit, o.Pfit + " fit", "L")
	L.Draw("same")
	tC.Print("plots/IndividualSlice_"+str(j)+("_v.png" if o.varbins else ".png"))
	ndof = fit.GetNDF()
	chi2 = fit.GetChisquare()
	c2.append(chi2/ndof)
	if (chi2/ndof) > float(o.reject): continue
	a[0].append(QCD.GetYaxis().GetBinCenter(j))
	a[1].append(0.)
	for p in range(nPar):
		Ps[p][0].append(fit.GetParameter(p))
		Ps[p][1].append(fit.GetParError(p))
# Show the Chi2/Ndof for the individual slices
gC2 = TGraph(len(a[0]), numpy.array(a[0]), numpy.array(c2))
gC2.SetLineWidth(2)
gC2.SetLineColor(kBlue)
gC2.SetTitle("#chi^{2}/Ndof from "+o.Pfit+" fit in individual slices of #alpha;#alpha bin center;#chi^{2}/Ndof")
C_chi2 = TCanvas()
C_chi2.cd()
gC2.Draw("APC")
C_chi2.Print("plots/SliceChi2Ndof.png")
# throw up some parameter plots
def MakeEPlot(V, E, N):
	g = TGraphErrors(len(a[0]), numpy.array(a[0]), numpy.array(V), numpy.array(a[1]), numpy.array(E))
	g.SetMarkerStyle(20)
	g.SetMarkerColor(1)
	g.SetLineColor(1)
	f = TF1("f", "pol"+o.Ofit, min(a[0]), max(a[0]))
	f.SetLineWidth(2)
	f.SetLineColor(kBlue)
	g.Fit(f)
	g.SetTitle(N+";#alpha;parameter value")
	C = TCanvas()
	C.cd()
	g.Draw("AP")
	f.Draw("same")
	C.Print("plots/"+N+".png")
	return f
fits = []
print Ps
for p in range(nPar):
	fits.append(MakeEPlot(Ps[p][0], Ps[p][1], "P"+str(p)))
	
BF = FUNC.BigFit(o.Pfit, o.Ofit)
print BF
F = TF2("hardFit", BF, 0., 4000., 0.15, 0.35)
p = 0
for P in range(nPar):
	for O in range(int(o.Ofit)):
		F.SetParameter(p, fits[P].GetParameter(O))
		p+=1

G = TGraph2DErrors(len(x[0]), numpy.array(x[0]), numpy.array(y[0]), numpy.array(z[0]), numpy.array(x[1]), numpy.array(y[1]), numpy.array(z[1]))
G.Fit(F, "EM0R")
ndof = F.GetNDF()
chi2 = F.GetChisquare()
print "Chi2Ndof = "  + str(chi2)+"/"+str(ndof) + " = " + str(chi2/ndof)


for n in range(len(z[0])):
		b = BKG.FindBin(x[0][n],y[0][n])
		est = F.Eval(x[0][n],y[0][n])
		FIT.SetBinContent(b, est)
		
oF = TFile("dists/PostFitFun"+o.Pfit+o.Ofit+".root","recreate")
oF.cd()
FIT.Write()
BKG.Write()
SLICE.Write()
oF.Write()
oF.Save()
oF.Close()
