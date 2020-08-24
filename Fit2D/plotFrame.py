import ROOT
import os
import numpy
import math
import sys

funcForms = {}
funcForms["P2"] = "TMath::Power(1-(x/13000.0),[0])/TMath::Power(x/13000.0,[1])"
funcForms["P3"] = "TMath::Power(1.-(x/13000.0),[0])/(TMath::Power(x/13000.0,[1]+([2]*TMath::Log(x/13000.0))))"
funcForms["P4"] = "TMath::Power(1.-(x/13000.0),[0])/(TMath::Power(x/13000.0,[1]+([2]*TMath::Log(x/13000.0))+([3]*TMath::Power(TMath::Log(x/13000.0),2))))"
funcForms["CDF,1"] = "(1/TMath::Power(x, [0]))*TMath::Power(1-(x/13000.0),[1])"
funcForms["CDF,2"] = "(1/TMath::Power(x, [0]))*TMath::Power(1-(x/13000.0)+TMath::Power(x/13000.0,2),[1])"
funcForms["CDF,3"] = "(1/TMath::Power(x, [0]))*TMath::Power(1-(x/13000.0)+TMath::Power(x/13000.0,2)-TMath::Power(x/13000.0,3),[1])"
funcForms["expoPoli,1"] = "TMath::Exp([0]+([1]*TMath::Log(x)))"
funcForms["expoPoli,2"] = "TMath::Exp([0]+([1]*TMath::Log(x))+([2]*TMath::Power(TMath::Log(x),2)))"
funcForms["expoPoli,3"] = "TMath::Exp([0]+([1]*TMath::Log(x))+([2]*TMath::Power(TMath::Log(x),2))+([3]*TMath::Power(TMath::Log(x),3)))"

sigDir = {}
sigDir["Diquark_chi500suu2000"] = "/home/th544/CMSSW_10_6_2/src/tests/AvJJ_alpha_test/Diquark_chi500suu2000_rebin.root"
sigDir["Diquark_chi500suu3000"] = "/home/th544/CMSSW_10_6_2/src/tests/AvJJ_alpha_test/Diquark_chi500suu3000_rebin.root"
sigDir["Diquark_chi1000suu4000"] = "/home/th544/CMSSW_10_6_2/src/tests/AvJJ_alpha_test/Diquark_chi1000suu4000_rebin.root"
sigDir["Diquark_chi1000suu5000"] = "/home/th544/CMSSW_10_6_2/src/tests/AvJJ_alpha_test/Diquark_chi1000suu5000_rebin.root"
sigDir["Diquark_chi1800suu8000"] = "/home/th544/CMSSW_10_6_2/src/tests/AvJJ_alpha_test/Diquark_chi1800suu8000_rebin.root"
sigDir["Diquark_chi2000suu8400"] = "/home/th544/CMSSW_10_6_2/src/tests/AvJJ_alpha_test/Diquark_chi2000suu8400_rebin.root"
sigDir["Diquark_chi2100suu9000"] = "/home/th544/CMSSW_10_6_2/src/tests/AvJJ_alpha_test/Diquark_chi2100suu9000_rebin.root"
sigDir["Diquark_chi3000suu8000"] = "/home/th544/CMSSW_10_6_2/src/tests/AvJJ_alpha_test/Diquark_chi3000suu8000_rebin.root"
sigDir["Stop_500"] = "/home/th544/CMSSW_10_6_2/src/tests/AvJJ_alpha_test/Stop_500_rebin.root"
sigDir["Stop_750"] = "/home/th544/CMSSW_10_6_2/src/tests/AvJJ_alpha_test/Stop_750_rebin.root"
sigDir["Stop_1000"] = "/home/th544/CMSSW_10_6_2/src/tests/AvJJ_alpha_test/Stop_1000_rebin.root"
sigDir["Stop_2000"] = "/home/th544/CMSSW_10_6_2/src/tests/AvJJ_alpha_test/Stop_2000_rebin.root"

lumifb = {}
lumifb[2016] = 35.9
lumifb[2017] = 49.9
lumifb[2018] = 59.7

xs = {}
xs["Diquark_chi500suu2000"] = 1.8374e-5
xs["Diquark_chi500suu3000"] = 3.4599e-6
xs["Diquark_chi1000suu4000"] = 6.7451e-7
xs["Diquark_chi1000suu5000"] = 1.39369e-7
xs["Diquark_chi1800suu8000"] = 0
xs["Diquark_chi2000suu8400"] = 2.26626e-10
xs["Diquark_chi2100suu9000"] = 6.31791e-11
xs["Diquark_chi3000suu8000"] = 4.08e-10
xs["Stop_500"] = 0
xs["Stop_750"] = 0
xs["Stop_1000"] = 0
xs["Stop_2000"] = 0


colors = [ROOT.kMagenta, ROOT.kRed, ROOT.kOrange-3, ROOT.kYellow-6, ROOT.kTeal-7, ROOT.kCyan+4, ROOT.kAzure-2, ROOT.kViolet+5]

# define run once decorator
def run_once(f):
    def wrapper(*args, **kwargs):
        if not wrapper.has_run:
            wrapper.has_run = True
            return f(*args, **kwargs)
    wrapper.has_run = False
    return wrapper

def AddCMSLumi(pad, fb, extra):
    cmsText     = "CMS" + (" " + str(extra) if extra is not None else "")
    cmsTextFont   = 61  
    lumiTextSize     = 0.45
    lumiTextOffset   = 0.15
    cmsTextSize      = 0.5
    cmsTextOffset    = 0.15
    H = pad.GetWh()
    W = pad.GetWw()
    l = pad.GetLeftMargin()
    t = pad.GetTopMargin()
    r = pad.GetRightMargin()
    b = pad.GetBottomMargin()
    e = 0.025
    pad.cd()
    lumiText = (str(fb)+" fb^{-1}" if fb is not None else "") + "(13 TeV)"
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)	
    extraTextSize = 0.76*cmsTextSize
    latex.SetTextFont(42)
    latex.SetTextAlign(31) 
    latex.SetTextSize(lumiTextSize*t)	
    latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText)
    pad.cd()
    latex.SetTextFont(cmsTextFont)
    latex.SetTextSize(cmsTextSize*t)
    latex.SetTextAlign(11)
    latex.DrawLatex(l, 1-t+cmsTextOffset*t, cmsText)
    pad.Update()

def importTH1(F, H, rebinX=0, png=False, name=None, logy=False):
    inFile = ROOT.TFile(F)
    inHist = inFile.Get(H)
    if name is not None: inHist.SetName(name)
    if rebinX != 0: QCD.RebinX(rebinX)
    C = ROOT.TCanvas()
    C.SetCanvasSize(600, 600)
    if logy: C.SetLogy()
    C.cd()
    inHist.Draw()
    if png: C.Print(H+".png")
    inHist.SetDirectory(0)
    return inHist

def importTH2(F, H, rebinX=0, rebinY=0, png=False, logz=False):
    inFile = ROOT.TFile(F)
    inHist = inFile.Get(H)
    if rebinX != 0: QCD.RebinX(rebinX)
    if rebinY != 0: QCD.RebinY(rebinY)
    C = ROOT.TCanvas()
    C.SetCanvasSize(600, 600)
    if logz: C.SetLogz()
    C.cd()
    inHist.Draw("COLZ")
    if png: C.Print(H+".png")
    inHist.SetDirectory(0)
    return inHist

def MakeAToy(Hi):
    R = numpy.random
    Ho = Hi.Clone("H"+str(R.randint(1)))
    Hw = Hi.Clone("Hw"+str(R.randint(1)))
    for j in range(1,Ho.GetNbinsY()):
        for i in range(1,Ho.GetNbinsX()):
            v = Ho.GetBinContent(i,j)
            e = Ho.GetBinError(i,j)
            n = max(0,numpy.rint(v + (e*R.normal())))
            Hw.SetBinContent(i,j,int(n))
            if Ho.GetBinContent(i,j) != 0: Hw.SetBinContent(i,j,Hw.GetBinContent(i,j)/Ho.GetBinContent(i,j))
            else: Hw.SetBinContent(i,j,0)
            Hw.SetBinError(i,j,0)
    Ho.Multiply(Hw)
    return Ho

def MakeNBinsFromMinToMax(N,Min,Max):
    BINS = []
    for i in range(N+1):
        BINS.append(Min+(i*(Max-Min)/N))
    return BINS
  
def MakeParStrings(npar, polN):
    Ps = []
    for n in range(npar):
        P = "([" + str(n*(polN + 1)) + "]"
        for m in range(polN):
            P = P + " + [" + str(n*(polN + 1) + m + 1) + "]"
            for o in range(m + 1): P = P + "*y"
        P = P + ")"
        Ps.append(P)
    return Ps

def MakeTG1(x, y, name, title=None, xTitle=None, yTitle=None, png=False):
    G = ROOT.TGraph(len(x), numpy.array(x), numpy.array(y))
    G.SetTitle(title if title is not None else name)
    if xTitle is not None: G.GetXaxis().SetTitle(xTitle)
    if yTitle is not None: G.GetYaxis().SetTitle(yTitle)
    if png:
        C = ROOT.TCanvas()
        C.SetCanvasSize(600, 600)
        C.cd()
        G.Draw("AL")
        C.Print(name +".png")
    return G

def MakeTG1E(x, y, ex, ey, name, title=None, xTitle=None, yTitle=None, png=False):
    GE = ROOT.TGraphErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(ex), numpy.array(ey))
    GE.SetTitle(title if title is not None else name)
    if xTitle is not None: GE.GetXaxis().SetTitle(xTitle)
    if yTitle is not None: GE.GetYaxis().SetTitle(yTitle)
    if png: 
        C = ROOT.TCanvas()
        C.SetCanvasSize(600, 600)
        C.cd()
        GE.Draw("AL")
        C.Print(name + ".png")
    return GE

def UnbiasPull(H):
    NRGBs = 8
    stops = [ 0.00, 0.15, 0.30, 0.45, 0.55, 0.70, 0.85, 1.00 ]
    #red = [ 1.00, 1.00, 1.00, 1.00, 1.00, 0.60, 0.20, 0.00 ]
    #green = [ 0.00, 0.20, 0.60, 1.00, 1.00, 0.60, 0.20, 0.00 ]
    NCont = 80

    red = [ 0.20, 0.40, 0.60, 1.00, 1.00, 0.60, 0.40, 0.20 ]
    green = [ 1.00, 1.00, 1.00, 1.00, 1.00, 0.60, 0.20, 0.10 ]
    blue = [ 0.10, 0.20, 0.60, 1.00, 1.00, 1.00, 1.00, 1.00 ]

    # red = [ 1.00, 1.00, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00 ]
    # green = [ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 ]
    # blue = [ 0.00, 0.00, 0.00, 0.00, 1.00, 1.00, 1.00, 1.00 ]

    ROOT.TColor.CreateGradientColorTable(NRGBs, numpy.array(stops), numpy.array(red), numpy.array(green), numpy.array(blue), NCont)
    H.UseCurrentStyle()
    return H

def UseWS(D):   
    if not os.path.exists(D): os.makedirs(D)
    os.chdir(D)

def TH2toTG2E(H):
    # H = ROOT.TH1
    x = []
    y = []
    z = []
    ex = []
    ey = []
    ez = []
    for j in range(1, H.GetNbinsY()):
        for i in range(1, H.GetNbinsX()):
            x.append(H.GetXaxis().GetBinCenter(i))
            y.append(H.GetYaxis().GetBinCenter(j))
            z.append(H.GetBinContent(i,j))
            ex.append(H.GetXaxis().GetBinWidth(i)/2.)
            ey.append(H.GetYaxis().GetBinWidth(j)/2.)
            ez.append(H.GetBinError(i,j))
    G = ROOT.TGraph2DErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(z), numpy.array(ex), numpy.array(ey), numpy.array(ez))
    G.SetTitle("TG2E")
    G.SetDirectory(0)
    return G

def styleH(H, E, name, xTitle=None, yTitle=None):
    # H = ROOT.TH1()
    # E = ROOT.TH1() or ROOT.TF1()
    # name = string

    C = ROOT.TCanvas()
    C.SetCanvasSize(600, 600)
    pad = ROOT.TPad(name, name, 0, 0, 1, 1)
    H.SetLineColor(ROOT.kBlack)
    if xTitle is not None: H.GetXaxis().SetTitle(xTitle)
    if yTitle is not None: H.GetYaxis().SetTitle(yTitle)
    H.SetLineWidth(1)
    E.SetLineColor(ROOT.kRed)
    E.SetLineWidth(1)
    pad.cd()
    if type(H) == ROOT.TH1F():
        H.SetMarkerStyle(20)
        H.SetStats(0)
        H.GetXaxis().SetRange(0, H.FindLastBinAbove(0)+5)
    H.Draw()
    E.Draw("same")

    C.cd()
    pad.Draw()
    C.Print(name + ".png")
    C.Close()
    ROOT.gSystem.ProcessEvents()

def stylePull(Hin, Ein, name, signals=None, title=None, xTitle=None, yTitle=None, CH=None, pullRange=None, displayCI=False, legTitle=None, lumi=None, CMSExtra=None):
    # H = ROOT.TH1()
    # E = ROOT.TH1() or ROOT.TF1()
    # name = string
    xbins = []
    MinNonZeroBins = []
    MaxNonZeroBins = []
    for n in range(Hin.GetXaxis().GetNbins()): xbins.append(Hin.GetXaxis().GetBinUpEdge(n))

    H = Hin.Clone(Hin.GetName()+"_clone")
    E = Ein.Clone(Ein.GetName()+"_clone")

    SignalErr = []
    if signals is not None:
        for S in signals: SignalErr.append(S.Clone())

    C = ROOT.TCanvas()
    C.SetCanvasSize(600, 600)
    pad = ROOT.TPad(name, name, 0, 0.25, 1, 1)
    pad.SetTopMargin(0.1)
    pad.SetBottomMargin(0.005)
    pad.SetLeftMargin(0.13)
    pad.SetRightMargin(0.15)
    pad.SetLogy()

    legendx = 0.80
    if displayCI: legendx = legendx - 0.03
    if signals is not None:
        for n in signals: legendx = legendx - 0.04
    legend = ROOT.TLegend(0.55, legendx, 0.84, 0.88)
    if legTitle is not None: legend.SetHeader(legTitle)
    legend.SetFillColor(0)
    legend.SetLineColor(0)

    for n in range(1, H.GetXaxis().GetNbins()+1):
        try:
            if H.GetBinError(n) == 0 and H.GetXaxis().GetBinUpEdge(n) <= E.GetHistogram().GetBinLowEdge(E.GetHistogram().FindFirstBinAbove(0)): H.SetBinError(n, 0)
        except:
            if H.GetXaxis().GetBinUpEdge(n) <= E.GetBinLowEdge(E.FindFirstBinAbove(0)): H.SetBinError(n, 0)
    H.SetLineColor(ROOT.kBlack)
    if yTitle is not None: H.GetYaxis().SetTitle(yTitle)
    H.SetMarkerStyle(20)
    H.SetMarkerSize(0.5)
    H.SetLineWidth(1)
    H.SetStats(0)
    if title is not None: H.SetTitle(title)
    E.SetLineWidth(1)
    pad.cd()
    E.SetLineColor(ROOT.kGreen+3)
    H.Draw("")
    E.SetFillColor(ROOT.kGreen-9)
    E.SetFillStyle(1001)
    E.Draw("same hist")
    try: MinNonZeroBins.append(H.FindFirstBinAbove(0)-5)
    except: MinNonZeroBins.append(H.FindFirstBinAbove(0))
    MaxNonZeroBins.append(H.FindLastBinAbove(0)+5)

    bins = []
    for n in range(H.GetNbinsX()): bins.append(H.GetXaxis().GetBinUpEdge(n))
    CI = CH
    CIpull = ROOT.TH1F()
    if displayCI:
        if CH is None:
            CI = ROOT.TH1F("CI_"+name, "CI", len(bins)-1, numpy.array(bins))
            CIpull = ROOT.TH1F("CIpull_"+name, "CIpull", len(bins)-1, numpy.array(bins))
            ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(CI)
            for n in range(1, CI.GetNbinsX()):
                err = 1.4
                if H.GetBinContent(n) != 0: err = H.GetBinError(n)
                CIpull.SetBinContent(n, (CI.GetBinContent(n)-E.Eval(H.GetBinCenter(n)))/err)
                CIpull.SetBinError(n, CI.GetBinError(n)/err)
                try:
                    if CI.GetXaxis().GetBinUpEdge(n) <= E.GetHistogram().GetBinLowEdge(E.GetHistogram().FindFirstBinAbove(0)):
                        CI.SetBinContent(n, 0)
                        CIpull.SetBinContent(n, 0)
                        CIpull.SetBinError(n, 0)
                except:
                    if CI.GetXaxis().GetBinUpEdge(n) <= E.GetBinLowEdge(E.FindFirstBinAbove(0)):
                        CI.SetBinContent(n, 0)
                        CIpull.SetBinContent(n, 0)
                        CIpull.SetBinError(n, 0)

        CI.SetFillStyle(3344)
        CI.SetFillColorAlpha(ROOT.kGreen-2, 0.8)
        CI.SetLineColor(ROOT.kGreen-2)
        CI.Draw("e2 e0 same")
        CIpull.SetFillStyle(3344)
        CIpull.SetFillColorAlpha(ROOT.kGreen-2, 0.8)
        CIpull.SetLineColor(ROOT.kGreen-2)

    legend.AddEntry(H, "QCD MC")
    legend.AddEntry(E, "Background Estimate", "f")
    if displayCI: legend.AddEntry(CI, "Uncertainty", "f")

    if signals is not None:
        for n in range(len(signals)):
            for x in range(1, signals[n].GetXaxis().GetNbins()):
                SignalErr[n].SetBinContent(x, signals[n].GetBinContent(x))
                SignalErr[n].SetBinError(x, signals[n].GetBinError(x))
                signals[n].SetBinError(x, 0)
            signals[n].SetMarkerColor(colors[n])
            signals[n].SetFillColor(ROOT.kWhite)
            signals[n].SetFillStyle(1001)
            signals[n].SetLineWidth(1)
            signals[n].SetLineColor(colors[n])
            SignalErr[n].SetMarkerColor(colors[n])
            SignalErr[n].SetFillColorAlpha(colors[n], 0.6)
            SignalErr[n].SetFillStyle(3344)
            SignalErr[n].SetLineColor(colors[n])
            # try:
            #     if signals[n].GetXaxis().GetBinUpEdge(n) <= E.GetHistogram().GetBinLowEdge(E.GetHistogram().FindFirstBinAbove(0)):
            #         SignalErr[n].SetBinContent(n, 0)
            #         signals[n].SetBinContent(n, 0)
            #         SignalErr[n].SetBinError(n, 0)
            # except:
            #     if signals[n].GetXaxis().GetBinUpEdge(n) <= E.GetBinLowEdge(E.FindFirstBinAbove(0)):
            #         SignalErr[n].SetBinContent(n, 0)
            #         signals[n].SetBinContent(n, 0)
            #         SignalErr[n].SetBinError(n, 0)
            signals[n].Draw("same")
            # SignalErr[n].Draw("e2 same")
            legend.AddEntry(signals[n], signals[n].GetName(), "f")
            try: MinNonZeroBins.append(signals[n].FindFirstBinAbove(0)-5)
            except: MinNonZeroBins.append(signals[n].FindFirstBinAbove(0))
            MaxNonZeroBins.append(signals[n].FindLastBinAbove(0)+5)

    H.Draw("same")
    H.GetXaxis().SetRange(min(*MinNonZeroBins), max(*MaxNonZeroBins))
    legend.Draw("same")

    pull = H.Clone("p_"+name)
    pull.Add(E, -1)
    for n in range(1, pull.GetXaxis().GetNbins()):
        pull.SetBinError(n, 1.)
        try:
            if pull.GetXaxis().GetBinUpEdge(n) <= E.GetHistogram().GetBinLowEdge(E.GetHistogram().FindFirstBinAbove(0)):
                pull.SetBinContent(n, 0)
                pull.SetBinError(n, 0)
            else:
                if H.GetBinError(n) != 0: pull.SetBinContent(n, pull.GetBinContent(n)/H.GetBinError(n))
                else: pull.SetBinContent(n, pull.GetBinContent(n)/1.4)
        except:
            if pull.GetXaxis().GetBinUpEdge(n) <= E.GetBinLowEdge(E.FindFirstBinAbove(0)):
                pull.SetBinContent(n, 0)
                pull.SetBinError(n, 0)
            else:
                if H.GetBinError(n) != 0: pull.SetBinContent(n, pull.GetBinContent(n)/H.GetBinError(n))
                else: pull.SetBinContent(n, pull.GetBinContent(n)/1.4)
        if H.GetBinContent(n) == 0:
            pull.SetBinContent(n, 0)
            pull.SetBinError(n, 0)

    pullpad = ROOT.TPad("pp_"+name, '', 0, 0, 1, 0.25)
    pullpad.SetLeftMargin(0.13)
    pullpad.SetRightMargin(0.15)
    pullpad.SetTopMargin(0.04)
    pullpad.SetBottomMargin(0.3)
    pullpad.SetGrid()
    pull.SetBit(ROOT.TH1.kNoTitle)
    pull.SetStats(0)
    pull.SetMarkerStyle(20)
    pull.SetMarkerColor(ROOT.kBlack)
    pull.SetMarkerSize(0.5)
    pull.SetLineWidth(1)
    pull.GetXaxis().SetLabelSize(0.1)
    if xTitle is not None: pull.GetXaxis().SetTitle(xTitle)
    pull.GetXaxis().SetTitleSize(0.1)
    pull.GetYaxis().SetNdivisions(505)
    pull.GetYaxis().SetLabelSize(0.1)
    pull.GetYaxis().SetTitle("#frac{Bkg - Est}{Uncertainty}")
    pull.GetYaxis().SetTitleSize(0.08)
    pull.GetYaxis().SetTitleOffset(0.5)
    pull.GetYaxis().CenterTitle()
    pull.SetMarkerStyle(20)
    indicator = ROOT.TF1('ind_'+name, '0')
    indicator.SetLineColor(ROOT.kBlack)
    indicator.SetLineWidth(1)

    pullpad.cd()
    pull.GetXaxis().SetRange(min(*MinNonZeroBins), max(*MaxNonZeroBins))
    pull.Draw("")
    if displayCI:
        pass
        # CIpull.Draw("same e2")
        # pull.Draw("same")
    maxpull = max(abs(pull.GetMaximum()),abs(pull.GetMinimum()))
    if pullRange is not None: pull.GetYaxis().SetRangeUser(pullRange[0], pullRange[1])
    else: pull.GetYaxis().SetRangeUser(-maxpull - 0.2, maxpull + 0.2) # asymmetric pulls
    indicator.Draw("same")

    AddCMSLumi(pad, lumi, CMSExtra)

    C.cd()
    pad.Draw()
    pullpad.Draw()
    C.Print(name + ".png")
    C.Close()
    ROOT.gSystem.ProcessEvents()

def styleRatio(H, E, name):
    # H = ROOT.TH1()
    # E = ROOT.TH1() or ROOT.TF1()
    # name = string
    C = ROOT.TCanvas()
    C.SetCanvasSize(600, 600)
    pad = ROOT.TPad(name, name, 0, 0.25, 1, 1)
    pad.SetTopMargin(0.1)
    pad.SetBottomMargin(0.005)
    pad.SetLeftMargin(0.13)
    pad.SetRightMargin(0.15)
    pad.SetLogy()
    H.SetLineColor(ROOT.kBlack)
    H.GetYaxis().SetTitle("Number of events")
    H.SetMarkerStyle(20)
    H.SetLineWidth(1)
    H.SetStats(0)
    E.SetLineColor(ROOT.kRed)
    E.SetLineWidth(1)
    # try: E.SetStats(0)
    pad.cd()
    H.GetXaxis().SetRange(0, H.FindLastBinAbove(0)+5)
    # H.GetXaxis().SetRangeUser(0, 1600)
    H.Draw()
    E.Draw("same H")

    ratio = H.Clone("r_"+name)
    ratio.Divide(E)
    for n in range(1, ratio.GetXaxis().GetNbins()):
        ratio.SetBinError(n, 0)
        try:
            if ratio.GetXaxis().GetBinUpEdge(n) <= E.GetHistogram().GetBinLowEdge(E.GetHistogram().FindFirstBinAbove(0)): ratio.SetBinContent(n, 1)
        except:
            if ratio.GetXaxis().GetBinUpEdge(n) <= E.GetBinLowEdge(E.FindFirstBinAbove(0)): ratio.SetBinContent(n, 1)

    ratiopad = ROOT.TPad("rp_"+name, '', 0, 0, 1, 0.25)
    ratiopad.SetLeftMargin(0.13)
    ratiopad.SetRightMargin(0.15)
    ratiopad.SetTopMargin(0.04)
    ratiopad.SetBottomMargin(0.3)
    ratiopad.SetGrid()
    ratio.SetBit(ROOT.TH1.kNoTitle)
    ratio.SetStats(0)
    ratio.SetLineColor(ROOT.kRed)
    ratio.SetFillColor(ROOT.kWhite)
    ratio.GetXaxis().SetLabelSize(0.1)
    ratio.GetXaxis().SetTitle('#bar{M}_{jj} [GeV]')
    ratio.GetXaxis().SetTitleSize(0.1)
    ratio.GetYaxis().SetNdivisions(505)
    ratio.GetYaxis().SetLabelSize(0.1)
    ratio.GetYaxis().SetTitle("#frac{Bkg}{E}")
    ratio.GetYaxis().SetTitleSize(0.08)
    ratio.GetYaxis().SetTitleOffset(0.5)
    ratio.GetYaxis().CenterTitle()
    ratio.SetMarkerStyle(20)
    indicator = ROOT.TF1('ind_'+name, '1')
    indicator.SetLineColor(ROOT.kBlack)
    indicator.SetLineWidth(1)

    ratiopad.cd()
    ratio.GetXaxis().SetRange(0, H.FindLastBinAbove(0)+5)
    ratio.Draw()
    maxpull = max(abs(ratio.GetMaximum()),abs(ratio.GetMinimum()))
    ratio.GetYaxis().SetRangeUser(0, maxpull + 0.2) # asymmetric pulls
    # ratio.GetYaxis().SetRangeUser(0, 5)
    indicator.Draw("same")
    C.cd()
    pad.Draw()
    ratiopad.SetLogy()
    ratiopad.Draw()
    C.Print(name + ".png")
    C.Close()
    ROOT.gSystem.ProcessEvents()

def styleRatioPull(H, E, name):
    # H = ROOT.TH1()
    # E = ROOT.TH1() or ROOT.TF1()  
    # name = string
    return 0

class Fit1D:
    def __init__(self, H, E, options="EM0R"):
        self.H = H
        self.E = E
        self.Pv = []
        self.Pe = []

        self.H.Fit(E, options)
        for n in range(E.GetNpar()):
            self.Pv.append(self.E.GetParameter(n))
            self.Pe.append(self.E.GetParError(n))
 
    def GetConfidenceIntervals(self, H):
        ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(H)
        for n in range(1, self.H.GetXaxis().GetNbins()):
            try:
                if H.GetXaxis().GetBinUpEdge(n) <= self.E.GetHistogram().GetBinLowEdge(self.E.GetHistogram().FindFirstBinAbove(0)): H.SetBinContent(n, 0)
            except:
                if H.GetXaxis().GetBinUpEdge(n) <= self.E.GetBinLowEdge(self.E.FindFirstBinAbove(0)): H.SetBinContent(n, 0)
        return 0

    def GetNorm(self):
        return self.H.Integral()

    def GetParameter(self, n):
        return self.Pv[n]

    def GetParameters(self):
        return self.Pv
    
    def GetParError(self, n):
        return self.Pe[n]

    def GetParErrors(self):
        return self.Pe

    def makeH(self, name, xTitle=None, yTitle=None):
        styleH(self.H, self.E, name, xTitle=xTitle, yTitle=yTitle)
        return 0

    def makePull(self, name, title=None, signals=None, xTitle=None, yTitle=None, displayCI=False, legTitle=None, lumi=None, CMSExtra=None):
        stylePull(self.H, self.E, name, title=title, signals=signals, xTitle=xTitle, yTitle=yTitle, displayCI=displayCI, legTitle=legTitle, lumi=lumi, CMSExtra=CMSExtra)
        return 0

class Fit2D:
    def __init__(self, H, SF=None, lumi=None, CMSExtra=None):
        self.H = H.Clone()
        self.HbinsX = []
        self.HbinsY = []
        for n in range(self.H.GetNbinsX()+1): self.HbinsX.append(self.H.GetXaxis().GetBinUpEdge(n))
        for n in range(self.H.GetNbinsY()+1): self.HbinsY.append(self.H.GetYaxis().GetBinUpEdge(n))
        self.signals=SF

        self.xTitle = ""
        self.yTitle = ""
        self.zTitle = ""
        self.lumi = None
        self.CMSExtra = None

        self.loXcut = 0
        self.hiXcut = len(self.HbinsX)
        self.loYcut = 0
        self.hiYcut = len(self.HbinsY)
        ROOT.gStyle.SetNumberContours(80)

    def GetBkg(self):
        out = self.GHB.Clone()
        out.SetDirectory(0)
        return out

    def GetEst(self):
        out = self.GHE.Clone()
        out.SetDirectory(0)
        return out

    def GetPull(self):
        out = self.GHP.Clone()
        out.SetDirectory(0)
        return out

    def MakeMJJM4J(self, name, title=None, options="COLZ", logz=True, zrange=None, zrangeP=None, png=True):
        m4j_binning = MakeNBinsFromMinToMax(60, 0., 10000.)
        self.TB = ROOT.TH2F("TB",title if title is not None else "Global background in #bar M_{jj} #times M_{4j} space;#bar M_{jj} [GeV];M_{4j} [GeV]", len(self.HbinsX)-1, numpy.array(self.HbinsX), len(m4j_binning)-1, numpy.array(m4j_binning))
        self.TE = ROOT.TH2F("TE",title if title is not None else "Global estimate in #bar M_{jj} #times M_{4j} space;#bar M_{jj} [GeV];M_{4j} [GeV]", len(self.HbinsX)-1, numpy.array(self.HbinsX), len(m4j_binning)-1, numpy.array(m4j_binning))
        self.TP = ROOT.TH2F("TP",title if title is not None else "Global pull in #bar M_{jj} #times M_{4j} space;#bar M_{jj} [GeV];M_{4j} [GeV]", len(self.HbinsX)-1, numpy.array(self.HbinsX), len(m4j_binning)-1, numpy.array(m4j_binning))
        self.Tsignals = []
        if self.signals is not None:
            for s in range(len(self.signals)): self.Tsignals.append(ROOT.TH2F("TS_" + str(s), self.signals[s].GetName() + " in #bar M_{jj} #times M_{4j} space; #bar{M}_{jj} [GeV]; M_{4j} [GeV]", len(self.HbinsX)-1, numpy.array(self.HbinsX), len(m4j_binning)-1, numpy.array(m4j_binning)))

        for i in range(1,self.H.GetNbinsX()):
            for j in range(1,len(m4j_binning)):
                M2J = self.TB.GetXaxis().GetBinCenter(i)
                M4J = self.TB.GetYaxis().GetBinCenter(j)
                a = M2J/M4J
                abin = self.H.GetYaxis().FindBin(a)
                
                TBval = self.GHB.GetBinContent(i, abin)
                TEval = self.GHE.GetBinContent(i, abin)
                TPval = self.GHP.GetBinContent(i, abin)

                self.TB.SetBinContent(i, j, TBval)
                self.TE.SetBinContent(i, j, TEval)
                self.TP.SetBinContent(i, j, TPval)
                if self.signals is not None:
                    for s in range(len(self.signals)):
                        TSval = self.signals[s].GetBinContent(i, abin)
                        self.Tsignals[s].SetBinContent(i, j, TSval)

        ROOT.gStyle.SetPalette(57)
        CB = ROOT.TCanvas()
        CB.SetCanvasSize(600,600)
        if logz: CB.SetLogz()
        CB.cd()
        if zrange is not None: self.TB.GetZaxis().SetRangeUser(zrange[0], zrange[1])
        self.TB.UseCurrentStyle()
        self.TB.SetStats(0)
        self.TB.GetXaxis().SetTitle(self.xTitle)
        self.TB.GetYaxis().SetTitle(self.yTitle)
        self.TB.Draw(options)
        if png: CB.Print(name+"_B.png")

        CE = ROOT.TCanvas()
        CE.SetCanvasSize(600,600)
        if logz: CE.SetLogz()
        CE.cd()
        if zrange is not None: self.TE.GetZaxis().SetRangeUser(zrange[0], zrange[1])
        self.TE.UseCurrentStyle()
        self.TE.SetStats(0)
        self.TE.GetXaxis().SetTitle(self.xTitle)
        self.TE.GetYaxis().SetTitle(self.yTitle)
        self.TE.Draw(options)
        if png: CE.Print(name+"_E.png")

        CP = ROOT.TCanvas()
        CP.SetCanvasSize(600,600)
        CP.cd()
        if zrangeP is not None: self.TP.GetZaxis().SetRangeUser(zrangeP[0], zrangeP[1])
        else: self.TP.GetZaxis().SetRangeUser(-5., 5.)
        UnbiasPull(self.TP)
        self.TP.SetStats(0)
        self.TP.GetXaxis().SetTitle(self.xTitle)
        self.TP.GetYaxis().SetTitle(self.yTitle)
        self.TP.Draw(options)
        if png: CP.Print(name+"_P.png")

    def Print2D(self, name, options="COLZ", logz=True, zrange=None, zrangeP=None):
        CB = ROOT.TCanvas()
        if logz: CB.SetLogz()
        CB.cd()
        if zrange is not None: self.GHB.GetZaxis().SetRangeUser(zrange[0], zrange[1])
        self.GHB.UseCurrentStyle()
        self.GHB.SetStats(0)
        self.GHB.GetXaxis().SetTitle(self.xTitle)
        self.GHB.GetYaxis().SetTitle(self.yTitle)
        self.GHB.Draw(options)
        CB.Print(name+"_B.png")

        CE = ROOT.TCanvas()
        if logz: CE.SetLogz()
        CE.cd()
        if zrange is not None: self.GHE.GetZaxis().SetRangeUser(zrange[0], zrange[1])
        self.GHE.UseCurrentStyle()
        self.GHE.SetStats(0)
        self.GHE.GetXaxis().SetTitle(self.xTitle)
        self.GHE.GetYaxis().SetTitle(self.yTitle)
        self.GHE.Draw(options)
        CE.Print(name+"_E.png")

        CP = ROOT.TCanvas()
        CP.cd()
        if zrangeP is not None: self.GHP.GetZaxis().SetRangeUser(zrangeP[0], zrangeP[1])
        # else: self.GHP.GetZaxis().SetRangeUser(-5., 5.)
        UnbiasPull(self.GHP)
        self.GHP.SetStats(0)
        self.GHP.GetXaxis().SetTitle(self.xTitle)
        self.GHP.GetYaxis().SetTitle(self.yTitle)
        self.GHP.Draw(options)
        CP.Print(name+"_P.png")
        ROOT.gStyle.SetPalette(57)
        
    def PrintChi2Ndof(self):
        print("|===> Chi2/Ndof = " + str(self.chi2) + "/" + str(self.ndof))

    def ProjectionX(self, name, title=None, legTitle=None, start=0, end=1, logy=True, pullRange=None):
        axisLength = self.hiYcut - self.loYcut
        startBin = self.loYcut + int(math.ceil(axisLength*start))
        endBin = self.loYcut + int(math.floor(axisLength*end))
        
        projB = self.GHB.ProjectionX(name+"_B", startBin, endBin, "e")
        projE = self.GHE.ProjectionX(name+"_E", startBin, endBin, "e")
        projSignals = []
        for S in range(len(self.signals)):
            projS = self.signals[S].ProjectionX(self.signals[S].GetName()+"_projX", startBin, endBin, "e")
            projSignals.append(projS)

        stylePull(projB, projE, name, signals=projSignals, legTitle=legTitle, title=title, xTitle=self.xTitle, yTitle=self.yTitle, pullRange=pullRange, lumi=self.lumi, CMSExtra=self.CMSExtra)

    def ProjectionY(self, name, title=None, legTitle=None, start=0, end=1, logy=True, pullRange=None):
        axisLength = self.hiXcut - self.loXcut
        startBin = self.loXcut + int(math.ceil(axisLength*start))
        endBin = self.loXcut + int(math.floor(axisLength*end))

        projB = self.GHB.ProjectionY(name+"_B", startBin, endBin, "e")
        projE = self.GHE.ProjectionY(name+"_E", startBin, endBin, "e")
        projSignals = []
        for S in range(len(self.signals)):
            projS = self.signals[S].ProjectionY(self.signals[S].GetName()+"_projX", startBin, endBin, "e")
            projSignals.append(projS)

        stylePull(projB, projE, name, signals=projSignals, legTitle=legTitle, title=title, xTitle=self.xTitle, yTitle=self.yTitle, pullRange=pullRange, lumi=self.lumi, CMSExtra=self.CMSExtra)

    def ProjectionM4J(self, name, title=None, legTitle=None, start=0, end=-1, logy=True, pullRange=None):
        axisLength = self.TB.GetNbinsY() - 1
        startBin = int(math.ceil(axisLength*start))
        endBin = int(math.floor(axisLength*end))

        projB = self.TB.ProjectionY(name+"_B", startBin, endBin, "e")
        projE = self.TE.ProjectionY(name+"_E", startBin, endBin, "e")
        projSignals = []
        for S in range(len(self.signals)):
            projS = self.Tsignals[S].ProjectionY(self.signals[S].GetName()+"_projX", startBin, endBin, "e")
            projSignals.append(projS)

        stylePull(projB, projE, name, signals=projSignals, legTitle=legTitle, title=title, xTitle=self.xTitle, yTitle=self.yTitle, pullRange=pullRange, lumi=self.lumi, CMSExtra=self.CMSExtra)

    def SetCMSText(self, lumi, CMSExtra):
        self.lumi = lumi
        self.CMSExtra = CMSExtra

    @run_once
    def SetLoYcut(self, cut):
        self.loYcut = cut

    def SetHiYcut(self, cut):
        self.hiYcut = cut

    def SetTitleX(self, title):
        self.xTitle = title

    def SetTitleY(self, title):
        self.yTitle = title

    def SetTitleZ(self, title):
        self.zTitle = title

    def SliceYFit2D(self, fname="P2", polN=3, options="EMR0", startbin=0, startcut=400):
        self.f = funcForms[fname]
        self.numpar = ROOT.TF1("t", self.f, 1, 2).GetNpar()
        self.polN = polN
        self.startcut = startcut
        self.slices = []

        self.pulls = [[] for i in range(self.H.GetNbinsX())]
        self.norms = []
        self.starts = []
        self.chi2Ndof = []
        self.yalph = []
        self.yalphe = []
        self.CIs = []
        self.ps = [[[],[]] for i in range(self.numpar)]

        for j in range(1,self.H.GetNbinsY()+1):
            print "|===> Slice "+ str(j)
            Hx = self.H.ProjectionX("q_"+str(j),j,j)
            # try:
            #     for n in self.signals: 
            # except: pass
            self.slices.append(Hx.Clone())
            self.SetHiYcut(j)
            if not (Hx.Integral() > 0): continue
            if self.H.GetYaxis().GetBinCenter(j) > 0.35 or self.H.GetYaxis().GetBinCenter(j) < 0.15: continue
            self.SetLoYcut(j)
            m = Hx.GetMaximumBin() + int(startbin)
            start = max(Hx.GetBinLowEdge(m), self.startcut)
            end = Hx.GetBinCenter(self.H.GetNbinsX())
            fit = ROOT.TF1("tF_sliceyfit2d_"+str(j), self.f, start, end)
            fFrame = Fit1D(Hx, fit)
            for i in range(1, self.H.GetNbinsX()+1):
                p = 0
                if i >= Hx.FindBin(start):
                    if Hx.GetBinContent(i) == 0: Hx.SetBinError(i, 1.4)
                    p = (Hx.GetBinContent(i)-fit.Eval(Hx.GetBinCenter(i)))/Hx.GetBinError(i)
                self.pulls[i-1].append(p)
            ndof = fit.GetNDF()
            chi2 = fit.GetChisquare()
            if (ndof < 5) or (chi2/ndof > 1.5): continue
            self.starts.append(start)
            self.norms.append(Hx.Integral())
            self.chi2Ndof.append(chi2/ndof)
            for n in range(self.numpar):
                self.ps[n][0].append(fit.GetParameter(n))
                self.ps[n][1].append(fit.GetParError(n))
            self.yalph.append(self.H.GetYaxis().GetBinCenter(j))
            self.yalphe.append(self.H.GetYaxis().GetBinWidth(j)/2.)

        oPs = []
        for n in range(self.numpar):
            GE = MakeTG1E(self.yalph, self.ps[n][0], self.yalphe, self.ps[n][1], "p"+str(n), "p"+str(n), png=False)
            f = ROOT.TF1("f", "pol"+str(self.polN), min(self.yalph), max(self.yalph))
            pFrame = Fit1D(GE, f, options="EMR0")
            pFrame.makeH("p"+str(n))
            oPs.append(f)
        
        self.G = TH2toTG2E(self.H)
        Ps = MakeParStrings(self.numpar, self.polN)
        self.f_aug = self.f
        for n in range(self.numpar): self.f_aug = self.f_aug.replace("[%d]" % (self.numpar-n-1), Ps[self.numpar-n-1]) # replacing real parameter with polN parametrization
        self.F2 = ROOT.TF2("hardFit", self.f_aug, self.HbinsX[0], self.HbinsX[-1], self.HbinsY[0], self.HbinsY[-1])
        for n in range(self.numpar):
            for m in range(self.polN + 1):
                self.F2.SetParameter(n*(self.polN + 1) + m, oPs[n].GetParameter(m))
        C = ROOT.TCanvas()
        C.cd()
        self.G.Draw("AP")
        self.G.Fit(self.F2, options)
        self.F2.Draw("same")
        self.ndof = self.F2.GetNDF()
        self.chi2 = self.F2.GetChisquare()
        self.GHB = ROOT.TH2F("GHB", "Global Background;%s;%s" % (self.xTitle, self.yTitle), len(self.HbinsX)-1, numpy.array(self.HbinsX), len(self.HbinsY)-1, numpy.array(self.HbinsY))
        self.GHE = ROOT.TH2F("GHE", "Global Estimate;%s;%s" % (self.xTitle, self.yTitle), len(self.HbinsX)-1, numpy.array(self.HbinsX), len(self.HbinsY)-1, numpy.array(self.HbinsY))
        self.GHP = self.H.Clone("GHP")
        self.GHP.SetTitle("Global Pull;%s;%s" % (self.xTitle, self.yTitle))
        self.GHP.Add(self.F2,- 1.)
        GHPmax = 0
        for j in range(1, self.H.GetNbinsY()+1):
            Hx = self.H.ProjectionX("q2_"+str(j),j,j)
            if (Hx.Integral() > 0) and (self.H.GetYaxis().GetBinCenter(j)) < 0.35 and (self.H.GetYaxis().GetBinCenter(j)) > 0.15:
                m = Hx.GetMaximumBin() + int(startbin)
                start = max(Hx.GetBinLowEdge(m), self.startcut)
                for i in range(m, self.H.GetNbinsX()+1):
                    if i >= Hx.FindBin(start):
                        q = self.H.GetBinContent(i,j)
                        qe = Hx.GetBinError(i)
                        if q == 0: qe = 1.4
                        est = self.F2.Eval(Hx.GetBinCenter(i), self.H.GetYaxis().GetBinCenter(j))
                        P = (q - est)/qe
                        self.GHP.SetBinContent(i, j, P)
                        self.GHE.SetBinContent(i, j, est)
                        self.GHB.SetBinContent(i, j, q)
                        GHPmax = max(abs(P), GHPmax)
                    else:
                        self.GHP.SetBinContent(i, j, 0)
                        self.GHE.SetBinContent(i, j, 0)
                        self.GHB.SetBinContent(i, j, 0)
                for i in range(1, m): self.GHP.SetBinContent(i, j, 0)
            else:
                for i in range(1, self.H.GetNbinsX()+1): self.GHP.SetBinContent(i, j, 0)
        self.GHP.GetZaxis().SetRangeUser(-GHPmax-0.19, GHPmax+0.2)
