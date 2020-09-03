#
import numpy 

def FindAndSetMax(*args):
	if len(args) == 1: args = args[0]
	maximum = 0.0
	for i in args:
		i.SetStats(0)
		t = i.GetMaximum()
		if t > maximum:
			maximum = t
	for j in args:
		j.GetYaxis().SetRangeUser(0.1,maximum*1.35)#should be 1.35 (below as well)
		j.SetLineWidth(2)
	return maximum*1.35

def MakeNBinsFromMinToMax(N,Min,Max):
    BINS = []
    for i in range(N+1):
        BINS.append(Min+(i*(Max-Min)/N))
    return numpy.array(BINS)
    
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
	
def WhichFit(F):
	if F == "P2": return ["TMath::Power(1-(x/13000.0),[0])/(TMath::Power(x/13000.0,[1]))", 2]
	if F == "P3": return ["TMath::Power(1-(x/13000.0),[0])/(TMath::Power(x/13000.0,[1]+([2]*TMath::Log(x/13000.0))))", 3]
	if F == "P4": return ["TMath::Power(1.-(x/13000.0),[0])/(TMath::Power(x/13000.0,[1]+([2]*TMath::Log(x/13000.0))+([3]*TMath::Power(TMath::Log(x/13000.0),2))))", 4]
	if F == "CDF_1": return ["(1/TMath::Power(x, [0]))*TMath::Power(1-(x/13000.0),[1])", 2]
	if F == "CDF_2": return ["(1/TMath::Power(x, [0]))*TMath::Power(1-(x/13000.0)+TMath::Power(x/13000.0,2),[1])", 2]
	if F == "CDF_3": return ["(1/TMath::Power(x, [0]))*TMath::Power(1-(x/13000.0)+TMath::Power(x/13000.0,2)-TMath::Power(x/13000.0,3),[1])", 2]
	if F == "expoPoli_1": return ["TMath::Exp([0]+([1]*TMath::Log(x)))", 2]
	if F == "expoPoli_2": return ["TMath::Exp([0]+([1]*TMath::Log(x))+([2]*TMath::Power(TMath::Log(x),2)))", 3]
	if F == "expoPoli_3": return ["TMath::Exp([0]+([1]*TMath::Log(x))+([2]*TMath::Power(TMath::Log(x),2))+([3]*TMath::Power(TMath::Log(x),3)))", 4]
	if F == "atlas_0": return["TMath::Power(1.0-TMath::Power(x/13000.0, 1/3.), [0])/TMath::Power(x/13000.0, [1])", 2]
	if F == "atlas_1": return["TMath::Power(1.0-TMath::Power(x/13000.0, 1/3.), [0])/TMath::Power(x/13000.0, [1] + [2]*TMath::Log(x))", 3]
	if F == "atlas_2": return["TMath::Power(1.0-TMath::Power(x/13000.0, 1/3.), [0])/TMath::Power(x/13000.0, [1] + [2]*TMath::Log(x) + [3]*TMath::Power(TMath::Log(x),2))", 4]
def BigFit(F, O):
	npar = WhichFit(F)[1]
	polN = int(O)
	Ps = []
	for n in range(npar):
		P = "([" + str(n*(polN + 1)) + "]"
		for m in range(polN):
			P = P + " + [" + str(n*(polN + 1) + m + 1) + "]"
			for o in range(m + 1): P = P + "*y"
		P = P + ")"
		Ps.append(P)
	FitFunc = WhichFit(F)[0]
	for n in range(npar): FitFunc = FitFunc.replace("[%d]" % (npar-n-1), Ps[npar-n-1])
	return FitFunc
	
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
		
