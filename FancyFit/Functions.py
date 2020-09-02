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
def BigFit(F, O):
	X = "x/13000.0"
	P0 = "("
	P1 = "("
	if F == "P3": P2 = "("
	for p in range(int(O)+1):
		P0 += ("["+str(p)+"]"+(p*"*y"))
		P1 += ("["+str(p+int(O)+1)+"]"+(p*"*y"))
		if F == "P3":
			P2 += ("["+str(p+(2*int(O))+2)+"]"+(p*"*y"))
		if p != int(O):
			P0+= " + "
			P1+= " + "
			if F == "P3": P2+= " + "
		else:
			P0+= ")"
			P1+= ")"
			if F == "P3": P2+= ")"
	print P0
	print P1
	if F == "P3": print P2
	if F == "P2": FitFunc = "TMath::Power((1.-"+X+"),"+P0+")/TMath::Power("+X+","+P1+")"
	if F == "P3": FitFunc = "TMath::Power((1.-"+X+"),"+P0+")/TMath::Power("+X+","+P1+"+("+P2+"*TMath::Log("+X+")))"
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
		
