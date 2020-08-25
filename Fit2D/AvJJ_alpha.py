import ROOT
import numpy
import os
import sys
RDF = ROOT.ROOT.RDataFrame

if not os.path.exists('AvJJ_alpha_test'):
    os.makedirs('AvJJ_alpha_test')
os.chdir('AvJJ_alpha_test')

resbins = [0.0, 30.0, 50.0, 70.0, 90.0, 110.0, 130.0, 150.0, 170.0, 190.0, 210.0, 230.0, 250.0, 270.0, 290.0, 310.0, 330.0, # uniform bins\ 
    350.0, 370.0, 390.0, 412.0, 433.0, 455.0, 478.0, 502.0, 526.0, 551.0, 577.0, 603.0, 631.0, 659.0, 688.0, 717.0, 748.0, 779.0, 812.0, 845.0, 879.0, 914.0, 950.0, 988.0, 1026.0, 1066.0, 1106.0, 1148.0, 1191.0, 1235.0, 1281.0, 1328.0, 1376.0, 1426.0, 1477.0, 1529.0, 1583.0, # resolution-based binning\
    1640, 1700, 1760, 1820, 1880, 1940, 2000, 2060, 2120, 2180, 2240, 2300, 2360, 2420, 2480, 2540, 2600, 2660, 2720, 2780, 2840, 2900, 2960, 3020, 3060, 3120, 3180, 3240, 3300, 3360, 3420, 3480, 3540, 3600, 3660, 3720, 3780, 3840, 3900, 3960, 4020, 4080, 4140, 4200, 4260, 4320, 4380, 4440, 4500, 4560, 4620, 4680, 4740, 4800, 4860, 4920, 4980, 5040] # uniform bins

def MakeNBinsFromMinToMax(N,Min,Max):
    BINS = []
    for i in range(N+1):
        BINS.append(Min+(i*(Max-Min)/N))
    return BINS

def RDFDeltaR(rdf, name, jet1, jet2):
    # ROOT.gInterpreter.ProcessLine("auto Use_Fit_"+newname+" = "+fit+";")
    usefit_code = 	'''
                    #include <TROOT.h>
                    float DeltaR(float pt1, float pt2, float eta1, float eta2, float phi1, float phi2, float m1, float m2)
                    {
                        TLorentzVector j1 = TLorentzVector();
                        j1.SetPtEtaPhiM(pt1, eta1, phi1, m1);
                        TLorentzVector j2 = TLorentzVector();
                        j2.SetPtEtaPhiM(pt2, eta2, phi2, m2);
                        return j1.DeltaR(j2);
                    }
                    '''
    ROOT.gInterpreter.Declare(usefit_code)
    new_rdf = rdf.Define(name, "DeltaR(j%d_pt, j%d_pt, j%d_eta, j%d_eta, j%d_phi, j%d_phi, j%d_m, j%d_m)" % (jet1,jet2,jet1,jet2,jet1,jet2,jet1,jet2))
    return new_rdf

# def MakeAToy(Hi):
#     R = numpy.random
#     Ho = Hi.Clone("H"+str(R.randint(1)))
#     Hw = Hi.Clone("Hw"+str(R.randint(1)))
#     for j in range(1,Ho.GetNbinsY()):
#         for i in range(1,Ho.GetNbinsX()):
#             v = Ho.GetBinContent(i,j)
#             e = Ho.GetBinError(i,j)
#             n = max(0,numpy.rint(v + (e*R.normal())))
#             Hw.SetBinContent(i,j,int(n))
#             if Ho.GetBinContent(i,j) != 0: Hw.SetBinContent(i,j,Hw.GetBinContent(i,j)/Ho.GetBinContent(i,j))
#             else: Hw.SetBinContent(i,j,0)
#             Hw.SetBinError(i,j,0)
#     Ho.Multiply(Hw)
#     return Ho

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

signal = str(sys.argv[1])
cut = str(sys.argv[2])
try: year = 2000+int(sys.argv[3])
except: year = None

C = ROOT.TChain("tree_nominal")
F = ""
if signal == "data": F = "/cms/xaastorage-2/PicoTrees/4JETS/" + str(year) + "/v6/data_SingleMuon_Run" + str(year) + "_all_" + str(year) + ".root"

if signal == "QCD": F = "/cms/xaastorage-2/PicoTrees/4JETS/" + str(year) + "/v6/QCD_HT_all_" + str(year) + ".root"
if signal == "QCD_v0": F = "/cms/xaastorage-2/PicoTrees/4JETS/" + str(year) + "/v6_v0/QCD_HT_all_" + str(year) + ".root"

if signal == "Diquark_chi500suu2000": F = "/cms/xaastorage-2/PicoTrees/4JETS/2016/v6/Diquark_S2000_chi500_2016.root"
if signal == "Diquark_chi500suu3000": F = "/cms/xaastorage-2/PicoTrees/4JETS/2016/v6/Diquark_S3000_chi500_2016.root"
if signal == "Diquark_chi1000suu4000": F = "/cms/xaastorage-2/PicoTrees/4JETS/2016/v6/Diquark_S4000_chi1000_2016.root"
if signal == "Diquark_chi1000suu5000": F = "/cms/xaastorage-2/PicoTrees/4JETS/2016/v6/Diquark_S5000_chi1000_2016.root"
if signal == "Diquark_chi1800suu8000": F = "/cms/evah/workspace/CMSSW_10_2_2/src/nanoTREE/trees/2017/Diquark_chi1800suu8000/tree_Diquark_chi1800suu8000.root"
# ask for v6 version of this tree
if signal == "Diquark_chi2000suu8400": F = "/cms/xaastorage-2/PicoTrees/4JETS/2016/v6/Diquark_S8400_chi2000_2016.root"
if signal == "Diquark_chi2100suu9000": F = "/cms/xaastorage-2/PicoTrees/4JETS/2016/v6/Diquark_S9000_chi2100_2016.root"
if signal == "Diquark_chi3000suu8000": F = "/cms/xaastorage-2/PicoTrees/4JETS/2016/v6/Diquark_S8000_chi3000_2016.root"

if signal == "Xaa_X2000a500": F = "/home/th544/CMSSW_10_2_2/src/picotreeMaker/trees/2016/Xaa_x2000a500/Xaa_x2000a500_2016_DR.root"

if signal == "Stop_500": F = "/cms/xaastorage-2/PicoTrees/4JETS/2017/v6/RPV_M500_2017.root"
if signal == "Stop_750": F = "/cms/xaastorage-2/PicoTrees/4JETS/2017/v6/RPV_M750_2017.root"
if signal == "Stop_1000": F = "/cms/xaastorage-2/PicoTrees/4JETS/2017/v6/RPV_M1000_2017.root"
if signal == "Stop_2000": F = "/cms/xaastorage-2/PicoTrees/4JETS/2017/v6/RPV_M2000_2017.root"

extractBins = True
if signal == "QCD" and year == 2018: extractBins = False

C.Add(F)
rdf = RDF(C)
htcut = 0

if signal == "data": rdf = rdf.Define("total_weight", "1.")

elif signal == "QCD":
    if year == 2016: rdf = rdf.Define("total_weight", "weight_xsN*weight_PU*35900.") # 2016
    elif year == 2017: rdf = rdf.Define("total_weight", "weight_xsN*weight_PU*49900.") # 2017
    else: #elif year == 2018:
        rdf = rdf.Define("total_weight", "weight_xsN*weight_PU*59740.") # 2018

elif signal == "Diquark_chi1800suu8000": rdf = rdf.Define("total_weight", "weight*59740.")

elif signal == "Xaa_X2000a500": rdf = rdf.Define("total_weight", "weight*0.0001228*36000./0.0000008862")

else: rdf = rdf.Define("total_weight", "weight_xsN*weight_PU*59740")

rdf = rdf.Define("testmass", "evt_2JetM")
rdf = rdf.Define("alpha", "evt_2JetM/evt_4JetM")
if not ("/v6/" in F):
    rdf = RDFDeltaR(rdf, "dj1_dR", 1, 2)
    rdf = RDFDeltaR(rdf, "dj2_dR", 3, 4)
rdf = rdf.Define("avgDR", "(dj1_dR + dj2_dR)/2.")

cuts = "evt_Masym < 0.1 && evt_Deta < 1.1"
if "HT" in cut: cuts = cuts + " && evt_HTAK4 > %f" % htcut
if "dR" in cut: cuts = cuts + " && dj1_dR < 2.0 && dj2_dR < 2.0"
if "M4J" in cut: cuts = cuts + " && evt_4JetM > 1530."

print("|===> Cuts: " + cuts)
finebins = MakeNBinsFromMinToMax(1000, 0.1, 0.4)
Mbins = MakeNBinsFromMinToMax(400, 0.0, 4000.0)

AvJJ_alpha_lazy = rdf.Filter(cuts).Histo2D(("AvJJ_alpha", signal + " #alpha v.s. #bar{M}_{jj};#bar{M}_{jj} [GeV];#alpha", len(resbins)-1, numpy.array(resbins), len(finebins)-1, numpy.array(finebins)), "evt_2JetM","alpha","total_weight")
AvJJ_alpha = AvJJ_alpha_lazy.GetValue()

binedges = []
# rebin into equal bins -- AvJJ_alpha, y axis
if extractBins is False:
    total = AvJJ_alpha.Integral()
    targetbins = 100
    limit = total/targetbins
    sum = 0
    sums = []
    for n in range(1, AvJJ_alpha.GetYaxis().GetNbins()):
        varInt = AvJJ_alpha.Integral(1, -1, n, n)
        sum = sum + varInt
        if sum >= limit:
            binedges.append(AvJJ_alpha.GetYaxis().GetBinUpEdge(n))
            sums.append(sum)
            sum = 0
    binedges.append(finebins[-1])
    binO = open("binEdges.txt", "w")
    binO.write(str(binedges))
    binO.close()
else:
    binF = open("binEdges.txt", "r").read()
    binedges = eval(binF)
AvJJ_alpha_rebin_lazy = rdf.Filter(cuts).Histo2D(("AvJJ_alpha_rebin_pre", signal + " #alpha v.s. #bar{M}_{jj}, equal #alpha bins;#bar{M}_{jj} [GeV];#alpha", len(Mbins)-1, numpy.array(Mbins), len(binedges)-1, numpy.array(binedges)), "evt_2JetM","alpha","total_weight")
AvJJ_alpha_rebin = AvJJ_alpha_rebin_lazy.GetValue()

HBW = AvJJ_alpha_rebin.Clone("widthH")
for i in range(1, AvJJ_alpha_rebin.GetNbinsX()+1):
    for j in range(1, AvJJ_alpha_rebin.GetNbinsY()+1):
        HBW.SetBinContent(i, j, AvJJ_alpha_rebin.GetYaxis().GetBinWidth(j)*500)
        HBW.SetBinError(i, j, 0)
AvJJ_alpha_rebin.Divide(HBW)
AvJJ_alpha_rebin.Sumw2()
AvJJ_alpha_rebin.SetStats(0)

# if signal == "QCD": AvJJ_alpha_rebin = MakeAToy(AvJJ_alpha_rebin)
AvJJ_alpha_rebin.SetName("AvJJ_alpha_rebin")

C2_2d = ROOT.TCanvas()
C2_2d.cd()
C2_2d.SetLogz()
AvJJ_alpha_rebin.Draw("colz")
C2_2d.Print(signal + ("_" + str(year) if year is not None else "") + "_MJJ_alpha" + ("_" + cut if cut != "" else "") + ".png")

outfile = ROOT.TFile(signal + ("_" + str(year) if year is not None else "") + "_rebin.root", "recreate")
outfile.cd()
AvJJ_alpha_rebin.Write()
outfile.Close()
