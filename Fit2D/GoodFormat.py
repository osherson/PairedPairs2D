from ROOT import gStyle
from ROOT import gROOT
from ROOT import TStyle
from ROOT import TColor
import numpy

def pullFormat(H):
    NRGBs = 8
    NCont = 80
    stops = [ 0.00, 0.15, 0.30, 0.45, 0.55, 0.70, 0.85, 1.00 ]
    #red = [ 1.00, 1.00, 1.00, 1.00, 1.00, 0.60, 0.20, 0.00 ]
    #green = [ 0.00, 0.20, 0.60, 1.00, 1.00, 0.60, 0.20, 0.00 ]
    red = [ 0.20, 0.40, 0.60, 1.00, 1.00, 0.60, 0.40, 0.20 ]
    green = [ 1.00, 1.00, 1.00, 1.00, 1.00, 0.60, 0.20, 0.10 ]
    blue = [ 0.10, 0.20, 0.60, 1.00, 1.00, 1.00, 1.00, 1.00 ]
    TColor.CreateGradientColorTable(NRGBs, numpy.array(stops), numpy.array(red), numpy.array(green), numpy.array(blue), NCont)
    gStyle.SetNumberContours(NCont)
    H.UseCurrentStyle()
    return H
