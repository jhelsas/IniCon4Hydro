#!/usr/bin/python

import sys
import os
import ROOT

from ROOT import gROOT, TCanvas, TGraph2D, TLegend

gROOT.Reset()

#if len(sys.argv)!=2:
#  print 'Please, provide a file containing the densities.'
#  exit()

fdens = sys.argv[1]

c1 = TCanvas('c1')

gr1 = TGraph2D(fdens, '%lg %lg %*s %lg')
gr1.SetName('SPHdens')
gr1.SetLineColor(2)
gr1.GetXaxis().SetTitle('x [fm]')
gr1.GetYaxis().SetTitle('y [fm]')
gr1.Draw('surf')

#gr2 = TGraph2D('ploting.dat', '%lg %lg %*s %lg')
#gr2.SetName('Oridens')
#gr2.SetLineColor(4)
#gr2.Draw('surfsame')

#c2 = TCanvas('c2')
#gr3 = TGraph2D('ploting.dat', '%lg %lg %*s %*s %lg')
#gr3.SetName('Densdiff')
#gr3.SetLineColor(1)
#gr3.GetXaxis().SetTitle('x [fm]')
#gr3.GetYaxis().SetTitle('y [fm]')
#gr3.Draw('surf')

# wait for input to keep the GUI (which lives on a ROOT event dispatcher) alive
if __name__ == '__main__':
    rep = ''
    while not rep in [ 'q', 'Q' ]:
        rep = raw_input( 'enter "q" to quit: ' )
        if 1 < len(rep):
            rep = rep[0]
