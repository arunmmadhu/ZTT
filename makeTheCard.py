#! /bin/env python


import ROOT
import os
from math import pi, sqrt
from glob import glob
from pdb import set_trace
from array import array 
from numpy import *
import math
import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--selection'      ,  help="Introduce your selection; [Default: %(default)s] "               , dest='selection'         , default='bdt_cv > -0.1')
parser.add_argument('--signalnorm'     ,  help="Signal Normalization; [Default: %(default)s] "                   , dest='signalnorm'        , type = float, default=0.005)
parser.add_argument('--category'       ,  help="Category; [Default: %(default)s] "                               , dest='category'          , default='taumu')
parser.add_argument('--bdt_point'      ,  help="Prefix_output; [Default: %(default)s] "                          , dest='bdt_point'         , default='0.0')
parser.add_argument('--datafile'       ,  help="Input Mini Tree; [Default: %(default)s] "                        , dest='datafile'          , default='Combine_Tree_ztau3mutau.root')
parser.add_argument('--nbins'          ,  help="Number of bins in the mass spectra; [Default: %(default)s] "     , dest='nbins'             , type = int  , default=30)
parser.add_argument('--signal_range_lo',  help="Signal mass window low edge; [Default: %(default)s] "            , dest='signal_range_lo'   , default=1.73)
parser.add_argument('--signal_range_hi',  help="Signal mass window high edge; [Default: %(default)s] "           , dest='signal_range_hi'   , default=1.82)
parser.add_argument('--fit_range_lo'   ,  help="Overal fit range, low edge; [Default: %(default)s] "             , dest='fit_range_lo'      , default=1.60)
parser.add_argument('--fit_range_hi'   ,  help="Overal fit range, high edge; [Default: %(default)s] "            , dest='fit_range_hi'      , default=2.02)
parser.add_argument('--blinded'        ,  help="Blind the signal range; [Default: %(default)s] "                 , dest='blinded'           , action='store_true', default = True )

parser.set_defaults(blinded=True)




args = parser.parse_args() 

blinded         = args.blinded
selection       = args.selection
nbins           = args.nbins
signal_range_lo = double(args.signal_range_lo)
signal_range_hi = double(args.signal_range_hi)
fit_range_lo    = double(args.fit_range_lo)
fit_range_hi    = double(args.fit_range_hi)
minitree        = args.datafile

print "selection: ", selection
print "bdtprefix: ", args.bdt_point

ROOT.gStyle.SetOptStat(True)
ROOT.TH1.SetDefaultSumw2()


#if len(args.category)>0:
#    args.category = '_'+args.category

gaus = ROOT.TF1('signalgaus', 'gaus', 1.7, 1.86)
MiniTreeFile = ROOT.TFile.Open(args.datafile)
MiniTreeFile.cd()

treeName=''
if(args.category == 'taue'):treeName  = 'ztau3mutaue'
if(args.category == 'taumu'):treeName = 'ztau3mutaumu'
if(args.category == 'tauhA'):treeName = 'ztau3mutauh_A'
if(args.category == 'tauhB'):treeName = 'ztau3mutauh_B'




tree = MiniTreeFile.Get(treeName)

mass_histo_mc = ROOT.TH1F('mass_histo_mc', 'mass_histo_mc', nbins, 1.6, 2.)
#tree.Draw('tripletMass>>mass_histo_mc', '(' + selection + '& dataMCtype !=1 ' + ') * event_weight * %f' %args.signalnorm)
tree.Draw('tripletMass>>mass_histo_mc', '(' + selection + '& isMC !=0 ' + ') * weight * %f' %args.signalnorm)
#print "INTEGRAL ",mass_histo_mc.Integral()



tripletMass          = ROOT.RooRealVar('tripletMass'                , '3#mu mass'           , fit_range_lo, fit_range_hi, 'GeV')
bdt_cv               = ROOT.RooRealVar('bdt_cv'                     , 'bdt_cv'              , -1 , 1)
dimu_OS1             = ROOT.RooRealVar('dimu_OS1'                   , 'dimu_OS1'            ,  0 , 2)
dimu_OS2             = ROOT.RooRealVar('dimu_OS2'                   , 'dimu_OS2'            ,  0 , 2)
event_weight         = ROOT.RooRealVar('weight'                     , 'event_weight'        ,  0,  5)  # this weight includes also the scale  mc signal scale
category             = ROOT.RooRealVar('category'                   , 'category'            ,  0,  5)
isMC                 = ROOT.RooRealVar('isMC'                       , 'isMC'                ,  0,  100)
scale                = ROOT.RooRealVar('scale'                      , 'scale'               ,  args.signalnorm)  

print isMC.getVal()

tripletMass.setRange('left' , fit_range_lo    , signal_range_lo)
tripletMass.setRange('right', signal_range_hi , fit_range_hi)
tripletMass.setRange('full' , fit_range_lo    , fit_range_hi)



slope = ROOT.RooRealVar('slope', 'slope', -0.001, -5, 5)
expo  = ROOT.RooExponential('bkg_expo', 'bkg_expo', tripletMass, slope)



nbkg = ROOT.RooRealVar('nbkg', 'nbkg', 1000, 0, 500000)
expomodel = ROOT.RooAddPdf('bkg_extended_expo', 'bkg_extended_expo', ROOT.RooArgList(expo), ROOT.RooArgList(nbkg))


mean  = ROOT.RooRealVar('mean' , 'mean' ,   1.78, 1.6, 1.9)
width = ROOT.RooRealVar('width', 'width',   0.02, 0.0, 0.1)
gaus  = ROOT.RooGaussian('sig_gaus', 'sig_gaus', tripletMass, mean, width)


cbwidth = ROOT.RooRealVar('cbwidth','cbwidth', 0.02, 0.0, 0.1)
cbalpha = ROOT.RooRealVar('cbalpha','cbalpha', 1.00, -20, 20 )
cbn     = ROOT.RooRealVar('cbn'    ,'cbn'    , 2,    0  , 5  )
cb      = ROOT.RooCBShape('cb'     ,'cb'     , tripletMass, mean, cbwidth, cbalpha, cbn)



fraction    = ROOT.RooRealVar('fraction' , 'fraction',   0.01, 0, 0.1)
#combined_model = ROOT.RooAddPdf("combined_model", "g+a", ROOT.RooArgList(gaus,cb), ROOT.RooArgList(gs_fraction,cb_fraction))
combined_model = ROOT.RooAddPdf("combined_model", "g+a", ROOT.RooArgList(gaus,cb), ROOT.RooArgList(fraction))




variables = ROOT.RooArgSet()
variables.add(tripletMass)
variables.add(bdt_cv)
variables.add(dimu_OS1)
variables.add(dimu_OS2)
variables.add(event_weight)
variables.add(category)
variables.add(isMC)
variables.add(scale)


MCSelector = ROOT.RooFormulaVar('MCSelector', 'MCSelector', selection + ' & isMC !=0 ', ROOT.RooArgList(variables))
fullmc = ROOT.RooDataSet('mc', 'mc', tree, variables, MCSelector,'scale')


frame = tripletMass.frame()
frame.SetTitle('')

fullmc.plotOn(frame, 
              ROOT.RooFit.Binning(nbins), 
              ROOT.RooFit.XErrorSize(0), 
              ROOT.RooFit.LineWidth(2),
              ROOT.RooFit.MarkerStyle(6),
              ROOT.RooFit.MarkerColor(ROOT.kCyan + 2),
              ROOT.RooFit.MarkerSize(0.75),
              ROOT.RooFit.FillColor(ROOT.kCyan  + 2),
)

#print " ", fullmc.sumEntries()


#####  tese several models; finally keep only gaus
results_gaus = gaus.fitTo(fullmc, ROOT.RooFit.Range(signal_range_lo, signal_range_hi), ROOT.RooFit.Save())
results_gaus.Print()


#results_cb   = cb.fitTo(fullmc, ROOT.RooFit.Range(signal_range_lo, signal_range_hi), ROOT.RooFit.Save())
#results_cb.Print()
 

#results_combined_pdf = combined_model.fitTo(fullmc, ROOT.RooFit.Range(signal_range_lo, signal_range_hi), ROOT.RooFit.Save())
#results_combined_pdf.Print()


#results_gaus = gaus.fitTo(fullmc, ROOT.RooFit.Save(), ROOT.RooFit.SumW2Error(True))
#results_gaus = gaus.chi2FitTo(fullmc, ROOT.RooFit.Save())

gaus.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed +2 ))
#cb.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kMagenta +2 ))
#combined_model.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kCyan +2 ))



frame.Draw()
cmd1 = 'mkdir plots; mkdir datacards; mkdir workspaces;'
cmd2 = 'mkdir plots/%s;'%args.category + 'mkdir datacards/%s;'%args.category + 'mkdir workspaces/%s;'%args.category
os.system(cmd1)
os.system(cmd2)
latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextAngle(0)
latex.SetTextFont(42)
latex.SetTextAlign(31)
latex.DrawLatex(0.57, 0.85, 'taushape%s_%dbins'%(args.category, nbins))

ROOT.gPad.SaveAs('plots/%s/taushape%s_%dbins_bdtcut%s.png'%(args.category,args.category, nbins, args.bdt_point))


DataSelector      = ROOT.RooFormulaVar('DataSelector', 'DataSelector', selection + ' & isMC == 0', ROOT.RooArgList(variables))
BlindDataSelector = ROOT.RooFormulaVar('DataSelector', 'DataSelector', selection + ' & isMC == 0 &  abs(tripletMass  - 1.776) > %s' %( (signal_range_hi -signal_range_lo)/2) , ROOT.RooArgList(variables))




if blinded:
    print 'BLIND'
    fulldata     = ROOT.RooDataSet('data', 'data', tree,  variables, BlindDataSelector)
else:
    fulldata     = ROOT.RooDataSet('data', 'data', tree,  variables, DataSelector)



fulldata.plotOn(frame, 
                ROOT.RooFit.Binning(nbins),
                ROOT.RooFit.MarkerStyle(21),
                ROOT.RooFit.MarkerColor(ROOT.kBlack), 
                ROOT.RooFit.MarkerSize(0.75))


if blinded:
    results_expo = expomodel.fitTo(fulldata, ROOT.RooFit.Range('left,right'), ROOT.RooFit.Save())
else:
    results_expo = expomodel.fitTo(fulldata, ROOT.RooFit.Range('full'), ROOT.RooFit.Save())




expomodel.plotOn(frame,  ROOT.RooFit.LineColor(ROOT.kBlack) )



ctmp_canvas = ROOT.TCanvas('Category %s' %args.category,"Categories",0,0,700,500);
ctmp_canvas.SetFrameLineWidth(3);
ctmp_canvas.SetTickx();
ctmp_canvas.SetTicky();




frame.Draw()
legend = ROOT.TLegend(0.12,0.70,0.43,0.86)
legend.AddEntry(frame.getObject(1),"Signal model gauss","L")
legend.AddEntry(frame.getObject(4),"data model exp","L")
#legend.Draw()


latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextAngle(0)
latex.SetTextFont(42)
latex.SetTextAlign(31)
latex.DrawLatex(0.57, 0.85, 'mass_fit%s_%dbins_bdtcut%s'%(args.category, nbins, args.bdt_point))


ROOT.gPad.SaveAs('plots/%s/massfit_%s_%dbins_bdtcut%s.png'%(args.category,args.category, nbins, args.bdt_point))


output = ROOT.TFile.Open('workspaces/%s/workspace%s_bdtcut%s.root' %(args.category,args.category, args.bdt_point), 'recreate')
seldata = ROOT.RooDataSet('data', 'data', tree, variables, DataSelector)
selsignal = ROOT.RooDataSet('mc', 'mc', tree, variables, MCSelector)



data =  ROOT.RooDataSet(
    'data_obs', 
    'data_obs',
    seldata, 
    ROOT.RooArgSet(tripletMass)
)

mc =  ROOT.RooDataSet(
    'mc', 
    'mc',
    selsignal, 
    ROOT.RooArgSet(tripletMass)
)


# create workspace

print 'creating workspace'
cb_fraction    = ROOT.RooRealVar('cb_fraction' , 'cbfraction' ,   1-fraction.getVal(), 1- fraction.getVal() , 1- fraction.getVal())
gs_fraction    = ROOT.RooRealVar('gs_fraction' , 'gs_fraction',   fraction.getVal(), fraction.getVal(), fraction.getVal())

workspace = ROOT.RooWorkspace('t3m_shapes')

workspace.factory('tripletMass[%f, %f]' % (fit_range_lo, fit_range_hi))
workspace.factory("Exponential::bkg(tripletMass, a0%s[%f,%f,%f])" %(args.category, slope.getVal(), slope.getError(), slope.getError()) )  


workspace.factory('cb_fraction[%f]'  % cb_fraction.getVal())
workspace.factory('gs_fraction[%f]'  % gs_fraction.getVal())
workspace.factory('mean[%f]'  % mean .getVal())
workspace.factory('sigma[%f]' % width.getVal())
workspace.factory('RooGaussian::sig(tripletMass, mean, sigma)')

workspace.factory('cbsigma[%f]' % cbwidth.getVal())
workspace.factory('cbalpha[%f]' % cbalpha.getVal())
workspace.factory('cbn[%f]' % cbn.getVal())
workspace.factory('RooCBShape::cbsig(tripletMass, mean, cbsigma, cbalpha, cbn)')






it = workspace.allVars().createIterator()
all_vars = [it.Next() for _ in range( workspace.allVars().getSize())]
for var in all_vars:
    if var.GetName() in ['mean', 'sigma','cb_fraction','gs_fraction','cbsigma','cbalpha','cbn']:
        var.setConstant()

combined_model_fit = ROOT.RooAddPdf("combined_model_fit", "g+a", ROOT.RooArgList(gaus,cb), ROOT.RooArgList(gs_fraction,cb_fraction), False)
#combined_model_fit = ROOT.RooAddPdf("combined_model_fit", "g+a", ROOT.RooArgList(gaus,cb), fraction, False)
getattr(workspace,'import')(combined_model_fit)

getattr(workspace,'import')(data)
getattr(workspace,'import')(mc)
workspace.Write()
output.Close()



# make  the datacard
with open('datacards/%s/ZTT_T3mu_%s_bdtcut%s.txt' %(args.category,args.category,args.bdt_point), 'w') as card:
   card.write(
'''
# ZTT
imax 1 number of bins
jmax * number of processes minus 1
kmax * number of nuisance parameters
--------------------------------------------------------------------------------
shapes background    category{cat}       ../../workspaces/{wdir}/workspace{cat}.root t3m_shapes:bkg
shapes signal        category{cat}       ../../workspaces/{wdir}/workspace{cat}.root t3m_shapes:sig_gaus
shapes data_obs      category{cat}       ../../workspaces/{wdir}/workspace{cat}.root t3m_shapes:data_obs
--------------------------------------------------------------------------------
bin               category{cat}
observation       {obs:d}
--------------------------------------------------------------------------------

bin                                     category{cat}       category{cat}
process                                 signal              background
process                                 0                   1
rate                                   {signal:.4f}        {bkg:.4f}
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
'''.format(
         cat      = args.category,
         wdir     = args.category,
         obs      = fulldata.numEntries() if blinded==False else -1,
         signal   = mass_histo_mc.Integral(),
#         signal   = fullmc.sumEntries(),         # use the ArgDataSet entries instead
         bkg      = nbkg.getVal(),
         )
)
