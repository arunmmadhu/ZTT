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
import CMS_lumi, tdrstyle
from CMSStyle import CMS_lumi


parser = argparse.ArgumentParser()

parser.add_argument('--selection'      ,  help="Introduce your selection; [Default: %(default)s] "               , dest='selection'         , default='bdt_cv > -0.1')
parser.add_argument('--signalnorm'     ,  help="Signal Normalization; [Default: %(default)s] "                   , dest='signalnorm'        , type = float, default=0.5)
parser.add_argument('--category'       ,  help="Category; [Default: %(default)s] "                               , dest='category'          , default='taumu')
parser.add_argument('--bdt_point'      ,  help="Prefix_output; [Default: %(default)s] "                          , dest='bdt_point'         , default='0.0')
parser.add_argument('--datafile'       ,  help="Input Mini Tree; [Default: %(default)s] "                        , dest='datafile'          , default='../../Combine_Tree_ztau3mutau.root')
parser.add_argument('--nbins'          ,  help="Number of bins in the mass spectra; [Default: %(default)s] "     , dest='nbins'             , type = int  , default=40)
parser.add_argument('--signal_range_lo',  help="Signal mass window low edge; [Default: %(default)s] "            , dest='signal_range_lo'   , default=1.74)
parser.add_argument('--signal_range_hi',  help="Signal mass window high edge; [Default: %(default)s] "           , dest='signal_range_hi'   , default=1.81)
parser.add_argument('--fit_range_lo'   ,  help="Overal fit range, low edge; [Default: %(default)s] "             , dest='fit_range_lo'      , default=1.60)
parser.add_argument('--fit_range_hi'   ,  help="Overal fit range, high edge; [Default: %(default)s] "            , dest='fit_range_hi'      , default=2.00)
parser.add_argument('--blinded'        ,  help="Blind the signal range; [Default: %(default)s] "                 , dest='blinded'           , action='store_true', default = False )
parser.add_argument('--alt_pdf'        ,  help="Whether to use an alternate PDF, True; [Default: %(default)s] "  , dest='alt_pdf', action='store_true', default = False )
parser.add_argument('--no-alt_pdf'     ,  help="Whether to use an alternate PDF, False; [Default: %(default)s] " , dest='alt_pdf', action='store_false' )
parser.add_argument('--pdf_switch_point' ,  help="Point where you switch to the alternate PDF; [Default: %(default)s] " , dest='pdf_switch_point'  , default = 0.0 )
parser.add_argument('--fixed_slope'      ,  help="Value of the sloped after the point where you switch to the alternate PDF; [Default: %(default)s] " , dest='fixed_slope'  , default = -0.001 )
parser.add_argument('--outdir'           ,  help="Output directory; [Default: %(default)s] "                     , dest='outdir'            , default='unfixed_slope')

parser.set_defaults(blinded=True)


args            = parser.parse_args() 
blinded         = args.blinded
selection       = args.selection
nbins           = args.nbins
signal_range_lo = double(args.signal_range_lo)
signal_range_hi = double(args.signal_range_hi)
fit_range_lo    = double(args.fit_range_lo)
fit_range_hi    = double(args.fit_range_hi)
minitree        = args.datafile
alt_pdf         = args.alt_pdf
pdf_switch_point= args.pdf_switch_point
output_dir      = args.outdir

print "selection: ", selection
print "bdtprefix: ", args.bdt_point

ROOT.gStyle.SetOptStat(True)
ROOT.TH1.SetDefaultSumw2()

# Enable batch mode
ROOT.gROOT.SetBatch(True)

# CMS style
CMS_lumi.cmsText = "CMS, work in progress"
CMS_lumi.extraText = ""
CMS_lumi.cmsTextSize = 0.65
CMS_lumi.outOfFrame = True
tdrstyle.setTDRStyle()

#if len(args.category)>0:
#    args.category = '_'+args.category

gaus = ROOT.TF1('signalgaus', 'gaus', 1.7, 1.86)
MiniTreeFile = ROOT.TFile.Open(args.datafile)
MiniTreeFile.cd()

treeName=''
if(args.category == 'all'):treeName   = 'ztautau'
if(args.category == 'taue'):treeName  = 'ztau3mutaue'
if(args.category == 'taumu'):treeName = 'ztau3mutaumu'
if(args.category == 'tauhA'):treeName = 'ztau3mutauh_A'
if(args.category == 'tauhB'):treeName = 'ztau3mutauh_B'




tree = MiniTreeFile.Get(treeName)

mass_histo_mc = ROOT.TH1F('mass_histo_mc', 'mass_histo_mc', nbins, fit_range_lo, fit_range_hi)
#tree.Draw('tripletMass>>mass_histo_mc', '(' + selection + '& isMC !=0 ' + ') * weight * %f' %args.signalnorm)
tree.Draw('tripletMass>>mass_histo_mc', '(' + selection + '& isMC !=0 & (isMC == 211 | isMC == 210231 | isMC == 210232 | isMC == 210233 ) ' + ') ' ) # weight is always one


signal_range = mass_histo_mc.Integral(mass_histo_mc.FindFixBin(args.signal_range_lo), mass_histo_mc.FindFixBin(args.signal_range_hi) )
full_range   = mass_histo_mc.Integral(mass_histo_mc.FindFixBin(args.fit_range_lo), mass_histo_mc.FindFixBin(args.fit_range_hi) )

ratioToSignal = signal_range/full_range

SignalIntegral = mass_histo_mc.GetEntries() * ratioToSignal *  args.signalnorm    


tripletMass          = ROOT.RooRealVar('tripletMass'                , '3#mu mass'           , fit_range_lo, fit_range_hi, 'GeV')
bdt_cv               = ROOT.RooRealVar('bdt_cv'                     , 'bdt_cv'              , -1 , 1)
dimu_OS1             = ROOT.RooRealVar('dimu_OS1'                   , 'dimu_OS1'            ,  0 , 2)
dimu_OS2             = ROOT.RooRealVar('dimu_OS2'                   , 'dimu_OS2'            ,  0 , 2)
event_weight         = ROOT.RooRealVar('weight'                     , 'event_weight'        ,  0,  5)  # this weight includes also the scale  mc signal scale
category             = ROOT.RooRealVar('category'                   , 'category'            ,  0,  5)
isMC                 = ROOT.RooRealVar('isMC'                       , 'isMC'                ,  0,  1000000)
scale                = ROOT.RooRealVar('scale'                      , 'scale'               ,  args.signalnorm)  



tripletMass.setRange('left' , fit_range_lo    , signal_range_lo)
tripletMass.setRange('right', signal_range_hi , fit_range_hi)
tripletMass.setRange('full' , fit_range_lo    , fit_range_hi)


tripletMass.setRange("SB1",fit_range_lo,1.75)
tripletMass.setRange("SB2",1.80,fit_range_hi)
tripletMass.setRange("fullRange",fit_range_lo,fit_range_hi);
tripletMass.setRange("SIG",signal_range_lo,signal_range_hi)




# Fixing a certain value of exponential slope after a certain value of the BDT cut
#if alt_pdf and args.bdt_point > pdf_switch_point:
#        myVar.setConstant(True)
        
if alt_pdf and (float(args.bdt_point) > float(pdf_switch_point)):
        slope = ROOT.RooRealVar('slope', 'slope', float(args.fixed_slope), float(args.fixed_slope), float(args.fixed_slope))
else:
        slope = ROOT.RooRealVar('slope', 'slope', float(args.fixed_slope), -100, 100)

expo  = ROOT.RooExponential('bkg_expo', 'bkg_expo', tripletMass, slope)

flat_poly  = ROOT.RooPolynomial('bkg_flat_poly', 'bkg_flat_poly', tripletMass)

nbkg = ROOT.RooRealVar('nbkg', 'nbkg', 1000, 0, 500000)

#Switching to flat polynomial after a certain point
#Redacted

#if alt_pdf and args.bdt_point > pdf_switch_point:
#        pdfmodel = ROOT.RooAddPdf('bkg_flat_expo', 'bkg_flat_expo', ROOT.RooArgList(flat_poly), ROOT.RooArgList(nbkg))
#else:
#        pdfmodel = ROOT.RooAddPdf('bkg_extended_expo', 'bkg_extended_expo', ROOT.RooArgList(expo), ROOT.RooArgList(nbkg))
pdfmodel = ROOT.RooAddPdf('bkg_extended_expo', 'bkg_extended_expo', ROOT.RooArgList(expo), ROOT.RooArgList(nbkg))


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


MCSelector = ROOT.RooFormulaVar('MCSelector', 'MCSelector', selection + ' & isMC !=0 & (isMC == 211 | isMC == 210231 | isMC == 210232 | isMC == 210233 ) ', ROOT.RooArgList(variables))


fullmc = ROOT.RooDataSet('mc', 'mc', tree, variables, MCSelector,'scale')


frame = tripletMass.frame()
frame.SetTitle('')

fullmc.plotOn(frame, 
              ROOT.RooFit.Binning(nbins), 
              ROOT.RooFit.XErrorSize(0), 
              ROOT.RooFit.LineWidth(2),
              ROOT.RooFit.MarkerStyle(6),
              ROOT.RooFit.MarkerColor(ROOT.kRed ),
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

gaus.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed ))
#cb.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kMagenta +2 ))
#combined_model.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kCyan +2 ))



frame.Draw()
cmd0 = 'mkdir '+output_dir+';'
cmd1 = 'mkdir '+output_dir+'/plots; mkdir '+output_dir+'/datacards; mkdir '+output_dir+'/workspaces;'
cmd2 = 'mkdir '+output_dir+'/plots/%s;'%args.category + 'mkdir '+output_dir+'/datacards/%s;'%args.category + 'mkdir '+output_dir+'/workspaces/%s;'%args.category
os.system(cmd0)
os.system(cmd1)
os.system(cmd2)
latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextAngle(0)
latex.SetTextFont(42)
latex.SetTextAlign(31)
latex.DrawLatex(0.57, 0.85, 'taushape%s_%dbins'%(args.category, nbins))

#ROOT.gPad.SaveAs(output_dir+'/plots/%s/taushape%s_%dbins_bdtcut%s.png'%(args.category,args.category, nbins, args.bdt_point))


DataSelector      = ROOT.RooFormulaVar('DataSelector', 'DataSelector', selection + ' & isMC == 0', ROOT.RooArgList(variables))
BlindDataSelector = ROOT.RooFormulaVar('DataSelector', 'DataSelector', selection + ' & isMC == 0 & (tripletMass<=%s || tripletMass>=%s) ' %(signal_range_lo,signal_range_hi) , ROOT.RooArgList(variables))




if blinded:
    print 'BLIND'
    fulldata     = ROOT.RooDataSet('data', 'data', tree,  variables, BlindDataSelector)
else:
    fulldata     = ROOT.RooDataSet('data', 'data', tree,  variables, DataSelector)



fulldata.plotOn(frame, 
                ROOT.RooFit.Binning(nbins),
                ROOT.RooFit.MarkerStyle(20),
                ROOT.RooFit.MarkerColor(ROOT.kBlack), 
                ROOT.RooFit.MarkerSize(0.75))


if blinded:
    results_pdf = pdfmodel.fitTo(fulldata, ROOT.RooFit.Range('left,right'), ROOT.RooFit.Save())
else:
    results_pdf = pdfmodel.fitTo(fulldata, ROOT.RooFit.Range('full'), ROOT.RooFit.Save())

#print("Normalization: ",nbkg.getVal())

SG_integral = pdfmodel.createIntegral(ROOT.RooArgSet(tripletMass), ROOT.RooArgSet(tripletMass), "SIG").getVal()
SB_integral = pdfmodel.createIntegral(ROOT.RooArgSet(tripletMass), ROOT.RooArgSet(tripletMass), "left,right").getVal()


if blinded:
    pdfmodel.plotOn(frame,  ROOT.RooFit.LineColor(ROOT.kBlue) , ROOT.RooFit.Normalization(nbkg.getVal()*SB_integral, ROOT.RooAbsReal.NumEvent), ROOT.RooFit.ProjectionRange('full') )
else:
    pdfmodel.plotOn(frame,  ROOT.RooFit.LineColor(ROOT.kBlue) , ROOT.RooFit.Normalization(nbkg.getVal(), ROOT.RooAbsReal.NumEvent), ROOT.RooFit.ProjectionRange('full') )



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

CMS_lumi(ROOT.gPad, 5, 0)
ROOT.gPad.Update()

ROOT.gPad.SaveAs(output_dir+'/plots/%s/massfit_%s_%dbins_bdtcut%s.png'%(args.category,args.category, nbins, args.bdt_point))


output = ROOT.TFile.Open(output_dir+'/workspaces/%s/workspace%s_bdtcut%s.root' %(args.category,args.category, args.bdt_point), 'recreate')
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
#cb_fraction    = ROOT.RooRealVar('cb_fraction' , 'cbfraction' ,   1-fraction.getVal(), 1- fraction.getVal() , 1- fraction.getVal())
#gs_fraction    = ROOT.RooRealVar('gs_fraction' , 'gs_fraction',   fraction.getVal(), fraction.getVal(), fraction.getVal())

workspace = ROOT.RooWorkspace('t3m_shapes')

workspace.factory('tripletMass[%f, %f]' % (fit_range_lo, fit_range_hi))

if alt_pdf and args.bdt_point > pdf_switch_point:
        workspace.factory("Polynomial::bkg(tripletMass, a0%s[%f])" %(args.category,1.0) ) 
else:
        workspace.factory("Exponential::bkg(tripletMass, a0%s[%f,%f,%f])" %(args.category, slope.getVal(), slope.getError(), slope.getError()) ) 


with open("Slopes_%s"%(args.category)+".txt", "a") as f:
   f.write("Cut: %s slopes: %s n_sideband: %s expected_bkg: %s \n"%(args.bdt_point,slope.getVal(),nbkg.getVal()*SB_integral,nbkg.getVal()*SG_integral if nbkg.getVal()*SG_integral > 0.001 else 0.001))


#workspace.factory('cb_fraction[%f]'  % cb_fraction.getVal())
#workspace.factory('gs_fraction[%f]'  % gs_fraction.getVal())
workspace.factory('mean[%f]'  % mean .getVal())
workspace.factory('sigma[%f]' % width.getVal())
workspace.factory('RooGaussian::sig(tripletMass, mean, sigma)')

#workspace.factory('cbsigma[%f]' % cbwidth.getVal())
#workspace.factory('cbalpha[%f]' % cbalpha.getVal())
#workspace.factory('cbn[%f]' % cbn.getVal())
#workspace.factory('RooCBShape::cbsig(tripletMass, mean, cbsigma, cbalpha, cbn)')

it = workspace.allVars().createIterator()
all_vars = [it.Next() for _ in range( workspace.allVars().getSize())]
for var in all_vars:
#    if var.GetName() in ['mean', 'sigma','cb_fraction','gs_fraction','cbsigma','cbalpha','cbn']:
    if var.GetName() in ['mean', 'sigma']:
        var.setConstant()

#combined_model_fit = ROOT.RooAddPdf("combined_model_fit", "g+a", ROOT.RooArgList(gaus,cb), ROOT.RooArgList(gs_fraction,cb_fraction), False)
#combined_model_fit = ROOT.RooAddPdf("combined_model_fit", "g+a", ROOT.RooArgList(gaus,cb), fraction, False)
#getattr(workspace,'import')(combined_model_fit)

getattr(workspace,'import')(data)
getattr(workspace,'import')(mc)
workspace.Write()
output.Close()

br = ''
if args.category=='taue':
        br = '022'#0.04/17.82

if args.category=='taumu':
        br = '023'#0.04/17.39

if args.category=='tauhA':
        br = '026'#0.16/61.48

if args.category=='tauhB':
        br = '026'#0.16/61.48
        
if args.category=='all':
        br = '000'

# make  the datacard
with open(output_dir+'/datacards/%s/ZTT_T3mu_%s_bdtcut%s.txt' %(args.category,args.category,args.bdt_point), 'w') as card:
   card.write(
'''
# ZTT
imax 1 number of bins
jmax * number of processes minus 1
kmax * number of nuisance parameters
--------------------------------------------------------------------------------
shapes background    category{cat}       ../../workspaces/{wdir}/workspace{cat}_bdtcut{bdt}.root t3m_shapes:bkg
shapes signal        category{cat}       ../../workspaces/{wdir}/workspace{cat}_bdtcut{bdt}.root t3m_shapes:sig
shapes data_obs      category{cat}       ../../workspaces/{wdir}/workspace{cat}_bdtcut{bdt}.root t3m_shapes:data_obs
--------------------------------------------------------------------------------
bin               category{cat}
observation       {obs:d}
--------------------------------------------------------------------------------

bin                                     category{cat}       category{cat}
process                                 signal              background
process                                 0                   1
rate                                   {signal:.4f}        {bkg:.4f}
--------------------------------------------------------------------------------
lumi              lnN                       1.025               -
Zxs               lnN                       1.0188              -
BrTauX            lnN                       1.0{br_taux}        -
bkgNorm_{cat}     rateParam  category{cat} background 1.
a0{cat}           param             {slopeval:.4f}        {slopeerr:.4f}
--------------------------------------------------------------------------------
'''.format(
         cat      = args.category,
         bdt      = args.bdt_point,
         wdir     = args.category,
         obs      = fulldata.numEntries() if blinded==False else -1,
         signal   = SignalIntegral,
         bkg      = nbkg.getVal()*SG_integral if nbkg.getVal()*SG_integral > 0.001 else 0.001,
         slopeval = slope.getVal(), 
         slopeerr = slope.getError(),
         br_taux  = br,
         )
)

#print "SG_integral: ", SG_integral
#print "nbkg.getVal(): ", nbkg.getVal()
#print "expo.getVal(): ", expo.getVal()
#print "nbkg.getVal()*SB_integral: ", nbkg.getVal()*SB_integral



