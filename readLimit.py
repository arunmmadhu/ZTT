#!/usr/bin/env python

import ROOT
from ROOT import TFile, TTree, TCanvas, TGraph, TMultiGraph, TGraphErrors, TLegend, TPaveLabel, TPaveText, TLatex
import CMS_lumi, tdrstyle
import subprocess # to execute shell command
import argparse
import numpy as np
import tdrstyle
from CMSStyle import CMS_lumi
import os


parser = argparse.ArgumentParser()
parser.add_argument('--category'       ,  help="Category; [Default: %(default)s] "                               , dest='category'          , default='taumu')
parser.add_argument('--inputdir'       ,  help="Output directory; [Default: %(default)s] "                       , dest='inputdir'          , default='flat_fit')
args = parser.parse_args()

input_dir  = args.inputdir

ROOT.gROOT.SetBatch(ROOT.kTRUE)
 
# CMS style
CMS_lumi.cmsText = "CMS, work in progress"
CMS_lumi.extraText = ""
CMS_lumi.cmsTextSize = 0.65
CMS_lumi.outOfFrame = True
tdrstyle.setTDRStyle()


#Function to generate BDT cuts
def generate_bdt_cuts(minimum, maximum, median):
    # Ensure minimum, median, and maximum are properly ordered
    if not (minimum < median < maximum):
        raise ValueError("Minimum, median, and maximum must be ordered as minimum < median < maximum.")
    
    nominal_step = 0.0001
    step_1 = 0.04
    step_2 = 0.06
    
    # How close to the median do we have granular steps (decrease to have a smaller region of granularity)
    how_close = 0.5
    
    # Generate more points around the median
    left_side_ll = np.arange(minimum, round((minimum*how_close+median*(1-how_close)), 2) , step_2)
    #left_side_l = np.arange(round((minimum*how_close+median*(1-how_close)), 2), median, step_1)
    #right_side_r = np.arange(median, round((maximum*how_close+median*(1-how_close)), 2), step_1)
    center = np.arange(round((minimum*how_close+median*(1-how_close)), 2), round((maximum*how_close+median*(1-how_close)), 2), step_1)
    right_side_rr = np.arange(round((maximum*how_close+median*(1-how_close)), 2), maximum + nominal_step, step_2)
    
    # Combine and round to two decimal places
    #bdt_cuts = np.concatenate((left_side_ll, left_side_l, right_side_r, right_side_rr))
    bdt_cuts = np.concatenate((left_side_ll, center, right_side_rr))
    bdt_cuts = np.round(bdt_cuts, 2)
    bdt_cuts = np.unique(bdt_cuts)
    
    return bdt_cuts.tolist()

 

 
# EXECUTE datacards
def executeDataCards(labels,values, category):
 
    for value in values:
        label = "%s" % (value)
        combine_command = "combineTool.py -M AsymptoticLimits  -n %s -d %s --cl 0.90 " % (category+label,input_dir+'/datacards/'+category+'/ZTT_T3mu_'+category+'_bdtcut'+label+'.txt')
#        combine_command = "combineTool.py -M BayesianSimple  -n %s -d %s --rMin 0 --rMax 50 --cl 0.90 " % (category+label,input_dir+'/datacards/'+category+'/ZTT_T3mu_'+category+'_bdtcut'+label+'.txt')
#        combine_command = "combineTool.py -M HybridNew --LHCmode LHC-limits  -n %s -d %s --rMin 0 --rMax 50 --cl 0.90 -t 5 --expectedFromGrid 0.5" % (category+label,input_dir+'/datacards/'+category+'/ZTT_T3mu_'+category+'_bdtcut'+label+'.txt')

#        combine_command = "combineTool.py -M HybridNew --LHCmode LHC-limits  -n %s -d %s --rMin 0 --rMax 50 --cl 0.90 -t 5 --expectedFromGrid 0.5" % (category+label,'datacards/'+category+'/ZTT_T3mu_'+category+'_bdtcut'+label+'.txt')

        print ""
        print ">>> " + combine_command
        os.system(combine_command)
        print ">>>   higgsCombine"+category+label+".Asymptotic.mH125.root created"



def executeDataCards_onCondor(labels,values, category):
 
    for value in values:
        label = "%s" % (value)
#        combine_command = "combineTool.py -M AsymptoticLimits  -n %s -d %s --cl 0.90  --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name HybridTest%s " % (category+label,input_dir+'/datacards/'+category+'/ZTT_T3mu_'+category+'_bdtcut'+label+'.txt',category+label)
#        combine_command = "combineTool.py -M BayesianSimple  -n %s -d %s --rMin 0 --rMax 50 --cl 0.90 " % (category+label,input_dir+'/datacards/'+category+'/ZTT_T3mu_'+category+'_bdtcut'+label+'.txt')
        combine_command = "combineTool.py -M HybridNew --LHCmode LHC-limits  -n %s -d %s --rMin 0 --rMax 50 --cl 0.90 -t 25 --expectedFromGrid 0.5 --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name HybridTest%s " % (category+label,input_dir+'/datacards/'+category+'/ZTT_T3mu_'+category+'_bdtcut'+label+'.txt',category+label)

        # this on needs to add to submit to condor:
        # --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name HybridTest
        print ""
        print ">>> " + combine_command
        os.system(combine_command)
        print ">>>   higgsCombine"+category+label+".Asymptotic.mH125.root created"
        
 
 
# GET limits from root file
def getLimits(file_name):
 
    file = TFile(file_name)
    tree = file.Get("limit")
 
    limits = [ ]
    for quantile in tree:
        limits.append(tree.limit)
        print ">>>   %.2f" % limits[-1]
 
    return limits[:6]
 
 
# PLOT upper limits
def plotUpperLimits(labels,values,prefix,outputLabel):
 
    N = len(labels)
    yellow = TGraph(2*N)   
    green = TGraph(2*N)    
    median = TGraph(N)     
 
    up2s = [ ]
    upm  = [ ]
    text_limits=open("TextLimits%s"%(prefix)+outputLabel+".txt","w")
    for i in range(N):
#        file_name = "higgsCombine"+prefix+labels[i]+".AsymptoticLimits.mH120.root"
#        file_name = "higgsCombine"+prefix+labels[i]+".HybridNew.mH120.quant0.500.root"
        file_name = "higgsCombine"+prefix+labels[i]+".HybridNew.mH120.123456.quant0.500.root"

        print "filename:  ", file_name
        limit = getLimits(file_name)
#        up2s.append(limit[4])
        upm.append(limit[2])
#        yellow.SetPoint(    i,    values[i], limit[4]) # + 2 sigma
#        green.SetPoint(     i,    values[i], limit[3]) # + 1 sigma
        median.SetPoint(    i,    values[i], limit[2]) #    median
#        green.SetPoint(  2*N-1-i, values[i], limit[1]) # - 1 sigma
#        yellow.SetPoint( 2*N-1-i, values[i], limit[0]) # - 2 sigma
        text_limits.write("bdt %.2f     median exp %.2f\n"%(values[i],limit[2]))

    W = 800
    H  = 600
    T = 0.08*H
    B = 0.12*H
    L = 0.12*W
    R = 0.04*W
    c = TCanvas("c","c",100,100,W,H)
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)
    c.SetLeftMargin( L/W )
    c.SetRightMargin( R/W )
    c.SetTopMargin( T/H )
    c.SetBottomMargin( B/H )
    c.SetGrid()
    c.SetFrameLineWidth(2);
    c.SetTickx();
    c.SetTicky();
    c.cd()

    frame = c.DrawFrame(1.4,0.001, 4.1, 1.2)
    frame.GetYaxis().CenterTitle()

    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleOffset(0.9)
    frame.GetXaxis().SetNdivisions(508)
    frame.GetYaxis().CenterTitle(False)
    frame.GetYaxis().SetTitle("B(#tau #rightarrow #mu#mu#mu) UL (10^{-7})")
    frame.GetXaxis().SetTitle("MVA cut value")




    #frame.SetMinimum(median.GetHistogram().GetMinimum()*0.8)
    #frame.SetMaximum(median.GetHistogram().GetMinimum()+(5.0/3.0)*(median.GetHistogram().GetMaximum()-median.GetHistogram().GetMinimum()))
    
    frame.SetMinimum(0.0)
    frame.SetMaximum(25.0)

#    frame.GetXaxis().SetLimits(min(values),max(values)*1.2)
    frame.GetXaxis().SetLimits(min(values) ,max(values)*1.2)


 
    yellow.SetFillColor(ROOT.kOrange)
    yellow.SetLineColor(ROOT.kOrange)
    yellow.SetFillStyle(1001)
#    yellow.Draw('F')
 
    green.SetFillColor(ROOT.kGreen+1)
    green.SetLineColor(ROOT.kGreen+1)
    green.SetFillStyle(1001)
#    green.Draw('Fsame')
 
    median.SetLineColor(1)
    median.SetLineWidth(2)
    median.SetLineStyle(2)
    median.Draw('Lsame')
    median.Draw()
 
#    CMS_lumi.CMS_lumi(c,13,11)
    ROOT.gPad.SetTicks(1,1)
    CMS_lumi(ROOT.gPad, 5, 0)
    ROOT.gPad.Update()
    frame.Draw('sameaxis')
 
    x1 = 0.15
    x2 = x1 + 0.24
    y2 = 0.86
    y1 = 0.70
    legend = TLegend(x1,y1,x2,y2)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.041)
    legend.SetTextFont(42)
#    legend.AddEntry(median, "Asymptotic CL_{s} expected upper limit",'L')
    legend.AddEntry(median, "HybridNew CL_{s} expected upper limit",'L')
#    legend.AddEntry(green, "#pm 1 std. deviation",'f')
#    legend.AddEntry(yellow,"#pm 2 std. deviation",'f')


    legend.Draw()
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    Text = ''
    if prefix=='taue':
        Text = 'Category: Z#rightarrow#tau_{e}#tau_{3#mu}'

    if prefix=='taumu':
        Text = 'Category: Z#rightarrow#tau_{#mu}#tau_{3#mu}'

    if prefix=='tauhA':
        Text = 'Category: Z#rightarrow#tau_{h,1-prong}#tau_{3#mu}'

    if prefix=='tauhB':
        Text = 'Category: Z#rightarrow#tau_{h,3-prong}#tau_{3#mu}'
        
    if prefix=='all':
        Text = 'Category: Z#rightarrow#tau#tau_{3#mu}'

    latex.SetTextAlign(1)
    latex.DrawLatex(0.15, 0.85, Text)
    latex.Draw('same') 
    print " "
    c.SaveAs("Limit_scan_Category_"+prefix+outputLabel+".png")
    c.Close()
 
 
# RANGE of floats
def frange(start, stop, step):
    i = start
    while i <= stop:
        yield i
        i += step
 
 
# MAIN
def main():
 

    #bdt_cuts = [-0.4,  -0.2,  0.00, 0.1, 0.2, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39,  0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.6, 0.7]      # a few entries for test
    #bdt_cuts = [-0.4,  -0.2,  0.00, 0.1, 0.2, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37]      # a few entries for test

#    categories = ['taue','taumu','tauhA','tauhB','all']
#    categories = ['taue']
    categories = ['taue','taumu','tauhA','tauhB','all']

    bdt_cuts_taue = generate_bdt_cuts(-0.2, 0.6, 0.2)
    bdt_cuts_taumu = generate_bdt_cuts(0.0, 0.7, 0.4)
    bdt_cuts_tauhA = generate_bdt_cuts(-0.2, 0.6, 0.32)
    bdt_cuts_tauhB = generate_bdt_cuts(-0.4, 0.55, 0.15)
    bdt_cuts_all = generate_bdt_cuts(0.15, 0.85, 0.55)
    
    #bdt_cuts_all = [0.33, 0.37, 0.41, 0.45, 0.49, 0.53, 0.57, 0.61, 0.65, 0.68, 0.74, 0.8]
    
    bdt_cut_dict = {
    'all': bdt_cuts_all,
    'taue': bdt_cuts_taue,
    'taumu': bdt_cuts_taumu,
    'tauhA': bdt_cuts_tauhA,
    'tauhB': bdt_cuts_tauhB
    }

    print "category", args.category 
    outputLabel = ''
    for cat in categories:
        labels = [ ]
        values = [ ]
        
        bdt_cuts = bdt_cut_dict[cat]

        for cl in bdt_cuts:
            values.append(cl)
            label = "%s" % (cl)
            labels.append(label)
        print "values", values
        print "labels", labels
        print "prefix", cat
#        executeDataCards(labels,values,cat)
#        executeDataCards_onCondor(labels,values,cat)

        plotUpperLimits(labels,values,cat,outputLabel)
 
 
 
if __name__ == '__main__':
    main()
