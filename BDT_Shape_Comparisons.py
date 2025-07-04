#!/usr/bin/env python3

import ROOT
from ROOT import TFile, TTree, TCanvas, TGraph, TMultiGraph, TGraphErrors, TLegend, TPaveLabel, TPaveText, TLatex
from ROOT import RooRealVar, RooFormulaVar, RooExponential, RooDataHist, RooArgList, RooAddPdf, RooFit, RooDataSet, RooGenericPdf, RooBifurGauss
import subprocess # to execute shell command
import argparse
import numpy as np
import tdrstyle
#from CMSStyle import CMS_lumi
import CMSStyle
import os
import re


# CMS style
CMSStyle.cmsText = "CMS"
CMSStyle.extraText = "Work in progress"
CMSStyle.relPosX = 0.070
CMSStyle.outOfFrame = False
CMSStyle.alignX_ = 1
CMSStyle.relPosX    = 0.04
#CMSStyle.relPosY    = 0.025
tdrstyle.setTDRStyle()

# references for T, B, L, R
H_ref = 600; 
W_ref = 800; 
W = W_ref
H  = H_ref
T = 0.08*H_ref
B = 0.18*H_ref 
L = 0.22*W_ref
R = 0.04*W_ref


class BDT_Shape_Comparisons:
        
        def __init__(self):
                
                self.bdt_cv = None
                
                self.bgausmeanMC = None
                self.bgaussigmaMC_a = None
                self.bgaussigmaMC_b = None
                self.bgaus_distMC = None
                self.BDTNorm_MC = None
                self.BDT_distribution_MC = None
                self.MCSelector = None
                self.fullmc = None
                
                self.a = None
                self.b = None
                self.c = None
                self.d = None
                self.quadratic = None
                self.expModel = None
                self.BDTNorm = None
                self.BDT_distribution_MC = None
        
        
        
        def Compare_BDT_Scores_MultipleDatasets(self, datafiles, categ, isMC=False):
            frame = ROOT.RooRealVar('bdt_cv', 'BDT Score', -1, 1).frame()
            colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta]
            
            labels = ["old", "new"]
            
            # Map categ to treename and LaTeX-style label
            if categ == 'taue':
                treename = 'ztau3mutaue'
                cat_label = r"#tau_{e}"
            elif categ == 'taumu':
                treename = 'ztau3mutaumu'
                cat_label = r"#tau_{#mu}"
            elif categ == 'tauhA':
                treename = 'ztau3mutauh_A'
                cat_label = r"#tau_{h,1-prong}"
            elif categ == 'tauhB':
                treename = 'ztau3mutauh_B'
                cat_label = r"#tau_{h,3-prong}"
            elif categ == 'all':
                treename = 'ztautau'
                cat_label = "Inclusive"
            else:
                treename = categ
                cat_label = categ
        
            legend = ROOT.TLegend(0.3, 0.3, 0.55, 0.55)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.SetHeader(cat_label, "C")
            
            for i, datafile in enumerate(datafiles):
                MiniTreeFile = ROOT.TFile.Open(datafile)
                tree = MiniTreeFile.Get(treename)  # e.g. 'ztau3mutauh_A'
                
                # Define all variables only once
                bdt_cv = ROOT.RooRealVar("bdt_cv", "bdt_cv", -1.01, 1.01)
                tripletMass = ROOT.RooRealVar('tripletMass', '3#mu mass', 0, 100)
                isMC = ROOT.RooRealVar("isMC", "isMC", 0, 1000000)
                weight = ROOT.RooRealVar("weight", "weight", 0, 5)
                dimu_OS1 = ROOT.RooRealVar('dimu_OS1', 'dimu_OS1', 0, 1000)
                dimu_OS2 = ROOT.RooRealVar('dimu_OS2', 'dimu_OS2', 0, 1000)
                
                # Choose data or mc without weights
                
                variables_noweight = ROOT.RooArgSet(tripletMass, bdt_cv, isMC, weight, dimu_OS1, dimu_OS2)
                selector_str = "isMC == 0" if not isMC else "( isMC == 211 || isMC == 210231 || isMC == 210232 || isMC == 210233 )"
                dataset_unweighted = ROOT.RooDataSet("data_{0}".format(i), "data_{0}".format(i), tree, variables_noweight, selector_str)
                
                
                # Normalizing the samples to 1.0
                norm_factor = 1.0 / dataset_unweighted.numEntries()
                print("dataset entries: ",dataset_unweighted.numEntries())
                scale = ROOT.RooRealVar("scale", "scale", norm_factor)
                dataset_vars = dataset_unweighted.get()
                dataset_vars.add(scale)
                
                dataset = ROOT.RooDataSet("scaled_ds", "scaled_ds",dataset_unweighted,dataset_vars,"","scale")
                
                curve_name = "curve_%d" % i
                dataset.plotOn(
                    frame,
                    ROOT.RooFit.Binning(100),
                    ROOT.RooFit.MarkerColor(colors[i % len(colors)]),
                    ROOT.RooFit.LineColor(colors[i % len(colors)]),
                    ROOT.RooFit.MarkerStyle(20 + i),
                    ROOT.RooFit.MarkerSize(0.75),
                    ROOT.RooFit.Name(curve_name)
                )
                
                legend.AddEntry(curve_name, labels[i], "lep")
            
            
            #Plotting
            bdt_canvas = ROOT.TCanvas("bdt_canvas", "bdt_canvas", 800, 800)
            bdt_canvas.cd()
            bdt_canvas.SetLeftMargin( L/W )
            bdt_canvas.SetRightMargin( R/W )
            bdt_canvas.SetTopMargin( T/H )
            bdt_canvas.SetBottomMargin( B/H )
            frame.GetXaxis().SetTitle("BDT score")
            frame.GetXaxis().SetNdivisions(505)
            frame.GetYaxis().SetTitleOffset(1.6)
            frame.GetYaxis().SetNdivisions(505)
            frame.SetTitle("BDT score comparison")
            frame.Draw()
            legend.Draw()
            CMSStyle.CMS_lumi(bdt_canvas, 5, 11)
            bdt_canvas.Update()
            
            tag = "mc" if isMC else "data"
            output_name = "bdt_score_comparison_multidatasets_%s_%s.png" % (tag, categ)
            bdt_canvas.SaveAs(output_name)
        
        
        
        
        
        def Compare_BDT_Scores_MultipleSelections(self, datafile, categ, selection_list, labels, isMC=False):
            frame = ROOT.RooRealVar('bdt_cv', 'BDT Score', -1, 1).frame()
            colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kOrange+7]
            
            # Map categ to tree name and label
            if categ == 'taue':
                treename = 'ztau3mutaue'
                cat_label = r"#tau_{e}"
            elif categ == 'taumu':
                treename = 'ztau3mutaumu'
                cat_label = r"#tau_{#mu}"
            elif categ == 'tauhA':
                treename = 'ztau3mutauh_A'
                cat_label = r"#tau_{h,1-prong}"
            elif categ == 'tauhB':
                treename = 'ztau3mutauh_B'
                cat_label = r"#tau_{h,3-prong}"
            elif categ == 'all':
                treename = 'ztautau'
                cat_label = "Inclusive"
            else:
                treename = categ
                cat_label = categ
                
            MiniTreeFile = ROOT.TFile.Open(datafile)
            tree = MiniTreeFile.Get(treename)
            
            bdt_cv = ROOT.RooRealVar("bdt_cv", "bdt_cv", -1.01, 1.01)
            tripletMass = ROOT.RooRealVar('tripletMass', '3#mu mass', 0, 100)
            isMCvar = ROOT.RooRealVar("isMC", "isMC", 0, 1000000)
            weight = ROOT.RooRealVar("weight", "weight", 0, 5)
            dimu_OS1 = ROOT.RooRealVar('dimu_OS1', 'dimu_OS1', 0, 1000)
            dimu_OS2 = ROOT.RooRealVar('dimu_OS2', 'dimu_OS2', 0, 1000)
            
            base_vars = ROOT.RooArgSet(tripletMass, bdt_cv, isMCvar, weight, dimu_OS1, dimu_OS2)
            
            legend = ROOT.TLegend(0.3, 0.3, 0.55, 0.55)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.SetHeader(cat_label, "C")
            
            for i, selection in enumerate(selection_list):
                selector_str = "isMC == 0" if not isMC else "( isMC == 211 || isMC == 210231 || isMC == 210232 || isMC == 210233 )"
                fullcut = "({}) && ({})".format(selector_str, selection)
                
                dataset_unweighted = ROOT.RooDataSet("sel_%d_raw" % i, "sel_%d_raw" % i, tree, base_vars, fullcut)
                
                entries = dataset_unweighted.numEntries()
                if entries == 0:
                    print("WARNING: No entries for selection:", selection)
                    continue
                    
                norm_factor = 1.0 / entries
                scale = ROOT.RooRealVar("scale", "scale", norm_factor)
                dataset_vars = dataset_unweighted.get()
                dataset_vars.add(scale)
                
                dataset = ROOT.RooDataSet("sel_%d" % i, "sel_%d" % i, dataset_unweighted, dataset_vars, "", "scale")
                
                dataset.plotOn(
                    frame,
                    ROOT.RooFit.Binning(100),
                    ROOT.RooFit.MarkerColor(colors[i % len(colors)]),
                    ROOT.RooFit.LineColor(colors[i % len(colors)]),
                    ROOT.RooFit.MarkerStyle(20 + i),
                    ROOT.RooFit.MarkerSize(0.75),
                    ROOT.RooFit.Name(labels[i])
                )
                legend.AddEntry(labels[i], labels[i], "lep")
                
            bdt_canvas = ROOT.TCanvas("bdt_canvas", "bdt_canvas", 800, 800)
            bdt_canvas.cd()
            bdt_canvas.SetLeftMargin( L/W )
            bdt_canvas.SetRightMargin( R/W )
            bdt_canvas.SetTopMargin( T/H )
            bdt_canvas.SetBottomMargin( B/H )
            frame.GetXaxis().SetTitle("BDT score")
            frame.GetXaxis().SetNdivisions(505)
            frame.GetYaxis().SetTitleOffset(1.6)
            frame.GetYaxis().SetNdivisions(505)
            frame.SetTitle("BDT score: Selection Comparisons")
            frame.Draw()
            legend.Draw()
            CMSStyle.CMS_lumi(bdt_canvas, 5, 11)
            bdt_canvas.Update()
        
            tag = "mc" if isMC else "data"
            output_name = "bdt_score_comparison_multiselections_%s_%s.png" % (tag, categ)
            ROOT.gPad.SaveAs(output_name)
        
        
        
        
        def Compare_BDT_Scores_MultipleDatasetsAndSelections(self, datafiles, datafile, categ, selection_list, labels, isMC=False):
            frame = ROOT.RooRealVar('bdt_cv', 'BDT Score', -1, 1).frame()
            colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kOrange+7]
            
            # Map categ to tree name and label
            if categ == 'taue':
                treename = 'ztau3mutaue'
                cat_label = r"#tau_{e}"
            elif categ == 'taumu':
                treename = 'ztau3mutaumu'
                cat_label = r"#tau_{#mu}"
            elif categ == 'tauhA':
                treename = 'ztau3mutauh_A'
                cat_label = r"#tau_{h,1-prong}"
            elif categ == 'tauhB':
                treename = 'ztau3mutauh_B'
                cat_label = r"#tau_{h,3-prong}"
            elif categ == 'all':
                treename = 'ztautau'
                cat_label = "Inclusive"
            else:
                treename = categ
                cat_label = categ
                
            legend = ROOT.TLegend(0.3, 0.3, 0.55, 0.55)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.SetHeader(cat_label, "C")
            
            # Common variables
            bdt_cv = ROOT.RooRealVar("bdt_cv", "bdt_cv", -1.01, 1.01)
            tripletMass = ROOT.RooRealVar('tripletMass', '3#mu mass', 0, 100)
            isMCvar = ROOT.RooRealVar("isMC", "isMC", 0, 1000000)
            weight = ROOT.RooRealVar("weight", "weight", 0, 5)
            dimu_OS1 = ROOT.RooRealVar('dimu_OS1', 'dimu_OS1', 0, 1000)
            dimu_OS2 = ROOT.RooRealVar('dimu_OS2', 'dimu_OS2', 0, 1000)
            
            base_vars = ROOT.RooArgSet(tripletMass, bdt_cv, isMCvar, weight, dimu_OS1, dimu_OS2)
            
            # -- First: compare multiple datasets (e.g., old/new)
            for i, datafile_in in enumerate(datafiles):
                MiniTreeFile = ROOT.TFile.Open(datafile_in)
                tree = MiniTreeFile.Get(treename)
                
                selector_str = "isMC == 0" if not isMC else "( isMC == 211 || isMC == 210231 || isMC == 210232 || isMC == 210233 )"
                dataset_unweighted = ROOT.RooDataSet("data_%d_raw" % i, "data_%d_raw" % i, tree, base_vars, selector_str)
                
                entries = dataset_unweighted.numEntries()
                if entries == 0:
                    print("WARNING: No entries in dataset:", datafile_in)
                    continue
                    
                norm_factor = 1.0 / entries
                scale = ROOT.RooRealVar("scale", "scale", norm_factor)
                dataset_vars = dataset_unweighted.get()
                dataset_vars.add(scale)
                
                dataset = ROOT.RooDataSet("scaled_ds_%d" % i, "scaled_ds_%d" % i, dataset_unweighted, dataset_vars, "", "scale")
                
                curve_name = "ds_curve_%d" % i
                dataset.plotOn(
                    frame,
                    ROOT.RooFit.Binning(100),
                    ROOT.RooFit.MarkerColor(colors[i % len(colors)]),
                    ROOT.RooFit.LineColor(colors[i % len(colors)]),
                    ROOT.RooFit.MarkerStyle(20 + i),
                    ROOT.RooFit.MarkerSize(0.75),
                    ROOT.RooFit.Name(curve_name)
                )
                legend.AddEntry(curve_name, "Dataset %s" % labels[i], "lep")
                
            # -- Second: compare multiple selections on a single dataset
            MiniTreeFile = ROOT.TFile.Open(datafile)
            tree = MiniTreeFile.Get(treename)
            
            for j, selection in enumerate(selection_list):
                selector_str = "isMC == 0" if not isMC else "( isMC == 211 || isMC == 210231 || isMC == 210232 || isMC == 210233 )"
                fullcut = "({}) && ({})".format(selector_str, selection)
                
                dataset_unweighted = ROOT.RooDataSet("sel_%d_raw" % j, "sel_%d_raw" % j, tree, base_vars, fullcut)
                
                entries = dataset_unweighted.numEntries()
                if entries == 0:
                    print("WARNING: No entries for selection:", selection)
                    continue
                    
                norm_factor = 1.0 / entries
                scale = ROOT.RooRealVar("scale", "scale", norm_factor)
                dataset_vars = dataset_unweighted.get()
                dataset_vars.add(scale)
                
                dataset = ROOT.RooDataSet("sel_%d" % j, "sel_%d" % j, dataset_unweighted, dataset_vars, "", "scale")
                
                curve_name = "sel_curve_%d" % j
                dataset.plotOn(
                    frame,
                    ROOT.RooFit.Binning(100),
                    ROOT.RooFit.MarkerColor(colors[(j + len(datafiles)) % len(colors)]),
                    ROOT.RooFit.LineColor(colors[(j + len(datafiles)) % len(colors)]),
                    ROOT.RooFit.MarkerStyle(24 + j),
                    ROOT.RooFit.MarkerSize(0.75),
                    ROOT.RooFit.Name(curve_name)
                )
                legend.AddEntry(curve_name, "Sel %s" % labels[j], "lep")
                
            bdt_canvas = ROOT.TCanvas("bdt_canvas", "bdt_canvas", 800, 800)
            bdt_canvas.cd()
            bdt_canvas.SetLeftMargin( L/W )
            bdt_canvas.SetRightMargin( R/W )
            bdt_canvas.SetTopMargin( T/H )
            bdt_canvas.SetBottomMargin( B/H )
            frame.GetXaxis().SetTitle("BDT score")
            frame.GetXaxis().SetNdivisions(505)
            frame.GetYaxis().SetTitleOffset(1.6)
            frame.GetYaxis().SetNdivisions(505)
            frame.SetTitle("BDT score: Combined Comparison")
            frame.Draw()
            legend.Draw()
            CMSStyle.CMS_lumi(bdt_canvas, 5, 11)
            bdt_canvas.Update()
            
            
        
            tag = "mc" if isMC else "data"
            output_name = "bdt_score_comparison_combined_%s_%s.png" % (tag, categ)
            ROOT.gPad.SaveAs(output_name)
            
        
        
        
        
        
        
        
        def Plot_Bdt_Symmetry(self, datafile, categ, isMC=False):
                
                Loose_BDT_min = -0.9
                Loose_BDT_max = 0.1
                Loose_BDT_N = 12
                Loose_BDT_step = (Loose_BDT_max - Loose_BDT_min) / Loose_BDT_N
                Loose_BDT_list = []
                Ratio_list = []
                
                
                fit_range_lo = 1.6
                fit_range_hi = 2.0
                
                signal_range_lo = 1.74
                signal_range_hi = 1.81
                
                MiniTreeFile = ROOT.TFile.Open(datafile)
                MiniTreeFile.cd()
                
                treeName=''
                signalnorm = 1.0
                cat_label = ""
                if(categ == 'taue'):
                        treeName  = 'ztau3mutaue'
                        signalnorm = 0.00000856928
                        cat_label = r"$\tau_{e}$"
                if(categ == 'taumu'):
                        treeName = 'ztau3mutaumu'
                        signalnorm = 0.00000822810
                        cat_label = r"$\tau_{\mu}$"
                if(categ == 'tauhA'):
                        treeName = 'ztau3mutauh_A'
                        signalnorm = 0.00000815958
                        cat_label = r"$\tau_{h,1-prong}$"
                if(categ == 'tauhB'):
                        treeName = 'ztau3mutauh_B'
                        signalnorm = 0.00000815958
                        cat_label = r"$\tau_{h,3-prong}$"
                if(categ == 'all'):
                        treeName   = 'ztautau'
                        signalnorm = 0.00000824176
                        cat_label = "Inclusive"
                
                tree = MiniTreeFile.Get(treeName)
                
                tripletMass          = ROOT.RooRealVar('tripletMass'                , '3#mu mass'           , fit_range_lo, fit_range_hi, 'GeV')
                bdt_cv               = ROOT.RooRealVar('bdt_cv'                     , 'bdt_cv'              , -1 , 1)
                dimu_OS1             = ROOT.RooRealVar('dimu_OS1'                   , 'dimu_OS1'            ,  0 , 10)
                dimu_OS2             = ROOT.RooRealVar('dimu_OS2'                   , 'dimu_OS2'            ,  0 , 10)
                event_weight         = ROOT.RooRealVar('weight'                     , 'event_weight'        ,  0,  5)  # this weight includes also the scale  mc signal scale
                category             = ROOT.RooRealVar('category'                   , 'category'            ,  0,  5)
                isMC                 = ROOT.RooRealVar('isMC'                       , 'isMC'                ,  0,  1000000)
                scale                = ROOT.RooRealVar('scale'                      , 'scale'               ,  signalnorm)  
                
                
                
                variables = ROOT.RooArgSet()
                variables.add(tripletMass)
                variables.add(bdt_cv)
                variables.add(dimu_OS1)
                variables.add(dimu_OS2)
                variables.add(event_weight)
                variables.add(category)
                variables.add(isMC)
                
                phivetoes="(fabs(dimu_OS1 - 1.020)>0.020)&(fabs(dimu_OS2 - 1.020)>0.020)&"
                omegavetoes="fabs(dimu_OS1 - 0.782)>0.020&fabs(dimu_OS2 - 0.782)>0.020&"
                
                graph_ratio_vs_bdt = ROOT.TGraphErrors()
                graph_ratio_vs_bdt.SetName("graph_ratio_vs_bdt")
                graph_ratio_vs_bdt.SetTitle("Sideband Symmetry vs BDT Cut;BDT Cut;Sideband Ratio")
                
                
                for i in range(Loose_BDT_N):
                        bdt_cut = Loose_BDT_min + i * Loose_BDT_step + Loose_BDT_step / 2.0
                        
                        # For fitting BDT Output in Data
                        
                        BDT_Score_Min=-0.3
                        
                        BlindDataSelector = RooFormulaVar('DataSelector', 'DataSelector',' bdt_cv > ' + str(bdt_cut)+' & '+ phivetoes+' isMC == 0 & (tripletMass<=%s || tripletMass>=%s) & (tripletMass>=%s & tripletMass<=%s) ' %(signal_range_lo,signal_range_hi,fit_range_lo,fit_range_hi) , RooArgList(variables))
                        
                        fulldata = RooDataSet('fulldata', 'fulldata', tree,  variables, BlindDataSelector)
                        
                        tripletMass.setRange('left' , fit_range_lo    , signal_range_lo)
                        tripletMass.setRange('right', signal_range_hi , fit_range_hi)
                        tripletMass.setRange('full' , fit_range_lo    , fit_range_hi)
                        
                        tripletMass.setRange("SB1",fit_range_lo,1.75)
                        tripletMass.setRange("SB2",1.80,fit_range_hi)
                        tripletMass.setRange("fullRange",fit_range_lo,fit_range_hi)
                        tripletMass.setRange("SIG",signal_range_lo,signal_range_hi)
                        
                        nbkg = ROOT.RooRealVar('nbkg', 'nbkg', 1000, 0, 500000)
                        slope = ROOT.RooRealVar('slope', 'slope', 1.0, -100, 100)
                        expo = ROOT.RooExponential('bkg_expo', 'bkg_expo', tripletMass, slope)
                        pdfmodel = ROOT.RooAddPdf('bkg_extended_expo', 'bkg_extended_expo', ROOT.RooArgList(expo), ROOT.RooArgList(nbkg))
                        results_pdf = pdfmodel.fitTo(fulldata, ROOT.RooFit.Range('left,right'), ROOT.RooFit.Save())
                        
                        
                        SG_abs = pdfmodel.createIntegral(ROOT.RooArgSet(tripletMass), ROOT.RooFit.Range("SIG"))
                        SB_abs = pdfmodel.createIntegral(ROOT.RooArgSet(tripletMass), ROOT.RooFit.Range("left,right"))
                        
                        # Get their values
                        SG_integral = SG_abs.getVal()
                        SB_integral = SB_abs.getVal()
                        Ratio_val = SG_integral / SB_integral if SB_integral > 0 else 0
                        
                        # Compute propagated errors from fit result
                        SG_error = SG_abs.getPropagatedError(results_pdf)
                        SB_error = SB_abs.getPropagatedError(results_pdf)
                        
                        # Propagate the error on the ratio manually
                        if SG_integral > 0 and SB_integral > 0:
                            Ratio_error = Ratio_val * ((SG_error / SG_integral)**2 + (SB_error / SB_integral)**2)**0.5
                        else:
                            Ratio_error = 0.0
                            
                        print(f"[{bdt_cut:.3f}] SG: {SG_integral:.5g} +/- {SG_error:.2g}, SB: {SB_integral:.5g} +/- {SB_error:.2g}")
                        
                        
                        
                        point_index = graph_ratio_vs_bdt.GetN()  # Current point index
                        graph_ratio_vs_bdt.SetPoint(point_index, bdt_cut, Ratio_val)
                        graph_ratio_vs_bdt.SetPointError(point_index, 0.0, Ratio_error)
                        
                c = ROOT.TCanvas("c", "BDT vs SB Symmetry", 800, 600)
                graph_ratio_vs_bdt.SetMarkerStyle(20)
                graph_ratio_vs_bdt.SetMarkerSize(1.0)
                graph_ratio_vs_bdt.SetLineWidth(2)
                graph_ratio_vs_bdt.Draw("AP")  # A = axis, P = points with errors
                
                c.SetGrid()
                
                output_name = "Ratio_vs_BDTCut_%s.png" % (categ)
                
                c.SaveAs(output_name)
                #c.SaveAs("Ratio_vs_BDTCut_%s.pdf")
                
        
        
        
        
        def get_signal_window(self, datafile, categ, isMC=False):

                    
                fit_range_lo = 1.6
                fit_range_hi = 2.0
                
                signal_range_lo = 1.74
                signal_range_hi = 1.81
                
                MiniTreeFile = ROOT.TFile.Open(datafile)
                MiniTreeFile.cd()
                
                treeName=''
                signalnorm = 1.0
                cat_label = ""
                bdt_cut = 0.1 #A stringent enough cut
                if(categ == 'taue'):
                        treeName  = 'ztau3mutaue'
                        signalnorm = 0.00000856928
                        cat_label = r"$\tau_{e}$"
                        bdt_cut = 0.28
                if(categ == 'taumu'):
                        treeName = 'ztau3mutaumu'
                        signalnorm = 0.00000822810
                        cat_label = r"$\tau_{\mu}$"
                        bdt_cut = 0.44
                if(categ == 'tauhA'):
                        treeName = 'ztau3mutauh_A'
                        signalnorm = 0.00000815958
                        cat_label = r"$\tau_{h,1-prong}$"
                        bdt_cut = 0.38
                if(categ == 'tauhB'):
                        treeName = 'ztau3mutauh_B'
                        signalnorm = 0.00000815958
                        cat_label = r"$\tau_{h,3-prong}$"
                        bdt_cut = 0.15
                if(categ == 'all'):
                        treeName   = 'ztautau'
                        signalnorm = 0.00000824176
                        cat_label = "Inclusive"
                        bdt_cut = 0.59
                
                tree = MiniTreeFile.Get(treeName)
                
                tripletMass          = ROOT.RooRealVar('tripletMass'                , '3#mu mass'           , fit_range_lo, fit_range_hi, 'GeV')
                bdt_cv               = ROOT.RooRealVar('bdt_cv'                     , 'bdt_cv'              , -1 , 1)
                dimu_OS1             = ROOT.RooRealVar('dimu_OS1'                   , 'dimu_OS1'            ,  0 , 2)
                dimu_OS2             = ROOT.RooRealVar('dimu_OS2'                   , 'dimu_OS2'            ,  0 , 2)
                event_weight         = ROOT.RooRealVar('weight'                     , 'event_weight'        ,  0,  5)  # this weight includes also the scale  mc signal scale
                category             = ROOT.RooRealVar('category'                   , 'category'            ,  0,  5)
                isMC                 = ROOT.RooRealVar('isMC'                       , 'isMC'                ,  0,  1000000)
                scale                = ROOT.RooRealVar('scale'                      , 'scale'               ,  signalnorm)  
                
                
                
                variables = ROOT.RooArgSet()
                variables.add(tripletMass)
                variables.add(bdt_cv)
                variables.add(dimu_OS1)
                variables.add(dimu_OS2)
                variables.add(event_weight)
                variables.add(category)
                variables.add(isMC)
                
                #Define mc and data which can be reduced later
                
                phivetoes="(fabs(dimu_OS1 - 1.020)>0.020)&(fabs(dimu_OS2 - 1.020)>0.020)&"
                omegavetoes="fabs(dimu_OS1 - 0.782)>0.020&fabs(dimu_OS2 - 0.782)>0.020&"
                
                MCSelector = RooFormulaVar('MCSelector', 'MCSelector', ' bdt_cv > ' + str(bdt_cut)+' & '+ phivetoes+' isMC !=0 & (isMC == 211 | isMC == 210231 | isMC == 210232 | isMC == 210233 ) & (tripletMass>=%s & tripletMass<=%s) ' %(signal_range_lo,signal_range_hi) , RooArgList(variables))
                
                fullmc_unweighted = RooDataSet('mc', 'mc', tree_norm, variables, MCSelector)
                dataset_vars = fullmc_unweighted.get()
                dataset_vars.add(scale)
                
                fullmc = RooDataSet('mc', 'mc', fullmc_unweighted, dataset_vars, "",'scale')
                
                
                BlindDataSelector = RooFormulaVar('DataSelector', 'DataSelector',' bdt_cv > ' + str(bdt_cut)+' & '+ phivetoes+' isMC == 0 & (tripletMass>=%s & tripletMass<=%s) ' %(fit_range_lo,fit_range_hi) , RooArgList(variables))
                        
                fulldata = RooDataSet('fulldata', 'fulldata', tree,  variables, BlindDataSelector)
                
                tripletMass.setRange('left' , fit_range_lo    , signal_range_lo)
                tripletMass.setRange('right', signal_range_hi , fit_range_hi)
                tripletMass.setRange('full' , fit_range_lo    , fit_range_hi)
                
                tripletMass.setRange("SB1",fit_range_lo,1.75)
                tripletMass.setRange("SB2",1.80,fit_range_hi)
                tripletMass.setRange("fullRange",fit_range_lo,fit_range_hi)
                tripletMass.setRange("SIG",signal_range_lo,signal_range_hi)
                
                mean = ROOT.RooRealVar("mean", "mean", 1.776, 0., 5.)
                sigma = ROOT.RooRealVar("sigma", "sigma", 0.5, 0.001, 10)
                Gauss = ROOT.RooGaussian("Gauss", "Gauss dist", InvMass, mean, sigma)
                GaussNorm = ROOT.RooRealVar("GaussNorm", "GaussNorm", 0.5, 0.001, 1.0)
                mc_pdf = ROOT.RooAddPdf("mc_pdf", "mc_pdf", ROOT.RooArgList(Gauss), ROOT.RooArgList(GaussNorm))
                mc_fitresult = mc_pdf.fitTo(fullmc, ROOT.RooFit.Range("R3"), ROOT.RooFit.Save())
                
                Gaussian_Sigma_From_Loose_BDT_Cut = sigma.getVal()
                
                sigma_scale = 0.0
                signal_peak_region_min = 0.0
                signal_peak_region_max = 0.0
                
                # Range of sigma factors
                Xa_min = 1.5
                Xa_max = 3.5
                N_a = 20
                
                triplet_mass_bins = 40
                
                for m in range(N_a):
                        step = (Xa_max - Xa_min) / N_a
                        sigma_scale = Xa_min + m * step
                
                        print(f"sigma_scale is: {sigma_scale}")
                
                        signal_peak_region_min = 1.77686 - sigma_scale * Gaussian_Sigma_From_Loose_BDT_Cut
                        signal_peak_region_max = 1.77686 + sigma_scale * Gaussian_Sigma_From_Loose_BDT_Cut
                
                        val_1 = int(((signal_peak_region_min - fit_range_lo) /
                                    ((fit_range_hi - fit_range_lo) / triplet_mass_bins)) + 0.5)
                        val_2 = int(((signal_peak_region_max - fit_range_lo) /
                                    ((fit_range_hi - fit_range_lo) / triplet_mass_bins)) + 0.5)
                
                        signal_peak_region_min = fit_range_lo + val_1 * ((fit_range_hi - fit_range_lo) / triplet_mass_bins)
                        signal_peak_region_max = fit_range_lo + val_2 * ((fit_range_hi - fit_range_lo) / triplet_mass_bins)
                
                        print(f"min: {signal_peak_region_min} max: {signal_peak_region_max}")
                        
                        MC_inside_peak = '(tripletMass>=%s & tripletMass<=%s)' %(signal_peak_region_min,signal_peak_region_max)
                        Data_outside_peak = '(tripletMass<=%s || tripletMass>=%s)' %(signal_peak_region_min,signal_peak_region_max)
                        
                        MC_reduced_dataset = fullmc.reduce(MC_inside_peak)
                        Data_reduced_dataset = fulldata.reduce(Data_outside_peak)
                        
                        
                        
                        nbkg = ROOT.RooRealVar('nbkg', 'nbkg', 1000, 0, 500000)
                        slope = ROOT.RooRealVar('slope', 'slope', 1.0, -100, 100)
                        expo = ROOT.RooExponential('bkg_expo', 'bkg_expo', tripletMass, slope)
                        pdfmodel = ROOT.RooAddPdf('bkg_extended_expo', 'bkg_extended_expo', ROOT.RooArgList(expo), ROOT.RooArgList(nbkg))
                        results_pdf = pdfmodel.fitTo(Data_reduced_dataset, ROOT.RooFit.Range('left,right'), ROOT.RooFit.Save())
                        
                        SG_integral = pdfmodel.createIntegral(ROOT.RooArgSet(tripletMass), ROOT.RooArgSet(tripletMass), "SIG").getVal()
                        SB_integral = pdfmodel.createIntegral(ROOT.RooArgSet(tripletMass), ROOT.RooArgSet(tripletMass), "left,right").getVal()
                        
                        N_s_1 = 0.0
                        N_b_1 = 0.0
                        
                        N_s_1 = MC_inside_peak.sumEntries()
                        N_b_1 = nbkg.getVal()*SG_integral
                
                        print(f"signal: {N_s_1} bkg: {N_b_1}")
                
                        if round(N_b_1 / 0.00001) * 0.00001 > 0.0:
                            cmd = f"python card_modifiers/card_mod.py --luminosity 59.0 --s {N_s_1} --b {N_b_1}"
                            os.system(cmd)
                
                            cmd_bayesian = f"combine -M BayesianSimple modified_simplest_card.txt --cl 0.9 -t 100 > Test_Bayesian_Output.txt"
                            os.system(cmd_bayesian)
                
                            with open(f"Test_Bayesian_Output.txt") as f:
                                for line in f:
                                    if "median expected limit" in line:
                                        tokens = line.split()
                                        print(f"Limit Bayesian: {tokens[5]}")
                
                            cmd_ul = f"python3 ../../workdirDataWithMCSkimmed_SeparateJul_08_2024/Projections/CLs_UL_Calculator.py {N_s_1} {N_b_1} > out_UL_Calc.txt"
                            os.system(cmd_ul)
                
                            with open(f"out_UL_Calc.txt") as f:
                                for line in f:
                                    if "upper limit" in line:
                                        print(f"Limit UL_Calc: {line.split()[-1]}")




if __name__ == "__main__":
    
        ROOT.gROOT.SetBatch(True)
        
        datafile = "../../Combine_Tree_ztau3mutau_PF_PostBDT.root"
        
        # Assume these methods belong to an object or class, e.g. BDTPlotter()
        BDTPlotter = BDT_Shape_Comparisons()  # Replace this with your actual class name or self instance
        
        # Categories to loop over
        #categories = ['taue', 'taumu', 'tauhA', 'tauhB','all']
        categories = ['taue']
        
        # Example BDT selections you want to compare
        bdt_selections = [
            "dimu_OS2 > 0.1",
            "dimu_OS2 > 0.5",
            "dimu_OS2 > 1.0",
            "dimu_OS2 > 1.5"
        ]
        selection_labels = ["dimu cut1", "dimu cut2", "dimu cut3", "dimu cut4"]
        
        # Sample datasets (for Compare_BDT_Scores_MultipleDatasets)
        dataset_files = [
            "../../Combine_Tree_ztau3mutau_orig_PostBDT.root",
            #"../../Combine_Tree_ztau3mutau_PFGL_PostBDT.root",
            "../../Combine_Tree_ztau3mutau_PF_PostBDT.root",
        ]
        
        # Loop over categories
        for categ in categories:
            print("Running for category: " + str(categ))
        
            # Construct correct tree name for this category
            tree = ''
            if categ == 'taue':
                tree = 'ztau3mutaue'
            elif categ == 'taumu':
                tree = 'ztau3mutaumu'
            elif categ == 'tauhA':
                tree = 'ztau3mutauh_A'
            elif categ == 'tauhB':
                tree = 'ztau3mutauh_B'
            elif categ == 'all':
                tree = 'ztautau'
        
            # 1. Compare different datasets (same selection)
            #BDTPlotter.Compare_BDT_Scores_MultipleDatasets(dataset_files, categ, True)
            
            # 2. Compare different selections on the same file. Pick any one file.
            #BDTPlotter.Compare_BDT_Scores_MultipleSelections(datafile, categ, bdt_selections, selection_labels, isMC=True)
            
            # 3. Full comparison: different files + different selections
            #BDTPlotter.Compare_BDT_Scores_MultipleDatasetsAndSelections(dataset_files, datafile, categ, bdt_selections, selection_labels, isMC=True)
            
            # 4. Constancy of 'r' plots
            #BDTPlotter.Plot_Bdt_Symmetry(datafile, categ, False)
            
            # 5. Get signal peak width
            BDTPlotter.get_signal_window(datafile, categ, False)