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
import array


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
H_ref = 800; 
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
        
        
        
        def Compare_BDT_Scores_MultipleDatasets(self, datafiles, categ, WhetherMC=False):
            #frame = ROOT.RooRealVar('bdt_cv', 'BDT Score', -1, 1).frame()
            colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta]
            
            labels = ["Relaxed cuts", "All cuts"]
            
            fit_range_lo = 1.6
            fit_range_hi = 2.0
            
            signal_range_lo = 1.74
            signal_range_hi = 1.81
            
            whichToFit = "Combine_Tree_ztau3mutau_PF_PostBDT_symmetric.root"
            which_col = ROOT.kRed
            
            # Map categ to treename and LaTeX-style label
            if categ == 'taue':
                treename = 'ztau3mutaue'
                cat_label = r"#tau_{e}"
                whichToFit = "Combine_Tree_ztau3mutau_orig_PostBDT.root"
                which_col = ROOT.kBlack
            elif categ == 'taumu':
                treename = 'ztau3mutaumu'
                cat_label = r"#tau_{#mu}"
                whichToFit = "Combine_Tree_ztau3mutau_PF_PostBDT_symmetric.root"
                which_col = ROOT.kRed
            elif categ == 'tauhA':
                treename = 'ztau3mutauh_A'
                cat_label = r"#tau_{h,1-prong}"
                whichToFit = "Combine_Tree_ztau3mutau_orig_PostBDT.root"
                which_col = ROOT.kBlack
            elif categ == 'tauhB':
                treename = 'ztau3mutauh_B'
                cat_label = r"#tau_{h,3-prong}"
                whichToFit = "Combine_Tree_ztau3mutau_orig_PostBDT.root"
                which_col = ROOT.kBlack
            elif categ == 'all':
                treename = 'ztautau'
                cat_label = "Inclusive"
                whichToFit = "Combine_Tree_ztau3mutau_PF_PostBDT_symmetric.root"
                which_col = ROOT.kRed
            else:
                treename = categ
                cat_label = categ
        
            if WhetherMC:
                    legend = ROOT.TLegend(0.3, 0.3, 0.55, 0.55)  # Left half
            else:
                    #legend = ROOT.TLegend(0.65, 0.3, 0.9, 0.55)  # Right half
                    legend = ROOT.TLegend(0.3, 0.3, 0.55, 0.55)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            # legend.SetHeader(cat_label, "C")
            
            # Define variables
            bdt_cv = ROOT.RooRealVar("bdt_cv", "bdt_cv", -1.0, 1.0)
            tripletMass = ROOT.RooRealVar('tripletMass', '3#mu mass', 0, 100)
            isMC = ROOT.RooRealVar("isMC", "isMC", 0, 1000000)
            weight = ROOT.RooRealVar("weight", "weight", 0, 5)
            dimu_OS1 = ROOT.RooRealVar('dimu_OS1', 'dimu_OS1', 0, 1000)
            dimu_OS2 = ROOT.RooRealVar('dimu_OS2', 'dimu_OS2', 0, 1000)
            
            frame = bdt_cv.frame()
            
            for i, datafile in enumerate(datafiles):
                    MiniTreeFile = ROOT.TFile.Open(datafile)
                    tree = MiniTreeFile.Get(treename)
                
                    variables = ROOT.RooArgSet(tripletMass, bdt_cv, isMC, weight, dimu_OS1, dimu_OS2)
                
                    phivetoes = "(fabs(dimu_OS1 - 1.020)>0.020)&&(fabs(dimu_OS2 - 1.020)>0.020)&&"
                
                    if not WhetherMC:
                        selector = ROOT.RooFormulaVar(
                            'DataSelector', 'DataSelector',
                            phivetoes + f"isMC == 0 && (tripletMass<={signal_range_lo} || tripletMass>={signal_range_hi}) && (tripletMass>={fit_range_lo} && tripletMass<={fit_range_hi})",
                            ROOT.RooArgList(variables)
                        )
                    else:
                        selector = ROOT.RooFormulaVar(
                            'MCSelector', 'MCSelector',
                            phivetoes + f"isMC != 0 && (isMC == 211 || isMC == 210231 || isMC == 210232 || isMC == 210233) && (tripletMass>={signal_range_lo} && tripletMass<={signal_range_hi})",
                            ROOT.RooArgList(variables)
                        )
                
                    dataset_unweighted = ROOT.RooDataSet(f"data_{i}", f"data_{i}", tree, variables, selector)
                
                    norm_factor = 1.0 / dataset_unweighted.numEntries()
                    print("dataset entries:", dataset_unweighted.numEntries())
                    scale = ROOT.RooRealVar("scale", "scale", norm_factor)
                    dataset_vars = dataset_unweighted.get()
                    dataset_vars.add(scale)
                
                    dataset = ROOT.RooDataSet("scaled_ds", "scaled_ds", dataset_unweighted, dataset_vars, "", "scale")
                    #dataset = dataset_unweighted
                
                    curve_name = f"curve_{i}"
                    dataset.plotOn(
                        frame,
                        ROOT.RooFit.Binning(100),
                        ROOT.RooFit.MarkerColor(colors[i % len(colors)]),
                        ROOT.RooFit.LineColor(colors[i % len(colors)]),
                        ROOT.RooFit.MarkerStyle(20 + i),
                        ROOT.RooFit.MarkerSize(0.75),
                        ROOT.RooFit.Name(curve_name)
                    )
                    curve = frame.findObject(curve_name)
                    legend.AddEntry(curve, labels[i], "lep")
                
                    # Only fit data from Combine_Tree_ztau3mutau_PF_PostBDT.root
                    if not WhetherMC and os.path.basename(datafile) == whichToFit:
                    #if True:
                        BDT_Score_Min = -0.3
                        bdt_cv.setRange("BDT_Fit_Range", BDT_Score_Min, 1.0)
                
                        a = ROOT.RooRealVar("a", "a", 1.0, 0.0, 10.0)
                        b = ROOT.RooRealVar("b", "b", 1.0, -10.0, 10.0)
                        c = ROOT.RooRealVar("c", "c", 1.0, -10.0, 10.0)
                        d = ROOT.RooRealVar("d", "d", 1.0, -10.0, 10.0)
                
                        quadratic = ROOT.RooFormulaVar("quadratic", "a + b*bdt_cv + c*bdt_cv*bdt_cv + d*bdt_cv*bdt_cv*bdt_cv",
                                                       ROOT.RooArgList(a, b, c, d, bdt_cv))
                        expModel = ROOT.RooGenericPdf("expModel", "exp(quadratic)", ROOT.RooArgList(quadratic))
                
                        BDTNorm = ROOT.RooRealVar("BDTNorm", "BDTNorm", 500.0, 0.1, 1000000000)
                        BDT_distribution = ROOT.RooAddPdf("BDT_distribution", "BDT_distribution",
                                                          ROOT.RooArgList(expModel), ROOT.RooArgList(BDTNorm))
                
                        results_pdf = BDT_distribution.fitTo(dataset_unweighted, ROOT.RooFit.Range("BDT_Fit_Range"), ROOT.RooFit.Save())
                        
                        BDT_distribution.plotOn(
                            frame,
                            ROOT.RooFit.Normalization(dataset.sumEntries("", "BDT_Fit_Range"), ROOT.RooAbsReal.NumEvent),
                            ROOT.RooFit.LineColor(which_col),
                            ROOT.RooFit.LineStyle(2),
                            ROOT.RooFit.Name("fit_curve")
                        )
                        legend.AddEntry(frame.findObject("fit_curve"), "Fit", "l")
            
            
            
            #Plotting
            bdt_canvas = ROOT.TCanvas("bdt_canvas", "bdt_canvas", 800, 800)
            bdt_canvas.cd()
            bdt_canvas.SetLeftMargin( L/W )
            bdt_canvas.SetRightMargin( R/W )
            bdt_canvas.SetTopMargin( T/H )
            bdt_canvas.SetBottomMargin( B/H )
            
            if not WhetherMC:
                    bdt_canvas.SetLogy()
                    frame.SetMinimum(1e-5)
            frame.GetXaxis().SetTitle("BDT score")
            frame.GetXaxis().SetNdivisions(505)
            frame.GetYaxis().SetTitleOffset(1.6)
            frame.GetYaxis().SetNdivisions(505)
            frame.SetTitle("BDT score comparison")
            frame.Draw('sameaxis')
            legend.Draw()
            if WhetherMC:
                    CMSStyle.extraText = "Work in progress"
                    CMSStyle.CMS_lumi(bdt_canvas, 5, 11)
            else:
                    CMSStyle.extraText = "        Work in progress"
                    CMSStyle.CMS_lumi(bdt_canvas, 5, 0)
            bdt_canvas.Update()
            
            tag = "mc" if WhetherMC else "data"
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
        
        
        
        
        def Compare_BDT_Scores_MultipleDatasetsAndSelections(self, datafiles, datafile, categ, selection_list, labels, WhetherMC=False):
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
                
                selector_str = "isMC == 0" if not WhetherMC else "( isMC == 211 || isMC == 210231 || isMC == 210232 || isMC == 210233 )"
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
                selector_str = "isMC == 0" if not WhetherMC else "( isMC == 211 || isMC == 210231 || isMC == 210232 || isMC == 210233 )"
                fullcut = "({}) ".format(selection)
                
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
            
        
        def Compare_BDT_Scores_MultipleDatasetsAndSelectionsZTTMass(self, datafiles, datafile, categ, selection_list, labels, WhetherMC=True):
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
            #legend.SetHeader(cat_label, "C")
            
            # Common variables
            bdt_cv = ROOT.RooRealVar("bdt_cv", "bdt_cv", -1.01, 1.01)
            tripletMass = ROOT.RooRealVar('tripletMass', '3#mu mass', 0, 100)
            isMC = ROOT.RooRealVar("isMC", "isMC", 0, 1000000)
            weight = ROOT.RooRealVar("weight", "weight", 0, 5)
            dimu_OS1 = ROOT.RooRealVar('dimu_OS1', 'dimu_OS1', 0, 1000)
            dimu_OS2 = ROOT.RooRealVar('dimu_OS2', 'dimu_OS2', 0, 1000)
            
            base_vars = ROOT.RooArgSet(tripletMass, bdt_cv, isMC, weight, dimu_OS1, dimu_OS2)
            
            # -- First: compare multiple datasets (e.g., old/new)
            for i, datafile_in in enumerate(datafiles):
                MiniTreeFile = ROOT.TFile.Open(datafile_in)
                tree = MiniTreeFile.Get(treename)
                
                selector_str = "isMC == 0" if not WhetherMC else "( isMC == 211 || isMC == 210231 || isMC == 210232 || isMC == 210233 )"
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
                legend.AddEntry(frame.findObject(curve_name), r"m_{#tau} = 1.777", "lep")
                
            # -- Second: compare multiple selections on a single dataset
            MiniTreeFile = ROOT.TFile.Open(datafile)
            tree = MiniTreeFile.Get(treename)
            
            for j, selection in enumerate(selection_list):
                selector_str = "isMC == 0" if not WhetherMC else "( isMC == 211 || isMC == 210231 || isMC == 210232 || isMC == 210233 )"
                fullcut = "({}) ".format(selection)
                
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
                legend.AddEntry(frame.findObject(curve_name), labels[j], "lep")
                
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
                
                
                fit_range_lo = 1.55
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
                        
                        SG_abs = pdfmodel.createIntegral(ROOT.RooArgSet(tripletMass), ROOT.RooArgSet(tripletMass), "SIG")
                        SB_abs = pdfmodel.createIntegral(ROOT.RooArgSet(tripletMass), ROOT.RooArgSet(tripletMass), "left,right")
                        
                        SG_abs_withNorm = ROOT.RooProduct("SG_abs_withNorm", "SG_abs_withNorm",ROOT.RooArgList(nbkg.getVal(), SG_abs))
                        SB_abs_withNorm = ROOT.RooProduct("SB_abs_withNorm", "SB_abs_withNorm",ROOT.RooArgList(nbkg.getVal(), SB_abs))
                        
                        # Get their values
                        SG_integral = SG_abs_withNorm.getVal()
                        SB_integral = SB_abs_withNorm.getVal()
                        
                        if SB_integral > 0:
                                inputs = ROOT.RooArgList(SG_abs_withNorm, SB_abs_withNorm)
                                ratio_formula = ROOT.RooFormulaVar("ratio", "(@0)/(@1)", inputs)
                                Ratio_val = ratio_formula.getVal()
                        else:
                                Ratio_val = 0
                        
                        # Propagate the error on the ratio
                        if SG_integral > 0 and SB_integral > 0:
                            #Ratio_error = ((SG_error / SB_integral)**2 + (SG_integral * SB_error / SB_integral**2)**2)**0.5
                            Ratio_error = ratio_formula.getPropagatedError(results_pdf)
                        else:
                            Ratio_error = 0.0
                            
                        #print(f"[{bdt_cut:.3f}] SG: {SG_integral:.5g} +/- {SG_error:.2g}, SB: {SB_integral:.5g} +/- {SB_error:.2g}, nbk: {nbkg.getVal():.5g}")
                        print(f"[{bdt_cut:.3f}] SG: {SG_integral:.5g}, SB: {SB_integral:.5g}, nbk: {nbkg.getVal():.5g}, ratio: {Ratio_val:.5g}, ratio error: {Ratio_error:.5g}")
                        
                        
                        
                        point_index = graph_ratio_vs_bdt.GetN()  # Current point index
                        graph_ratio_vs_bdt.SetPoint(point_index, bdt_cut, Ratio_val)
                        graph_ratio_vs_bdt.SetPointError(point_index, 0.0, Ratio_error)
                        
                
                
                
                bdt_canvas = ROOT.TCanvas("bdt_canvas", "BDT vs SB Symmetry", 800, 800)
                bdt_canvas.cd()
                bdt_canvas.SetLeftMargin( L/W )
                bdt_canvas.SetRightMargin( R/W )
                bdt_canvas.SetTopMargin( T/H )
                bdt_canvas.SetBottomMargin( B/H )
                graph_ratio_vs_bdt.SetMarkerStyle(20)
                graph_ratio_vs_bdt.SetMarkerSize(1.0)
                graph_ratio_vs_bdt.SetLineWidth(2)
                
                # Set the y-axis range
                graph_ratio_vs_bdt.SetMinimum(0.1)
                graph_ratio_vs_bdt.SetMaximum(0.3)
                
                graph_ratio_vs_bdt.Draw("AP")  # A = axis, P = points with errors
                
                graph_ratio_vs_bdt.GetXaxis().SetTitle("BDT cut")
                graph_ratio_vs_bdt.GetXaxis().SetNdivisions(505)
                #graph_ratio_vs_bdt.GetXaxis().SetTitleSize(0.045)
                graph_ratio_vs_bdt.GetYaxis().SetTitle(r"B_{S}/B_{SB}")
                graph_ratio_vs_bdt.GetYaxis().SetTitleOffset(1.6)
                graph_ratio_vs_bdt.GetYaxis().SetNdivisions(505)
                #graph_ratio_vs_bdt.GetYaxis().SetTitleSize(0.045)
                
                bdt_canvas.SetGrid()
                
                CMSStyle.CMS_lumi(bdt_canvas, 5, 11)
                bdt_canvas.Update()
                
                output_name = "Ratio_vs_BDTCut_%s.png" % (categ)
                
                bdt_canvas.SaveAs(output_name)
                #bdt_canvas.SaveAs("Ratio_vs_BDTCut_%s.pdf")
                
        
        
        
        
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
                        #bdt_cut = 0.24
                if(categ == 'taumu'):
                        treeName = 'ztau3mutaumu'
                        signalnorm = 0.00000822810
                        cat_label = r"$\tau_{\mu}$"
                        #bdt_cut = 0.36
                if(categ == 'tauhA'):
                        treeName = 'ztau3mutauh_A'
                        signalnorm = 0.00000815958
                        cat_label = r"$\tau_{h,1-prong}$"
                        #bdt_cut = 0.28
                if(categ == 'tauhB'):
                        treeName = 'ztau3mutauh_B'
                        signalnorm = 0.00000815958
                        cat_label = r"$\tau_{h,3-prong}$"
                        #bdt_cut = 0.20
                if(categ == 'all'):
                        treeName   = 'ztautau'
                        signalnorm = 0.00000824176
                        cat_label = "Inclusive"
                        #bdt_cut = 0.59
                
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
                
                MCSelector = RooFormulaVar('MCSelector', 'MCSelector', ' bdt_cv > ' + str(bdt_cut)+' & '+ phivetoes+' isMC !=0 & (isMC == 211 | isMC == 210231 | isMC == 210232 | isMC == 210233 ) & (tripletMass>=%s & tripletMass<=%s) ' %(fit_range_lo,fit_range_hi) , RooArgList(variables))
                
                fullmc_unweighted = RooDataSet('mc', 'mc', tree, variables, MCSelector)
                dataset_vars = fullmc_unweighted.get()
                dataset_vars.add(scale)
                
                fullmc = RooDataSet('mc', 'mc', fullmc_unweighted, dataset_vars, "",'scale')
                
                
                BlindDataSelector = RooFormulaVar('DataSelector', 'DataSelector',' bdt_cv > ' + str(bdt_cut)+' & '+ phivetoes+' isMC == 0 & (tripletMass>=%s & tripletMass<=%s) ' %(fit_range_lo,fit_range_hi) , RooArgList(variables))
                        
                fulldata = RooDataSet('fulldata', 'fulldata', tree,  variables, BlindDataSelector)
                
                tripletMass.setRange("fullRange",1.73,1.82)
                
                # Crystal Ball #1 (left tail)
                mean = ROOT.RooRealVar("mean", "mean", 1.776, 1.75, 1.80)
                sigma = ROOT.RooRealVar("sigma", "sigma", 0.02, 0.001, 0.1)
                alpha1 = ROOT.RooRealVar("alpha1", "alpha1", 1.5, 0.5, 5.0)
                n1 = ROOT.RooRealVar("n1", "n1", 2.0, 0.5, 1000.0)
                cb1 = ROOT.RooCrystalBall("cb1", "cb1", tripletMass, mean, sigma, alpha1, n1)
                
                # Crystal Ball #2 (right tail)
                alpha2 = ROOT.RooRealVar("alpha2", "alpha2", -1.5, -5.0, -0.5)
                n2 = ROOT.RooRealVar("n2", "n2", 2.0, 0.5, 1000.0)
                cb2 = ROOT.RooCrystalBall("cb2", "cb2", tripletMass, mean, sigma, alpha2, n2)
                
                # Combine both CBs into a double CB
                frac = ROOT.RooRealVar("frac", "frac", 0.5, 0.0, 1.0)
                mc_pdf = ROOT.RooAddPdf("mc_pdf", "Double Crystal Ball", cb1, cb2, frac)
                
                fitresult = mc_pdf.fitTo(fullmc_unweighted, ROOT.RooFit.Range("fullRange"), ROOT.RooFit.Save())
                
                Gaussian_Sigma_From_Loose_BDT_Cut = sigma.getVal()
                
                sigma_scale = 0.0
                signal_peak_region_min = 0.0
                signal_peak_region_max = 0.0
                
                # Range of sigma factors
                Xa_min = 1.0
                Xa_max = 3.5
                N_a = 25
                
                triplet_mass_bins = 40
                
                output_name = "limits_summary_%s.txt" % (categ)
                
                with open(output_name, "w") as out_file:
                        out_file.write("sigma_scale\tmin\tmax\tsignal\tbackground\tLimit_UL_Calc\tLimit_Bayesian\n")
                
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
                        
                        tripletMass.setRange('left' , fit_range_lo    , signal_peak_region_min)
                        tripletMass.setRange('right', signal_peak_region_max , fit_range_hi)
                        tripletMass.setRange("SIG",signal_peak_region_min,signal_peak_region_max)
                        
                        #nbkg = ROOT.RooRealVar('nbkg', 'nbkg', 1000, 0, 500000)
                        #slope = ROOT.RooRealVar('slope', 'slope', 1.0, -100, 100)
                        #expo = ROOT.RooExponential('bkg_expo', 'bkg_expo', tripletMass, slope)
                        #pdfmodel = ROOT.RooAddPdf('bkg_extended_expo', 'bkg_extended_expo', ROOT.RooArgList(expo), ROOT.RooArgList(nbkg))
                        #results_pdf = pdfmodel.fitTo(Data_reduced_dataset, ROOT.RooFit.Range('left,right'), ROOT.RooFit.Save())
                        
                        nbkg = ROOT.RooRealVar('nbkg', 'nbkg', 1000, 0, 500000)
                        flat = ROOT.RooPolynomial('bkg_flat', 'bkg_flat', tripletMass)
                        pdfmodel = ROOT.RooAddPdf('bkg_extended_flat', 'bkg_extended_flat', ROOT.RooArgList(flat), ROOT.RooArgList(nbkg))
                        results_pdf = pdfmodel.fitTo(Data_reduced_dataset, ROOT.RooFit.Range('left,right'), ROOT.RooFit.Save())
                        
                        SG_integral = pdfmodel.createIntegral(ROOT.RooArgSet(tripletMass), ROOT.RooArgSet(tripletMass), "SIG").getVal()
                        SB_integral = pdfmodel.createIntegral(ROOT.RooArgSet(tripletMass), ROOT.RooArgSet(tripletMass), "left,right").getVal()
                        
                        N_s_1 = 0.0
                        N_b_1 = 0.0
                        
                        N_s_1 = MC_reduced_dataset.sumEntries()
                        N_b_1 = nbkg.getVal()*SG_integral
                
                        print(f"signal: {N_s_1} bkg: {N_b_1}")
                
                        if round(N_b_1 / 0.00001) * 0.00001 > 0.0:
                            cmd = f"python3 card_modifiers/card_mod.py --luminosity 59.0 --s {N_s_1} --b {N_b_1}"
                            os.system(cmd)
                
                            #cmd_bayesian = f"combine -M BayesianSimple modified_simplest_card.txt --cl 0.9 -t 100 > Test_Bayesian_Output.txt"
                            #os.system(cmd_bayesian)
                            #
                            #with open(f"Test_Bayesian_Output.txt") as f:
                            #    for line in f:
                            #        if "median expected limit" in line:
                            #            tokens = line.split()
                            #            print(f"Limit Bayesian: {tokens[5]}")
                
                            cmd_ul = f"python3 ../Projections/CLs_UL_Calculator_Integral.py {N_s_1} {N_b_1} > out_UL_Calc.txt"
                            os.system(cmd_ul)
                            
                            with open("out_UL_Calc.txt") as f:
                                for line in f:
                                    if "r_val" in line:
                                        limit_val = line.split()[-1]
                                        print(f"Limit UL_Calc: {limit_val}")
                                        
                                        with open(output_name, "a") as out_file:
                                            out_file.write(f"{sigma_scale:.4f}\t{signal_peak_region_min:.5f}\t{signal_peak_region_max:.5f}\t{N_s_1:.3f}\t{N_b_1:.2f}\t{limit_val}\n")
                                            #out_file.write(f"{sigma_scale:.4f}\t{signal_peak_region_min:.5f}\t{signal_peak_region_max:.5f}\t{N_s_1:.3f}\t{N_b_1:.2f}\t{limit_val}\t{tokens[5]}\n")
                                            
                                            
                                            
        def get_signal_fit(self, datafile, categ, WhetherMC=False):
                
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
                        bdt_cut = 0.24
                if(categ == 'taumu'):
                        treeName = 'ztau3mutaumu'
                        signalnorm = 0.00000822810
                        cat_label = r"$\tau_{\mu}$"
                        bdt_cut = 0.36
                if(categ == 'tauhA'):
                        treeName = 'ztau3mutauh_A'
                        signalnorm = 0.00000815958
                        cat_label = r"$\tau_{h,1-prong}$"
                        bdt_cut = 0.28
                if(categ == 'tauhB'):
                        treeName = 'ztau3mutauh_B'
                        signalnorm = 0.00000815958
                        cat_label = r"$\tau_{h,3-prong}$"
                        bdt_cut = 0.20
                if(categ == 'all'):
                        treeName   = 'ztautau'
                        signalnorm = 0.00000824176
                        cat_label = "Inclusive"
                        bdt_cut = 0.59
                
                bdt_cut = 0.1
                bdt_cut = -1.0
                
                tree = MiniTreeFile.Get(treeName)
                
                tripletMass          = ROOT.RooRealVar('tripletMass'                , '3#mu mass'           , fit_range_lo, fit_range_hi, 'GeV')
                bdt_cv               = ROOT.RooRealVar('bdt_cv'                     , 'bdt_cv'              , -1 , 1)
                dimu_OS1             = ROOT.RooRealVar('dimu_OS1'                   , 'dimu_OS1'            ,  0 , 100)
                dimu_OS2             = ROOT.RooRealVar('dimu_OS2'                   , 'dimu_OS2'            ,  0 , 100)
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
                
                MCSelector = RooFormulaVar('MCSelector', 'MCSelector', ' bdt_cv > ' + str(bdt_cut)+' & '+ phivetoes+' isMC !=0 & (isMC == 211 | isMC == 210231 | isMC == 210232 | isMC == 210233 ) & (tripletMass>=%s & tripletMass<=%s) ' %(fit_range_lo,fit_range_hi) , RooArgList(variables))
                
                fullmc_unweighted = RooDataSet('mc', 'mc', tree, variables, MCSelector)
                dataset_vars = fullmc_unweighted.get()
                dataset_vars.add(scale)
                
                fullmc = RooDataSet('mc', 'mc', fullmc_unweighted, dataset_vars, "",'scale')
                
                tripletMass.setRange("fullRange",1.6,2.0)
                
                # Crystal Ball #1 (left tail)
                mean = ROOT.RooRealVar("mean", "mean", 1.776, 1.75, 1.80)
                sigma = ROOT.RooRealVar("sigma", "sigma", 0.02, 0.001, 0.1)
                alpha1 = ROOT.RooRealVar("alpha1", "alpha1", 1.5, 0.5, 5.0)
                n1 = ROOT.RooRealVar("n1", "n1", 2.0, 0.5, 1000.0)
                cb1 = ROOT.RooCrystalBall("cb1", "cb1", tripletMass, mean, sigma, alpha1, n1)
                
                # Crystal Ball #2 (right tail)
                alpha2 = ROOT.RooRealVar("alpha2", "alpha2", -1.5, -5.0, -0.5)
                n2 = ROOT.RooRealVar("n2", "n2", 2.0, 0.5, 1000.0)
                cb2 = ROOT.RooCrystalBall("cb2", "cb2", tripletMass, mean, sigma, alpha2, n2)
                
                # Combine both CBs into a double CB
                frac = ROOT.RooRealVar("frac", "frac", 0.5, 0.0, 1.0)
                mc_pdf = ROOT.RooAddPdf("mc_pdf", "Double Crystal Ball", cb1, cb2, frac)
                
                fitresult = mc_pdf.fitTo(fullmc_unweighted, ROOT.RooFit.Range("fullRange"), ROOT.RooFit.Save())
                
                
                
                
                # Create canvas
                canvas = ROOT.TCanvas("canvas", "canvas", 800, 800)
                canvas.SetLeftMargin(L/W)
                canvas.SetRightMargin(R/W)
                canvas.SetTopMargin(T/H)
                canvas.SetBottomMargin(B/H)
                canvas.cd()
                
                # Create frame for tripletMass
                frame = tripletMass.frame(ROOT.RooFit.Title(f"Signal DoubleCB Fit: {categ}"))
                
                # Plot the dataset (scaled signal)
                fullmc.plotOn(
                    frame,
                    ROOT.RooFit.Binning(40),
                    ROOT.RooFit.MarkerStyle(20),
                    ROOT.RooFit.MarkerColor(ROOT.kBlack),
                    ROOT.RooFit.Name("sig_data")
                )
                
                # Plot the fitted Double Crystal Ball
                mc_pdf.plotOn(
                    frame,
                    ROOT.RooFit.LineColor(ROOT.kRed),
                    ROOT.RooFit.Name("fit_curve")
                )
                
                # Extract fit values
                mean_val = mean.getVal()
                mean_err = mean.getError()
                sigma_val = sigma.getVal()
                sigma_err = sigma.getError()
                alphaL_val = alpha1.getVal()
                alphaR_val = alpha2.getVal()
                nL_val = n1.getVal()
                nR_val = n2.getVal()
                
                frame.GetXaxis().SetNdivisions(505)
                frame.GetYaxis().SetTitleOffset(1.6)
                frame.GetXaxis().SetTitle("m_{3#mu} (GeV)")
                frame.Draw()
                
                # Add legend
                legend = ROOT.TLegend(0.60, 0.65, 0.88, 0.88)
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                legend.AddEntry(frame.findObject("sig_data"), "Signal MC", "lep")
                legend.AddEntry(frame.findObject("fit_curve"), "Double CB Fit", "l")
                legend.Draw()
                
                # Add TLatex for fit parameters
                latex = ROOT.TLatex()
                latex.SetNDC()
                latex.SetTextSize(0.035)
                latex.DrawLatex(0.60, 0.60, f"#mu = {mean_val:.5f} #pm {mean_err:.5f}")
                latex.DrawLatex(0.60, 0.55, f"#sigma = {sigma_val:.5f} #pm {sigma_err:.5f}")
                latex.DrawLatex(0.60, 0.50, f"#alpha_{{L}} = {alphaL_val:.3f}, n_{{L}} = {nL_val:.3f}")
                latex.DrawLatex(0.60, 0.45, f"#alpha_{{R}} = {alphaR_val:.3f}, n_{{R}} = {nR_val:.3f}")
                
                # Save the plot
                CMSStyle.CMS_lumi(canvas, 5, 11)
                canvas.Update()
                canvas.SaveAs(f"signal_fit_{categ}.png")
                print(f"Saved: signal_fit_{categ}.png")

                
                
                
        
        
        def get_signal_fit_norefit(self, datafile, categ, WhetherMC=False):
                
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
                        bdt_cut = 0.24
                if(categ == 'taumu'):
                        treeName = 'ztau3mutaumu'
                        signalnorm = 0.00000822810
                        cat_label = r"$\tau_{\mu}$"
                        bdt_cut = 0.36
                if(categ == 'tauhA'):
                        treeName = 'ztau3mutauh_A'
                        signalnorm = 0.00000815958
                        cat_label = r"$\tau_{h,1-prong}$"
                        bdt_cut = 0.28
                if(categ == 'tauhB'):
                        treeName = 'ztau3mutauh_B'
                        signalnorm = 0.00000815958
                        cat_label = r"$\tau_{h,3-prong}$"
                        bdt_cut = 0.20
                if(categ == 'all'):
                        treeName   = 'ztautau'
                        signalnorm = 0.00000824176
                        cat_label = "Inclusive"
                        bdt_cut = 0.59
                
                bdt_cut = 0.1
                bdt_cut = -1.0
                
                tree = MiniTreeFile.Get(treeName)
                
                tripletMassNonRefit  = ROOT.RooRealVar('tripletMassNonRefit'        , '3#mu mass'           , fit_range_lo, fit_range_hi, 'GeV')
                bdt_cv               = ROOT.RooRealVar('bdt_cv'                     , 'bdt_cv'              , -1 , 1)
                dimu_OS1             = ROOT.RooRealVar('dimu_OS1'                   , 'dimu_OS1'            ,  0 , 100)
                dimu_OS2             = ROOT.RooRealVar('dimu_OS2'                   , 'dimu_OS2'            ,  0 , 100)
                event_weight         = ROOT.RooRealVar('weight'                     , 'event_weight'        ,  0,  5)  # this weight includes also the scale  mc signal scale
                category             = ROOT.RooRealVar('category'                   , 'category'            ,  0,  5)
                isMC                 = ROOT.RooRealVar('isMC'                       , 'isMC'                ,  0,  1000000)
                scale                = ROOT.RooRealVar('scale'                      , 'scale'               ,  signalnorm)  
                
                
                
                variables = ROOT.RooArgSet()
                variables.add(tripletMassNonRefit)
                variables.add(bdt_cv)
                variables.add(dimu_OS1)
                variables.add(dimu_OS2)
                variables.add(event_weight)
                variables.add(category)
                variables.add(isMC)
                
                #Define mc and data which can be reduced later
                
                phivetoes="(fabs(dimu_OS1 - 1.020)>0.020)&(fabs(dimu_OS2 - 1.020)>0.020)&"
                omegavetoes="fabs(dimu_OS1 - 0.782)>0.020&fabs(dimu_OS2 - 0.782)>0.020&"
                
                MCSelector = RooFormulaVar('MCSelector', 'MCSelector', ' bdt_cv > ' + str(bdt_cut)+' & '+ phivetoes+' isMC !=0 & (isMC == 211 | isMC == 210231 | isMC == 210232 | isMC == 210233 ) & (tripletMassNonRefit>=%s & tripletMassNonRefit<=%s) ' %(fit_range_lo,fit_range_hi) , RooArgList(variables))
                
                fullmc_unweighted = RooDataSet('mc', 'mc', tree, variables, MCSelector)
                dataset_vars = fullmc_unweighted.get()
                dataset_vars.add(scale)
                
                fullmc = RooDataSet('mc', 'mc', fullmc_unweighted, dataset_vars, "",'scale')
                
                tripletMassNonRefit.setRange("fullRange",1.6,2.0)
                
                # Crystal Ball #1 (left tail)
                mean = ROOT.RooRealVar("mean", "mean", 1.776, 1.75, 1.80)
                sigma = ROOT.RooRealVar("sigma", "sigma", 0.02, 0.001, 0.1)
                alpha1 = ROOT.RooRealVar("alpha1", "alpha1", 1.5, 0.5, 5.0)
                n1 = ROOT.RooRealVar("n1", "n1", 2.0, 0.5, 1000.0)
                cb1 = ROOT.RooCrystalBall("cb1", "cb1", tripletMassNonRefit, mean, sigma, alpha1, n1)
                
                # Crystal Ball #2 (right tail)
                alpha2 = ROOT.RooRealVar("alpha2", "alpha2", -1.5, -5.0, -0.5)
                n2 = ROOT.RooRealVar("n2", "n2", 2.0, 0.5, 1000.0)
                cb2 = ROOT.RooCrystalBall("cb2", "cb2", tripletMassNonRefit, mean, sigma, alpha2, n2)
                
                # Combine both CBs into a double CB
                frac = ROOT.RooRealVar("frac", "frac", 0.5, 0.0, 1.0)
                mc_pdf = ROOT.RooAddPdf("mc_pdf", "Double Crystal Ball", cb1, cb2, frac)
                
                fitresult = mc_pdf.fitTo(fullmc_unweighted, ROOT.RooFit.Range("fullRange"), ROOT.RooFit.Save())
                
                
                
                
                # Create canvas
                canvas = ROOT.TCanvas("canvas", "canvas", 800, 800)
                canvas.SetLeftMargin(L/W)
                canvas.SetRightMargin(R/W)
                canvas.SetTopMargin(T/H)
                canvas.SetBottomMargin(B/H)
                canvas.cd()
                
                # Create frame for tripletMassNonRefit
                frame = tripletMassNonRefit.frame(ROOT.RooFit.Title(f"Signal DoubleCB Fit: {categ}"))
                
                # Plot the dataset (scaled signal)
                fullmc.plotOn(
                    frame,
                    ROOT.RooFit.Binning(40),
                    ROOT.RooFit.MarkerStyle(20),
                    ROOT.RooFit.MarkerColor(ROOT.kBlack),
                    ROOT.RooFit.Name("sig_data")
                )
                
                # Plot the fitted Double Crystal Ball
                mc_pdf.plotOn(
                    frame,
                    ROOT.RooFit.LineColor(ROOT.kRed),
                    ROOT.RooFit.Name("fit_curve")
                )
                
                # Extract fit values
                mean_val = mean.getVal()
                mean_err = mean.getError()
                sigma_val = sigma.getVal()
                sigma_err = sigma.getError()
                alphaL_val = alpha1.getVal()
                alphaR_val = alpha2.getVal()
                nL_val = n1.getVal()
                nR_val = n2.getVal()
                
                frame.GetXaxis().SetNdivisions(505)
                frame.GetYaxis().SetTitleOffset(1.6)
                frame.GetXaxis().SetTitle("m_{3#mu} (GeV)")
                frame.Draw()
                
                # Add legend
                legend = ROOT.TLegend(0.60, 0.65, 0.88, 0.88)
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                legend.AddEntry(frame.findObject("sig_data"), "Signal MC", "lep")
                legend.AddEntry(frame.findObject("fit_curve"), "Double CB Fit", "l")
                legend.Draw()
                
                # Add TLatex for fit parameters
                latex = ROOT.TLatex()
                latex.SetNDC()
                latex.SetTextSize(0.035)
                latex.DrawLatex(0.60, 0.60, f"#mu = {mean_val:.5f} #pm {mean_err:.5f}")
                latex.DrawLatex(0.60, 0.55, f"#sigma = {sigma_val:.5f} #pm {sigma_err:.5f}")
                latex.DrawLatex(0.60, 0.50, f"#alpha_{{L}} = {alphaL_val:.3f}, n_{{L}} = {nL_val:.3f}")
                latex.DrawLatex(0.60, 0.45, f"#alpha_{{R}} = {alphaR_val:.3f}, n_{{R}} = {nR_val:.3f}")
                
                # Save the plot
                CMSStyle.CMS_lumi(canvas, 5, 11)
                canvas.Update()
                canvas.SaveAs(f"signal_fit_norefit_{categ}.png")
                print(f"Saved: signal_fit_norefit_{categ}.png")
                
                
                
        def get_sig_bkg_efficiencies(self, datafile, categ, WhetherMC=False):
                
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
                filename = ''
                if(categ == 'taue'):
                        treeName  = 'ztau3mutaue'
                        signalnorm = 0.00000856928
                        cat_label = r"$\tau_{e}$"
                        bdt_cut = 0.24
                        filename = '../Projections/makeYield/ZTT/TextLimits_taue.txt'
                if(categ == 'taumu'):
                        treeName = 'ztau3mutaumu'
                        signalnorm = 0.00000822810
                        cat_label = r"$\tau_{\mu}$"
                        bdt_cut = 0.36
                        filename = '../Projections/makeYield/ZTT/TextLimits_taumu.txt'
                if(categ == 'tauhA'):
                        treeName = 'ztau3mutauh_A'
                        signalnorm = 0.00000815958
                        cat_label = r"$\tau_{h,1-prong}$"
                        bdt_cut = 0.28
                        filename = '../Projections/makeYield/ZTT/TextLimits_tauhA.txt'
                if(categ == 'tauhB'):
                        treeName = 'ztau3mutauh_B'
                        signalnorm = 0.00000815958
                        cat_label = r"$\tau_{h,3-prong}$"
                        bdt_cut = 0.20
                        filename = '../Projections/makeYield/ZTT/TextLimits_tauhB.txt'
                if(categ == 'all'):
                        treeName   = 'ztautau'
                        signalnorm = 0.00000824176
                        cat_label = "Inclusive"
                        bdt_cut = 0.59
                        filename = '../Projections/makeYield/ZTT/TextLimits_tauhB.txt'
                
                bdt_cut = 0.1
                bdt_cut = -1.0
                
                tree = MiniTreeFile.Get(treeName)
                
                tripletMass          = ROOT.RooRealVar('tripletMass'                , '3#mu mass'           , fit_range_lo, fit_range_hi, 'GeV')
                bdt_cv               = ROOT.RooRealVar('bdt_cv'                     , 'bdt_cv'              , -1 , 1)
                dimu_OS1             = ROOT.RooRealVar('dimu_OS1'                   , 'dimu_OS1'            ,  0 , 100)
                dimu_OS2             = ROOT.RooRealVar('dimu_OS2'                   , 'dimu_OS2'            ,  0 , 100)
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
                
                MCSelector = RooFormulaVar('MCSelector', 'MCSelector', ' bdt_cv > ' + str(bdt_cut)+' & '+ phivetoes+' isMC !=0 & (isMC == 211 | isMC == 210231 | isMC == 210232 | isMC == 210233 ) & (tripletMass>=%s & tripletMass<=%s) ' %(fit_range_lo,fit_range_hi) , RooArgList(variables))
                
                fullmc_unweighted = RooDataSet('mc', 'mc', tree, variables, MCSelector)
                dataset_vars = fullmc_unweighted.get()
                dataset_vars.add(scale)
                
                fullmc = RooDataSet('mc', 'mc', fullmc_unweighted, dataset_vars, "",'scale')
                
                total_mc_scaled_signal_events = fullmc.sumEntries()
                
                tripletMass.setRange("fullRange",1.73,1.82)
                
                mean = ROOT.RooRealVar("mean", "mean", 1.776, 0., 5.)
                sigma = ROOT.RooRealVar("sigma", "sigma", 0.5, 0.001, 10)
                Gauss = ROOT.RooGaussian("Gauss", "Gauss dist", tripletMass, mean, sigma)
                GaussNorm = ROOT.RooRealVar("GaussNorm", "GaussNorm", 0.5, 0.001, 1.0)
                mc_pdf = ROOT.RooAddPdf("mc_pdf", "mc_pdf", ROOT.RooArgList(Gauss), ROOT.RooArgList(GaussNorm))
                mc_fitresult = mc_pdf.fitTo(fullmc_unweighted, ROOT.RooFit.Range("fullRange"), ROOT.RooFit.Save())
                
                Gaussian_Sigma_From_Loose_BDT_Cut = sigma.getVal()
                
                BlindDataSelector = RooFormulaVar('DataSelector', 'DataSelector', phivetoes+' isMC == 0 & (tripletMass<=%s || tripletMass>=%s) & (tripletMass>=%s & tripletMass<=%s) ' %(signal_range_lo,signal_range_hi,fit_range_lo,fit_range_hi) , RooArgList(variables))
                
                fulldata_shape = RooDataSet('fulldata_shape', 'fulldata_shape', tree,  variables, BlindDataSelector)
                
                fulldata = fulldata_shape
                
                total_data_sideband_events = fulldata.sumEntries()
                
                bdt = []
                lumi = []
                sig = []
                bkg = []
                bkg_err = []
                
                # Read the file and extract the data
                with open( filename, "r") as file:
                        for line in file:
                            parts = line.split()
                            bdt.append(float(parts[1]))       # bdt
                            lumi.append(float(parts[3]))      # lumi
                            sig.append(float(parts[7]))       # sig
                            bkg.append(float(parts[9]))       # bkg
                            bkg_err.append(float(parts[11]))  # bkg_err
                            
                # Compute efficiencies
                sig_eff = []
                bkg_eff = []
                
                for i in range(len(bdt)):
                    lumi_scale = lumi[i] / 59.83

                    # Scale total expected MC/data yields for that lumi
                    expected_signal_events = total_mc_scaled_signal_events * lumi_scale
                    expected_background_events = total_data_sideband_events * lumi_scale

                    # Compute efficiencies
                    sig_eff.append(sig[i] / expected_signal_events if expected_signal_events > 0 else 0)
                    bkg_eff.append( (bkg[i] * (1.0/0.21212121) ) / expected_background_events if expected_background_events > 0 else 0)

                # Optional: print or store
                for i in range(len(bdt)):
                    print(f"BDT: {bdt[i]:.2f}, SigEff: {sig_eff[i]:.4f}, BkgEff: {bkg_eff[i]:.4f}")
                
                # === Plot 1: Signal Efficiency vs BDT Cut ===
                canvas_sig = ROOT.TCanvas("canvas_sig", "Signal Efficiency", 800, 800)
                canvas_sig.SetLeftMargin( L/W )
                canvas_sig.SetRightMargin( R/W )
                canvas_sig.SetTopMargin( T/H )
                canvas_sig.SetBottomMargin( B/H )
                canvas_sig.cd()
                
                graph_sig = ROOT.TGraph(len(bdt), array.array('d', bdt), array.array('d', sig_eff))
                graph_sig.SetTitle(f"Signal Efficiency: {categ}")
                graph_sig.SetLineColor(ROOT.kRed+1)
                graph_sig.SetLineWidth(2)
                graph_sig.SetMarkerColor(ROOT.kRed+1)
                graph_sig.SetMarkerStyle(20)
                graph_sig.Draw("ALP")
                graph_sig.GetXaxis().SetTitle("BDT Cut")
                graph_sig.GetYaxis().SetTitle("Signal Efficiency")
                graph_sig.GetYaxis().SetTitleOffset(1.4)
                graph_sig.GetXaxis().SetNdivisions(505)
                graph_sig.GetYaxis().SetNdivisions(505)
                graph_sig.GetYaxis().SetRangeUser(0, 1.05)
                
                CMSStyle.CMS_lumi(canvas_sig, 5, 11)
                canvas_sig.Update()
                canvas_sig.SaveAs(f"signal_efficiency_vs_bdt_{categ}.png")
                print(f"Saved: signal_efficiency_vs_bdt_{categ}.png")
                
                
                # === Plot 2: Background Efficiency vs BDT Cut ===
                canvas_bkg = ROOT.TCanvas("canvas_bkg", "Background Efficiency", 800, 800)
                canvas_bkg.SetLeftMargin( L/W )
                canvas_bkg.SetRightMargin( R/W )
                canvas_bkg.SetTopMargin( T/H )
                canvas_bkg.SetBottomMargin( B/H )
                canvas_bkg.cd()
                
                graph_bkg = ROOT.TGraph(len(bdt), array.array('d', bdt), array.array('d', bkg_eff))
                graph_bkg.SetTitle(f"Background Efficiency: {categ}")
                graph_bkg.SetLineColor(ROOT.kBlue+1)
                graph_bkg.SetLineWidth(2)
                graph_bkg.SetMarkerColor(ROOT.kBlue+1)
                graph_bkg.SetMarkerStyle(21)
                graph_bkg.Draw("ALP")
                graph_bkg.GetXaxis().SetTitle("BDT Cut")
                graph_bkg.GetYaxis().SetTitle("Background Efficiency")
                graph_bkg.GetYaxis().SetTitleOffset(1.8)
                graph_bkg.GetXaxis().SetNdivisions(505)
                graph_bkg.GetYaxis().SetNdivisions(505)
                graph_bkg.GetYaxis().SetRangeUser(0, 0.006)
                
                CMSStyle.CMS_lumi(canvas_bkg, 5, 11)
                canvas_bkg.Update()
                canvas_bkg.SaveAs(f"background_efficiency_vs_bdt_{categ}.png")
                print(f"Saved: background_efficiency_vs_bdt_{categ}.png")


                
                
if __name__ == "__main__":
    
        ROOT.gROOT.SetBatch(True)
        
        #datafile = "../../Combine_Tree_ztau3mutau_orig_PostBDT.root"
        
        #datafile = "../../Combine_Tree_ztau3mutau_PF_PostBDT.root"
        
        datafile = "../../Combine_Tree_ztau3mutau_PF_PostBDT_symmetric.root"
        
        datafile_norefit = "../../Combine_Tree_ztau3mutau_PF_PostBDT_unrefit_mass.root"
        
        #datafile_ZTTmass = "../../Combine_Tree_ztau3mutau_ZTTMass_PF_PostBDT.root"
        datafile_ZTTmass = "../../Combine_Tree_ztau3mutau_ZTTMass_origTracker_PostBDT.root"
        
        # Categories to loop over
        #categories = ['taue', 'taumu', 'tauhA', 'tauhB','all']
        categories = ['taue', 'taumu', 'tauhA', 'tauhB']
        #categories = ['tauhB']
        
        # Example BDT selections you want to compare
        bdt_selections = [
            "( isMC == 251 || isMC == 251231 || isMC == 251232 || isMC == 251233 )",
            "( isMC == 252 || isMC == 252231 || isMC == 252232 || isMC == 252233 )",
            "( isMC == 253 || isMC == 253231 || isMC == 253232 || isMC == 253233 )",
            "( isMC == 254 || isMC == 254231 || isMC == 254232 || isMC == 254233 )",
            "( isMC == 255 || isMC == 255231 || isMC == 255232 || isMC == 255233 )"
        ]
        
        r"#tau_{e}"
        selection_labels = [r"m_{#tau} = 1.65", r"m_{#tau} = 1.70", r"m_{#tau} = 1.85", r"m_{#tau} = 1.90", r"m_{#tau} = 1.95"]
        
        # Sample datasets (for Compare_BDT_Scores_MultipleDatasets)
        dataset_files = [
            "../../Combine_Tree_ztau3mutau_orig_PostBDT.root",
            #"../../Combine_Tree_ztau3mutau_PFGL_PostBDT.root",
            #"../../Combine_Tree_ztau3mutau_PF_PostBDT.root",
            "../../Combine_Tree_ztau3mutau_PF_PostBDT_symmetric.root",
        ]
        
        # Loop over categories
        for categ in categories:
            print("Running for category: " + str(categ))
            
            # These methods belong to a class, BDTPlotter()
            BDTPlotter = BDT_Shape_Comparisons()
        
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
            BDTPlotter.Compare_BDT_Scores_MultipleDatasets(dataset_files, categ, False)
            
            # 2. Compare different selections on the same file. Pick any one file.
            #BDTPlotter.Compare_BDT_Scores_MultipleSelections(datafile, categ, bdt_selections, selection_labels, isMC=True)
            
            # 3. Full comparison: different files + different selections
            #BDTPlotter.Compare_BDT_Scores_MultipleDatasetsAndSelections(dataset_files, datafile, categ, bdt_selections, selection_labels, isMC=True)
            
            # 4. Full comparison: different files + different selections on same file
            #BDTPlotter.Compare_BDT_Scores_MultipleDatasetsAndSelectionsZTTMass(dataset_files, datafile_ZTTmass, categ, bdt_selections, selection_labels, True)
            
            # 5. Constancy of 'r' plots
            #BDTPlotter.Plot_Bdt_Symmetry(datafile, categ, False)
            
            # 6. Get signal peak width
            #BDTPlotter.get_signal_window(datafile, categ, False)
            
            # 7. Get signal peak fit
            #BDTPlotter.get_signal_fit(datafile, categ, False)
            
            # 7. Get signal peak fit
            #BDTPlotter.get_signal_fit_norefit(datafile_norefit, categ, False)
            
            # 8. Get signal sig/bkg efficiencies
            #BDTPlotter.get_sig_bkg_efficiencies(datafile, categ, False)