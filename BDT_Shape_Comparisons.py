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
            
        
        
        def BDT_Shape_Comparisons(self,datafile,categ):
                
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
                self.bdt_cv          = ROOT.RooRealVar('bdt_cv'                     , 'bdt_cv'              , -1 , 1)
                dimu_OS1             = ROOT.RooRealVar('dimu_OS1'                   , 'dimu_OS1'            ,  0 , 2)
                dimu_OS2             = ROOT.RooRealVar('dimu_OS2'                   , 'dimu_OS2'            ,  0 , 2)
                event_weight         = ROOT.RooRealVar('weight'                     , 'event_weight'        ,  0,  5)  # this weight includes also the scale  mc signal scale
                category             = ROOT.RooRealVar('category'                   , 'category'            ,  0,  5)
                isMC                 = ROOT.RooRealVar('isMC'                       , 'isMC'                ,  0,  1000000)
                scale                = ROOT.RooRealVar('scale'                      , 'scale'               ,  signalnorm)  
                
                
                
                variables = ROOT.RooArgSet()
                variables.add(tripletMass)
                variables.add(self.bdt_cv)
                variables.add(dimu_OS1)
                variables.add(dimu_OS2)
                variables.add(event_weight)
                variables.add(category)
                variables.add(isMC)
                variables.add(scale)
                
                phivetoes="(fabs(dimu_OS1 - 1.020)>0.020)&(fabs(dimu_OS2 - 1.020)>0.020)&"
                omegavetoes="fabs(dimu_OS1 - 0.782)>0.020&fabs(dimu_OS2 - 0.782)>0.020&"
                
                
                # For fitting BDT Output in Data
                
                BDT_Score_Min=-0.3
                
                BlindDataSelector = ROOT.RooFormulaVar('DataSelector', 'DataSelector', phivetoes+' isMC == 0 & (tripletMass<=%s || tripletMass>=%s) & (tripletMass>=%s & tripletMass<=%s) ' %(signal_range_lo,signal_range_hi,fit_range_lo,fit_range_hi) , ROOT.RooArgList(variables))
                
                fulldata = ROOT.RooDataSet('data', 'data', tree,  variables, BlindDataSelector)
                
                self.bdt_cv.setRange("BDT_Fit_Range", BDT_Score_Min, 1.0);
                
                self.a = ROOT.RooRealVar("a", "a", 1.0, 0.0, 10.0)
                self.b = ROOT.RooRealVar("b", "b", 1.0, -10.0, 10.0)
                self.c = ROOT.RooRealVar("c", "c", 1.0, -10.0, 10.0)
                self.d = ROOT.RooRealVar("d", "d", 1.0, -10.0, 10.0)
                
                #quadratic = ROOT.RooFormulaVar("quadratic", "a + b*self.bdt_cv + c*self.bdt_cv*self.bdt_cv", ROOT.RooArgList(a, b, c, self.bdt_cv))
                #expModel = ROOT.RooGenericPdf("expModel", "exp(quadratic)", ROOT.RooArgList(quadratic)) #Exponential of the quadratic polynomial
                
                self.quadratic = ROOT.RooFormulaVar("quadratic", "a + b*bdt_cv + c*bdt_cv*bdt_cv + d*bdt_cv*bdt_cv*bdt_cv", ROOT.RooArgList(self.a, self.b, self.c, self.d, self.bdt_cv))
                self.expModel = ROOT.RooGenericPdf("expModel", "exp(quadratic)", ROOT.RooArgList(self.quadratic)) #Exponential of the cubic polynomial
                
                self.BDTNorm = ROOT.RooRealVar("BDTNorm", "BDTNorm", 500.0, 0.1, 1000000000)
                self.BDT_distribution = ROOT.RooAddPdf("BDT_distribution", "BDT_distribution",ROOT.RooArgList(self.expModel), ROOT.RooArgList(self.BDTNorm))
                #BDT_distribution = ROOT.RooAddPdf("BDT_distribution", "BDT_distribution",ROOT.RooArgList(quadratic), ROOT.RooArgList(BDTNorm))
                
                results_pdf = self.BDT_distribution.fitTo(fulldata, ROOT.RooFit.Range('BDT_Fit_Range'), ROOT.RooFit.Save())
                results_pdf.Print()
                
                
                
                
                # For fitting BDT Output in Signal
                
                self.MCSelector = ROOT.RooFormulaVar('MCSelector', 'MCSelector', phivetoes+' isMC !=0 & (isMC == 211 | isMC == 210231 | isMC == 210232 | isMC == 210233 ) & (tripletMass>=%s & tripletMass<=%s) ' %(fit_range_lo,fit_range_hi) , ROOT.RooArgList(variables))
                
                self.fullmc = ROOT.RooDataSet('mc', 'mc', tree, variables, self.MCSelector,'scale')
                
                self.bdt_cv.setRange("BDT_MC_Fit_Range", -1.0, 1.0);
                
                self.bgausmeanMC = ROOT.RooRealVar("bgausmeanMC", "bgausmeanMC", 0.5, 0.0, 0.9)
                self.bgaussigmaMC_a = ROOT.RooRealVar("bgaussigmaMC_a", "bgaussigmaMC_a", 0.2, 0.000001, 1.0)
                self.bgaussigmaMC_b = ROOT.RooRealVar("bgaussigmaMC_b", "bgaussigmaMC_b", 0.2, 0.000001, 1.0)
                
                self.bgaus_distMC = ROOT.RooBifurGauss("bgaus_distMC", "bgaus dist MC", self.bdt_cv, self.bgausmeanMC, self.bgaussigmaMC_a, self.bgaussigmaMC_b)
                
                self.BDTNorm_MC = ROOT.RooRealVar("BDTNorm_MC", "BDTNorm_MC", 500.0, 0.1, 50000)
                self.BDT_distribution_MC = ROOT.RooAddPdf("BDT_distribution", "BDT_distribution",ROOT.RooArgList(self.bgaus_distMC), ROOT.RooArgList(self.BDTNorm_MC))
                
                results_mcpdf = self.BDT_distribution_MC.fitTo(self.fullmc, ROOT.RooFit.Range('BDT_MC_Fit_Range'), ROOT.RooFit.Save())
                results_mcpdf.Print()
                
                
                # Plot BDT for data and MC
                
                frame1 = self.bdt_cv.frame()
                frame1.SetTitle('')
                frame1.GetXaxis().SetTitle("BDT score, "+cat_label)
                
                frame2 = self.bdt_cv.frame()
                frame2.SetTitle('')
                frame2.GetXaxis().SetTitle("BDT score, "+cat_label)
                
                nbins = 100
                
                self.fullmc.plotOn(frame1, 
                              ROOT.RooFit.Binning(nbins), 
                              ROOT.RooFit.XErrorSize(0), 
                              ROOT.RooFit.LineWidth(2),
                              ROOT.RooFit.MarkerStyle(6),
                              ROOT.RooFit.MarkerColor(ROOT.kRed ),
                              ROOT.RooFit.MarkerSize(0.75),
                              ROOT.RooFit.FillColor(ROOT.kCyan  + 2)
                )
                #self.BDT_distribution_MC.plotOn(frame1, ROOT.RooFit.LineColor(ROOT.kRed ))
                
                
                fulldata.plotOn(frame2, 
                        ROOT.RooFit.Binning(nbins),
                        ROOT.RooFit.MarkerStyle(20),
                        ROOT.RooFit.MarkerColor(ROOT.kBlack), 
                        ROOT.RooFit.MarkerSize(0.75))
                        
                #BDT_distribution.plotOn(frame2,  ROOT.RooFit.LineColor(ROOT.kBlue) , ROOT.RooFit.Normalization(fulldata.sumEntries("1", "BDT_Fit_Range"), ROOT.RooAbsReal.NumEvent), ROOT.RooFit.ProjectionRange('BDT_Fit_Range') )
                #BDT_distribution.plotOn(frame2,  ROOT.RooFit.LineColor(ROOT.kBlue) , ROOT.RooFit.Normalization(BDTNorm.getVal()*BDT_distribution.createIntegral(ROOT.RooArgSet(self.bdt_cv), ROOT.RooArgSet(self.bdt_cv), "BDT_Fit_Range").getVal(), ROOT.RooAbsReal.NumEvent), ROOT.RooFit.ProjectionRange('BDT_Fit_Range') )
                
                self.BDT_distribution.plotOn(frame2,  ROOT.RooFit.LineColor(ROOT.kBlue) )
                
                frame1.Draw()
                ROOT.gPad.SetTicks(1,1)
                CMSStyle.CMS_lumi(ROOT.gPad, 5, 0)
                ROOT.gPad.Update()
                frame1.Draw('sameaxis')
                ROOT.gPad.SaveAs('bdt_fit_mc_'+categ+'.png')
                ROOT.gPad.Clear()
                
                ROOT.gPad.SetLogy()
                frame2.Draw()
                ROOT.gPad.SetTicks(1,1)
                CMSStyle.CMS_lumi(ROOT.gPad, 5, 0)
                ROOT.gPad.Update()
                frame2.Draw('sameaxis')
                ROOT.gPad.SaveAs('bdt_fit_bkg_'+categ+'.png')
                ROOT.gPad.SetLogy(0)
                
                
                #print("Certain BDT cut: ", self.fullmc.reduce('bdt_cv > 0.5').sumEntries())
                
                
                return frame1, frame2, self.fullmc, fulldata
                
                
        
        



if __name__ == "__main__":
    
        ROOT.gROOT.SetBatch(True)
        
        datafile = "../../Combine_Tree_ztau3mutau_orig_PostBDT.root"
        
        # Assume these methods belong to an object or class, e.g. BDTPlotter()
        BDTPlotter = BDT_Shape_Comparisons()  # Replace this with your actual class name or self instance
        
        # Categories to loop over
        categories = ['taue', 'taumu', 'tauhA', 'tauhB','all']
        
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
            BDTPlotter.Compare_BDT_Scores_MultipleDatasetsAndSelections(dataset_files, datafile, categ, bdt_selections, selection_labels, isMC=True)
            