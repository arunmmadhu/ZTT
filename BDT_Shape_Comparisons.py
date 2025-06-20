#!/usr/bin/env python2

import ROOT
from ROOT import TFile, TTree, TCanvas, TGraph, TMultiGraph, TGraphErrors, TLegend, TPaveLabel, TPaveText, TLatex
from ROOT import RooRealVar, RooFormulaVar, RooExponential, RooDataHist, RooArgList, RooAddPdf, RooFit, RooDataSet, RooGenericPdf, RooBifurGauss
import subprocess # to execute shell command
import argparse
import numpy as np
import CMS_lumi, tdrstyle
from CMSStyle import CMS_lumi
import os
import re


# CMS style
CMS_lumi.cmsText = "CMS, work in progress"
CMS_lumi.extraText = ""
CMS_lumi.cmsTextSize = 0.65
CMS_lumi.outOfFrame = True
tdrstyle.setTDRStyle()


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
        
        
        
        
        
        def Compare_BDT_Scores_MultipleDatasets(self, datafiles, categ):
            frame = ROOT.RooRealVar('bdt_cv', 'BDT Score', -1, 1).frame()
            colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta]
            
            for i, datafile in enumerate(datafiles):
                MiniTreeFile = ROOT.TFile.Open(datafile)
                tree = MiniTreeFile.Get(categ)  # e.g. 'ztau3mutauh_A'
                
                # Define all variables only once
                bdt_cv = ROOT.RooRealVar("bdt_cv", "bdt_cv", -1.0, 1.0)
                tripletMass = ROOT.RooRealVar('tripletMass', '3#mu mass', 1.6, 2.0)
                isMC = ROOT.RooRealVar("isMC", "isMC", 0, 1000000)
                weight = ROOT.RooRealVar("weight", "weight", 0, 5)
                dimu_OS1 = ROOT.RooRealVar('dimu_OS1', 'dimu_OS1', 0, 2)
                dimu_OS2 = ROOT.RooRealVar('dimu_OS2', 'dimu_OS2', 0, 2)
                
                variables = ROOT.RooArgSet(tripletMass, bdt_cv, isMC, weight, dimu_OS1, dimu_OS2)
                
                selector_str = "(fabs(dimu_OS1 - 1.020)>0.020)&&(fabs(dimu_OS2 - 1.020)>0.020)&&(fabs(dimu_OS1 - 0.782)>0.020)&&(fabs(dimu_OS2 - 0.782)>0.020)&&(tripletMass>1.6)&&(tripletMass<2.0)"
                
                dataset = ROOT.RooDataSet("data_{0}".format(i), "data_{0}".format(i), tree, variables, selector_str, "weight")
                
                dataset.plotOn(
                    frame,
                    ROOT.RooFit.Binning(100),
                    ROOT.RooFit.MarkerColor(colors[i % len(colors)]),
                    ROOT.RooFit.LineColor(colors[i % len(colors)]),
                    ROOT.RooFit.MarkerStyle(20 + i),
                    ROOT.RooFit.MarkerSize(0.75)
                )
            
            frame.GetXaxis().SetTitle("BDT score")
            frame.SetTitle("BDT score comparison")
            frame.Draw()
            ROOT.gPad.SaveAs("bdt_score_comparison_datasets.png")
        
        
        
        def Compare_BDT_Scores_MultipleSelections(self, datafile, categ, selection_list, labels):
            frame = ROOT.RooRealVar('bdt_cv', 'BDT Score', -1, 1).frame()
            colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kOrange+7]
            
            MiniTreeFile = ROOT.TFile.Open(datafile)
            tree = MiniTreeFile.Get(categ)
            
            bdt_cv = ROOT.RooRealVar("bdt_cv", "bdt_cv", -1.0, 1.0)
            tripletMass = ROOT.RooRealVar('tripletMass', '3#mu mass', 1.6, 2.0)
            isMC = ROOT.RooRealVar("isMC", "isMC", 0, 1000000)
            weight = ROOT.RooRealVar("weight", "weight", 0, 5)
            dimu_OS1 = ROOT.RooRealVar('dimu_OS1', 'dimu_OS1', 0, 5)
            dimu_OS2 = ROOT.RooRealVar('dimu_OS2', 'dimu_OS2', 0, 5)
            
            variables = ROOT.RooArgSet(tripletMass, bdt_cv, isMC, weight, dimu_OS1, dimu_OS2)
            
            for i, selection in enumerate(selection_list):
                #fullcut = "(fabs(dimu_OS1 - 1.020)>0.020)&&(fabs(dimu_OS2 - 1.020)>0.020)&&(fabs(dimu_OS1 - 0.782)>0.020)&&(fabs(dimu_OS2 - 0.782)>0.020)&&" + selection
                fullcut = "{}".format(selection)
                
                dataset = ROOT.RooDataSet("sel_{0}".format(i), "sel_{0}".format(i), tree, variables, fullcut, "weight")
                
                dataset.plotOn(
                    frame,
                    ROOT.RooFit.Binning(100),
                    ROOT.RooFit.MarkerColor(colors[i % len(colors)]),
                    ROOT.RooFit.LineColor(colors[i % len(colors)]),
                    ROOT.RooFit.MarkerStyle(20 + i),
                    ROOT.RooFit.MarkerSize(0.75),
                    ROOT.RooFit.Name(labels[i])
                )
            
            frame.GetXaxis().SetTitle("BDT score")
            frame.SetTitle("BDT score: Selection Comparisons")
            frame.Draw()
            ROOT.gPad.BuildLegend()
            ROOT.gPad.SaveAs("bdt_score_comparison_selections.png")
        
        
        
        def Plot_Multiple_BDT_Shapes(self, root_files, tree_name, selections, labels, signalnorms):
            """
            Plot multiple BDT score distributions (from multiple datasets and selection cuts) on the same frame.
            
            Args:
                root_files   (list[str])   : Paths to ROOT files
                tree_name    (str)         : TTree name inside each ROOT file
                selections   (list[str])   : List of selection strings (one per curve)
                labels       (list[str])   : Labels for legend
                signalnorms  (list[float]) : Normalization scale factors for each curve
            """
            assert len(root_files) == len(selections) == len(labels) == len(signalnorms)
            
            # Variables
            bdt_cv = ROOT.RooRealVar("bdt_cv", "BDT Score", -1.0, 1.0)
            tripletMass = ROOT.RooRealVar('tripletMass', '3#mu mass', 1.6, 2.0)
            isMC = ROOT.RooRealVar("isMC", "isMC", 0, 1000000)
            weight = ROOT.RooRealVar("weight", "weight", 0, 5)
            scale = ROOT.RooRealVar("scale", "scale", 1.0)
            dimu_OS1 = ROOT.RooRealVar('dimu_OS1', 'dimu_OS1', 0, 2)
            dimu_OS2 = ROOT.RooRealVar('dimu_OS2', 'dimu_OS2', 0, 2)
            
            variables = ROOT.RooArgSet(bdt_cv, tripletMass, isMC, weight, scale, dimu_OS1, dimu_OS2)
            
            # Frame
            frame = bdt_cv.frame()
            frame.SetTitle("BDT Score Comparison")
            frame.GetXaxis().SetTitle("BDT score")
            
            colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kOrange+7]
            markers = [20, 21, 22, 23, 33, 34]
            
            legend = ROOT.TLegend(0.6, 0.65, 0.88, 0.88)
            
            for i in range(len(root_files)):
                f = ROOT.TFile.Open(root_files[i])
                tree = f.Get(tree_name)
                
                # Apply fixed vetoes + custom selection
                phiveto = "(fabs(dimu_OS1 - 1.020) > 0.020) && (fabs(dimu_OS2 - 1.020) > 0.020)"
                omegaveto = "(fabs(dimu_OS1 - 0.782) > 0.020) && (fabs(dimu_OS2 - 0.782) > 0.020)"
                basecut = "{} && {} && (tripletMass > 1.6 && tripletMass < 2.0)".format(phiveto, omegaveto)
                fullcut = "{} && ({})".format(basecut, selections[i])
                
                scale.setVal(signalnorms[i])
                
                dataset = ROOT.RooDataSet(
                    "ds_{0}".format(i), "ds_{0}".format(i),
                    tree, variables,
                    fullcut, "scale"
                )
                
                dataset.plotOn(
                    frame,
                    ROOT.RooFit.Binning(100),
                    ROOT.RooFit.MarkerColor(colors[i % len(colors)]),
                    ROOT.RooFit.LineColor(colors[i % len(colors)]),
                    ROOT.RooFit.MarkerStyle(markers[i % len(markers)]),
                    ROOT.RooFit.MarkerSize(0.75),
                    ROOT.RooFit.Name("curve_{0}".format(i))
                )
                
                legend.AddEntry("curve_{0}".format(i), labels[i], "lep")
                
            frame.Draw()
            legend.Draw()
            ROOT.gPad.SetTicks(1, 1)
            CMS_lumi(ROOT.gPad, 5, 0)
            ROOT.gPad.SaveAs("bdt_comparison_all.png")
            
            
        
        
        
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
                
                phivetoes="(fabs(dimu_OS1 - 1.020)>0.020)&(fabs(dimu_OS2 - 1.020)>0.020)"
                omegavetoes="&fabs(dimu_OS1 - 0.782)>0.020&fabs(dimu_OS2 - 0.782)>0.020&"
                
                
                # For fitting BDT Output in Data
                
                BDT_Score_Min=-0.3
                
                BlindDataSelector = RooFormulaVar('DataSelector', 'DataSelector', phivetoes+omegavetoes+' isMC == 0 & (tripletMass<=%s || tripletMass>=%s) & (tripletMass>=%s & tripletMass<=%s) ' %(signal_range_lo,signal_range_hi,fit_range_lo,fit_range_hi) , RooArgList(variables))
                
                fulldata = RooDataSet('data', 'data', tree,  variables, BlindDataSelector)
                
                self.bdt_cv.setRange("BDT_Fit_Range", BDT_Score_Min, 1.0);
                
                self.a = RooRealVar("a", "a", 1.0, 0.0, 10.0)
                self.b = RooRealVar("b", "b", 1.0, -10.0, 10.0)
                self.c = RooRealVar("c", "c", 1.0, -10.0, 10.0)
                self.d = RooRealVar("d", "d", 1.0, -10.0, 10.0)
                
                #quadratic = RooFormulaVar("quadratic", "a + b*self.bdt_cv + c*self.bdt_cv*self.bdt_cv", RooArgList(a, b, c, self.bdt_cv))
                #expModel = RooGenericPdf("expModel", "exp(quadratic)", RooArgList(quadratic)) #Exponential of the quadratic polynomial
                
                self.quadratic = RooFormulaVar("quadratic", "a + b*bdt_cv + c*bdt_cv*bdt_cv + d*bdt_cv*bdt_cv*bdt_cv", RooArgList(self.a, self.b, self.c, self.d, self.bdt_cv))
                self.expModel = RooGenericPdf("expModel", "exp(quadratic)", RooArgList(self.quadratic)) #Exponential of the cubic polynomial
                
                self.BDTNorm = RooRealVar("BDTNorm", "BDTNorm", 500.0, 0.1, 1000000000)
                self.BDT_distribution = RooAddPdf("BDT_distribution", "BDT_distribution",RooArgList(self.expModel), RooArgList(self.BDTNorm))
                #BDT_distribution = RooAddPdf("BDT_distribution", "BDT_distribution",RooArgList(quadratic), RooArgList(BDTNorm))
                
                results_pdf = self.BDT_distribution.fitTo(fulldata, RooFit.Range('BDT_Fit_Range'), RooFit.Save())
                results_pdf.Print()
                
                
                
                
                # For fitting BDT Output in Signal
                
                self.MCSelector = RooFormulaVar('MCSelector', 'MCSelector', phivetoes+omegavetoes+' isMC !=0 & (isMC == 211 | isMC == 210231 | isMC == 210232 | isMC == 210233 ) & (tripletMass>=%s & tripletMass<=%s) ' %(fit_range_lo,fit_range_hi) , RooArgList(variables))
                
                self.fullmc = RooDataSet('mc', 'mc', tree, variables, self.MCSelector,'scale')
                
                self.bdt_cv.setRange("BDT_MC_Fit_Range", -1.0, 1.0);
                
                self.bgausmeanMC = RooRealVar("bgausmeanMC", "bgausmeanMC", 0.5, 0.0, 0.9)
                self.bgaussigmaMC_a = RooRealVar("bgaussigmaMC_a", "bgaussigmaMC_a", 0.2, 0.000001, 1.0)
                self.bgaussigmaMC_b = RooRealVar("bgaussigmaMC_b", "bgaussigmaMC_b", 0.2, 0.000001, 1.0)
                
                self.bgaus_distMC = RooBifurGauss("bgaus_distMC", "bgaus dist MC", self.bdt_cv, self.bgausmeanMC, self.bgaussigmaMC_a, self.bgaussigmaMC_b)
                
                self.BDTNorm_MC = RooRealVar("BDTNorm_MC", "BDTNorm_MC", 500.0, 0.1, 50000)
                self.BDT_distribution_MC = RooAddPdf("BDT_distribution", "BDT_distribution",RooArgList(self.bgaus_distMC), RooArgList(self.BDTNorm_MC))
                
                results_mcpdf = self.BDT_distribution_MC.fitTo(self.fullmc, RooFit.Range('BDT_MC_Fit_Range'), RooFit.Save())
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
                CMS_lumi(ROOT.gPad, 5, 0)
                ROOT.gPad.Update()
                frame1.Draw('sameaxis')
                ROOT.gPad.SaveAs('bdt_fit_mc_'+categ+'.png')
                ROOT.gPad.Clear()
                
                ROOT.gPad.SetLogy()
                frame2.Draw()
                ROOT.gPad.SetTicks(1,1)
                CMS_lumi(ROOT.gPad, 5, 0)
                ROOT.gPad.Update()
                frame2.Draw('sameaxis')
                ROOT.gPad.SaveAs('bdt_fit_bkg_'+categ+'.png')
                ROOT.gPad.SetLogy(0)
                
                
                #print "Certain BDT cut: ", self.fullmc.reduce('bdt_cv > 0.5').sumEntries()
                
                
                return frame1, frame2, self.fullmc, fulldata
                
                
        
        



                        
if __name__ == "__main__":
    
        ROOT.gROOT.SetBatch(True)
        
        datafile = "../Combine_Tree_ztau3mutau.root"
        
        # Assume these methods belong to an object or class, e.g. BDTPlotter()
        BDTPlotter = BDT_Shape_Comparisons()  # Replace this with your actual class name or self instance
        
        # Categories to loop over
        categories = ['taue', 'taumu', 'tauhA', 'tauhB','all']
        #categories = ['taumu']
        #categories = ['taue','taumu','tauhA','tauhB','all']
        #categories = ['tauhA','tauhB','all']
        
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
            "ztautau_A.root",
            "ztautau_B.root",
            "ztautau_C.root"
        ]
        
        # Common tree name in all datasets
        tree_name = "ztau3mutauh_A"
        
        # For full combination (dataset + selection), define:
        combined_files = [
            "ztautau_A.root",
            "ztautau_B.root",
            "ztautau_A.root",
            "ztautau_C.root"
        ]
        combined_selections = [
            "bdt_cv > -0.3",
            "bdt_cv > 0.0",
            "bdt_cv > 0.2",
            "bdt_cv > 0.4"
        ]
        combined_labels = [
            "A: bdt > -0.3",
            "B: bdt > 0.0",
            "A: bdt > 0.2",
            "C: bdt > 0.4"
        ]
        combined_scales = [0.0000081] * 4  # Replace with correct normalizations
        
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
            # BDTPlotter.Compare_BDT_Scores_MultipleDatasets(dataset_files, tree)
        
            # 2. Compare different selections on the same file
            # Pick any one file as reference (e.g., A)
            BDTPlotter.Compare_BDT_Scores_MultipleSelections(datafile, tree, bdt_selections, selection_labels)
        
            # 3. Full comparison: different files + different selections
            #BDTPlotter.Plot_Multiple_BDT_Shapes(
            #    combined_files,
            #    tree,
            #    combined_selections,
            #    combined_labels,
            #    combined_scales
            #)
            