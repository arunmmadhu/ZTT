#!/usr/bin/env python

import os
import argparse
import ROOT
from ROOT import TFile
import math
import array
import numpy as np

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


if __name__ == "__main__":


    phivetoes="\'(fabs(dimu_OS1 - 1.020)>0.020)&(fabs(dimu_OS2 - 1.020)>0.020)"
    omegavetoes="&fabs(dimu_OS1 - 0.782)>0.020&fabs(dimu_OS2 - 0.782)>0.020"




    bdt_cuts_all = generate_bdt_cuts(0.1, 0.8, 0.55)
    bdt_cuts_taue = generate_bdt_cuts(-0.2, 0.6, 0.2)
    bdt_cuts_taumu = generate_bdt_cuts(-0.1, 0.6, 0.4)
    bdt_cuts_tauhA = generate_bdt_cuts(-0.2, 0.6, 0.2)
    bdt_cuts_tauhB = generate_bdt_cuts(-0.4, 0.55, 0.2)
    
    open('Slopes_all.txt', 'w').close()
    open('Slopes_taue.txt', 'w').close()
    open('Slopes_taumu.txt', 'w').close()
    open('Slopes_tauhA.txt', 'w').close()
    open('Slopes_tauhB.txt', 'w').close()

    command=""
    
    for bdt_cut in bdt_cuts_taue:
        command+="./makeTheCard.py --selection=" + phivetoes + omegavetoes + " &bdt_cv > %s\'"%bdt_cut   + " --category=\"taue\""  + "  --signalnorm=0.00000856928"  + "  --bdt_point=%s"%bdt_cut + "  --outdir=unfixed_slope " + ";"
        
        #command+="./makeTheCard.py --selection=" + phivetoes + omegavetoes + " &bdt_cv > %s\'"%bdt_cut   + " --category=\"taue\""  + "  --signalnorm=0.00000856928"  + "  --bdt_point=%s"%bdt_cut + " --alt_pdf "  + "  --pdf_switch_point=-0.02 " + "  --fixed_slope=-0.29 " + "  --outdir=fixed_slope " + ";"
        #pass
       
    for bdt_cut in bdt_cuts_taumu:
        command+="./makeTheCard.py --selection=" + phivetoes + omegavetoes + " &bdt_cv > %s\'"%bdt_cut   + " --category=\"taumu\"" + "  --signalnorm=0.00000822810"  + "  --bdt_point=%s"%bdt_cut + "  --outdir=unfixed_slope " + ";"
        
        #command+="./makeTheCard.py --selection=" + phivetoes + omegavetoes + " &bdt_cv > %s\'"%bdt_cut   + " --category=\"taumu\"" + "  --signalnorm=0.00000822810"  + "  --bdt_point=%s"%bdt_cut + " --alt_pdf "  + "  --pdf_switch_point=0.19 " + "  --fixed_slope=0.86 " + "  --outdir=fixed_slope " + ";"
        #pass
    
    for bdt_cut in bdt_cuts_tauhA:
        command+="./makeTheCard.py --selection=" + phivetoes + omegavetoes + " &bdt_cv > %s\'"%bdt_cut   + " --category=\"tauhA\"" + "  --signalnorm=0.00000815958"  + "  --bdt_point=%s"%bdt_cut + "  --outdir=unfixed_slope " + ";"
        
        #command+="./makeTheCard.py --selection=" + phivetoes + omegavetoes + " &bdt_cv > %s\'"%bdt_cut   + " --category=\"tauhA\"" + "  --signalnorm=0.00000815958"  + "  --bdt_point=%s"%bdt_cut + " --alt_pdf "  + "  --pdf_switch_point=0.12 " + "  --fixed_slope=1.86 " + "  --outdir=fixed_slope " + ";"
        #pass
        
    for bdt_cut in bdt_cuts_tauhB:
        command+="./makeTheCard.py --selection=" + phivetoes + omegavetoes + " &bdt_cv > %s\'"%bdt_cut   + " --category=\"tauhB\"" + "  --signalnorm=0.00000815958"  + "  --bdt_point=%s"%bdt_cut + "  --outdir=unfixed_slope " + ";"
        
        #command+="./makeTheCard.py --selection=" + phivetoes + omegavetoes + " &bdt_cv > %s\'"%bdt_cut   + " --category=\"tauhB\"" + "  --signalnorm=0.00000815958"  + "  --bdt_point=%s"%bdt_cut + " --alt_pdf "  + "  --pdf_switch_point=-0.16 " + "  --fixed_slope=0.39 " + "  --outdir=fixed_slope " + ";"
        #pass
        
    for bdt_cut in bdt_cuts_all:
        command+="./makeTheCard.py --selection=" + phivetoes + omegavetoes + " &bdt_cv > %s\'"%bdt_cut   + " --category=\"all\""  + "  --signalnorm=0.00000824176"  + "  --bdt_point=%s"%bdt_cut + "  --outdir=unfixed_slope " + ";"
        
        #command+="./makeTheCard.py --selection=" + phivetoes + omegavetoes + " &bdt_cv > %s\'"%bdt_cut   + " --category=\"all\""  + "  --signalnorm=0.00000824176"  + "  --bdt_point=%s"%bdt_cut + " --alt_pdf "  + "  --pdf_switch_point=0.57 " + "  --fixed_slope=-0.11 " + "  --outdir=fixed_slope " + ";"
        #pass
    
    print command
    os.system(command)
