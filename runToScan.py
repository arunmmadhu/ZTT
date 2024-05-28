#!/usr/bin/env python

import os
import argparse
import ROOT
from ROOT import TFile
import math
import array





if __name__ == "__main__":


    phivetoes="\'(abs(dimu_OS1 - 1.020)>0.025)&(abs(dimu_OS2 - 1.020)>0.025)"
    omegavetoes="&abs(dimu_OS1 - 0.782)>0.04&fabs(dimu_OS2 - 1.020)>0.04"



#    bdt_cuts = [-0.4, -0.3, -0.2, -0.15, -0.05, 0.00, 0.05, 0.15,  0.2, 0.25, 0.30, 0.35, 0.40, 0.45]    
    bdt_cuts = [-0.4,  -0.2, -0.05,  0.00, 0.05,   0.2,  0.4]      # a few entries for test
    command=""
    for bdt_cut in bdt_cuts:

        command+="./makeTheCard.py --selection=" + phivetoes + omegavetoes + " &bdt_cv > %s\'"%bdt_cut   + " --category=\"taue\""  + "  --signalnorm=0.001"  + "  --bdt_point=%s"%bdt_cut + ";"

        command+="./makeTheCard.py --selection=" + phivetoes + omegavetoes + " &bdt_cv > %s\'"%bdt_cut   + " --category=\"taumu\"" + "  --signalnorm=0.001"  + "  --bdt_point=%s"%bdt_cut + ";"

        command+="./makeTheCard.py --selection=" + phivetoes + omegavetoes + " &bdt_cv > %s\'"%bdt_cut   + " --category=\"tauhA\"" + "  --signalnorm=0.001"  + "  --bdt_point=%s"%bdt_cut + ";"

        command+="./makeTheCard.py --selection=" + phivetoes + omegavetoes + " &bdt_cv > %s\'"%bdt_cut   + " --category=\"tauhB\"" + "  --signalnorm=0.001"  + "  --bdt_point=%s"%bdt_cut + ";"

    print command
#    os.system(command)
