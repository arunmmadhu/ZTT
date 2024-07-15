#!/usr/bin/env python

import os
import argparse
import ROOT
from ROOT import TFile
import math
import array





if __name__ == "__main__":


    phivetoes="\'(fabs(dimu_OS1 - 1.020)>0.020)&(fabs(dimu_OS2 - 1.020)>0.020)"
    omegavetoes="&fabs(dimu_OS1 - 0.782)>0.020&fabs(dimu_OS2 - 0.782)>0.020"




    bdt_cuts = [-0.4,  -0.2,  0.00, 0.1, 0.2, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39,  0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.6, 0.7]   

    command=""
    for bdt_cut in bdt_cuts:

        command+="./makeTheCard.py --selection=" + phivetoes + omegavetoes + " &bdt_cv > %s\'"%bdt_cut   + " --category=\"all\""  + "  --signalnorm=0.00000712772"  + "  --bdt_point=%s"%bdt_cut + ";"
        
        command+="./makeTheCard.py --selection=" + phivetoes + omegavetoes + " &bdt_cv > %s\'"%bdt_cut   + " --category=\"taue\""  + "  --signalnorm=0.00000741099"  + "  --bdt_point=%s"%bdt_cut + ";"

        command+="./makeTheCard.py --selection=" + phivetoes + omegavetoes + " &bdt_cv > %s\'"%bdt_cut   + " --category=\"taumu\"" + "  --signalnorm=0.00000711593"  + "  --bdt_point=%s"%bdt_cut + ";"

        command+="./makeTheCard.py --selection=" + phivetoes + omegavetoes + " &bdt_cv > %s\'"%bdt_cut   + " --category=\"tauhA\"" + "  --signalnorm=0.00000705558"  + "  --bdt_point=%s"%bdt_cut + ";"

        command+="./makeTheCard.py --selection=" + phivetoes + omegavetoes + " &bdt_cv > %s\'"%bdt_cut   + " --category=\"tauhB\"" + "  --signalnorm=0.00000705558"  + "  --bdt_point=%s"%bdt_cut + ";"

    print command
    os.system(command)
