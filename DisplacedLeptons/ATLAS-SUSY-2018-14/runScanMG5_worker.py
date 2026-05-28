#!/usr/bin/env python3

# 1) Run MadGraph using the options set in the input file 
# (the proc_card.dat, parameter_card.dat and run_card.dat...).

from __future__ import print_function
import sys,os
from runScanMG5_helper import generateEvents, moveFolders
import subprocess
import time
import pickle


if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser( description=
            "Runs a single job defined by a config file for running MadGraph scans with condor." )
    ap.add_argument('-c', '--configfile', required=True,
            help='path to the config file for this job, i.e. config.pkl')
    
    t0 = time.time()
    args = ap.parse_args()
    configFile = args.configfile

    with open(configFile, "rb") as f:
        parserDict = pickle.load(f)

    output = generateEvents(parserDict)

    moveFolders(output)

    print("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
