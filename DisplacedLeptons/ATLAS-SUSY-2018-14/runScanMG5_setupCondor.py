#!/usr/bin/env python3

# 1) Run MadGraph using the options set in the input file 
# (the proc_card.dat, parameter_card.dat and run_card.dat...).

from __future__ import print_function
import sys,os
from runScanMG5_helper import generateInputFiles,generateCondorFile,logger
import logging
import subprocess
import time
 

if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser( description=
            "Creates a series of config files for running MadGraph scans with condor." )
    ap.add_argument('-o', '--outputFile', default='runScanMG5_condor.sub',
            help='name for the the condor submit file [runScanMG5_condor.sub].')
    ap.add_argument('-w', '--workerFile', default='runScanMG5_worker.py',
            help='path to the the worker file [runScanMG5_worker.py].')
    ap.add_argument('-p', '--parfile', default='scan_parameters.ini',
            help='path to the parameters file [scan_parameters.ini].')
    ap.add_argument('-v', '--verbose', default='info',
            help='verbose level (debug, info, warning or error). Default is info')

    # First make sure the correct env variables have been set:
    LDPATH = subprocess.check_output('echo $LD_LIBRARY_PATH',shell=True,text=True)
    ROOTINC = subprocess.check_output('echo $ROOT_INCLUDE_PATH',shell=True,text=True)
    pythiaDir = os.path.abspath('./MG5/HEPTools/pythia8/lib')
    delphesDir = os.path.abspath('./DelphesLLP/external')
    if pythiaDir not in LDPATH or delphesDir not in ROOTINC:
        print('Enviroment variables not properly set. Run source setenv.sh first.')
        sys.exit()
    

    t0 = time.time()

    args = ap.parse_args()

    level = args.verbose
    levels = { "debug": logging.DEBUG, "info": logging.INFO,
               "warn": logging.WARNING,
               "warning": logging.WARNING, "error": logging.ERROR }
    if level in levels:       
        logger.setLevel(level = levels[level])

    scanFolders = generateInputFiles(args.parfile)
    if len(scanFolders) == 0:
        logger.error("No input files created. Exiting.")
        sys.exit()
    elif len(scanFolders) > 1:
        logger.warning(f"Multiple scan folders created: {scanFolders}. Make sure to submit condor jobs for all folders.")
        
    for folderTuple in scanFolders:
        inputFolder = folderTuple.inputFolder
        resultsFolder = folderTuple.outputFolder
        subFile = generateCondorFile(inputFolder,resultsFolder,args.outputFile,args.workerFile,args.verbose)
        logger.info(f"Submit file {subFile} created for config folder {inputFolder}")
            
    print("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
