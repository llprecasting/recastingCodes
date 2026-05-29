#!/usr/bin/env python3

# 1) Run MadGraph using the options set in the input file 
# (the proc_card.dat, parameter_card.dat and run_card.dat...).

from __future__ import print_function
import os
from runScanMG5_helper import generateEvents,logger,moveFolders
import time, logging
import shutil
from typing import Dict

def compressOutputFolder(output : Dict) -> None:
    """
    Compresses the output folder of a run and removes the original folder.
    """
    outputFolder = os.path.abspath(output['runFolder'])

    runFolder = output['runFolder']
    shutil.rmtree(runFolder)
    
    compactOutputFolder = outputFolder + '.tar.gz'
    shutil.make_archive(outputFolder, 'gztar', outputFolder)
    logger.info(f"Output folder {runFolder} compressed to {compactOutputFolder}")

def runSingleJob(configFile : str,verbose : str) -> None:
    
    t0 = time.time()
    levels = { "debug": logging.DEBUG, "info": logging.INFO,
               "warn": logging.WARNING,
               "warning": logging.WARNING, "error": logging.ERROR }
    if verbose in levels:       
        logger.setLevel(level = levels[verbose])

    logger.info(f"Running job with config file: {configFile}")
    output = generateEvents(configFile)
    try:
        moveFolders(output)
    except Exception as e:
        compressOutputFolder(output)
    logger.info(f"Done in {(time.time()-t0)/60.0:.2f} minutes.")
    
if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser( description=
            "Runs a single job defined by a config file for running MadGraph scans with condor. After completion, the output folder will be compressed and removed." )
    ap.add_argument('-c', '--configfile', required=True,
            help='path to the config file for this job, i.e. config.ini')
    ap.add_argument('-v', '--verbose', required=False, default='info',
            help='verbose level [info]')
    
    args = ap.parse_args()
    runSingleJob(args.configfile,args.verbose)

    
