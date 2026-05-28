#!/usr/bin/env python3

# 1) Run MadGraph using the options set in the input file 
# (the proc_card.dat, parameter_card.dat and run_card.dat...).

from __future__ import print_function
import sys,os
from configParserWrapper import ConfigParserExt
from runScanMG5_helper import generateEvents,logger
import time, logging


def compressOutputFolder(outputFolder):
    import shutil
    compactOutputFolder = outputFolder + '.tar.gz'
    shutil.make_archive(outputFolder, 'gztar', outputFolder)
    return compactOutputFolder

if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser( description=
            "Runs a single job defined by a config file for running MadGraph scans with condor." )
    ap.add_argument('-c', '--configfile', required=True,
            help='path to the config file for this job, i.e. config.ini')
    ap.add_argument('-v', '--verbose', required=False, default='info',
            help='verbose level [info]')
    
    t0 = time.time()
    args = ap.parse_args()
    configFile = args.configfile

    level = args.verbose
    levels = { "debug": logging.DEBUG, "info": logging.INFO,
               "warn": logging.WARNING,
               "warning": logging.WARNING, "error": logging.ERROR }
    if level in levels:       
        logger.setLevel(level = levels[level])

    parser = ConfigParserExt(inline_comment_prefixes="#")
    ret = parser.read(configFile)
    if ret == []:
        print(f"Could not read config file: {configFile}")
        sys.exit(1)
    parserDict = parser.toDict(raw=False,abspath_existing=True)

    logger.info(f"Running job with config file: {configFile}")
    output = generateEvents(parserDict)
    runFolder = os.path.abspath(output['runFolder'])
    compactOutputFolder = compressOutputFolder(runFolder)
    shutil.rmtree(runFolder)
    logger.info(f"Output folder {runFolder} compressed to {compactOutputFolder}")
    logger.info(f"Done in {(time.time()-t0)/60.0:.2f} minutes.")
