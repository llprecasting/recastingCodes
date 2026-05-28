#!/usr/bin/env python3

# 1) Run MadGraph using the options set in the input file 
# (the proc_card.dat, parameter_card.dat and run_card.dat...).

from __future__ import print_function
import sys
from configParserWrapper import ConfigParserExt
from runScanMG5_helper import generateEvents, moveFolders
import time


if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser( description=
            "Runs a single job defined by a config file for running MadGraph scans with condor." )
    ap.add_argument('-c', '--configfile', required=True,
            help='path to the config file for this job, i.e. config.ini')
    
    t0 = time.time()
    args = ap.parse_args()
    configFile = args.configfile

    parser = ConfigParserExt(inline_comment_prefixes="#")
    ret = parser.read(configFile)
    if ret == []:
        print(f"Could not read config file: {configFile}")
        sys.exit(1)
    parserDict = parser.toDict(raw=False,abspath_existing=True)

    output = generateEvents(parserDict)

    moveFolders(output)

    print("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
