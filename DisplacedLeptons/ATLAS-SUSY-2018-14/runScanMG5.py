#!/usr/bin/env python3

from __future__ import print_function
import sys,os,glob
from configParserWrapper import ConfigParserExt
from runScanMG5_helper import generateInputFiles,logger
from runScanMG5_worker import runSingleJob
import logging
import subprocess
import multiprocessing
import time,datetime

FORMAT = '%(levelname)s in %(module)s.%(funcName)s(): %(message)s at %(asctime)s'
logging.basicConfig(format=FORMAT,datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("MG5Scan")

def main(parfile,verbose):

    level = verbose
    levels = { "debug": logging.DEBUG, "info": logging.INFO,
               "warn": logging.WARNING,
               "warning": logging.WARNING, "error": logging.ERROR }
    if level in levels:       
        logger.setLevel(level = levels[level])
   
    parser = ConfigParserExt(inline_comment_prefixes="#")   
    ret = parser.read(parfile)
    if ret == []:
        logger.error( f"No such file or directory: {parfile}")
        sys.exit(1)

    scanFolders = generateInputFiles(parfile)
    if len(scanFolders) == 0:
        logger.error( f"No valid input files generated from {parfile}.")
        sys.exit(1)

    # Get a list of all generated input config files
    allInputFiles = []
    for folderTuple in scanFolders:
        allInputFiles += list(glob.glob(os.path.join(folderTuple.inputFolder, '*.ini')))

    # Start multiprocessing pool
    ncpus = -1
    if parser.has_option("options","ncpu"):
        ncpus = int(parser.get("options","ncpu"))
    if ncpus  < 0:
        ncpus =  multiprocessing.cpu_count()
    ncpus = min(ncpus,len(allInputFiles))
    pool = multiprocessing.Pool(processes=ncpus)
    if ncpus > 1:
        logger.info('Running %i jobs in parallel with %i processes' %(len(allInputFiles),ncpus))
    else:
        logger.info('Running %i jobs in series with a single process' %(len(allInputFiles)))

    now = datetime.datetime.now()
    children = []
    for inputFile in allInputFiles:
        logger.debug(f'submitting with input file: {inputFile}')
        p = pool.apply_async(runSingleJob, args=(inputFile,verbose,))
        children.append(p)

#     Wait for jobs to finish:
    output = [p.get() for p in children]
    logger.info("Finished all runs (%i) at %s" %(len(allInputFiles),now.strftime("%Y-%m-%d %H:%M")))

    return output

if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser( description=
            "Run a (serial) MadGraph scan for the parameters defined in the parameters file." )
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
        logger.error('Environment variables not properly set. Run source setenv.sh first.')
        sys.exit()

    t0 = time.time()

    args = ap.parse_args()
    output = main(args.parfile,args.verbose)
            
    logger.info("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
