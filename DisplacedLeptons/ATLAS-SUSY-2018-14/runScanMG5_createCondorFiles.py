#!/usr/bin/env python3

# 1) Run MadGraph using the options set in the input file 
# (the proc_card.dat, parameter_card.dat and run_card.dat...).

from __future__ import print_function
from copyreg import pickle
import sys,os,glob
from configParserWrapper import ConfigParserExt
from runScanMG5_helper import generateProcess
import logging
import subprocess
import tempfile
import time,datetime
import pickle
from typing import Set

FORMAT = '%(levelname)s in %(module)s.%(funcName)s(): %(message)s at %(asctime)s'
logging.basicConfig(format=FORMAT,datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("MG5Scan")
    
 
def generate_configFiles(parfile,verbose) -> Set[str]:
   
    level = verbose
    levels = { "debug": logging.DEBUG, "info": logging.INFO,
               "warn": logging.WARNING,
               "warning": logging.WARNING, "error": logging.ERROR }
    if not level in levels:
        logger.error ( "Unknown log level ``%s'' supplied!" % level )
        sys.exit()
    logger.setLevel(level = levels[level])    

    parser = ConfigParserExt(inline_comment_prefixes="#")   
    ret = parser.read(parfile)
    if ret == []:
        logger.error( "No such file or directory: '%s'" % args.parfile)
        sys.exit()
            
    #Get a list of parsers (in case loops have been defined)    
    parserList = parser.expandLoops()
    
    now = datetime.datetime.now()
    configFolders = set([])
    for irun,newParser in enumerate(parserList):
        processFolder = newParser.get('MadGraphPars','processFolder')
        processFolder = os.path.abspath(processFolder)
        configFolder = os.path.join(processFolder,'config')
        if processFolder[-1] == '/':
            processFolder = processFolder[:-1]
        if not os.path.isdir(processFolder):
            logger.info('Folder %s not found. Running MG5 to create folder.' %processFolder)
            generateProcess(newParser)
        if not os.path.isdir(configFolder):
            os.makedirs(configFolder, exist_ok=True)


        # Get largest existing events folder:
        run0 = 1
        eventsFolder = os.path.join(processFolder,'Events')
        if os.path.isdir(eventsFolder):
            for runF in glob.glob(os.path.join(eventsFolder,'run*')):
                run0 = max(run0,int(os.path.basename(runF).replace('run_',''))+1)

        # Create temporary folder names if running in parallel
        if len(parserList) > 1:
            # Create temporary folders
            runFolder = tempfile.mkdtemp(prefix='%s_'%(processFolder),suffix='_run_%02d' %(run0+irun))
            os.removedirs(runFolder)
        else:
            runFolder = processFolder

        newParser.set('MadGraphPars','runFolder',runFolder)
        newParser.set('MadGraphPars','runNumber','%02d' %(run0+irun))

        parserDict = newParser.toDict(raw=False)
        # Create config file to store the dictionary of parameters used for this run:
        outfile = f"{configFolder}/job_{irun:05d}.pkl"
        with open(outfile, "wb") as f:
            pickle.dump(parserDict, f)
            
        configFolders.add(configFolder)
    logger.info(f"Created {len(parserList)} config files at {now.strftime('%Y-%m-%d %H:%M')}")

    return configFolders
    

def generate_condorSubmitFile(configFolder,subFile,worker_file='runScanMG5_worker.py'):

    worker = os.path.abspath(worker_file)
    if not os.path.isfile(worker):
        logger.error(f"Worker file {worker} not found. Make sure the file exists and the path is correct.")
        sys.exit()
    
    configFolder = os.path.abspath(configFolder)
    submitFile = os.path.abspath(os.path.join(configFolder,subFile))
    with open(submitFile, 'w') as f:
        f.write(f"executable = python3\n")
        f.write(f"arguments =  {worker} -c $(config)\n")
        f.write("getenv = True\n")
        f.write("request_memory = 2GB\n")
#        f.write("request_cpus = 1\n")
        f.write("output = condor_output/job.$(Cluster).$(Process).out\n")
        f.write("error = condor_output/job.$(Cluster).$(Process).err\n")
        f.write("log = condor_output/job.$(Cluster).$(Process).log\n")
        f.write(f"queue config matching {configFolder}/*pkl'\n")

    return submitFile


if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser( description=
            "Creates a series of config files for running MadGraph scans with condor." )
    ap.add_argument('-o', '--outputFile', default='runScanMG5_condor.sub',
            help='name for the the condor submit file [runScanMG5_condor.sub].')
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
    configFolders = generate_configFiles(args.parfile,args.verbose)
    if len(configFolders) == 0:
        logger.error("No config files created. Exiting.")
        sys.exit()
    elif len(configFolders) > 1:
        logger.warning(f"Multiple config folders created: {configFolders}. Make sure to submit condor jobs for all folders.")
        
    for configFolder in configFolders:
        subFile = generate_condorSubmitFile(configFolder,args.outputFile)
        logger.info(f"Submit file {subFile} created for config folder {configFolder}")
            
    print("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
