#!/usr/bin/env python3

# 1) Run MadGraph using the options set in the input file 
# (the proc_card.dat, parameter_card.dat and run_card.dat...).

from __future__ import print_function
import sys,os,glob
from configParserWrapper import ConfigParserExt
import logging,shutil
import subprocess
import multiprocessing
import tempfile
import time,datetime

FORMAT = '%(levelname)s in %(module)s.%(funcName)s(): %(message)s at %(asctime)s'
logging.basicConfig(format=FORMAT,datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("MG5Scan")

def generateProcess(parser):
    """
    Runs the madgraph process generation.
    This step just need to be performed once for a given
    model, since it is independent of the 
    numerical values of the model parameters.
    
    :param parser: Dictionary with parser sections.
    
    :return: True if successful. Otherwise False.
    """
    
    
    #Get run folder:    
    pars = parser["MadGraphPars"]
    processCard = os.path.abspath(pars["proccard"])    
    if not os.path.isfile(processCard):
        logger.error("Process card %s not found" %processCard)
        raise ValueError()

    processFolder = os.path.abspath(pars["processFolder"])
    if os.path.isdir(processFolder):
        logger.warning("Process folder %s found. Skipping process generation." %processFolder)
        return False

    logger.info('Generating process using %s' %processCard)

    # Create copy of process card to replace output folder
    procCard = tempfile.mkstemp(suffix='.dat', prefix='procCard_')
    os.close(procCard[0])
    procCard = procCard[1]
    shutil.copy(processCard,procCard)
    with open(procCard,'r') as f:
        lines = f.readlines()
    lines = [l for l in lines[:] if l.strip()[:6] != 'output']
    lines.append('output %s\n' %processFolder)
    with open(procCard,'w') as f:
        for l in lines:
            f.write(l)
    
    #Generate process
    mg5Folder = os.path.abspath('./MG5')
    run = subprocess.Popen('./bin/mg5_aMC -f %s' %procCard,shell=True,
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                                cwd=mg5Folder)
         
    output,errorMsg = run.communicate()
    logger.debug('MG5 process error:\n %s \n' %errorMsg)
    logger.debug('MG5 process output:\n %s \n' %output)
    logger.info("Finished process generation")

    os.remove(procCard)
        
    return True

def getInfoFromOutput(outputStr):
    """
    Try to fetch the run name, the number of events
    and the total cross-section from the MG5 output.

    :param outputStr: String containing the MG5 output

    :return: Dictionary with run summary
    """

    outputStr = str(outputStr)
    # Get summary block
    summaryBlock = outputStr.split("Results Summary for run")[-1]
    summaryBlock = summaryBlock.split("Done")[0]
    # Split block into lines
    summaryLines = [l for l in summaryBlock.split('\\n') if l.strip()]
    # Get info from lines:
    xsec, runNumb, runTag, nevts = None,None,None,None
    for l in summaryLines:
        if 'tag' in l:
            runNumb,runTag = l.split('tag:')
            runNumb = runNumb.replace(':','').replace('=','').strip()
            runTag = runTag.replace(':','').replace('=','').strip()
        elif 'cross-section' in l.lower():
            xsec = eval(l.lower().split(':')[1].split('+-')[0])
        elif 'nb of events'in l.lower():
            nevts = eval(l.lower().split(':')[1])

    runInfo = {'run number' : runNumb, 'run tag' : runTag, 
               'cross-section (pb)' : xsec, 'Number of events' : int(nevts)}
    return runInfo

def runMG5(parser):
    
    """
    Runs the madgraph event generation.
    
    :param parser: Dictionary with parser sections.
    
    :return: Dictionary with run info. False if failed.
    """

    pars = parser["MadGraphPars"]
    if not 'runFolder' in pars:
        logger.error('Run folder not defined.')
        return False        
    else:
        runFolder = pars['runFolder']

    if not 'processFolder' in pars:
        logger.error('Process folder not defined.')
        return False        
    else:
        processFolder = pars['processFolder']
        if not os.path.isdir(processFolder):
            logger.error('Process folder %s not found.' %processFolder)
            return False            

    # If run folder does not exist, create it using processFolder as a template:
    if not os.path.isdir(runFolder):
        runFolder = shutil.copytree(processFolder,runFolder,
                                    ignore=shutil.ignore_patterns('Events','*.lhe'),
                                    symlinks=True)
        os.makedirs(os.path.join(runFolder,'Events'))       
        logger.info("Created temporary folder %s" %runFolder) 
            
    if not os.path.isdir(runFolder):
        logger.error('Run folder %s not found.' %runFolder)
        return False

    if 'runcard' in pars:
        if os.path.isfile(pars['runcard']):    
            shutil.copyfile(pars['runcard'],os.path.join(runFolder,'Cards/run_card.dat'))
        else:
            raise ValueError("Run card %s not found" %pars['runcard'])
    if 'paramcard' in pars:
        if os.path.isfile(pars['paramcard']):
            shutil.copyfile(pars['paramcard'],os.path.join(runFolder,'Cards/param_card.dat'))    
        else:
            raise ValueError("Param card %s not found" %pars['paramcard'])
    
    #Generate commands file:       
    commandsFile = tempfile.mkstemp(suffix='.txt', prefix='MG5_commands_', dir=runFolder)
    os.close(commandsFile[0])
    commandsFileF = open(commandsFile[1],'w')
    commandsFileF.write('shower=OFF\n')
    commandsFileF.write('detector=OFF\n')
    commandsFileF.write('done\n')
    comms = parser["MadGraphSet"]
    # Set the MadGraph parameters defined in the ini file
    for key,val in comms.items():
        commandsFileF.write('set %s %s\n' %(key,val))

    commandsFileF.write('set iseed %s\n' %(parser["MadGraphPars"]["runNumber"]))
    #Done setting up options
    commandsFileF.write('done\n')
    commandsFileF.write('done\n')

    commandsFileF.close()
    commandsFile = commandsFile[1]      

    ncore = parser['options']['ncore']
    
    logger.debug("Generating MG5 events with command file %s" %commandsFile)
    run = subprocess.Popen('./bin/generate_events --multicore --nb_core=%i < %s' %(ncore,commandsFile),
                           shell=True,stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,cwd=runFolder)
      
    output,errorMsg= run.communicate()
    runInfo = {}
    runInfo.update(pars)

    logger.debug('MG5 event error:\n %s \n' %errorMsg)
    logger.debug('MG5 event output:\n %s \n' %output)    
    # Get info from output
    runInfo.update(getInfoFromOutput(output))

    os.remove(commandsFile)

    return runInfo

def runDelphesPythia8(parser,runInfo):
    """
    Runs the Pythia8 and Delphes directly (bypassing hepmc generation)
    on the LHE file.
    
    :param parser: Dictionary with parser sections.
    :param runInfo: Dictionary with MadGraph run information.
    
    :return: Dictionary with run info. False if failed.
    """

    cleanOutput = parser['options']['cleanOutput']

    pars = parser["DelphesPars"]
    delphesDir = os.path.abspath(pars['delphesDir'])
    delphescard = os.path.abspath(pars['delphescard'])
    pythiaCard = os.path.abspath(pars['pythia8card'])
    
    runFolder = parser['MadGraphPars']['runFolder']
    nevts = parser["MadGraphSet"]["nevents"]
    lheFile = os.path.join(runFolder,'Events',runInfo['run number'],'unweighted_events.lhe.gz')
    rootFile = os.path.join(runFolder,'Events',runInfo['run number'], '%s_delphes_events.root'  %runInfo['run tag'])

    #Generate pythia card with additional info
    pythia_tmp = tempfile.mkstemp(suffix='.dat', prefix='pythia_card_', dir=delphesDir)
    os.close(pythia_tmp[0])
    pythia_tmp = pythia_tmp[1]
    shutil.copyfile(pythiaCard,pythia_tmp)
    with open(pythia_tmp,'a') as f:
        f.write('\n\n')
        f.write('### Added lines:\n')
        f.write('Beams:frameType = 4\n')
        f.write('Beams:LHEF = %s\n' %lheFile)
        f.write('Main:numberOfEvents   = %i\n' %nevts)

    logger.debug("Running DelphesPythia8 with config files %s and %s" %(pythia_tmp,delphescard))
    run = subprocess.Popen('./DelphesPythia8 %s %s %s' %(delphescard,pythia_tmp,rootFile),shell=True,
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,
                                cwd=delphesDir)
    output,errorMsg = run.communicate()

    runInfo.update({'DelphesOutput' : output, 'DelphesError' : errorMsg})

    os.remove(pythia_tmp)

    #Generate pythia log file
    pythia_log = rootFile.replace('_delphes_events.root','_pythia_delphes.log')
    with open(pythia_log,'w') as f:
        f.write(output)

    if cleanOutput:
        logger.debug('Removing  %s' %lheFile)
        if os.path.isfile(lheFile):
            os.remove(lheFile)
        

    return runInfo

def generateEvents(parser):
    """
    Run MadGraph5 to generate LHE events and then DelphesPythia8.

    :param parser: Dictionary with parser sections.
    
    :return: Dictionary with run info. False if failed.
    """

    t0 = time.time()
    if parser["options"]["runMadGraph"]:
        runInfo = runMG5(parser)
        if parser["options"]["runPythiaDelphes"]:
            runInfo = runDelphesPythia8(parser,runInfo)
    elif parser["options"]["runPythiaDelphes"]:
        logger.error("Pythia and Delphes can only run if runMadGraph = True")
        return {}
    
    runInfo.update({'time (s)' : time.time()-t0})
    
    return runInfo

def moveFolders(runInfo):
    """
    Move the run folders from the temporary running folder
    to the process folder.
    """

    logger.info('Finished event generation for run %i in %1.2f min' 
                %(int(runInfo['runNumber']),runInfo['time (s)']/60.))

    # Get run folder:
    runFolder = os.path.abspath(runInfo['runFolder'])
    runNumber = int(runInfo['runNumber'])
    processFolder = os.path.abspath(runInfo['processFolder'])
    # If run folder and process folder are the same, do nothing
    if runFolder == processFolder:
        return
    
    # Move run folder results to process folder
    eventFolder = os.path.join(processFolder,'Events')
    runDirs = list(glob.glob(os.path.join(runFolder,'Events','run_*')))
    if len(runDirs) != 1:
        logger.error('Something went wrong. Found %i run folders in %s' %(len(runDirs),runFolder))
        return False
    runDir = runDirs[0]
    finalRunDir = os.path.join(eventFolder,'run_%02d' %runNumber)
    logger.info('Moving %s to %s' %(runDir,finalRunDir))
    shutil.move(runDir,finalRunDir)
    logger.info('Deleting temporary folder %s' %(runFolder))
    shutil.rmtree(runFolder)


def main(parfile,verbose):
   
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

    # Start multiprocessing pool
    ncpus = -1
    if parser.has_option("options","ncpu"):
        ncpus = int(parser.get("options","ncpu"))
    if ncpus  < 0:
        ncpus =  multiprocessing.cpu_count()
    ncpus = min(ncpus,len(parserList))
    pool = multiprocessing.Pool(processes=ncpus)
    if ncpus > 1:
        logger.info('Running %i jobs in parallel with %i processes' %(len(parserList),ncpus))
    else:
        logger.info('Running %i jobs in series with a single process' %(len(parserList)))

    now = datetime.datetime.now()
    children = []
    for irun,newParser in enumerate(parserList):
        processFolder = newParser.get('MadGraphPars','processFolder')
        processFolder = os.path.abspath(processFolder)
        if processFolder[-1] == '/':
            processFolder = processFolder[:-1]
        if not os.path.isdir(processFolder):
            logger.info('Folder %s not found. Running MG5 to create folder.' %processFolder)
            generateProcess(newParser)

        # Get largest existing events folder:
        run0 = 1
        eventsFolder = os.path.join(processFolder,'Events')
        if os.path.isdir(eventsFolder):
            for runF in glob.glob(os.path.join(eventsFolder,'run*')):
                run0 = max(run0,int(os.path.basename(runF).replace('run_',''))+1)

        # Create temporary folder names if running in parallel
        if ncpus > 1:
            # Create temporary folders
            runFolder = tempfile.mkdtemp(prefix='%s_'%(processFolder),suffix='_run_%02d' %(run0+irun))
            os.removedirs(runFolder)
        else:
            runFolder = processFolder

        newParser.set('MadGraphPars','runFolder',runFolder)
        newParser.set('MadGraphPars','runNumber','%02d' %(run0+irun))

        parserDict = newParser.toDict(raw=False)
        logger.debug('submitting with pars:\n %s \n' %parserDict)
        p = pool.apply_async(generateEvents, args=(parserDict,), callback=moveFolders)                       
        children.append(p)

#     Wait for jobs to finish:
    output = [p.get() for p in children]
    logger.info("Finished all runs (%i) at %s" %(len(parserList),now.strftime("%Y-%m-%d %H:%M")))

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
        print('Enviroment variables not properly set. Run source setenv.sh first.')
        sys.exit()

    t0 = time.time()

    args = ap.parse_args()
    output = main(args.parfile,args.verbose)
            
    print("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
