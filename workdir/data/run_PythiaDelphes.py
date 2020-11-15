#!/usr/bin/env python3

import os,time,subprocess
import logging

FORMAT = '%(levelname)s in %(module)s: %(message)s at %(asctime)s'
logging.basicConfig(format=FORMAT,datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger(__name__)


def runAll(args):

    #Run pythia
    logger.info('Running main_pythia2HepMC.exe')
    outputfile = os.path.splitext(args.inputfile)[0]+'.hep'
    logger.debug('./main_pythia2HepMC.exe -f %s -c %s -o %s -n %i' %(args.inputfile,args.pythiacfg,
                                                          outputfile,args.nevents))
    run = subprocess.Popen('./main_pythia2HepMC.exe -f %s -c %s -o %s -n %i' %(args.inputfile,args.pythiacfg,
                                                          outputfile,args.nevents)
                           ,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)

    output,errorMsg= run.communicate()
    logger.debug('Pythia error:\n %s \n' %errorMsg)
    logger.debug('Pythia output:\n %s \n' %output)

    logger.info('Done.')

    #Run Delphes
    logger.info('Running DelphesHepMC')
    inputfile = outputfile
    outputfile = os.path.splitext(inputfile)[0]+'.root'
    if os.path.isfile(outputfile):
        os.remove(outputfile)
    logger.debug('../Delphes/DelphesHepMC %s %s %s' %(args.delphescard,
                                                          outputfile,inputfile))
    run = subprocess.Popen('../Delphes/DelphesHepMC %s %s %s' %(args.delphescard,
                                                          outputfile,inputfile)
                           ,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)


    output,errorMsg= run.communicate()
    logger.debug('Delphes error:\n %s \n' %errorMsg)
    logger.debug('Delphes output:\n %s \n' %output)

    logger.info('Done.')
    return


if __name__ == "__main__":

    import argparse
    ap = argparse.ArgumentParser( description=
            "Runs Pythia to generate HepMC events and then Delphes." )
    ap.add_argument('-f', '--inputfile', required=True,
            help='path to the Pythia input file (SLHA or LHE).')
    ap.add_argument('-n', '--nevents', default = 100, type = int,
            help='Number of events to be generated [100].')
    ap.add_argument('-c', '--pythiacfg', default='pythia8.cfg',
            help='path to the Pythia8 config file [pythia8_chargino.cfg].')
    ap.add_argument('-C', '--delphescard', default='gen_card.tcl',
            help='path to the Delphes card file [gen_card.tcl].')
    ap.add_argument('-v', '--verbose', default='info',
            help='verbose level (debug, info, warning or error) [info]. Default is error')



    t0 = time.time()

    args = ap.parse_args()

    level = args.verbose.lower()
    levels = { "debug": logging.DEBUG, "info": logging.INFO,
               "warn": logging.WARNING,
               "warning": logging.WARNING, "error": logging.ERROR }
    if not level in levels:
        logger.error ( "Unknown log level ``%s'' supplied!" % level )
        sys.exit()
    logger.setLevel(level = levels[level])


    t0 = time.time()

    args = ap.parse_args()
    r = runAll(args)

    logger.info('Done in %1.1f min' %((time.time()-t0)/60.))
