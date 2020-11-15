#!/usr/bin/env python3

from __future__ import print_function
import time,os,sys
try:
    from ConfigParser import SafeConfigParser as ConfigParser
except ImportError:
    from configparser import ConfigParser


def runDT(inputFile,code_dir):
    #import random
    #random.seed(10) #Set the random seeed to generate reproduciable results

    code_dir = os.path.abspath('../../DisappearingTracks/ATLAS-SUSY-2016-06')


    start_time = time.time()


    #Set default parameters:
    parser = ConfigParser({'kfactor': 1, 'PIDchargino' : 1000024, 'PIDneutralino' : 1000022,
                            'nloop' : 1000, 'Weight' : 'root',
                            'effFile' : './Modules/DisappearingTrack2016-TrackAcceptanceEfficiency.root'})
    #Get user-defined parameters
    inputFile = 'input_param_DT.dat'
    if not os.path.isfile(inputFile):
        print('Could not find %s file' %inputFile)
        sys.exit()

    parser.read(inputFile)
    if not parser.has_option('options','DelphesPath'):
        print('The DelphePath variable must be set in %s' %inputFile)
        sys.exit()

    delphesPath = os.path.abspath(os.path.expanduser(parser.get('options','DelphesPath')))
    if not os.path.isdir(delphesPath):
        print('Delphes folder %s not found' %delphesPath)
        sys.exit()

    current_dir = os.getcwd()
    Modules_dir = os.path.abspath(os.path.join(code_dir,'Modules'))

    sys.path.append(Modules_dir)
    from eventSelector import EventSelector

    #Create an event selector:
    evtSelector = EventSelector(parser)

    #Set Delphes path and load required libraries:
    evtSelector.loadDelphesLib(delphesPath)

    #Load ROOT/Delphes events file:
    evtSelector.loadRootFile()

    #Select events and compute weights for the lifetimes in tau_array:
    evtSelector.selectEvents()

    #Save output to file:
    evtSelector.writeResults()


    final_time = time.time() - start_time
    if (final_time < 60.0):
        print("time final =",final_time , "seconds")
    elif (final_time >= 60.0 and final_time < 3600):
        print("time final =",final_time/60.0 , "minutes")
    else:
        print("time final =",final_time/3600.0 , "hours")


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
