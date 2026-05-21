#!/usr/bin/env python3
import numpy as np
import os,sys
from pathlib import Path
import tqdm
import logging
from recastCode.helper import (filterObjects,getLLPLifetime, \
                    overlapRemoval, minDphilist, eff_trigger, \
                    getLLPDecayRadius,getLLPDecayTime,electronPtSmear,\
                    eff_track_EWK,eff_track_Strong, cutFlow)
from recastCode.getEffsFromROOT import getObjects,preSelection,getEfficiencies,saveOutput,main
FORMAT = '%(levelname)s: %(message)s'
logging.basicConfig(format=FORMAT,datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger()   

# Fix seed so results are reproducible!
np.random.seed(seed=123)

DelphesLLP_path = Path(os.path.abspath("./DelphesLLP"))
os.environ['ROOT_INCLUDE_PATH'] = os.path.join(DelphesLLP_path,"external")

import ROOT


ROOT.gSystem.Load(os.path.join(DelphesLLP_path,"libDelphes.so"))

ROOT.gInterpreter.Declare('#include "classes/SortableObject.h"')
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')


# Define SRs and Cutflow
ewk_cutflow = cutFlow(name='EWK_cutflow',levels=['All', 'GRL and Cleaning', 'MET Trigger', 'Lepton Veto', 
                    'MET > 200 GeV', 'Jet pT > 100 GeV', 'min(DeltaPhi(JetMET)) > 1.0'])
strong_cutflow = cutFlow(name='Strong_cutflow',levels=['All', 'GRL and Cleaning', 'MET Trigger', 'Lepton Veto',
                      'MET > 250 GeV', 'Jet pT > 100,20,20 GeV', 'min(DeltaPhi(JetMET)) > 0.4'])
ewk_SR = cutFlow(name='EWK_SR',levels=['All', 'Kinematic', 'Tracklet Emulation', 'Leading tracklet',
                                'DeltaR(jet) > 0.4', 'DeltaR(electron) > 0.4', 'DeltaR(muon) > 0.4',
                                 '0.1 < Eta < 1.9'])
strong_SR = cutFlow(name='Strong_SR',levels=['All', 'Kinematic', 'Tracklet Emulation', 'Leading tracklet',
                                'DeltaR(jet) > 0.4', 'DeltaR(electron) > 0.4', 'DeltaR(muon) > 0.4',
                                 '0.1 < Eta < 1.9'])

if __name__ == "__main__":
      
    import argparse
    parser = argparse.ArgumentParser(description='Analyse delphesLLP output to produce efficiencies for ATLAS-SUSY-2019-18 DT search')
    parser.add_argument('-f','--inputfile', help='Path to the Delphes root file with the event sample to be analysed.')
    parser.add_argument('-tau0','--tau0',metavar='tau0', help='Proper lifetime (in ns) used for event generation',type=float, required=False, default=0.0)
    parser.add_argument('-tauF','--tau_file',metavar='tau_file', help='CSV file containing the lifetime values (in ns) used for reweighting',
                        type=str, required=False, default=None)
    parser.add_argument('-v', '--verbose', default='info',
                        help='verbose level (debug, info, warning or error). Default is info')

    args = parser.parse_args()



    import subprocess
    # First make sure the correct env variables have been set:
    LDPATH = subprocess.check_output('echo $LD_LIBRARY_PATH',shell=True,text=True)
    ROOTINC = subprocess.check_output('echo $ROOT_INCLUDE_PATH',shell=True,text=True)
    pythiaDir = os.path.abspath('./MG5/HEPTools/pythia8/lib')
    delphesDir = os.path.abspath('./DelphesLLP/external')
    if pythiaDir not in LDPATH or delphesDir not in ROOTINC:
        print('Enviroment variables not properly set. Run source setenv.sh first.')
        sys.exit()

    level = args.verbose
    levels = { "debug": logging.DEBUG, "info": logging.INFO,
                "warn": logging.WARNING,
                "warning": logging.WARNING, "error": logging.ERROR }
    if level in levels:       
        logger.setLevel(level = levels[level])



    inputfile = args.inputfile
    tau0 = args.tau0
    tau_file = args.tau_file
    main(inputfile,tau0,tau_file)
