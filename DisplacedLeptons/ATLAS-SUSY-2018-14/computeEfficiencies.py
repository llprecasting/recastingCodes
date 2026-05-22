#!/usr/bin/env python3
import numpy as np
import os,sys,glob
from pathlib import Path
import logging
from helper import (filterObjects,getModelInfo,saveOutput, \
                    electron_reco, muon_reco, deltaR, cutFlow)
from numpy import ndarray
from typing import Any, Dict, List, Tuple, Union
import multiprocessing
import subprocess

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
from ROOT import TFile,Electron, Jet, MissingET, Muon, TTree


# Define SRs and Cutflow
ee_cutflow = cutFlow(name='ee_cutflow',levels=['All', 'PreSelection', 'Trigger', '2e', 
                    'pT > 65 GeV', '3 mm < d0', 'DeltaRll > 0.2'])
mm_cutflow = cutFlow(name='mm_cutflow',levels=['All', 'PreSelection', 'Trigger', '2mu', 
                    'pT > 65 GeV', '3 mm < d0', 'DeltaRll > 0.2'])
em_cutflow = cutFlow(name='em_cutflow',levels=['All', 'PreSelection', 'Trigger', 'emu', 
                    'pT > 65 GeV', '3 mm < d0', 'DeltaRll > 0.2'])

ee_SR = cutFlow(name='ee_SR',levels=['All', 'Final Selection' ])
mm_SR = cutFlow(name='mm_SR',levels=['All', 'Final Selection' ])

em_SR = cutFlow(name='em_SR',levels=['All', 'Final Selection' ])



def getObjects(DelphesTree: TTree) -> Any:
    """
    Returns the needed objects from the tree after applying
    minimum requirements on pT, eta and overlap.
    """

    llps = DelphesTree.bsmMothers
    # muons = DelphesTree.MuonNonIso
    # electrons = DelphesTree.ElectronNonIso
    # muons = DelphesTree.MuonSmear
    # electrons = DelphesTree.ElectronSmear
    muons = DelphesTree.Muon
    electrons = DelphesTree.Electron
    
    llps = filterObjects(llps,pTmin=0.0,etaMax=5.0)
    muons = filterObjects(muons,pTmin=20.0,etaMax=2.5)
    electrons = filterObjects(electrons, pTmin=20.0, etaMax=2.5)
    for el in electrons:
        el.PID = -11*el.Charge
    for mu in muons:
        mu.PID = -13*mu.Charge

    return llps,muons,electrons

def getSR(leptons: List[Electron | Muon]) -> str:
    """
    Returns the signal region for a given list of leptons.
    """

    if len(leptons) != 2:
        raise ValueError(f"Error getting SR, expected 2 leptons (found {len(leptons)})")

    leptonIDs = [abs(lep.PID) for lep in leptons]
    if leptonIDs == [13, 13]:
        return "mm"
    elif leptonIDs == [11, 11]:
        return "ee"
    elif leptonIDs == [11, 13]:
        return "em"
    elif leptonIDs == [13, 11]:
        return "me"
    else:
        raise ValueError(f"Invalid lepton IDs: {leptonIDs}!")

def preSelection(muons: List[Union[Any, Muon]],electrons: List[Union[Electron, Any]]) -> Union[None,List[Electron | Muon]]:
    """
    Applies the pre-selection requirements for the different SRs and computes the trigger efficiency.
    """


    if len(muons) + len(electrons) < 2:
        # logger.debug(f"Event failed pre-selection: less than 2 leptons (muons: {len(muons)}, electrons: {len(electrons)})")
        return None
    
    allLeptons = sorted(muons + electrons, key=lambda lep: lep.PT,reverse=True)
    leptons = allLeptons[:2]

    sr = getSR(leptons)
    
     # Different eta cuts for electrons and muons -- based on atlas recommendations
    if (sr == "ee" or sr == "em") and abs(leptons[0].Eta) > 2.47:
        return None
    if (sr == "ee" or sr == "me") and abs(leptons[1].Eta) > 2.47:
        return None
    if (sr == "mm" or sr == "me") and abs(leptons[0].Eta) > 2.5: 
        return None
    if (sr == "mm" or sr == "em") and abs(leptons[1].Eta) > 2.5: 
        return None
    
    # logger.debug(f"Event passed pre-selection:  (muons: {len(muons)}, electrons: {len(electrons)})")
    return leptons


def passTrigger(leptons : List[Electron | Muon]) -> bool:
    """
    Applies the trigger requirements for the different SRs.
    """

    
    sr = getSR(leptons)

    pass_trigger = False
    if (sr == "ee" or sr == "em") and leptons[0].PT > 160: pass_trigger = True
    if (sr == "ee" or sr == "me") and leptons[1].PT > 160: pass_trigger = True
    if  sr == "ee" and (leptons[0].PT > 60 and leptons[1].PT > 60): pass_trigger = True
    if (sr == "mm" or sr == "me") and (leptons[0].PT > 60 and abs(leptons[0].Eta) < 1.07): pass_trigger = True
    if (sr == "mm" or sr == "em") and (leptons[1].PT > 60 and abs(leptons[1].Eta) < 1.07): pass_trigger = True
    
    return pass_trigger

def getEfficiencies(inputFile: str) -> Dict[str, Any]:

    
    f = TFile(inputFile,'read')
    DelphesTree = f.Get('Delphes')
    nevts = DelphesTree.GetEntries()
    
    totalweight = 0
    ct=0

    for entry in range(nevts):
        DelphesTree.GetEntry(entry)
        ct+=1

        weight = float(DelphesTree.Weight.At(0).Weight)    
        totalweight += weight
        _,muons,electrons = getObjects(DelphesTree)

        # Reset cutflows to beginning
        ee_cutflow.reset()
        mm_cutflow.reset()
        em_cutflow.reset()
        
        ee_SR.reset()
        mm_SR.reset()
        em_SR.reset()

        # Fill first key (All)
        ee_cutflow.fill(1.0)
        em_cutflow.fill(1.0)
        mm_cutflow.fill(1.0)

        ee_SR.fill(weight)
        mm_SR.fill(weight)
        em_SR.fill(weight)

        leptons_preSel = preSelection(muons,electrons)
        if leptons_preSel is None:
            continue

        # Fill Pre-selection pass
        ee_cutflow.fill_next(1.0)
        em_cutflow.fill_next(1.0)
        mm_cutflow.fill_next(1.0)

        if not passTrigger(leptons_preSel):
            continue
        # Fill trigger pass
        ee_cutflow.fill_next(1.0)
        em_cutflow.fill_next(1.0)
        mm_cutflow.fill_next(1.0)

        signal_region = getSR(leptons_preSel)
        
        # Selectt cutflow to fill based on SR
        if signal_region == "ee":
            cutflow = ee_cutflow
            eff_SR = ee_SR
        elif signal_region == "mm":
            cutflow = mm_cutflow
            eff_SR = mm_SR
        else:
            cutflow = em_cutflow
            eff_SR = em_SR
        
        cutflow.fill_next(1.0)


        ## signal pt and d0 cuts 
        if leptons_preSel[0].PT < 65 or leptons_preSel[1].PT < 65:
            continue
        cutflow.fill_next(1.0)

        if abs(leptons_preSel[0].D0) < 3 or abs(leptons_preSel[1].D0) < 3:
            continue
        # if abs(leptons_preSel[0].D0) > 300 or abs(leptons_preSel[1].D0) > 300:
            # continue
        cutflow.fill_next(1.0)
    
        if deltaR(leptons_preSel[0], leptons_preSel[1]) < 0.2:
            continue
        cutflow.fill_next(1.0)

        # Get lepton reconstruction efficiencies
        lepton_effs = []
        for lep in leptons_preSel:
            if abs(lep.PID) == 11:
                lepton_effs.append(electron_reco.efficiency(Lepton_p_textT_GeV=lep.PT, 
                                                            Lepton_d_0_mm=abs(lep.D0)))
            elif abs(lep.PID) == 13:
                lepton_effs.append(muon_reco.efficiency(Lepton_p_textT_GeV=lep.PT, 
                                                        Lepton_d_0_mm=abs(lep.D0)))
            
        recoEff = float(np.prod(lepton_effs))
        if recoEff == 0:
            continue
        eff_SR.fill_next(weight*recoEff)


    #End of loop
    logger.info(f"Loop Ended! Evts analysed: {ct}")
    for eff_SR in [ee_SR,mm_SR,em_SR]:
        eff_SR.divide(totalweight)
    eff_dict = {}
    eff_dict['Eff SR(ee)'] = ee_SR
    eff_dict['Eff SR(mm)'] = mm_SR
    eff_dict['Eff SR(em)'] = em_SR
    eff_dict['totalweight'] = totalweight
    eff_dict['Nevents'] = ct
    eff_dict['inputFile'] = inputFile

    for cutflow in [ee_cutflow,mm_cutflow,em_cutflow]:
        cutflow.divide(cutflow.weights[0]) # Normalize to total number of events (first level of cutflow)
        logger.debug(f"{cutflow.to_string()}\n\n")

    for eff_SR in [ee_SR,mm_SR,em_SR]:
        logger.debug(f"{eff_SR.to_string()}\n\n")

    
    return eff_dict


def main(inputfile: str,llpPDG :int = 1000011) -> None:

    # Read banner file to extract information about LLP mass, LLP lifetime and total cross-section
    bannerFile = None
    d = os.path.dirname(inputfile)
    b_files = list(glob.glob(os.path.join(d,'*_banner.txt')))
    if not b_files:
        logger.error(f"No banner files found in {d}!")
        raise ValueError()
    
    bannerFile = b_files[0]
    modelDict = getModelInfo(bannerFile,llpPDG)
    tau0 = modelDict['tau0_ns']

    resDict = getEfficiencies(inputfile)
    resDict.update(modelDict)
    if 'totalweight' in resDict and "Number of Events" in resDict:
        resDict['Cross-Section (pb)'] = resDict['totalweight']/resDict["Number of Events"]
    else:
        resDict['Cross-Section (pb)'] = None

    outFile = inputfile.split('.root')[0].split('.hepmc')[0]
    # outFile = outFile +'_effs.json'

    # effs = resDict.pop('Eff SR').to_dict()
    # tauList = resDict.pop("tau_ns")
    # effsList = []
    # for itau,tau in enumerate(tauList):
    #     effDict = {'tau_ns' : tau}
    #     for sr in effs:
    #         effDict[sr] = effs[sr][0][itau]
    #         effDict[sr+' Error'] = effs[sr][1][itau]
    #     effsList.append(effDict)

    # i, = np.where(np.isclose(tauList, tau0,rtol=1e-3))
    # i = i[0]
    # logger.info(f'tau(ns) = {tauList[i]:1.3g}:')
    # for sr in effs:
    #     logger.info(f'  {sr} = {effs[sr][0][i]:1.3e} +- {effs[sr][1][i]:1.3e}')
    
    
    # resDict['Efficiencies'] = effsList
    # saveOutput(resDict,outFile)


    


if __name__ == "__main__":
      
    import argparse
    parser = argparse.ArgumentParser(description='Analyse the Delphes output to produce efficiencies for the ATLAS-SUSY-2019-18 DT search')
    parser.add_argument('-i','--input', help='Path to  Delphes ROOT file or to a folder containing Delphes ROOT files with the event samples to be analysed.')
    parser.add_argument('-l','--llpPDG',help='LLP PDG [1000011]',type=int, required=False, default=1000011)
    parser.add_argument('-n', '--ncpus',type=int,default=1,help='number of parallel jobs to run when running over multiple files [default=1].')
    parser.add_argument('-v', '--verbose', default='info',
                        help='verbose level (debug, info, warning or error). Default is warning')

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


    inputF = args.input
    llpPDG = args.llpPDG

    if os.path.isfile(inputF):
        inputFiles = [os.path.abspath(inputF)]
    elif os.path.isdir(args.input):
        # Find root files:
        pattern = os.path.join(args.input, "**", f"*.root")
        inputFiles = list(glob.glob(pattern, recursive=True))
        if not inputFiles:
            logger.error(f"No .root files found in {args.input}!")
            raise ValueError()
            
    else:
        logger.error(f"File/Folder {args.input} not found!")
        raise ValueError()
    
    logger.info(f"Running over {len(inputFiles)} files")
    ncpus = min(len(inputFiles),args.ncpus)
    pool = multiprocessing.Pool(processes=ncpus)
    children = []
    if ncpus > 1:
        ijob = -1
    else:
        ijob = 0
    for rootFile in inputFiles:
        p = pool.apply_async(main, args=(rootFile,llpPDG,))
        children.append(p)

    logger.info(f'Running {len(inputFiles)} jobs in {ncpus} instances')
    for ichildren in range(len(children)):
        children[ichildren].get()

