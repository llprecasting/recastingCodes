#!/usr/bin/env python3
import numpy as np
import os,sys,glob
from pathlib import Path
import logging
from helper import (filterObjects,getModelInfo,saveOutput, \
                    effMap, deltaR, cutFlow, getD0, getZ0, getR, count_tracker_layer_crossings)
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
from ROOT import TFile,Electron, Muon, TTree


#Initialize efficiency maps
# Use smoothed 2D binned data for electrons from https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-14/figaux_19a.png
electron_reco = effMap(filepath="./ATLAS_data/HEPData-ins1831504-v2-csv/pt-d0electronefficiency.csv")

# Use smoothed 2D binned data for muons from https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-14/figaux_19b.png
muon_reco = effMap(filepath="./ATLAS_data/HEPData-ins1831504-v2-csv/pt-d0muonefficiency.csv")
ee_acceptance = effMap(filepath="./ATLAS_data/HEPData-ins1831504-v2-csv/pt-ptselectronacceptance.csv")
mm_acceptance = effMap(filepath="./ATLAS_data/HEPData-ins1831504-v2-csv/pt-ptsmuonacceptance.csv")
em_acceptance = effMap(filepath="./ATLAS_data/HEPData-ins1831504-v2-csv/pt-ptstauacceptance.csv")




def getObjects(DelphesTree: TTree) -> Any:
    """
    Returns the needed objects from the tree after applying
    minimum requirements on pT, eta and overlap.
    """

    llps = DelphesTree.bsmMothers
    
    # muons = list(DelphesTree.MuonNonIso)
    # electrons = list(DelphesTree.ElectronNonIso)
    # muons = list(DelphesTree.MuonSmear)
    # electrons = list(DelphesTree.ElectronSmear)
    # muons = list(DelphesTree.Muon)
    # electrons = list(DelphesTree.Electron)

    daughters = DelphesTree.bsmFinalDaughters
    electrons = [ptc for ptc in daughters if abs(ptc.PID) == 11]
    muons = [ptc for ptc in daughters if abs(ptc.PID) == 13]
    
    # llps = filterObjects(llps,pTmin=0.0,etaMax=5.0)
    muons = filterObjects(muons,pTmin=20.0,etaMax=2.5)
    electrons = filterObjects(electrons, pTmin=20.0, etaMax=2.5)
    for el in electrons:
        if not hasattr(el,'PID'):
            el.PID = -11*el.Charge
        # if not hasattr(el,'D0'):
        el.D0 = getD0(el)
        el.Z0 = getZ0(el)
        el.R = getR(el)
    for mu in muons:
        if not hasattr(mu,'PID'):
            mu.PID = -13*mu.Charge
        # if not hasattr(mu,'D0'):
        mu.D0 = getD0(mu)
        mu.Z0 = getZ0(mu)
        mu.R = getR(mu)

    return llps,muons,electrons

def getSR(leptons: List[Union[Electron, Muon]]) -> str:
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

def preSelection(muons: List[Union[Any, Muon]],electrons: List[Union[Electron, Any]]) -> Union[None, List[Union[Electron, Muon]]]:
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


def passTrigger(leptons : List[Union[Electron, Muon]]) -> bool:
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

def createdBeforeECAL(ptc)-> bool:
    """
    Returns True if the particle has been created before the ECAL, False otherwise. 
    Since the triggers require the electrons to deposit their energy in the ECAL, we must impose
    these requirement for electrons.
    """

    zmax = 3700.0
    rhomax = 1400.0
    rho = np.sqrt(ptc.X**2 + ptc.Y**2)
    z = np.abs(ptc.Z)
    if rho > rhomax:
        return False
    elif z > zmax:
        return False
    else:
        return True
    
def createdBeforeSCT(ptc)-> bool:
    """
    Returns True if the particle has been created before the SCT, False otherwise. 
    Since the tracking reconstruction algorithm requires several hits in the SCT, we must impose
    these requirement for electrons.
    """

    zmax = 750.0
    rhomax = 500.0
    rho = np.sqrt(ptc.X**2 + ptc.Y**2)
    z = np.abs(ptc.Z)
    if rho > rhomax:
        return False
    elif z > zmax:
        return False
    else:
        return True
    
def numberOfHits(ptc) -> int:
    """
    Returns the number of hits for a given particle in the pixel and/or SCT layers of the tracker. Since the tracking algorithm requires hits in these layers, we must impose these requirement for electrons.
    """

    # Production vertex coordinates
    R = np.sqrt(ptc.X**2 + ptc.Y**2)
    Z = ptc.Z
    # Velocity components in cylindrical coordinates
    # (only the direction is relevant)
    vR = np.sqrt(ptc.Px**2 + ptc.Py**2)
    vz = ptc.Pz
    vtot = np.sqrt(ptc.Px**2 + ptc.Py**2 + ptc.Pz**2)
    vR = vR/vtot
    vz = vz/vtot
    hits = count_tracker_layer_crossings(R,Z,vR,vz)
    n_total = hits['total']

    return n_total



def getEfficiencies(inputFile: str) -> Dict[str, Any]:


    # Define SRs and Cutflow
    ee_cutflow = cutFlow(name='ee_cutflow',levels=['All','PreSelection',
                                                    'Trigger', '2e', 
                        'pT > 65 GeV', '3 mm < d0', 'DeltaRll > 0.2','ee SR'])
    mm_cutflow = cutFlow(name='mm_cutflow',levels=['All','PreSelection',
                                                'Trigger', '2mu', 
                        'pT > 65 GeV', '3 mm < d0', 'DeltaRll > 0.2','mm SR'])
    em_cutflow = cutFlow(name='em_cutflow',levels=['All','PreSelection',
                                                    'Trigger', 'emu', 
                        'pT > 65 GeV', '3 mm < d0', 'DeltaRll > 0.2','em SR'])

    eff_SRs = cutFlow(name='eff_SRs',levels=['All',
                                            'Acceptance_ee',
                                            'AcceptanceCuts_ee',
                                            'AccEff_ee',
                                            'AccEffCuts_ee',
                                            'Acceptance_mm',
                                            'AcceptanceCuts_mm',
                                            'AccEff_mm',
                                            'AccEffCuts_mm',
                                            'Acceptance_em',
                                            'AcceptanceCuts_em',
                                            'AccEff_em',
                                            'AccEffCuts_em'])

    
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
        
        # Fill first key (All)
        ee_cutflow.fill(weight)
        em_cutflow.fill(weight)
        mm_cutflow.fill(weight)

        eff_SRs.fill_level('All', weight)

        leptons_preSel = preSelection(muons,electrons)
        if leptons_preSel is None:
            continue

        # Fill pre-Selection pass
        ee_cutflow.fill_next(weight)
        em_cutflow.fill_next(weight)
        mm_cutflow.fill_next(weight)

        signal_region = getSR(leptons_preSel)
        if signal_region == "me": 
            signal_region = "em" # collapse me and em for cutflow filling
        # Get acceptance for event
        acc = 0.0
        if signal_region == "ee":
            acc = ee_acceptance.efficiency(leading_lepton_p_textT_GeV=leptons_preSel[0].PT, 
                                           subleading_lepton_p_textT_GeV=leptons_preSel[1].PT)
            sr_cutflow = ee_cutflow
        elif signal_region == "mm":
            acc = mm_acceptance.efficiency(leading_lepton_p_textT_GeV=leptons_preSel[0].PT, 
                                           subleading_lepton_p_textT_GeV=leptons_preSel[1].PT)
            sr_cutflow = mm_cutflow
        else:
            acc = em_acceptance.efficiency(leading_lepton_p_textT_GeV=leptons_preSel[0].PT, 
                                           subleading_lepton_p_textT_GeV=leptons_preSel[1].PT)
            sr_cutflow = em_cutflow
        
        eff_SRs.fill_level(f'Acceptance_{signal_region}', acc*weight)
        # Get efficiency for event
        recoEff = 0.0
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
        eff_SRs.fill_level(f'AccEff_{signal_region}', recoEff*acc*weight)

        if not passTrigger(leptons_preSel):
            continue
        # Fill trigger pass
        ee_cutflow.fill_next(weight)
        em_cutflow.fill_next(weight)
        mm_cutflow.fill_next(weight)

        
        # Fille 2 leptons pass
        sr_cutflow.fill_next(weight)


        ## signal pt and d0 cuts 
        if leptons_preSel[0].PT < 65 or leptons_preSel[1].PT < 65:
            continue
        sr_cutflow.fill_next(weight)

        if abs(leptons_preSel[0].D0) < 3 or abs(leptons_preSel[1].D0) < 3:
            continue
        # if abs(leptons_preSel[0].D0) > 300 or abs(leptons_preSel[1].D0) > 300:
            # continue
        sr_cutflow.fill_next(weight)
    
        if deltaR(leptons_preSel[0], leptons_preSel[1]) < 0.2:
            continue
        sr_cutflow.fill_next(weight)

        eff_SRs.fill_level(f'AcceptanceCuts_{signal_region}', weight)

        if recoEff == 0.0:
            continue

        # Apply large Z0 cut imposed by the large radius tracking
        # (see Table 1 in https://cds.cern.ch/record/2275635/files/ATL-PHYS-PUB-2017-014.pdf)
        if any(abs(lep.Z0) > 1500. for lep in leptons_preSel):
            continue
        # Since the trigger requires the electrons to deposit their energy in the ECAL, we must impose that they are created before the ECAL.
        # if any((not createdBeforeECAL(lep) and abs(lep.PID) == 11) for lep in leptons_preSel):
            # continue
        # In order to have enough hits in the inner detector for the large radius tracking, we must impose that the leptons are created within R ~ 440 mm from the beamline (see Table 3 in https://cds.cern.ch/record/2275635/files/ATL-PHYS-PUB-2017-014.pdf)
        # if any(lep.R > 300. for lep in leptons_preSel):
            # continue
        nhits = [numberOfHits(lep) for lep in leptons_preSel]
        if any(nh < 3 for nh in nhits):
            continue

        sr_cutflow.fill_next(weight*recoEff)
        eff_SRs.fill_level(f'AccEffCuts_{signal_region}', recoEff*weight)

        


    #End of loop
    f.Close()
    logger.info(f"Loop Ended! Evts analysed: {ct}")
    for cutflow in [ee_cutflow,mm_cutflow,em_cutflow]:
        cutflow.divide(cutflow.weights[0]) # Normalize to total number of events (first level of cutflow)
        logger.debug(f"{cutflow.to_string()}\n\n")

    eff_SRs.divide(totalweight)
    eff_dict = {}
    eff_dict['Eff SR'] = eff_SRs.to_dict()
    eff_dict['totalweight'] = totalweight
    eff_dict['Nevents'] = ct
    eff_dict['inputFile'] = inputFile
    
    return eff_dict


def main(inputfile: str,llpPDG :int = 1000011) -> None:

    # Read banner file to extract information about LLP mass, LLP lifetime and total cross-section
    bannerFile = None
    d = os.path.dirname(inputfile)
    b_files = list(glob.glob(os.path.join(d,'*_banner.txt')))
    if not b_files:
        logger.error(f"No banner files found in {d}!")
        raise ValueError()

    logger.debug(f'\nRunning over input file {inputfile}')
    
    bannerFile = b_files[0]
    modelDict = getModelInfo(bannerFile,llpPDG)

    resDict = getEfficiencies(inputfile)
    resDict.update(modelDict)
    if 'totalweight' in resDict and "Number of Events" in resDict:
        resDict['Cross-Section (pb)'] = resDict['totalweight']/resDict["Number of Events"]
    else:
        resDict['Cross-Section (pb)'] = None

    outFile = inputfile.split('.root')[0].split('.hepmc')[0]
    outFile = outFile +'_effs.json'
    effsDict = resDict.pop('Eff SR')
    for key,(val,val_err) in list(effsDict.items()):
        resDict[key] = val
        resDict[key+'_err'] = val_err

    logger.info(f'tau(ns) = {modelDict["tau0_ns"]:1.3g}, mLLP(GeV) = {modelDict["mLLP"]:1.3g}')
    saveOutput(resDict,outFile)


    


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
    
    inputFiles = sorted(inputFiles)
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

