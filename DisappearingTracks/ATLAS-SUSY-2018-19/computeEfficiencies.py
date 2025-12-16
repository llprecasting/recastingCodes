#!/usr/bin/env python3
import numpy as np
import os,sys,glob
from pathlib import Path
import tqdm
import logging
from helper import (filterObjects,getModelInfo,saveOutput, \
                    overlapRemoval, minDphilist, eff_trigger, \
                    getLLPDecayRadius,getLLPDecayTime,electronPtSmear,\
                    eff_track_EWK,eff_track_Strong, cutFlow)
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



def getObjects(DelphesTree: TTree) -> Any:
    """
    Returns the needed objects from the tree after applying
    minimum requirements on pT, eta and overlap.
    """

    llps = DelphesTree.bsm
    jets = DelphesTree.Jet
    bmet = DelphesTree.MissingET
    muons = DelphesTree.Muon
    electrons = DelphesTree.Electron
    
    llps = filterObjects(llps,pTmin=0.0,etaMax=5.0)
    met = bmet.At(0)
    jets = filterObjects(jets,pTmin=20.0,etaMax=2.8)
    muons = filterObjects(muons,pTmin=10.0,etaMax=2.7)
    electrons = filterObjects(electrons, pTmin=10.0, etaMax=2.47)
    
    #Overlap Removal
    electrons = overlapRemoval(electrons, muons, 0.05)
    electrons = overlapRemoval(electrons, jets, 0.4)
    muons = overlapRemoval(muons, jets, 0.4)
    jets.sort(key=lambda j: j.PT,reverse=True)

    return llps,muons,electrons,jets,met


def preSelection(muons: List[Union[Any, Muon]],electrons: List[Union[Electron, Any]],jets: List[Union[Any, Jet]],met: MissingET,weight: float,
                 ewk_cutflow: cutFlow,strong_cutflow: cutFlow) -> Tuple[float, float]:
    """
    Apply the pre-selection including the trigger efficiency.
    Returns the pre-selection efficiency for each signal region 
    (eff = 0 if the event did not pass the cuts)
    """
 
    preSel_eff_EWK = 0.0
    preSel_eff_Strong = 0.0
    
    #Event Cleaning
    #Clean Bad Jets
    if any(abs(jet.Eta) > 2.4 for jet in jets):
        return preSel_eff_EWK,preSel_eff_Strong

    ewk_cutflow.fill_next(weight)
    strong_cutflow.fill_next(weight)

    if not jets:
        return preSel_eff_EWK,preSel_eff_Strong

    jet1pt = jets[0].PT

    if met.MET>500:
        preSel_eff = 1.0
    else:
        preSel_eff = eff_trigger.efficiency(met.MET,min(499.0,jet1pt))
        if np.isnan(preSel_eff) or preSel_eff==0:
            return preSel_eff_EWK,preSel_eff_Strong

    ewk_cutflow.fill_next(weight)
    strong_cutflow.fill_next(weight)

    if len(electrons)>0 or len(muons)>0:
        return preSel_eff_EWK,preSel_eff_Strong

    ewk_cutflow.fill_next(weight)
    strong_cutflow.fill_next(weight)

    #EWK-specific pre-selection
    passedKinEWK = True

    #MET Cut
    if met.MET < 200:
        passedKinEWK=False
    else:
        ewk_cutflow.fill_next(weight)

    if passedKinEWK:
        if jet1pt < 100:
            passedKinEWK=False
        else:
            ewk_cutflow.fill_next(weight)

    if passedKinEWK:
        if minDphilist(met,jets,4,50.0) > 1.0:
            ewk_cutflow.fill_next(weight)
        else:
            passedKinEWK=False

    #Strong-specific Pre-selection
    passedKinStrong=True

    #MET Cut
    if met.MET < 250:
        passedKinStrong=False
    else:
        strong_cutflow.fill_next(weight)

    #Jet PT Cuts
    if passedKinStrong:
        if jet1pt < 100:
            passedKinStrong=False
        elif len(jets)<3:
            passedKinStrong=False
        elif jets[1].PT<20 or jets[2].PT < 20:
            passedKinStrong=False
        else:
            strong_cutflow.fill_next(weight)

    #min(DeltaPhi(Jet,MET)) cut
    if passedKinStrong:
        if minDphilist(met,jets,4,50.0) > 0.4:
            strong_cutflow.fill_next(weight)
        else:
            passedKinStrong=False

    #Event must pass at least one preselection
    preSel_eff_EWK = preSel_eff*float(passedKinEWK)
    preSel_eff_Strong = preSel_eff*float(passedKinStrong)

    return preSel_eff_EWK,preSel_eff_Strong

def getEfficiencies(inputFile: str,tau0: float,tauList: ndarray,ijob: int=0) -> Dict[str, Any]:

    tauList = np.array(tauList)
    eff_SR = cutFlow(name="Efficiencies",levels=['EWK SR', 'Strong SR'],
                      zero_weight=np.zeros(len(tauList)))

    f = TFile(inputFile,'read')
    DelphesTree = f.Get('Delphes')
    nevts = DelphesTree.GetEntries()
    
    totalweight = 0
    ct=0


    for entry in tqdm.tqdm(range(nevts),position=ijob,
                            desc=inputFile,
                            leave=False):
        DelphesTree.GetEntry(entry)
        ct+=1
        # weights = float(DelphesTree.Weight.At(0).Weight)
        weight = 1.0
        totalweight += weight
        llps,muons,electrons,jets,met = getObjects(DelphesTree)

        # Reset cutflows to beginning
        ewk_cutflow.reset()
        strong_cutflow.reset()
        ewk_SR.reset()
        strong_SR.reset()

        # Fill first key (All)
        ewk_cutflow.fill(weight)
        strong_cutflow.fill(weight)
        ewk_SR.fill(weight)
        strong_SR.fill(weight)

        preSel_eff_EWK,preSel_eff_Strong = preSelection(muons,electrons,jets,met,weight,
                                                        ewk_cutflow,strong_cutflow)
    
        if (not preSel_eff_EWK) and (not preSel_eff_Strong):
            continue

        if not llps:
            continue

        # Compute relevant LLP variables
        for llp in llps:
            llp.daughter = DelphesTree.bsmDirectDaughters.At(llp.D1)
            llp.decayR = getLLPDecayRadius(llp)
            llp.decayT = getLLPDecayTime(llp)
            llp.smearedPt = electronPtSmear(llp.PT, llp.Charge)
            track_eff_EWK =  eff_track_EWK.efficiency(llp.Eta,llp.decayR)
            if np.isnan(track_eff_EWK):
                track_eff_EWK = 0.0
            track_eff_Strong =  eff_track_Strong.efficiency(llp.Eta,llp.decayR)
            if np.isnan(track_eff_Strong):
                track_eff_Strong = 0.0
                
            llp.tracklet_eff_EWK = track_eff_EWK
            llp.tracklet_eff_Strong = track_eff_Strong
            # Lifetime reweighting:
            if tau0 > 0.0:
                # llp.lifetime_reweight = np.exp(llp.decayT0/tau0-llp.decayT0/tauList)
                gamma = llp.P4().Gamma()
                llp.lifetime_reweight = (tau0/tauList)*np.exp(-(llp.decayT/gamma)*(1/tauList-1/tau0))
            else:
                llp.lifetime_reweight = np.ones(tauList.shape)
    

        # Sort by smearedPt:
        llps = sorted(llps, key=lambda llp: llp.smearedPt,reverse=True)

        # Add one entry for each llp
        ewk_SR.fill_next(weight*preSel_eff_EWK*len(llps))
        strong_SR.fill_next(weight*preSel_eff_Strong*len(llps))
        
        # Selected llps with track effciency > 0:
        llps_EWK = [llp for llp in llps if llp.tracklet_eff_EWK > 0.0]
        fill_EWK = 0.0
        if llps_EWK:
            fill_EWK = weight*preSel_eff_EWK*llps_EWK[0].tracklet_eff_EWK
        ######### Use only leading LLP!!
            # llps_EWK = llps_EWK[:1]
        ewk_SR.fill_next(fill_EWK)

        llps_Strong = [llp for llp in llps if llp.tracklet_eff_Strong > 0.0]
        fill_Strong = 0.0
        if llps_Strong:
            fill_Strong = weight*preSel_eff_Strong*llps_Strong[0].tracklet_eff_Strong
        ######### Use only leading LLP!!
            # llps_Strong = llps_Strong[:1]
        strong_SR.fill_next(fill_Strong)

        # Select llps with smearedPt > minPT:
        # minPT = 60 (20) GeV for the model-independent (model-dependent) search strategy
        minPT = 60.0
        llps_EWK = [llp for llp in llps_EWK[:] if llp.smearedPt > minPT]
        if llps_EWK:
            ewk_SR.fill_next(fill_EWK)

        llps_Strong = [llp for llp in llps_Strong[:] if llp.smearedPt > minPT]
        if llps_Strong:
            strong_SR.fill_next(fill_Strong)

        # Remove LLPs with overlap to jets, electrons and muons:
        for objList in [jets,electrons,muons]:
            llps_EWK = overlapRemoval(llps_EWK,objList,0.4)
            if llps_EWK:
                ewk_SR.fill_next(fill_EWK)
            llps_Strong = overlapRemoval(llps_Strong,objList,0.4)
            if llps_Strong:
                strong_SR.fill_next(fill_Strong)
        # Apply eta cut: 0.1 < abs(eta) < 1.9
        llps_EWK = [llp for llp in llps_EWK if 0.1 < abs(llp.Eta) < 1.9]
        if llps_EWK:
            ewk_SR.fill_next(fill_EWK)
        llps_Strong = [llp for llp in llps_Strong if 0.1 < abs(llp.Eta) < 1.9]
        if llps_Strong:
            strong_SR.fill_next(fill_Strong)
        
        # Finally compute event weight:
        # evt_weight = weight*preSelectionEff*llp_eff    

        evt_weight_EWK = np.zeros(len(tauList))
        evt_weight_Strong = np.zeros(len(tauList))

        if llps_EWK:
        # Require at least one LLP to be reconstructed and isolated
            llp_eff = 1.0 - np.prod([(1.0-llp.tracklet_eff_EWK*llp.lifetime_reweight)
                                    for llp in llps_EWK],axis=0)
            evt_weight_EWK = weight*preSel_eff_EWK*llp_eff
        
        if llps_Strong:
        # Require at least one LLP to be reconstructed and isolated
            llp_eff = 1.0 - np.prod([(1.0-llp.tracklet_eff_Strong*llp.lifetime_reweight)
                                    for llp in llps_Strong],axis=0)   
            evt_weight_Strong = weight*preSel_eff_Strong*llp_eff


        eff_SR.fill_level('EWK SR',evt_weight_EWK)
        eff_SR.fill_level('Strong SR',evt_weight_Strong)

    #End of loop
    logger.info(f"Loop Ended! Evts analysed: {ct}")
    eff_SR.divide(totalweight)
    eff_dict = {}
    eff_dict['Eff SR'] = eff_SR
    eff_dict['totalweight'] = totalweight
    eff_dict['Nevents'] = ct
    eff_dict['inputFile'] = inputFile
    eff_dict['tau_ns'] = tauList
    eff_dict['tau0_ns'] = tau0

    logger.debug(f"{ewk_cutflow.to_string()}\n\n")
    logger.debug(f"{strong_cutflow.to_string()}\n\n")

    logger.debug(f"{ewk_SR.to_string()}\n\n")
    logger.debug(f"{strong_SR.to_string()}\n\n")

    
    return eff_dict



def main(inputfile: str,llpPDG :int, tau_file: Union[str,None],ijob: int=0):

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

    tauList = [float(tau0)]
    if tau_file is not None:
        if not os.path.isfile(tau_file):
            raise ValueError(f"Reweighting file {tau_file} not found!")
        try:
            import csv
            with open(tau_file, mode='r', newline='') as file:
                csv_reader = csv.reader(l for l in file.readlines() 
                                        if not l.strip().startswith('#'))
                tauList += [float(row[0]) for row in csv_reader if row]
            
        except Exception as e:
            logger.error(str(e))
            logger.error(f"Error reding {tau_file}. Reweighting will not be applied.")

    tauList = np.array([float(f"{tau:1.4e}") for tau in tauList[:]])
    tauList = np.sort(np.unique(tauList))
    resDict = getEfficiencies(inputfile,tau0,tauList,ijob)
    resDict.update(modelDict)

    outFile = inputfile.split('.root')[0].split('.hepmc')[0]
    outFile = outFile +'_effs.json'

    effs = resDict.pop('Eff SR').to_dict()
    tauList = resDict.pop("tau_ns")
    effsList = []
    for itau,tau in enumerate(tauList):
        effDict = {'tau_ns' : tau}
        for sr in effs:
            effDict[sr] = effs[sr][0][itau]
            effDict[sr+' Error'] = effs[sr][1][itau]
        effsList.append(effDict)

    i, = np.where(np.isclose(tauList, tau0,rtol=1e-3))
    i = i[0]
    logger.info(f'tau(ns) = {tauList[i]:1.3g}:')
    for sr in effs:
        logger.info(f'  {sr} = {effs[sr][0][i]:1.3e} +- {effs[sr][1][i]:1.3e}')
    
    
    resDict['Efficiencies'] = effsList
    saveOutput(resDict,outFile)


    


if __name__ == "__main__":
      
    import argparse
    parser = argparse.ArgumentParser(description='Analyse the Delphes output to produce efficiencies for the ATLAS-SUSY-2019-18 DT search')
    parser.add_argument('-i','--input', help='Path to  Delphes ROOT file or to a folder containing Delphes ROOT files with the event samples to be analysed.')
    parser.add_argument('-l','--llpPDG',help='LLP PDG [1000024]',type=int, required=False, default=1000024)
    parser.add_argument('-tauF','--tau_file',metavar='tau_file', help='CSV file containing the lifetime values (in ns) used for reweighting. If empty or file not found, it will not apply reweighting [tau_list.csv].',
                        type=str, required=False, default='tau_list.csv')
    parser.add_argument('-n', '--ncpus',type=int,default=1,help='number of parallel jobs to run when running over multiple files [default=1].')
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


    inputF = args.input
    tau_file = args.tau_file
    if not tau_file or not os.path.isfile(tau_file):
        tau_file = None
    llpPDG = args.llpPDG

    if os.path.isfile(inputF):
        inputFiles = [os.path.abspath(inputF)]
    elif os.path.isdir(args.folder):
        # Find root files:
        pattern = os.path.join(args.folder, "**", f"*.root")
        inputFiles = list(glob.glob(pattern, recursive=True))
        if not inputFiles:
            logger.error(f"No .root files found in {args.folder}!")
            raise ValueError()
            
    else:
        logger.error(f"File/Folder {args.input} not found!")
        raise ValueError()
    
    logger.info(f"Running over {len(inputFiles)} files")
    ncpus = min(len(inputFiles),args.ncpus)
    pool = multiprocessing.Pool(processes=ncpus)
    children = []
    ijob = 0
    for rootFile in inputFiles:
        p = pool.apply_async(main, args=(rootFile,args.llpPDG,args.tau_file,ijob,))
        ijob += 1
        children.append(p)

    logger.info(f'Running {ijob} jobs in {ncpus} instances')
    for p in children: 
        p.get()
