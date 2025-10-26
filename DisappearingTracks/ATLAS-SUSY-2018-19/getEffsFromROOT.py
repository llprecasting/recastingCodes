#!/usr/bin/env python3
import numpy as np
import os,sys
from pathlib import Path
import tqdm
import logging
from helper import (DeltaPhi, filterObjects, \
                    overlapRemoval, minDphilist, eff_trigger, \
                    getLLPDecayRadius,getLLPLifetime,electronPtSmear,\
                    eff_track_EWK,eff_track_Strong, cutFlow)
FORMAT = '%(levelname)s: %(message)s at %(asctime)s'
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


def getObjects(DelphesTree):

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


def preSelection(muons,electrons,jets,met,weight,
                 ewk_cutflow,strong_cutflow):
 

  weight_EWK = 0.0
  weight_Strong = 0.0
  
  #Event Cleaning
  #Clean Bad Jets
  passedbadJets=True  ### ????
  for jet in jets:
    if jet.Eta>2.4:
      passedbadJets=False
      break
  if not passedbadJets:
    return weight_EWK,weight_Strong

  ewk_cutflow.fill_next(weight)
  strong_cutflow.fill_next(weight)

  if not jets:
    return weight_EWK,weight_Strong

  jet1pt = jets[0].PT

  temp = weight

  if met.MET>500:
    pass
  else:
    weight = weight*eff_trigger.reweight(met.MET,min(499.0,jet1pt))
    if np.isnan(weight) or weight==0:
      return weight_EWK,weight_Strong

  if np.isnan(weight) or weight < 0:
    print("Previous weight: ",temp," / Current Weight: ",weight)
    raise ValueError

  ewk_cutflow.fill_next(weight)
  strong_cutflow.fill_next(weight)

  if len(electrons)>0 or len(muons)>0:
    return weight_EWK,weight_Strong

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
  weight_EWK = weight*float(passedKinEWK)
  weight_Strong = weight*float(passedKinStrong)

  return weight_EWK,weight_Strong

def getEfficiencies(inputFile,tau0,tauList):

  tauList = np.array(tauList)
  eff_dict = {'EWK' : np.zeros((len(tauList),2)),
              'Strong' : np.zeros((len(tauList),2)),}

  f = ROOT.TFile(inputFile,'read')
  DelphesTree = f.Get('Delphes')
  nevts = DelphesTree.GetEntries()
 
  totalweight = 0
  ct=0

  for entry in tqdm.tqdm(range(0,nevts)):
    DelphesTree.GetEntry(entry)
    ct+=1
    weights = DelphesTree.Weight.At(0).Weight
    totalweight += weights
    # weight = float(weights)
    weight = 1.0
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

    weight_EWK,weight_Strong = preSelection(muons,electrons,jets,met,weight,
                                            ewk_cutflow,strong_cutflow)
    
    if (not weight_EWK) and (not weight_Strong):
      continue

    if not llps:
      continue

    for illp,llp in enumerate(llps):
      llp.daughter = DelphesTree.bsmDirectDaughters.At(llp.D1)
      llp.decayR = getLLPDecayRadius(llp)
      llp.decayT0 = getLLPLifetime(llp)
      llp.smearedPt = electronPtSmear(llp.PT, llp.Charge)
      llp.tracklet_weight = {'EWK' : np.zeros(len(tauList)),
                             'Strong' : np.zeros(len(tauList))}

      if weight_EWK:
        if illp == 0:
          ewk_SR.fill_next(weight_EWK)
        else:
          ewk_SR.fill(weight_EWK)
        tracklet_weight = eff_track_EWK.reweight(llp.Eta,llp.decayR)
        if np.isnan(tracklet_weight):
          print("NAN weight detected: ",weight_EWK," entry number: ",entry)
          continue
        llp.tracklet_weight['EWK'] = tracklet_weight

      if weight_Strong:
        if illp == 0:
          strong_SR.fill_next(weight_Strong)
        else:
          strong_SR.fill(weight_Strong)
        tracklet_weight = eff_track_Strong.reweight(llp.Eta,llp.decayR)
        if np.isnan(tracklet_weight):
          print("NAN weight detected: ",weight," entry number: ",entry)
          continue
        llp.tracklet_weight['Strong'] = tracklet_weight


    # Consider only the LLP with highest smeared PT:
    sr_weight_EWK = np.zeros(len(tauList))
    sr_weight_Strong = np.zeros(len(tauList))
    clevel_EWK = ewk_SR.current_level
    clevel_Strong = strong_SR.current_level
    for illp,llp in enumerate(llps):
      llp.weight_EWK = 0.0
      llp.weight_Strong = 0.0
      # Lifetime reweighting:
      llp.lifetime_reweight = np.exp(llp.decayT0/tau0-llp.decayT0/tauList)

      # Reset cutflows to correct level
      ewk_SR.reset(to_level=clevel_EWK)
      strong_SR.reset(to_level=clevel_Strong)

      track_weight_EWK = weight_EWK*llp.tracklet_weight['EWK']
      track_weight_Strong = weight_Strong*llp.tracklet_weight['Strong']

      # Skip LLPs with zero total weight for both SRs
      if (not track_weight_EWK) and (not track_weight_Strong):
        continue
      
      ewk_SR.fill_next(track_weight_EWK)
      strong_SR.fill_next(track_weight_Strong)
      
      if llp.smearedPt < 20: #Since we've sorted by descending smeared pT, no further charginos to consider once one hits the threshold
        continue

      #The following are included in the efficiency map:
      # 4 pixel laters and nSCT Hits == 0
      # nGangedFlaggedFake == 0
      # Pixel spoilt hits == 0
      # nPixel outliers == 0
      # |d0significance| < 1.5
      # |z0sin(theta)| < 0.5
      # ptcone40pT<0.04 (isolated)

      ewk_SR.fill_next(track_weight_EWK)
      strong_SR.fill_next(track_weight_Strong)

      #DeltaR(jets) > 0.4
      if any(llp.P4().DeltaR(jet.P4())<0.4 for jet in jets):
        continue
      ewk_SR.fill_next(track_weight_EWK)
      strong_SR.fill_next(track_weight_Strong)
      #DeltaR(electron) > 0.4
      if any(llp.P4().DeltaR(elec.P4())<0.4 for elec in electrons):
        continue
      ewk_SR.fill_next(track_weight_EWK)
      strong_SR.fill_next(track_weight_Strong)

      #DeltaR(muon) > 0.4
      if any(llp.P4().DeltaR(muon.P4())<0.4 for muon in muons):
        continue
      ewk_SR.fill_next(track_weight_EWK)
      strong_SR.fill_next(track_weight_Strong)

      # 0.1 < abs(eta) < 1.9
      trackletEta = abs(llp.Eta)
      if not (0.1 < trackletEta < 1.9):
        continue
      ewk_SR.fill_next(track_weight_EWK)
      strong_SR.fill_next(track_weight_Strong)

      # If passed all the cuts set the LLP total weight,
      # (weight*tracklet_weight) for the corresponding SR
      llp.weight_EWK = track_weight_EWK
      llp.weight_Strong = track_weight_Strong
     
    
    
    # Reweight by lifetime:
    
    # Compute the event efficiency:
    # eff = 1 - prod_i (1-llp[i].eff*reweight)
    evt_eff_EWK = 1.0 - np.prod([(1.0-llp.weight_EWK*llp.lifetime_reweight) 
                                 for llp in llps],axis=0)
    evt_eff_Strong = 1.0 - np.prod([(1.0-llp.weight_Strong*llp.lifetime_reweight) 
                                    for llp in llps],axis=0)

    eff_dict['EWK'] += evt_eff_EWK
    eff_dict['Strong'] += evt_eff_Strong
        
    #Calo-veto would enter here, but also folded in the efficiency map
    #fill("hist_Tracklet_Pt",chargino.PT)

  #End of loop
  print("Loop Ended! Evts analysed: ",ct,'\n')
  eff_dict['EWK'] = eff_dict['EWK']/totalweight
  eff_dict['Strong'] = eff_dict['Strong']/totalweight
  eff_dict['totalweight'] = totalweight
  eff_dict['Nevents'] = ct

  logger.debug(f"{ewk_cutflow.to_string()}\n\n")
  logger.debug(f"{strong_cutflow.to_string()}\n\n")

  logger.info(f"{ewk_SR.to_string()}\n\n")
  logger.info(f"{strong_SR.to_string()}\n\n")

  
  return eff_dict

if __name__ == "__main__":
      
  #Process events
  import argparse
  parser = argparse.ArgumentParser(description='Analyse delphesLLP output to produce efficiencies for ATLAS-ANA-SUSY-2019-018 DT search')
  parser.add_argument('inputfile', metavar='inputfile_path', help='Path to the delphesLLP root file with the event sample to be analysed.')
  parser.add_argument('-tau0','--tau0',metavar='tau0', help='Proper lifetime (in ns) used for event generation',type=float, required=True)
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
  tauList = [float(tau0)]
  if tau_file:
    if not os.path.isfile(tau_file):
      raise ValueError(f"Reweighting file {tau_file} not found!")
    try:
      import csv
      with open(tau_file, mode='r', newline='') as file:
        csv_reader = csv.reader(file)
        tauList = [float(row[0]) for row in csv_reader if row[0]]
    except:
      logger.error(f"Error reding {tau_file}. Reweighting will not be applied.")
    tauList = np.sort(np.unique(tauList))

  effsDict = getEfficiencies(inputfile,tau0,tauList)

  print(effsDict)

#print("Final SR error: %1.6e\n"%(np.sqrt(err)))
#print(err_sqr)
