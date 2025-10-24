#!/usr/bin/env python3
import numpy as np
import os,sys
from pathlib import Path
import tqdm
import logging
from helper import (DeltaPhi, filterParticles, filterJets, \
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
ewk_cutflow = cutFlow(['All', 'GRL and Cleaning', 'MET Trigger', 'Lepton Veto', 
                    'MET > 200 GeV', 'Jet pT > 100 GeV', 'min(DeltaPhi(JetMET)) > 1.0'])
strong_cutflow = cutFlow(['All', 'GRL and Cleaning', 'MET Trigger', 'Lepton Veto',
                      'MET > 250 GeV', 'Jet pT > 100,20,20 GeV', 'min(DeltaPhi(JetMET)) > 0.4'])
ewk_SR = cutFlow(['All', 'Kinematic', 'Tracklet Emulation', 'Leading tracklet',
                                'DeltaR(jet) > 0.4', 'DeltaR(electron) > 0.4', 'DeltaR(muon) > 0.4',
                                 '0.1 < Eta < 1.9'])
strong_SR = cutFlow(['All', 'Kinematic', 'Tracklet Emulation', 'Leading tracklet',
                                'DeltaR(jet) > 0.4', 'DeltaR(electron) > 0.4', 'DeltaR(muon) > 0.4',
                                 '0.1 < Eta < 1.9'])


def getObjects(DelphesTree):

  llps = DelphesTree.bsm
  jets = DelphesTree.Jet
  bmet = DelphesTree.MissingET
  muons = DelphesTree.Muon
  electrons = DelphesTree.Electron
  
  llps = filterParticles(llps,pTmin=0.0,etaMax=5.0)
  met = bmet.At(0).MET
  jets = filterJets(jets,pTmin=20.0,etaMax=2.8)
  muons = filterParticles(muons,pTmin=10.0,etaMax=2.7)
  electrons = filterParticles(electrons, pTmin=10.0, etaMax=2.47)
  
  #Overlap Removal
  electrons = overlapRemoval(electrons, muons, 0.05)
  electrons = overlapRemoval(electrons, jets, 0.4)
  muons = overlapRemoval(muons, jets, 0.4)
  jets.sort(key=lambda j: j.PT,reverse=True)

  return llps,muons,electrons,jets,met


def preSelection(muons,electrons,jets,met,weight,
                 ewk_cutflow,strong_cutflow):
 
  
  #Event Cleaning
  #Clean Bad Jets
  passedbadJets=True  ### ????
  for jet in jets:
    if jet.Eta>2.4:
      passedbadJets=False
      break
  if not passedbadJets:
    return 0.0

  ewk_cutflow.fill(weight)
  strong_cutflow.fill(weight)

  weight_EWK = 0.0
  weight_Strong = 0.0

  if not jets:
    return 0.0

  jet1pt = jets[0].PT

  temp = weight

  if met>500:
    pass
  else:
    weight = weight*eff_trigger.reweight(met,min(499.0,jet1pt))
    if np.isnan(weight) or weight==0:
      return weight_EWK,weight_Strong

  if np.isnan(weight) or weight < 0:
    print("Previous weight: ",temp," / Current Weight: ",weight)
    raise ValueError

  ewk_cutflow.fill(weight)
  strong_cutflow.fill(weight)

  if len(electrons)>0 or len(muons)>0:
    return weight_EWK,weight_Strong

  ewk_cutflow.fill(weight)
  strong_cutflow.fill(weight)

  #EWK-specific pre-selection
  passedKinEWK = True

  #MET Cut
  if met < 200:
    passedKinEWK=False
  else:
    ewk_cutflow.fill(weight)

  if passedKinEWK:
    if jet1pt < 100:
      passedKinEWK=False
    else:
      ewk_cutflow.fill(weight)

  if passedKinEWK:
    if minDphilist(met,jets,4,50.0) > 1.0:
      ewk_cutflow.fill(weight)
    else:
      passedKinEWK=False

  #Strong-specific Pre-selection
  passedKinStrong=True

  #MET Cut
  if met < 250:
    passedKinStrong=False
  else:
    strong_cutflow.fill(weight)

  #Jet PT Cuts
  if passedKinStrong:
    if jet1pt < 100:
      passedKinStrong=False
    elif len(jets)<3:
      passedKinStrong=False
    elif jets[1].PT<20 or jets[2].PT < 20:
      passedKinStrong=False
    else:
      strong_cutflow.fill(weight)

  #min(DeltaPhi(Jet,MET)) cut
  if passedKinStrong:
    if minDphilist(met,jets,4,50.0) > 0.4:
      strong_cutflow.fill(weight)
    else:
      passedKinStrong=False

  #Event must pass at least one preselection
  weight_EWK = weight*float(passedKinEWK)
  weight_Strong = weight*float(passedKinStrong)

  return weight_EWK,weight_Strong

def getEfficiencies(inputFile,tauList):

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
    weight = float(weights)
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

    #Signal Regions#
    #Emulate track and tracklet reco efficiency. From the SimpleAnalysis code, there seems to be folded in acceptances here.
    EWKLLPs = []
    StrongLLPs = []

    for llp in llps:
      llp.daughter = DelphesTree.bsmDirectDaughters.At(llp.D1)
      llp.decayR = getLLPDecayRadius(llp)
      llp.lifetime0 = getLLPLifetime(llp)
      llp.smearedPt = electronPtSmear(llp.PT, llp.Charge)

      if weight_EWK:
        ewk_SR.fill(weight)
        tracklet_weight = eff_track_EWK.reweight(llp.Eta,llp.decayR)
        if np.isnan(tracklet_weight):
          print("NAN weight detected: ",weight_EWK," entry number: ",entry)
          continue
        if len(tauList) > 0:
          tracklet_weight = tracklet_weight*np.exp(llp.lifetime0/tau_base-llp.lifetime0/args.tau_reweighting)
        if tracklet_weight > 0:
          EWKLLPs.append([llp,tracklet_weight])

      if weight_Strong:
        strong_SR.fill(weight)
        tracklet_weight = eff_track_Strong.reweight(llp.Eta,llp.decayR)
        if np.isnan(tracklet_weight):
          print("NAN weight detected: ",weight," entry number: ",entry)
          continue
        if args.tau_reweighting > 0:
          tracklet_weight = tracklet_weight*np.exp(llp.lifetime0/tau_base-llp.lifetime0/args.tau_reweighting)
        if tracklet_weight > 0:
          StrongLLPs.append([llp,tracklet_weight])

    if len(EWKLLPs) + len(StrongLLPs) == 0:
      continue
    if EWKLLPs:
      EWKLLPs.sort(key=lambda x: x[1],reverse=True)
    if StrongLLPs:
      StrongLLPs.sort(key=lambda x: x[1],reverse=True)

    #EWK tracklet selection
    
    loop_over_SR = [(ewk_SR,EWKLLPs),
                 (strong_SR,StrongLLPs)]
    for cutflow,llpList in loop_over_SR:
      first = True
      for llp,reweight in llpList:
        cutflow.fill(weight)
        
        if llp.smearedPt < 20: #Since we've sorted by descending smeared pT, no further charginos to consider once one hits the threshold
          break

        #The following are included in the efficiency map:
        # 4 pixel laters and nSCT Hits == 0
        # nGangedFlaggedFake == 0
        # Pixel spoilt hits == 0
        # nPixel outliers == 0
        # |d0significance| < 1.5
        # |z0sin(theta)| < 0.5
        # ptcone40pT<0.04 (isolated)

        #Isolated leading
        if not first: # sorted by pt - first=leading - only consider the first chargino past this point
          break
        first = False

        cutflow.fill(weight)

        #DeltaR(jets) > 0.4
        if any(llp.P4().DeltaR(jet.P4())<0.4 for jet in jets):
          continue
        cutflow.fill(weight)

        #DeltaR(electron) > 0.4
        if any(llp.P4().DeltaR(elec.P4())<0.4 for elec in electrons):
          continue
        cutflow.fill(weight)

        #DeltaR(muon) > 0.4
        if any(llp.P4().DeltaR(muon.P4())<0.4 for muon in muons):
          continue
        cutflow.fill(weight)

        # 0.1 < abs(eta) < 1.9
        trackletEta = abs(llp.Eta)
        if not (0.1 < trackletEta < 1.9):
          continue
        cutflow.fill(weight)
        
        #Calo-veto would enter here, but also folded in the efficiency map
        #fill("hist_Tracklet_Pt",chargino.PT)

  #End of loop
  print("Loop Ended! Evts analysed: ",ct,'\n')
  return ewk_cutflow,strong_cutflow,ewk_SR,strong_SR

if __name__ == "__main__":
      
  #Process events
  import argparse
  parser = argparse.ArgumentParser(description='Analyse delphesLLP output to produce efficiencies for ATLAS-ANA-SUSY-2019-018 DT search')
  parser.add_argument('inputfile', metavar='inputfile_path', help='Path to the delphesLLP root file with the event sample to be analysed.')
  parser.add_argument('-rw','--tau_reweighting',metavar='tau_reweighting_mval', help='Turn on lifetime reweighting for computation; default is 0 (no reweighting).',type=float,default=0)
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



  if args.tau_reweighting < 0:
    raise TypeError('Invalid value for lifetime, check input again.')

  inputfile = args.inputfile
  inputpath = Path(inputfile)
  run_tag = inputpath.parent.name

  masstaupairs=pd.read_pickle('run_taupairs.pcl')
  tau_base = masstaupairs.loc[run_tag]['tau_ns']

  effsDict = getEfficiencies(inputfile,masstaupairs)





#outfile = 'wino400GeV_0p2ns_hist.pkl'
#print("Saving histograms to file",outfile)
#pickle.dump(histlist, open(outfile,'wb'))
#outdf = pd.DataFrame(bkeeping)
#outdf.to_pickle(outfile)

#norm = 136*(9.12+19.37)/nevts
norm=1
histnames = ['hist_cutflow_kin_EWK','hist_cutflow_SR_EWK']
kin_cuts = ['Underflow','All,','GRL & Cleaning', 'MET Trigger', 'Lepton Veto', 'MET > 200 GeV', '1st Jet pT > 100 GeV', 'min(Dphi(Jet, MET)) > 1.0','Overflow']
SR_cuts = ['Underflow','All','Kinematic','Tracklet Emulation', 'Leading tracklet', 'DeltaR(jets) > 0.4',  'DeltaR(electron) > 0.4', 'DeltaR(muon) > 0.4', '0.1 < |eta| < 1.9','Overflow']
cutnames = [kin_cuts,SR_cuts]
#err=0
#for val in err_sqr:
#    err+=val*val

print("Processed Event Number:",nevts)
print('\n')
for hist,cuts in zip(histnames,cutnames):
    ahist = histlist[hist]
    print("Histogram: ",hist)
    for cut,entry in zip(cuts,range(ahist.GetXaxis().GetNbins()+2)):
        print(f'{cut} : {ahist.GetBinContent(entry)*norm:1.6e} +- {np.sqrt(ahist.GetSumw2().GetAt(entry))*norm:1.6e}')
    print('\n')
#print("Final SR error: %1.6e\n"%(np.sqrt(err)))
#print(err_sqr)
