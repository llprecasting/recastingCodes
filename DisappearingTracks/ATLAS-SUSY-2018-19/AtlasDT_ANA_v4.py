#!/usr/bin/env python3
import numpy as np
import sys,os
from scipy.stats import rv_continuous,uniform
from scipy.special import erf
from pathlib import Path
import tqdm
import math
import pickle
import pandas as pd
import re

# Fix seed so results are reproducible!
np.random.seed(seed=123)

DelphesLLP_path = os.path.abspath("./DelphesLLP")
os.environ['ROOT_INCLUDE_PATH'] = os.path.join(DelphesLLP_path,"external")

import ROOT


ROOT.gSystem.Load(os.path.join(DelphesLLP_path,"libDelphes.so"))

ROOT.gInterpreter.Declare('#include "classes/SortableObject.h"')
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')


class effMap:
  def __init__(self, mapname, filepath="DisappearingTrack2018-EfficiencyMaps.root"):
    self.fh = ROOT.TFile(filepath)
    try:
      self.h_eff = self.fh.Get(mapname)
    except:
      raise ValueError("No TH2D efficiency map with tag ", mapname," found in file ", filepath)

  def getEff(self,x,y):
    return self.h_eff.GetBinContent(self.h_eff.FindBin(x,y))

  def passes(self, x, y):
    eff = self.getEff(x,y)
    rnd = uniform.rvs()
    if rnd < eff:
      return True
    return False

  def reweight(self, x, y):
    return self.getEff(x,y)

#Initialize efficiency maps

eff_trigger = effMap('eff_trigger_average',filepath='DisappearingTrack2018-EfficiencyMaps.root')
eff_track_EWK = effMap('h_effmap_average_EWK',filepath='DisappearingTrack2018-EfficiencyMaps.root')
eff_track_Strong = effMap('h_effmap_average_Strong',filepath='DisappearingTrack2018-EfficiencyMaps.root')


#Define class for smearing functions
class smearingFunction(rv_continuous):
    """
    It is extremely slow, unless the _ppf is defined).
    It has the following useful methods: .pdf (probability density function),
    .cdf (cumulative density function), .ppf (percent point function or inverse of cdf),
    .rvs (random number generator).
    """

    def setPars(self,alpha, sigma, mean0=0.0):
      self.alpha = alpha
      self.sigma = sigma
      self.mean0 = mean0
      zmin = (self.a-mean0)/sigma
      zmax = (self.b-mean0)/sigma      
      self.norm = np.exp(-alpha**2/2)*(1-np.exp(alpha*(alpha+zmin)))/alpha
      self.norm += np.exp(-alpha**2/2)*(1-np.exp(alpha*(alpha-zmax)))/alpha
      self.norm += np.sqrt(2*np.pi)*erf(alpha/np.sqrt(2))
      self.norm = float(self.norm*sigma)


    def _pdf(self, x):
        sigma = self.sigma
        alpha = self.alpha
        mean0 = self.mean0
        z = (x-mean0)/sigma
        if z < -alpha:
          return np.exp(alpha*(z + (alpha/2.0)))/self.norm
        elif z > alpha:
          return np.exp(-alpha*(z - (alpha/2.0)))/self.norm
        else:
          return np.exp(-z**2/2.0)/self.norm
    
    def _cdf(self,x):
        sigma = self.sigma
        alpha = self.alpha
        mean0 = self.mean0
        a = (self.a-mean0)/sigma
        b = (self.b-mean0)/sigma
        z = (x[0]-mean0)/sigma
        z = min(b,z)
        cdf_ret = 0.0
        if z < a:
           return cdf_ret
        if z >= a:
           i_a = np.exp(alpha**2/2)*np.exp(a*alpha)/alpha
           i_z1 = np.exp(alpha**2/2)*np.exp(z*alpha)/alpha
           cdf_ret = i_z1 - i_a           
        if z > -alpha:
           i_z1 = np.exp(alpha**2/2)*np.exp(-alpha*alpha)/alpha
           i_malpha = -np.sqrt(np.pi/2)*erf(alpha/np.sqrt(2))
           i_z2 = np.sqrt(np.pi/2)*erf(z/np.sqrt(2))
           cdf_ret = (i_z1 - i_a) + (i_z2 - i_malpha)
        if z > alpha:
           i_z2 = np.sqrt(np.pi/2)*erf(alpha/np.sqrt(2))
           i_alpha = np.exp(-alpha**2/2)*(-1.0)/alpha
           i_z3 = np.exp(-alpha**2/2)*(-np.exp(alpha*(alpha-z)))
           cdf_ret = (i_z1 - i_a) + (i_z2 - i_malpha) + (i_z3 - i_alpha)
        
        cdf_ret = sigma*cdf_ret/self.norm
        return cdf_ret

# Create smearing functions for each pT range
a,b = -800.0,800.0
pTalphaSigmaPairs = [(10.0,1.86, 20.94),(15.0,1.86, 19.54),(20.0,1.86, 18.33),(25.0,1.86, 17.01),(35.0,1.82, 15.42),(45.0,1.66, 14.49),(60.0,1.54, 13.90),(100.0,1.64, 14.03)]
electronSmearList = []
for pT,alpha,sigma in pTalphaSigmaPairs:
  electronSmearF = smearingFunction(a=a,b=b,momtype=0)
  electronSmearF.setPars(alpha=alpha,sigma=sigma)
  electronSmearList.append((pT,electronSmearF))

def smear(inputPt, charge, sf):
    QoverPt = charge / (inputPt*1e-3) # [TeV^-1]
    QoverPtSmeared = abs(QoverPt + sf.rvs())
    PtSmeared = (1 / QoverPtSmeared) * 1e+3
    return PtSmeared

def electronSmear(inputPt, charge):
    if (inputPt < 10.0):
      return -1.0
    for pT,sf in electronSmearList:
      if inputPt > pT:
        return smear(inputPt, charge, sf)
    

#Define auxiliary functions for event cleaning and analysis

#Object readers
def getJets(branch, pTmin, etaMax, ContainsSUSY=False):
  outputjets = []  
  for ijet in range(branch.GetEntries()):
    hasSUSY=False
    jet = branch.At(ijet)
    if jet.PT < pTmin:
      # print(ijet,'pT')
      continue
    if abs(jet.Eta)>etaMax:
      # print(ijet, 'Eta', jet.Eta)
      continue
    if not ContainsSUSY:
      for ptc in jet.Particles:
        PID = 0
        try:
          PID = ptc.PID
        except:
          if isinstance(ptc,ROOT.Muon):
            # print(ijet,'mu')
            continue
          elif isinstance(ptc,ROOT.Electron):
            # print(ijet,'el')
            continue
          elif isinstance(ptc,ROOT.Photon):
            # print(ijet,'pho')
            continue
        if PID >= 1000000 and PID < 3000000:
          hasSUSY = True
          break
    if hasSUSY:
      # print(ijet)
      continue
    outputjets.append(jet)
  return outputjets

def getMuons(branch, pTmin, etaMax, mode='MuMedium'):
  outputmuons = []
  for imuon in range(branch.GetEntries()):
    muon = branch.At(imuon)
    if muon.PT < pTmin:
      continue
    if abs(muon.Eta) > etaMax:
      continue
    if mode=='MuMedium':
      pass
    outputmuons.append(muon)
  return outputmuons

def getCharginos(branch, pTmin, etaMax):
  outputchg=[]
  for ichg in range(branch.GetEntries()):
    chargino = branch.At(ichg)
    if chargino.PT < pTmin:
      continue
    if abs(chargino.Eta) > etaMax:
      continue
    outputchg.append(chargino)
  return outputchg

def getElectrons(branch, pTmin, etaMax, mode='ELooseBLLH'):
  outputelec = []
  for iel in range(branch.GetEntries()):
    elec = branch.At(iel)
    if elec.PT < pTmin:
      continue
    if abs(elec.Eta)<etaMax:
      continue
    if mode=='ELooseBLLH':
      pass
    outputelec.append(elec)
  return outputelec


#Get Kinematic variables from objects

def deltaR(ptc1,ptc2):
  lv1 = ptc1.P4() #Check if this is the proper way to read 4vector and use as input for DeltaR
  lv2 = ptc2.P4()
  return lv1.DeltaR(lv2)

def minDphi(ptc1, ptc2):
  # dphi = abs(ptc1.Phi - ptc2.Phi)
  # if dphi <= np.pi:
  #   return dphi
  # else:
  #   return 2*np.pi - dphi
  return abs(ptc1.P4().DeltaPhi(ptc2.P4()))

def minDphilist(ptc1, listptc2, length, cut):
  if len(listptc2)==0:
    return 0
  infDphi = 99999999
  for iptc,ptc2 in enumerate(listptc2):
    if iptc>=length:
      break
    if ptc2.PT<cut:
      continue
    infDphi=min(infDphi,minDphi(ptc1,ptc2))
  return infDphi

def getCharginoDecayRadius(chargino, branch):
  dau1 = branch.At(chargino.D1)
  # print(dau1.PID)
  return np.sqrt(dau1.X**2 + dau1.Y**2)

def getCharginoPLifetime(chargino, branch):
  dau1 = branch.At(chargino.D1)
  chP4 = chargino.P4()
  return 1e09*(dau1.T - chargino.T)/chP4.Gamma(),1e09*(dau1.T-chargino.T) #Assume genpart T is in sec, convert to ns

def getCharginoDecayLength_REST(chargino, branch):
  dau1 = branch.At(chargino.D1)
  return np.sqrt(dau1.X**2 + dau1.Y**2 + dau1.Z**2)/chargino.P4().Gamma()

def nmuon(muon):
  '''The missing transverse momentum is reconstructed as the negative vector sum of the transverse momenta
of photons, electrons, muons and jets, and a soft term. The soft term is reconstructed from tracks that are
associated with the hard-scatter vertex but not with any object already counted'''
  return 0,0


#Array Filters and Overlap removal
def filterMuons(DelphesTree, muons, mode='MuCaloTaggedOnly'):
  #Probably easier to define a new calomuon branch in Delphes? Not sure how to connect the Tower and Muon infos...
  return []

def overlapRemoval(input,filter,dR=0.05,mode='None'):
  if len(input)==0 or len(filter)==0:
    return input
  output=[]
  for ptc2 in filter:
    # print("Ops!")
    for ptc1 in input:
      if mode=='LessThan3Trakcs':
        ctTracks=0
        for track in ptc1.Constituents:
          if track.PT > 500:
            ctTracks+=1
        if ctTracks>2:
          continue
      if deltaR(ptc1,ptc2)>dR:
        output.append(ptc1)
  return output

def passed(ptc, condition): # Not sure how to implement most of the conditions brought onto the code
  return True

#Define SRs and Cutflow

histlist={}

def addHistogram(histlist, name, nbins, xlow, xup):
  histlist[name]=ROOT.TH1F(name,name,nbins,xlow,xup)
  histlist[name].Sumw2()

def addHistogram2D(histlist,name,nbinsx,xlow,xup,nbinsy,ylow,yup):
  histlist[name]=ROOT.TH2F(name,name,nbinsx,xlow,xup,nbinsy,ylow,yup)

def gen_fill(histlist, weight):
  def fill(histname, entry):
    histlist[histname].Fill(entry,weight)
  def fill2D(histname,entryx,entryy):
    histlist[histname].Fill(entryx,entryy,weight)
  return fill,fill2D

# Histograms
#addHistogram(histlist,"hist_MET",100,0,2000);
#addHistogram(histlist,"hist_Chargino_Pt",100,0,2000);
#addHistogram(histlist,"hist_Chargino_Pt_smeared",100,0,2000);

# Tracklet = Chargino surviving selection
#addHistogram(histlist,"hist_Tracklet_Pt",100,0,1000);
addHistogram(histlist,"hist_cutflow_kin_EWK",7,0,7); # All, GRL & Cleaning, MET Trigger, Lepton Veto, MET > 200 GeV, 1st Jet pT > 100 GeV, min(∆𝜙(Jet, MET)) > 1.0
addHistogram(histlist,"hist_cutflow_kin_Strong",7,0,7); # All, GRL & Cleaning, MET Trigger, Lepton Veto, MET > 250 GeV, 1st Jet pT > 100 GeV & 2nd & 3rd > 20 GeV, min(∆𝜙(Jet, MET)) > 0.4
addHistogram(histlist,"hist_cutflow_SR_EWK",8,0,8); # All, Kinematic, Tracklet Emulation, Leading tracklet, ∆R(jets) > 0.4,  ∆R(electron) > 0.4, ∆R(muon) > 0.4, 0.1 < |η| < 1.9
addHistogram(histlist,"hist_cutflow_SR_Strong",8,0,8); # All, Kinematic, Tracklet Emulation, Leading tracklet, ∆R(jets) > 0.4,  ∆R(electron) > 0.4, ∆R(muon) > 0.4, 0.1 < |η| < 1.9
#addHistogram(histlist,"hist_jet1pt",100,0,2000);
#addHistogram(histlist,"hist_jet2pt",100,0,2000);
#addHistogram(histlist,"hist_jet3pt",100,0,2000);
#addHistogram(histlist,"hist_jet4pt",100,0,2000);
#addHistogram(histlist,"hist_mDp_met_jets",50,0,3.15);
#addHistogram(histlist,"hist_mDp_met_jet1",50,0,3.15);
#addHistogram(histlist,"hist_mDp_met_jet2",50,0,3.15);
#addHistogram(histlist,"hist_mDp_met_jet3",50,0,3.15);
#addHistogram(histlist,"hist_mDp_met_jet4",50,0,3.15);
#addHistogram(histlist,"hist_mDp_C1_jets",50,0,3.15);
#addHistogram(histlist,"hist_mDp_C1_jet1",50,0,3.15);
#addHistogram(histlist,"hist_dR_C1_jets",100,0,5);
#addHistogram(histlist,"hist_dR_C1_jet1",100,0,5);
#addHistogram(histlist,"hist_charginoEta",50,-2.5,2.5);
#addHistogram(histlist,"hist_charginoDecayRad",50,0,500);
#addHistogram(histlist,"hist_charginoDecayLength",50,0,500);
#addHistogram(histlist,"hist_charginosPerEvent",5,0,5);
#addHistogram(histlist,"hist_chargino_smearedPT",100,0,2000);

addHistogram2D(histlist,"hist_jet1pt_met",50,0,500,50,0,500)

#Process events
import argparse
parser = argparse.ArgumentParser(description='Analyse delphesLLP output to produce efficiencies for ATLAS-ANA-SUSY-2019-018 DT search')
parser.add_argument('inputfile', metavar='inputfile_path', help='Path to the delphesLLP root file with the event sample to be analysed.')
parser.add_argument('-rw','--tau_reweighting',metavar='tau_reweighting_mval', help='Turn on lifetime reweighting for computation; default is 0 (no reweighting).',type=float,default=0)
args = parser.parse_args()
#if len(sys.argv)<1:
#   raise TypeError("No input file provided, please use 'python AtlasDT_ANA.py inputfile_path [outputfile_path]'")

#inputfile = sys.argv[1]
#if args.tau_reweighting==0:
#  def tau_reweighting(x):
#    return x
#elif args.tau_reweighting>0:
#  def tau_reweighting(x):
#    return np.exp()
if args.tau_reweighting < 0:
  raise TypeError('Invalid value for lifetime, check input again.')

inputfile = args.inputfile
inputpath = Path(inputfile)
run_tag = inputpath.parent.name

masstaupairs=pd.read_pickle('run_taupairs.pcl')
tau_base = masstaupairs.loc[run_tag]['tau_ns']

f = ROOT.TFile(inputfile,'read')
DelphesTree = f.Get('Delphes')
nevts = DelphesTree.GetEntries()

DelphesTree.GetEntry(0)

bchg = DelphesTree.bsm
bjets = DelphesTree.Jet
bmet = DelphesTree.MissingET
bmuons = DelphesTree.Muon
belec = DelphesTree.Electron
bweight = DelphesTree.Weight
bdd = DelphesTree.bsmDirectDaughters

totalweight = 0

ct=0

bkeeping=[]
err_sqr=[]

for entry in tqdm.tqdm(range(0,nevts)):

  # print("Scanning entry: ", entry)

  DelphesTree.GetEntry(entry)
  ct+=1
  charginos = getCharginos(bchg,0,5)
  met = bmet.At(0)
  metpt = met.MET
  jets = getJets(bjets,20.0,2.8,ContainsSUSY=False)

  muons = getMuons(bmuons,10.0,2.7,mode='MuMedium') # For now, didn't implement MuMedium criteria, only placeholder
  calomuons = [] # filterMuons(DelphesTree, muons, mode='MuCaloTaggedOnly') # See function definition to get issues
  notcalomuons = [muon for muon in muons if muon not in calomuons]

  electrons = getElectrons(belec, 10.0, 2.47, mode='ELooseBLLH')

  weights = bweight.At(0).Weight
  totalweight += weights

  # print("Reading complete")

  # fill = gen_fill(histlist,weights)
  weight = 1
  fill,fill2D = gen_fill(histlist,weight)

  # print("Fill method successfully defined")

  #Overlap Removal
  calomuons = overlapRemoval(calomuons, electrons, 0.05)
  electrons = overlapRemoval(electrons, notcalomuons, 0.05)
  muons = calomuons+notcalomuons

  # jets = overlapRemoval(jets, electrons, 0.2)

  electrons = overlapRemoval(electrons, jets, 0.4)

  # jets = overlapRemoval(jets, muons, 0.2, 'LessThan3Tracks')

  muons = overlapRemoval(muons, jets, 0.4)

  # print("Overlap removal successfull")

  #Sort jet list
  jets.sort(key=lambda j: j.PT,reverse=True)

  # print("Hist filling successfull")

  nCharginos = 0.5

  ##Preselections

  fill("hist_cutflow_kin_EWK",0.5);
  fill("hist_cutflow_kin_Strong",0.5);
  fill("hist_cutflow_SR_EWK",0.5);
  fill("hist_cutflow_SR_Strong",0.5);

  # print("Preselection initializing")

  #Event Cleaning

  #Clean Bad Jets
  passedbadJets=True
  if len(jets)>0:
    for jet in jets:
    #   if not passed(jet,'LooseBadJet'): # What is the definition?
    #     passedbadJets=False
    #     break
      if jet.Eta>2.4:
        passedbadJets=False
        break
    # if not passed(jets[0],'TightBadJet'):
    #   passedbadJets=False
    #   break
  if not passedbadJets:
    continue

  #Clean Bad Muons
  # passedbadMuons=True
  # if len(muons)>0:
  #   for muon in muons:
  #     if not passed(muon,'MuQoPSignificance'): # What is the definition?
  #       passedbadMuons=False
  #       break
  # if not passedbadMuons:
  #   continue

  #Clean Bad MET
  # passedbadMET=True
  # if len(muons)>0:# Analysis specifies that MET is obtained from the negative vector sum of mu, elec, jet, gamma and a soft term from tracks.
  #   for muon in muons: #I guess they wanted to check if this applied properly here? But the implementation doesn't make sense to me.
  #     nmuonPt,nmuonPhi = negmuon(muon) #Maybe they should have been summing over these?
  #     if nmuonPt/metpt * np.cos(nmuonPhi - met.Phi()) > 0.5: # And this should be outside the loop? But still, why only muons?
  #       passedbadMET=False
  #       break
  # if not passedbadMET:
  #   continue

  # print("Event cleaning successful")

  fill("hist_cutflow_kin_EWK",1.5);
  fill("hist_cutflow_kin_Strong",1.5);

  if len(jets)>0:
    jet1pt = jets[0].PT
  else:
    continue

  fill2D("hist_jet1pt_met",metpt,jet1pt)

  temp = weight

  if metpt>500:
    pass
  elif jet1pt>500:
    weight = weight * eff_trigger.reweight(metpt,499.0)
    if math.isnan(weight) or weight==0:
      continue
    fill,fill2D=gen_fill(histlist,weight)
  else:
    weight = weight * eff_trigger.reweight(metpt,jet1pt)
    if math.isnan(weight) or weight==0:
      continue
    fill,fill2D=gen_fill(histlist,weight)

  if math.isnan(weight) or weight < 0:
    print("Previous weight: ",temp," / Current Weight: ",weight)
    raise ValueError

  fill('hist_cutflow_kin_EWK',2.5)
  fill('hist_cutflow_kin_Strong',2.5)

  if len(electrons)>0 or len(muons)>0:
    continue

  fill('hist_cutflow_kin_EWK',3.5)
  fill('hist_cutflow_kin_Strong',3.5)

  #EWK-specific pre-selection
  passedKinEWK = True

  #MET Cut
  if metpt < 200:
    passedKinEWK=False
  else:
    fill('hist_cutflow_kin_EWK',4.5)

  if passedKinEWK:
    if jet1pt < 100:
      passedKinEWK=False
    else:
      fill('hist_cutflow_kin_EWK',5.5)

  if passedKinEWK:
    if (minDphiJetMET:=minDphilist(met,jets,4,50.0))>1.0:
      fill('hist_cutflow_kin_EWK',6.5)
    else:
      passedKinEWK=False

  #Strong-specific Pre-selection
  passedKinStrong=True

  #MET Cut
  if metpt < 250:
    passedKinStrong=False
  else:
    fill('hist_cutflow_kin_Strong',4.5)

  #Jet PT Cuts
  if passedKinStrong:
    if jet1pt < 100:
      passedKinStrong=False
    elif len(jets)<3:
      passedKinStrong=False
    elif jets[1].PT<20 or jets[2].PT < 20:
      passedKinStrong=False
    else:
      fill('hist_cutflow_kin_Strong',5.5)

  #min(DeltaPhi(Jet,MET)) cut
  if passedKinStrong:
    if (minDphiJetMET:=minDphilist(met,jets,4,50.0))>0.4:
      fill('hist_cutflow_kin_Strong',6.5)
    else:
      passedKinStrong=False

  #Event must pass at least one preselection
  if not (passedKinEWK or passedKinStrong):
    continue

  #Signal Regions#

  #Emulate track and tracklet reco efficiency. From the SimpleAnalysis code, there seems to be folded in acceptances here.
  EWKCharginos = []
  StrongCharginos = []

  for chargino in charginos:
    charginoEta = chargino.Eta
    charginoDecayRad = getCharginoDecayRadius(chargino,bdd)
    charginoPLifetime,charginoLabLifetime = getCharginoPLifetime(chargino,bdd)
    smearedPt = electronSmear(chargino.PT, chargino.Charge)

    if passedKinEWK:
      fill('hist_cutflow_SR_EWK',1.5)
    # if eff_track_EWK.passes(charginoEta,charginoDecayRad):
    #   EWKCharginos.append([chargino,smearedPt])
      if (tracklet_weight:=eff_track_EWK.reweight(charginoEta,charginoDecayRad))==0:
        continue
      if math.isnan(tracklet_weight):
        print("NAN weight detected: ",weight," entry number: ",entry)
        continue
      if args.tau_reweighting > 0:
        #print('Proper Lifetime: %2.6e'%charginoPLifetime)
        #print('LAB lifetime: %2.6e'%charginoLabLifetime)
        #print('Tau Base: %f\tTau rw: %f'%(tau_base,args.tau_reweighting))
        #print('Before tau reweight: %2.6e'%tracklet_weight)
        tracklet_weight = tracklet_weight*np.exp(charginoPLifetime/tau_base-charginoPLifetime/args.tau_reweighting)
        #print('After tau reweight: %2.6e\n'%tracklet_weight)
      EWKCharginos.append([chargino,smearedPt,tracklet_weight])

    if passedKinStrong:
      fill('hist_cutflow_SR_Strong',1.5)
      if eff_track_Strong.passes(charginoEta,charginoDecayRad):
        StrongCharginos.append([chargino,smearedPt])

  if len(EWKCharginos) + len(StrongCharginos) == 0:
    continue
  if len(EWKCharginos):
    EWKCharginos.sort(key=lambda x: x[1],reverse=True)
  if len(StrongCharginos):
    StrongCharginos.sort(key=lambda x: x[1],reverse=True)

  #EWK tracklet selection
  firstEWK = True
  for smearedChargino in EWKCharginos:
    chargino,smearedPt,reweight = smearedChargino[0],smearedChargino[1],smearedChargino[2]
    fill,fill2D = gen_fill(histlist,weight*reweight)
    bkeeping.append({'entry':entry,'weight':weight,'reweight':reweight,'chgEta':charginoEta,'chgDR':charginoDecayRad})
    fill("hist_cutflow_SR_EWK",2.5)
    #fill("hist_chargino_smearedPT",smearedPt)

    #20 GeV < track(let)pT
    if smearedPt < 20: #Since we've sorted by descending smeared pT, no further charginos to consider once one hits the threshold
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
    if not firstEWK: # sorted by pt - first=leading - only consider the first chargino past this point
      break
    firstEWK = False

    fill("hist_cutflow_SR_EWK",3.5)

    #DeltaR(jets) > 0.4
    char_iso_jets=True
    for jet in jets:
      if chargino.P4().DeltaR(jet.P4())<0.4:
        char_iso_jets=False
        break
    if not char_iso_jets:
      continue
    fill("hist_cutflow_SR_EWK",4.5)

    #DeltaR(electron) > 0.4
    char_iso_elec=True
    for electron in electrons:
      if chargino.P4().DeltaR(electron.P4())<0.4:
        char_iso_elec=False
        break
    if not char_iso_elec:
      continue
    fill("hist_cutflow_SR_EWK",5.5)

    #DeltaR(muon) > 0.4
    char_iso_mu = True
    for muon in muons:
      if chargino.P4().DeltaR(muon.P4())<0.4:
        char_iso_mu=False
        break
    if not char_iso_mu:
      continue
    fill("hist_cutflow_SR_EWK",6.5)

    # 0.1 < abs(eta) < 0.9
    trackletEta = abs(chargino.Eta)
    if trackletEta < 0.1:
      continue
    if trackletEta > 1.9:
      continue
    fill("hist_cutflow_SR_EWK",7.5)
    #err_sqr.append(weight*reweight)

    #Calo-veto would enter here, but also folded in the efficiency map
    #fill("hist_Tracklet_Pt",chargino.PT)

  #Strong tracklet selection
  firstStrong = True
  for smearedChargino in StrongCharginos:
    chargino,smearedPt = smearedChargino[0],smearedChargino[1]
    fill("hist_cutflow_SR_Strong",2.5)

    #20 GeV < track(let)pT
    if smearedPt < 20: #Since we've sorted by descending smeared pT, no further charginos to consider once one hits the threshold
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
    if not firstStrong: # sorted by pt - first=leading - only consider the first chargino past this point
      continue
    firstStrong = False

    fill("hist_cutflow_SR_Strong",3.5)

    #DeltaR(jets) > 0.4
    char_iso_jets=True
    for jet in jets:
      if chargino.P4().DeltaR(jet.P4())<0.4:
        char_iso_jets=False
        break
    if not char_iso_jets:
      continue
    fill("hist_cutflow_SR_Strong",4.5)

    #DeltaR(electron) > 0.4
    char_iso_elec=True
    for electron in electrons:
      if chargino.P4().DeltaR(electron.P4())<0.4:
        char_iso_elec=False
        break
    if not char_iso_elec:
      continue
    fill("hist_cutflow_SR_Strong",5.5)

    #DeltaR(muon) > 0.4
    char_iso_mu=True
    for muon in muons:
      if chargino.P4().DeltaR(muon.P4())<0.4:
        char_iso_mu=False
        break
    if not char_iso_mu:
      continue
    fill("hist_cutflow_SR_Strong",6.5)

    # 0.1 < abs(eta) < 0.9
    trackletEta = abs(chargino.Eta)
    if trackletEta < 0.1:
      continue
    if trackletEta > 1.9:
      continue
    fill("hist_cutflow_SR_Strong",7.5)

    #Calo-veto would enter here, but also folded in the efficiency map

#End of loop
print("Loop Ended! Evts analysed: ",ct,'\n')


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