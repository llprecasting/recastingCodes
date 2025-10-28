#!/usr/bin/env python3
import numpy as np
import os
from scipy.stats import rv_continuous,uniform
from scipy.special import erf
from pathlib import Path
import tqdm
import math
import pandas as pd

# Fix seed so results are reproducible!
np.random.seed(seed=123)

import ROOT

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

  def efficiency(self, x, y):
    return self.getEff(x,y)


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
    
    def smear(self, pT, charge):
      QoverPt = charge / (pT*1e-3) # [TeV^-1]
      QoverPtSmeared = abs(QoverPt + self.rvs())
      PtSmeared = (1 / QoverPtSmeared) * 1e+3
      return PtSmeared

class cutFlow(object):

  def __init__(self,name,levels,zero_weight = 0.0) -> None:
    self.name = name
    self.keys = levels[:]
    self.weights = np.array([zero_weight for _ in levels])
    self.weightsErr = np.array([zero_weight for _ in levels])
    self._current_level = 0

  def __repr__(self) -> str:
    return str(self)
  
  def __str__(self) -> str:
    return self.name

  def reset(self,to_level=0):
    self._current_level = to_level

  def fill_next(self,weight):
    """
    Add weight to the next level in the cutflow
    """
    self._current_level += 1
    self.fill(weight)

  def fill_level(self,level,weight):
    """
    Fill a specific level with the weight
    """
    if not (level in self.keys):
      raise ValueError(f"Level {level} not defined for {self}")
    clevel = self.keys.index(level)
    self.weights[clevel] += weight
    self.weightsErr[clevel] = np.sqrt(self.weightsErr[clevel]**2 + weight**2)
    
  def fill(self,weight):
    """
    Add weight to the current level in the cutflow
    """
    
    clevel = self.current_level
    if not (0 <= clevel < len(self.weights)):
      msg = f"Trying to fill the cutflow at level {clevel+1}, but it only has {len(self.weights)} levels."
      msg += " Check your cutflow definitions and if it is being properly reset."
      raise ValueError(msg)
    self.weights[clevel] += weight
    self.weightsErr[clevel] = np.sqrt(self.weightsErr[clevel]**2 + weight**2)

  def divide(self,factor):
    """
    Divide cutflow levels by factor.
    """
    for iw,w in enumerate(self.weights):
      self.weights[iw] = w/factor
      self.weightsErr[iw] = self.weightsErr[iw]/factor

  @property
  def current_level(self):
    """
    Simple method for getting the current level of the cutflow
    """

    return self._current_level

  def to_dict(self):

    cDict = {k : (w,wErr) for k,w,wErr in zip(self.keys,self.weights,self.weightsErr)}

    return cDict
  
  def to_string(self):

    d = self.to_dict()
    lines = [f"==== {self.name} ==="]
    for k,(w,wErr) in d.items():
      lines.append(f"{k} = {w:1.4e} +- {wErr:1.4e}")
    lines.append(f"===" + "="*len(self.name) + "===")
    return '\n'.join(lines)

#Initialize efficiency maps

eff_trigger = effMap('eff_trigger_average',filepath='DisappearingTrack2018-EfficiencyMaps.root')
eff_track_EWK = effMap('h_effmap_average_EWK',filepath='DisappearingTrack2018-EfficiencyMaps.root')
eff_track_Strong = effMap('h_effmap_average_Strong',filepath='DisappearingTrack2018-EfficiencyMaps.root')


# Create smearing functions for each pT range
a,b = -800.0,800.0
pTalphaSigmaPairs = [(10.0,1.86, 20.94),(15.0,1.86, 19.54),(20.0,1.86, 18.33),(25.0,1.86, 17.01),(35.0,1.82, 15.42),(45.0,1.66, 14.49),(60.0,1.54, 13.90),(100.0,1.64, 14.03)]
electronSmearList = []
for pT,alpha,sigma in pTalphaSigmaPairs:
    electronSmearF = smearingFunction(a=a,b=b,momtype=0)
    electronSmearF.setPars(alpha=alpha,sigma=sigma)
    electronSmearList.append((pT,electronSmearF))

def electronPtSmear(pT, charge):
  if (pT < 10.0):
      return -1.0
  for pT_bin,smearFunc in electronSmearList:
    if pT > pT_bin:
      return smearFunc.smear(pT, charge)
  return -1.0

#Object readers
def filterObjects(particleList, pTmin, etaMax):
  filteredParticles = []
  for ptc in particleList:
    if ptc.PT < pTmin:
      continue
    if abs(ptc.Eta) > etaMax:
      continue

    filteredParticles.append(ptc)
  
  return filteredParticles


#Get Kinematic variables from objects
def deltaR(ptc1,ptc2):
  lv1 = ptc1.P4() #Check if this is the proper way to read 4vector and use as input for DeltaR
  lv2 = ptc2.P4()
  return lv1.DeltaR(lv2)

def DeltaPhi(ptc1, ptc2):
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
    infDphi=min(infDphi,DeltaPhi(ptc1,ptc2))
  return infDphi

def overlapRemoval(input,filter,dR=0.4):
  
  if len(input)==0 or len(filter)==0:
    return input[:]
  
  output=[]
  for ptc1 in input:
    if any(deltaR(ptc1,ptc2)<dR for ptc2 in filter):
      continue
    output.append(ptc1)

  return output

def getLLPDecayRadius(llp):
  return np.sqrt(llp.daughter.X**2 + llp.daughter.Y**2)

def getLLPLifetime(llp):
  p4 = llp.P4()
  return 1e9*(llp.daughter.T - llp.T)/p4.Gamma() #Assume genpart T is in sec, convert to ns

def getLLPDecayTime(llp):
  return 1e9*(llp.daughter.T - llp.T) #Assume genpart T is in sec, convert to ns

def getLLPDecayLength_REST(llp, llp_daughter):
  return np.sqrt(llp_daughter.X**2 + llp_daughter.Y**2 + llp_daughter.Z**2)/llp.P4().Gamma()
