#!/usr/bin/env python3
import numpy as np
from scipy.interpolate import NearestNDInterpolator
from numpy import float64, ndarray
from typing import Any, Dict, List, Tuple, Union
import pyslha
import json
# Fix seed so results are reproducible!
np.random.seed(seed=123)


class effMap:
  def __init__(self, filepath: str) -> None:
    
    data = self.load_efficiency_map(filepath)
    self.setInterp(data)
    
  def load_efficiency_map(self, csv_path: str) -> Any:
    """
    Load (pT, d0) -> efficiency points from the HEPData CSV and return
    an evaluator function: efficiency(pt, d0, clip=True).
    """

    with open(csv_path, 'r') as f:
      lines = [l for l in f.readlines() if not l.startswith('#')]
    data = np.genfromtxt(lines, delimiter=",", comments="#",
                        names=True,   dtype=float)
    return data
  
  def setInterp(self, data: Any) -> Any:
    
    self.vars_limits = {var: (data[var].min(), data[var].max()) for var in data.dtype.names[:-1]}
    self.eff_label = data.dtype.names[-1]
    self.interp = NearestNDInterpolator(list(zip(data[list(self.vars_limits.keys())[0]],
                                                 data[list(self.vars_limits.keys())[1]])),
                                                 data[self.eff_label])
    
  def efficiency(self,**kwargs) -> float:

    var_values = [kwargs.get(var,None) for var in self.vars_limits.keys()]
    if any(v is None for v in var_values):
      raise ValueError(f"Missing variable(s) for efficiency map: {self.vars_limits.keys()}. Got {kwargs.keys()}")
    if any(v < self.vars_limits[var][0] or v > self.vars_limits[var][1] for var,v in zip(self.vars_limits.keys(),var_values)):
      return 0.0

    return float(self.interp(*var_values))


class cutFlow(object):

  def __init__(self,name: str,levels: List[str],zero_weight : Union[float,ndarray] = 0.0) -> None:
    self.name = name
    self.keys = levels[:]
    self.weights = np.array([zero_weight for _ in levels])
    self.weightsErr = np.array([zero_weight for _ in levels])
    self._current_level = 0

  def __repr__(self) -> str:
    return str(self)
  
  def __str__(self) -> str:
    return self.name

  def reset(self,to_level: int=0):
    self._current_level = to_level

  def fill_next(self,weight: Union[float,ndarray]):
    """
    Add weight to the next level in the cutflow
    """
    self._current_level += 1
    self.fill(weight)

  def fill_level(self,level: str,weight: Union[float,ndarray]):
    """
    Fill a specific level with the weight
    """
    if not (level in self.keys):
      raise ValueError(f"Level {level} not defined for {self}")
    clevel = self.keys.index(level)
    self.weights[clevel] += weight
    self.weightsErr[clevel] = np.sqrt(self.weightsErr[clevel]**2 + weight**2)
    
  def fill(self,weight: Union[float,ndarray]):
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

  def divide(self,factor: float):
    """
    Divide cutflow levels by factor.
    """
    for iw,w in enumerate(self.weights):
      self.weights[iw] = w/factor
      self.weightsErr[iw] = self.weightsErr[iw]/factor

  @property
  def current_level(self) -> int:
    """
    Simple method for getting the current level of the cutflow
    """

    return self._current_level

  def to_dict(self) -> Dict[str, Union[Tuple[ndarray, ndarray], Tuple[float, float]]]:

    cDict = {k : (w,wErr) for k,w,wErr in zip(self.keys,self.weights,self.weightsErr)}

    return cDict
  
  def to_string(self) -> str:

    d = self.to_dict()
    lines = [f"==== {self.name} ==="]
    for k,(w,wErr) in d.items():
      if isinstance(w,(float,int)):
        lines.append(f"{k} = {w:1.4e} +- {wErr:1.4e}")
      elif isinstance(w,(list,ndarray)):
        l = f"{k} = "
        l += ' / '.join([f"{wx:1.4e} +- {wErrx:1.4e}" for wx,wErrx in zip(w,wErr)])
        lines.append(l)
    lines.append(f"===" + "="*len(self.name) + "===")
    return '\n'.join(lines)

#Initialize efficiency maps

electron_reco = effMap(filepath="./ATLAS_data/HEPData-ins1831504-v2-csv/pt-d0electronefficiency.csv")
muon_reco = effMap(filepath="./ATLAS_data/HEPData-ins1831504-v2-csv/pt-d0muonefficiency.csv")


#Object readers
def filterObjects(particleList: Any, 
                  pTmin: float, etaMax: float) -> List[Any]:
  filteredParticles = []
  for ptc in particleList:
    if ptc.PT < pTmin:
      continue
    if abs(ptc.Eta) > etaMax:
      continue

    filteredParticles.append(ptc)
  
  return filteredParticles


#Get Kinematic variables from objects
def deltaR(ptc1,ptc2) -> float:
  lv1 = ptc1.P4() #Check if this is the proper way to read 4vector and use as input for DeltaR
  lv2 = ptc2.P4()
  return lv1.DeltaR(lv2)

def DeltaPhi(ptc1, ptc2) -> float:
  return abs(ptc1.P4().DeltaPhi(ptc2.P4()))

def minDphilist(ptc1, listptc2, length, cut) -> float:
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

def overlapRemoval(input: List[Any],filter: List[Any],dR: float=0.4) -> List[Any]:
  
  if len(input)==0 or len(filter)==0:
    return input[:]
  
  output=[]
  for ptc1 in input:
    if any(deltaR(ptc1,ptc2)<dR for ptc2 in filter):
      continue
    output.append(ptc1)

  return output

def getLLPDecayRadius(llp) -> float:
  return np.sqrt(llp.daughter.X**2 + llp.daughter.Y**2)

def getLLPDecayTime(llp) -> float:
  return 1e9*(llp.daughter.T - llp.T) #Assume genpart T is in sec, convert to ns

def getModelInfo(bannerFile : str, llpPDG : int) -> Dict[str, Union[float,int]]:

    modelInfoDict : Dict[str, Union[float,int]] = {'llpPDG' : llpPDG}
    slhaData = None
    with open(bannerFile,'r') as ff:
        data = ff.read()
        if '<slha>' in data:
            slhaData = data.split('<slha>')[1].split('</slha>')[0]
            slhaData = pyslha.readSLHA(slhaData)
        if '<MGGenerationInfo>' in data:
            mgInfo = data.split('<MGGenerationInfo>')[1].split('</MGGenerationInfo>')[0]
            for l in mgInfo.split('\n'):
                l = l.strip()
                if not  l: continue
                k,v = l.split(':')
                modelInfoDict[k.replace('#','').strip()] = float(v)
        if '<MGRunCard>' in data:
            runInfo = data.split('<MGRunCard>')[1].split('</MGRunCard>')[0]
            fields = ['custom_fcts','pt_bias_target',
                      'pt_bias_enhancement_power', 'pt_bias_min']
            for l in runInfo.split('\n'):
                l = l.strip()
                if not  l: continue
                if l[0] == '#': #skip comments
                  continue
                for field in fields:
                  if field in l:
                    value = l.split('=')[0].strip()
                    if not value:
                      value = None
                    if field != 'custom_fcts' and value is not None:
                       value = float(value)
                    modelInfoDict[field] = value
                  
    if slhaData is None:
       raise ValueError(f'Error reading banner file {bannerFile}')
    
    modelInfoDict['mLLP'] = slhaData.blocks['MASS'][llpPDG]
    if slhaData.decays:
        modelInfoDict['tau0_ns'] = (6.58212e-16/slhaData.decays[llpPDG].totalwidth)
    
    return modelInfoDict

def saveOutput(resultsDict: Dict[str, Any],outputFile: str):

  saveDict = {}
  for k,v in resultsDict.items():
    if isinstance(v,np.ndarray):
      saveDict[k] = v.tolist()
    else:
       saveDict[k] = v
  # Sort results by type of key:
  saveDict = {k : v for k,v in sorted(list(saveDict.items()), key = lambda x: isinstance(x[1],list),
                                           reverse=False)}

        
  with open(outputFile, 'w') as f:
    json.dump(saveDict,f,indent=4)

