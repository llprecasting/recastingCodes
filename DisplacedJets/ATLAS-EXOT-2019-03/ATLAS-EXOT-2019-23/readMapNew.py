#!/usr/bin/env python3
import yaml
from yaml.loader import SafeLoader
import numpy as np

# to use: pip install pyyaml numpy
# and configure the map names you downlaoded from HEPdata below 
mapFileNames={
  "high-ET": "HEPData-ins2043503-v3-Figure_11_of_Aux._Mat..yaml",
  "low-ET": "HEPData-ins2043503-v3-Figure_12_of_Aux._Mat..yaml",
  }

mapFileIndexRefs = {} #new!
maps = {}

# load the maps
for selection, filename in mapFileNames.items():
  fIn = open(filename)
  maps[selection] = yaml.load(fIn, Loader=SafeLoader)
  mapFileIndexRefs[selection] = []
  for i in range(len(maps[selection]['independent_variables'][0]['values'])): # new!
    mapFileIndexRefs[selection].append((maps[selection]['independent_variables'][0]['values'][i]["value"], maps[selection]['independent_variables'][1]['values'][i]["value"])) # new!

# Define the bin limits
varBin = {
       "pT" : np.array( [0, 50, 100, 200, 400, 1600], 'd'),
       "Lxy" : np.array( [0, 1.5, 2, 2.5, 3, 3.5, 3.9, np.inf],'d'),
       "Lz" : np.array( [0, 3.6, 4.2, 4.8, 5.5, 6, np.inf], 'd'),
       "decay":  { 4:0, 5:1, 6:2, 15:3} ,
       }

# and the number of bins in each relevant dimension
_nbinsLxy = len(varBin["Lxy"])-1
_nbinsPt = len(varBin["pT"]) -1
_nbinsDecay = len(varBin["decay"])

# get the global bin index
def getGlobalBinIndex(pT, eta, Lxy, Lz, decay):
  if abs(eta) <  1.5 :
    l_bin  =  np.digitize(Lxy, varBin["Lxy"])
  else:
    l_bin  = _nbinsLxy +  np.digitize(Lz, varBin["Lz"])
  l_bin -= 1

  pT_bin = np.digitize(pT, varBin["pT"])-1

  if decay not in varBin["decay"].keys():
    return -1
  else:
    decay_bin = varBin["decay"][decay]

  globalBinIndex =  l_bin* (_nbinsPt*_nbinsDecay) + pT_bin*(_nbinsDecay) + decay_bin 
  return globalBinIndex

# get a text-label from the bin index summarising the characteristics of LLPs in the bin
def getBinLabelFromIndex(globalBinIndex):
  ## # Lxyz |  pT | decay
  l_bin = globalBinIndex // (_nbinsPt*_nbinsDecay)
  globalBinIndex  = globalBinIndex % (_nbinsPt*_nbinsDecay)
  pT_bin = globalBinIndex // (_nbinsDecay)
  globalBinIndex  = globalBinIndex % (_nbinsDecay)
  decay_bin = globalBinIndex

  decayType = None
  for  k, v in varBin["decay"].items():
   if decay_bin==v: decayType=k
  decay = "decay=%d" % decayType

  pT = "pT in [%d--%d]" % (varBin["pT"][pT_bin], varBin["pT"][pT_bin+1])
  if l_bin >= _nbinsLxy:
    l ="|eta| > 1.5"
    l += ", Lz in [%.1f--%.1f]" % (varBin["Lz"][l_bin-_nbinsLxy], varBin["Lz"][l_bin-_nbinsLxy+1])
  else:
    l ="|eta| < 1.5"
    l += ", Lxy in [%.1f--%.1f]" % (varBin["Lxy"][l_bin], varBin["Lxy"][l_bin+1])

  return ", ".join ([decay, pT, l])

def queryMapFromIndices(index_1, index_2, selection="high-ET", return_labels=False):
  result = 0
  # independent variables: LLP1 index, LLP2 index.
  # The dict is organised as a list of all possible (index_1, index_2) pairs with nonzero probabilities (those with prob 0 are skipped)
  # the map is symmetric: you could interchange llp1/llp2. To avoid duplicated info to store, we only have one side of the map here,
  # so test both combinations when looking for a match
  index = -1
  if (index_1, index_2) in mapFileIndexRefs[selection]: index = mapFileIndexRefs[selection].index((index_1, index_2))
  elif (index_2, index_1) in mapFileIndexRefs[selection]: index = mapFileIndexRefs[selection].index((index_2, index_1))
  if index == -1 :
    return result
  else:
    prob = maps[selection]['dependent_variables'][0]['values'][index]["value"]
    result = prob
    label_1 = maps[selection]['dependent_variables'][1]['values'][index]["value"]
    label_2 = maps[selection]['dependent_variables'][2]['values'][index]["value"]
    if return_labels:
      result = prob, label_1, label_2
  return result

# query the map directly from the kinematics
def queryMapFromKinematics(pT1, eta1, Lxy1, Lz1, decay1, pT2, eta2, Lxy2, Lz2, decay2, selection="high-ET"):
  index_1 = getGlobalBinIndex(pT1, eta1, Lxy1, Lz1, decay1) 
  index_2 = getGlobalBinIndex(pT2, eta2, Lxy2, Lz2, decay2)
  if index_1 < 0 or index_2 < 0 :
    return 0

  else:
    return queryMapFromIndices(index_1, index_2, selection)

