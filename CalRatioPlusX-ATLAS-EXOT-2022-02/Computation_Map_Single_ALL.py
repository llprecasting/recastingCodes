#!/usr/bin/env python
import sys
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import scipy
import time
import readMapNew as rmN
import tqdm
import mplhep as hep
import lhe_parser as lhe
import hepmc_parser as hepmc
import Computation_Functions as cmfp
import random
import math
import os
import pandas as pd
import copy
import glob

#mass_mediator = int(sys.argv[1])
#mass_LLP = int(sys.argv[2])
mode  = sys.argv[1]
params = sys.argv[2]
mass_mediator = -1
mass_LLP = -1
mass_LLP2 = None
LLP2=None
model = sys.argv[3]
genMode="py8"

if mode=="W": vectorBoson=24
elif mode=="Z": vectorBoson=23
elif mode=="gg": vectorBoson=-1
else: 
  print ("UNKNOWN MODE ", mode)
  exit(1)

if model=="ALP":
  mass_LLP = float(params)
  params=f"{float(params):.6g}"
  slug = f'mALP{params}_{mode}' 
  mediator=-999 
  LLP=9000005
elif model=="HSS":
  mass_mediator = int(params.split("_")[0])
  mass_LLP = int(params.split("_")[1])
  slug = f'{mode}H{mass_mediator}_S{mass_LLP}_{mode}' 
  mediator=25 
  LLP=35
elif model=="HS1S2":
  mass_mediator = int(params.split("_")[0])
  mass_LLP = int(params.split("_")[1])
  mass_LLP2 = int(params.split("_")[2])
  slug = f'{mode}H{mass_mediator}_S1{mass_LLP}_S2{mass_LLP2}_{mode}' 
  mediator=25 
  LLP=35
  LLP2=36
elif model=="ZZd":
  mass_mediator = int(params.split("_")[0])
  mass_LLP = int(params.split("_")[1])
  slug = f'mHZZd_{mass_mediator}_{mass_LLP}_{mode}' 
  mediator=25 
  LLP=1023
else:
  print("unknown model", model, "... exit")
  exit(1)

InDir = sys.argv[4]
OutDir = sys.argv[5]

hasHEPData=0

random.seed(124)
hep.style.use("ATLAS") # Define a style for the plots

#Path Pythia8 file
file_selection = f"{OutDir}/Script_{slug}/Events/run_01/tag_1_pythia8_events.hepmc.gz"
MG_file_selection = f"{OutDir}/Script_{slug}/Events/run_01/unweighted_events.lhe.gz"
#factor = 0.048 if mass_mediator==125 else 1

#Constant
c = 3e8# Light velocity in m/s

os.system("mkdir -p Plots")

#tauN=np.geomspace(0.1,1e2, 30)
tauN=np.geomspace(1e-4,1e2, 20)
#tauN=[0.309]
iTauN = list(tauN).index(tauN[tauN > 1][0])

for gen in [genMode]:
  if gen=='py8':
    events = hepmc.HEPMC_EventFile(file_selection) # Open HEPMC file
    res = cmfp.parsing_hepmc_generic(events, mediator=mediator, LLP=LLP, LLP2=LLP2, vectorBoson=vectorBoson) # Parsing the HEPMC file
  else:
    events = lhe.EventFile(MG_file_selection)
    res = cmfp.parsing_LHE(events, mediator=mediator, LLP=LLP, vectorBoson=vectorBoson) # Parsing the HEPMC file
  keys = list(res.keys())
  Lxy = {}
  Lz = {}
  for k in keys:
    df = res[k]
    print (f"---- {k}{gen} ----", len(df))
    print(df.head())
    res[k]=cmfp.kinematics(df)
    if 'llp' in k: 
      thisLxy, thisLz = cmfp.decayLength(df, tauN)
      Lxy[k] = thisLxy
      Lz[k] = thisLz
      df[f'{k}_Lxy']=Lxy[k][0]
      df[f'{k}_Lz']=Lz[k][0]
  
  if len(res['llp2'])==0: 
    res['llp2']=copy.deepcopy(res['llp1'] )
    Lxy['llp2']=copy.deepcopy(Lxy['llp1'] )
    Lz['llp2']=copy.deepcopy(Lz['llp1'] )
    print ("LC DEBUG SETTING LLP1=LLP2")
  if len(res['V'])==0: 
    res['V']["V_pt"]=np.zeros(len(res['llp2']))
    res['V']["V_eta"]=np.zeros(len(res['llp2']))
    print ("LC DEBUG SETTING V")
  
  values={
  "llp1_eta" : res['llp1']['eta'],
  "llp1_pT" : res['llp1']['pT'] ,
  "llp1_ET" : np.sqrt(res['llp1']['pT']**2 + res['llp1']['mass']**2), 
  "llp1_Lxy" : Lxy['llp1'],
  "llp1_Lz" : Lz['llp1'], 
  'llp1_child_pdgId' : abs(res['llp1']['decay']), 
  "llp2_eta" : res['llp2']['eta'],
  "llp2_pT" : res['llp2']['pT'] ,
  "llp2_ET" : np.sqrt(res['llp2']['pT']**2 + res['llp2']['mass']**2), 
  "llp2_Lxy" : Lxy['llp2'], 
  "llp2_Lz" : Lz['llp2'],
  'llp2_child_pdgId' : abs(res['llp2']['decay']), 
  "V_pt" : res['V']['pT'],
  "V_eta" : res['V']['eta'],
  "Phi_eta" : res['Phi']['eta'],
  "Phi_pT" : res['Phi']['pT'] ,
  }
  keys_to_fix = []
  for k, v in  values.items():
    print (k, len(v))
    if len(v) == 0 :
      keys_to_fix+=[k]
  for k in keys_to_fix:
    values[k] = np.ones(len(res['llp1']['eta']))
    
  for k, v in  values.items():
    print (k, len(v))
    if len(v) == 0 :
      keys_to_fix+=[k]
  
  #values_specific = copy.deepcopy(values)
  #values_specific['llp1_Lxy'] = values['llp1_Lxy'][iTauN]
  #values_specific['llp2_Lxy'] = values['llp2_Lxy'][iTauN]
  #values_specific['llp1_Lz'] = values['llp1_Lz'][iTauN]
  #values_specific['llp2_Lz'] = values['llp2_Lz'][iTauN]
  #pd.DataFrame(values_specific).to_csv(f"Plots/{slug}_{gen}_ct{tauN[iTauN]}.csv")
  #if gen=='mg5': continue  
  if mode=="W":
    sels = ["WHS_highET","WHS_lowET", "WALP"]
  elif mode=="Z":
    sels = ["ZHS_highET","ZHS_lowET"]
  else:
    sels = ["CR+2J"]
  
  for sel in sels:
    eff = cmfp.eff_bdt_tauN(values,tauN, sel=sel) # Compute the efficiency from Pythia
    cmfp.plt_eff( eff , tauN, mass_mediator, mass_LLP, model=model, sel=sel, mass_s2=mass_LLP2) # Ploting and saving a comparison of all the results of efficiencies
    eff = cmfp.eff_bdt_tauN(values,tauN, sel=sel, do2D=True) # Compute the efficiency from Pythia
    cmfp.plt_eff2D( eff , tauN, mass_mediator, mass_LLP, model=model, sel=sel, mass_s2=mass_LLP2) # Ploting and saving a comparison of all the results of efficiencies
  
