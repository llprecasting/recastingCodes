#!/usr/bin/env python3
import os

# The number of events you want to generate for each sample
nevent = 50

#The Phi and S particles masses models you want to generate, i.e [m_Phi, m_S]
#WARNING: If you plan to use your own machine, run the different masses one by one
#or consider using the "Writing_Scripts_MG+P8.py" code (that loop over the chosen masses)
masses = [                         
          #["W", "ALP", 0.1],
          #["W", "ALP", 0.4],
          ["W", "ALP", 1],
          #["W", "ALP", 4],
#          ["W", "ALP", 10],
#          ["W", "ALP", 40],
#          #["W", "ALP", 100],
#          #["Z", "ALP", 0.1],
#          #["Z", "ALP", 0.4],
#          ["Z", "ALP", 1],
#          #["Z", "ALP", 4],
#          ["Z", "ALP", 10],
#          ["Z", "ALP", 40],
#          #["Z", "ALP", 100],
#          ["Z", "HSS", "1000_475"], 
#          ["Z", "HSS", "1000_275"], 
#          ["Z", "HSS", "1000_50"], 
#          ["Z", "HSS", "600_275"], 
#          ["Z", "HSS", "600_150"], 
#          ["Z", "HSS", "600_50"], 
#          ["Z", "HSS", "400_175"], 
##          ["Z", "HSS", "400_100"], 
#          ["Z", "HSS", "125_55"], 
#          ["Z", "HSS", "125_35"], 
#          ["Z", "HSS", "125_16"], 
#          #["Z", "HSS", "125_5"], 
#          ["Z", "HSS", "60_15"], 
#          #["Z", "HSS", "60_5"], 
#          ["Z", "ZZd", "600_400"], 
#          ["Z", "ZZd", "600_150"], 
#          ["Z", "ZZd", "400_200"], 
#          ["Z", "ZZd", "400_100"], 
#          ["Z", "ZZd", "250_100"], 
#          ["Z", "ZZd", "250_50"], 
#          ["W", "HSS", "1000_475"], 
#          ["W", "HSS", "1000_275"], 
#          ["W", "HSS", "1000_50"], 
#          ["W", "HSS", "600_275"], 
#          ["W", "HSS", "600_150"], 
#          ["W", "HSS", "600_50"], 
#          ["W", "HSS", "400_175"], 
#          ["W", "HSS", "400_100"], 
#           ["W", "HSS", "200_50"], 
#          ["W", "HSS", "125_55"], 
#          ["W", "HSS", "125_35"], 
#          ["W", "HSS", "125_16"], 
#          #["W", "HSS", "125_5"], 
#          ["W", "HSS", "60_15"], 
#          #["W", "HSS", "60_5"], 
#          ["W", "HSS", "200_50"], 
#          ["W", "HSS", "200_80"], 
#          ['gg', 'HSS', '200_50'],
##          ['gg', 'HS1S2', '200_50_50'],
#          ['gg', 'HSS', '400_100'],
##          ['gg', 'HS1S2', '400_100_100'],
#          ['gg', 'HSS', '600_150'],
#          ['gg', 'HSS', '125_55'],
##          ['gg', 'HS1S2', '125_55_55'],
##          ['gg', 'HSS', '125_16'],
##          ['gg', 'HS1S2', '125_16_16'],
##          ['gg', 'HSS', '60_15'],
##          ['gg', 'HS1S2', '60_15_15'],
          ]

#for m1 in range(0,350, 50):
#  for m2 in range(0,350, 50):
#    masses.append(['gg', 'HS1S2', f'600_{m1}_{m2}'])
#for m1 in range(0,150, 30):
#  for m2 in range(0,150, 30):
#    masses.append(['gg', 'HS1S2', f'200_{m1}_{m2}'])
#for m1 in range(0,80, 10):
#  for m2 in range(0,80, 10):
#    masses.append(['gg', 'HS1S2', f'125_{m1}_{m2}'])



#The full path to where the MadGraph folder is!
#InDir = "/users/divers/atlas/haddad/home2/MG5_aMC_v3_4_2"  
InDir = "/AtlasDisk/user/corpe/LLPRecasting2024/MG5_aMC_v3_4_2"  
workDir = "/AtlasDisk/user/corpe/LLPRecasting2024/MG5_aMC_v3_4_2/new/recastingCodes/CalRatioPlusX"  
mgDir =  InDir #"/AtlasDisk/user/corpe/LLPRecasting2024/MG5_aMC_v2_9_19"  

#The full path to where you want all your outputs to be saved 
#OutDir = "/users/divers/atlas/haddad/scratch/Recasting"  
OutDir = "/users/divers/atlas/corpe/scratch/RecastingScratch"  


for mode, model, params in masses:
     
     
    if model=="ALP":
      params=f"{params:.6g}"
      slug = f'mALP{params}_{mode}' 
    elif model=="HSS":
      mass_Phi = int(params.split("_")[0])
      mass_S = int(params.split("_")[1])
      slug = f'{mode}H{mass_Phi}_S{mass_S}_{mode}' 
    elif model=="HS1S2":
      mass_Phi = int(params.split("_")[0])
      mass_S1 = int(params.split("_")[1])
      mass_S2 = int(params.split("_")[2])
      slug = f'{mode}H{mass_Phi}_S1{mass_S1}_S2{mass_S2}_{mode}' 
    elif model=="ZZd":
      mass_Phi = int(params.split("_")[0])
      mass_Zd = int(params.split("_")[1])
      slug = f'mHZZd_{mass_Phi}_{mass_Zd}_{mode}' 
    else:
      print("unknown model", model, "... exit")
      exit(1)
  
    f = open(f"{OutDir}/Job_{slug}.sh",'w')
    
    #These two lines will setup ATLAS and Athena in order to use an updated verision of python (>= 3.7)
    #To be removed if you haven't access to cvmfs, or simply if you are using your up-to-date python  
    f.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
    f.write("asetup AnalysisBase,24.2.35\n")
    #f.write(f"cd {InDir}\n")
    #f.write(f"pip install numpy -t $PWD\n")
    f.write(f"export PYTHONPATH={InDir}:$PYTHONPATH\n")
    f.write(f"export TQDM_DISABLE=1\n")
    #Or avtivate your generated python environment
    ##f.write(f"{InDir}/env/Scripts/activate")
    
    f.write(f"python3 {workDir}/Writing_Scripts_MG+P8_Single_{model}.py {mode} {params} {nevent} {InDir} {OutDir}\n")
    f.write(f"cd {OutDir}\n")
    f.write(f"{mgDir}/bin/mg5_aMC -f {OutDir}/script_{slug}.txt\n")
    f.write(f"cd {workDir}/\n")
    f.write(f"python3 {workDir}//Computation_Map_Single_ALL.py {mode} {params} {model} {InDir} {OutDir} \n")
    f.close()
    
    os.system(f"chmod +x {OutDir}/Job_{slug}.sh")
    
    #Depending on where you want to run MadGraph, use the first command line to run locally
    #Or the second one to run in your lab servers (to be updated accordingly, here is for LPC)
    ##os.system(f"{InDir}/bin/mg5_aMC -f {OutDir}/script_mH{mass_Phi}_mS{mass_S}.txt")
    #os.system(f"qsub -q prod2C7@clratlserv04 -o {OutDir} -e {OutDir} {OutDir}/Job_{slug}.sh")
    print(f"{OutDir}/Job_{slug}.sh")
    
    
    
