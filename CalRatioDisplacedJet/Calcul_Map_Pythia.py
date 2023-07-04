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
import Calcul_Map_function_Pythia as cmfp
import random
import math
import uproot
import yaml
from yaml.loader import SafeLoader
random.seed(123)
hep.style.use("ATLAS") # Define a style for the plots


File_selection = ["./Script_mH1000_mS275/Events/run_01/tag_1_pythia8_events.hepmc.gz", "./Script_mH600_mS150/Events/run_01/tag_1_pythia8_events.hepmc.gz","./Script_mH400_mS100/Events/run_01/tag_1_pythia8_events.hepmc.gz", "./Script_mH200_mS50/Events/run_01/tag_1_pythia8_events.hepmc.gz", "./Script_mH125_mS55/Events/run_01/tag_1_pythia8_events.hepmc.gz"] # file that Pythia outputs

MG_File_selection = ['./Script_mH1000_mS275/Events/run_01/unweighted_events.lhe.gz', './Script_mH600_mS150/Events/run_01/unweighted_events.lhe.gz', './Script_mH400_mS100/Events/run_01/unweighted_events.lhe.gz', './Script_mH200_mS50/Events/run_01/unweighted_events.lhe.gz' ,'./Script_mH125_mS55/Events/run_01/unweighted_events.lhe.gz' ] # file that MadGraph outputs

File_HEP = ["./ATLAS_data/HEPData-ins2043503-v3-Figure_2e_of_Aux._Mat._1000_275.root","./ATLAS_data/HEPData-ins2043503-v3-Figure_2d_of_Aux._Mat._600_150.root","./ATLAS_data/HEPData-ins2043503-v3-Figure_2c_of_Aux._Mat._400_100.root", "./ATLAS_data/HEPData-ins2043503-v3-Figure_2c_of_Aux._Mat._200_50.root", "./ATLAS_data/HEPData-ins2043503-v3-Figure_2b_of_Aux._Mat._125_55.root"] #HEP data files

Branch_HEP = ["Figure 2e of Aux. Mat. 1000_275/Graph1D_y1;1",'Figure 2d of Aux. Mat. 600_150/Graph1D_y1;1', 'Figure 2c of Aux. Mat. 400_100/Graph1D_y1;1', 'Figure 2c of Aux. Mat. 200_50/Graph1D_y1;1','Figure 2b of Aux. Mat. 125_55/Graph1D_y1;1'] # Branch from HEP data

#Constante
c = 3e8# Light velocity in m/s
tauN=np.geomspace(0.1,1e2,50) # New lifetime range

Mass_phi = [1000,600,400, 200, 125] # Heavy Neutral Boson mass
Mass_s = [275,150,100, 50, 55] # LLP mass (DH = Dark Higgs)

count=0
for fichier_selection, MG_fichier_selection,mass_phi,mass_s in zip(Fichier_selection, MG_Fichier_selection, Mass_H,Mass_s):

    #Pythia
    events = hepmc.HEPMC_EventFile(fichier_selection) # Open HEPMC file
    px_TOT, py_TOT, pz_TOT, E_TOT, mass_TOT,pdg_TOT = cmfp.parsing_hepmc(events) # Parsing the HEPMC file
    px_tot, py_tot, pz_tot, E_tot, mass_tot, pdg_tot = cmfp.conversion_one_list(px_TOT, py_TOT, pz_TOT, E_TOT, mass_TOT, pdg_TOT) # Obtaining data in one list
    px_DH1, px_DH2, py_DH1, py_DH2, pz_DH1, pz_DH2, pdg_tot_DH1, pdg_tot_DH2, E_DH1, E_DH2, mass_DH1, mass_DH2 = cmfp.recover(px_tot, py_tot, pz_tot, E_tot, mass_tot, pdg_tot) # Separate data from DH1 and DH2
    beta_DH1, gamma_DH1, pT_DH1, eta_DH1 = cmfp.kinematics_DH1(px_DH1, py_DH1, pz_DH1, E_DH1) # Computing kinematics for DH1
    beta_DH2, gamma_DH2, pT_DH2, eta_DH2 = cmfp.kinemamtics_DH2(px_DH2, py_DH2, pz_DH2, E_DH2) # Computing kinematics for DH2
    Lxy_tot_DH1, Lz_tot_DH1 = cmfp.decaylenghtDH1(px_DH1, py_DH1, pz_DH1, E_DH1, gamma_DH1, tauN) # Computing the decay lenght for DH1
    Lxy_tot_DH2, Lz_tot_DH2 = cmfp.decaylenghtDH2(px_DH2, py_DH2, pz_DH2, E_DH2, gamma_DH2, tauN) # Computing the decay lenght for DH2

    #MG
    MG_events = lhe.EventFile(MG_fichier_selection) # Open LHE file
    px, py, pz, pdg, E, MASS = cmfp.parsing_LHE(MG_events) #Parsing the LHE file
    MG_px_DH1, MG_py_DH1,MG_pz_DH1,MG_E_DH1,MG_mass_DH1,MG_pdg_DH1_1 = cmfp.recover_MG_DH1(px, py, pz, E, MASS, pdg) # Separate data from DH1 and DH2
    MG_pT_DH1,MG_eta_DH1, MG_gamma_DH1 = cmfp.kinematics_MG_DH1(MG_px_DH1,MG_py_DH1,MG_pz_DH1,MG_E_DH1) # Computing kinematics for DH1
    MG_px_DH2, MG_py_DH2,MG_pz_DH2,MG_E_DH2,MG_mass_DH2,MG_pdg_DH2_1 = cmfp.recover_MG_DH2(px, py, pz, E, MASS, pdg) # Separate data from DH1 and DH2
    MG_pT_DH2,MG_eta_DH2, MG_gamma_DH2 = cmfp.kinemamtics_MG_DH2(MG_px_DH2,MG_py_DH2,MG_pz_DH2,MG_E_DH2) # Computing kinematics for DH2
    MG_Lxy_tot_DH1, MG_Lz_tot_DH1 = cmfp.decaylenght_MG_DH1(MG_px_DH1, MG_py_DH1, MG_pz_DH1, E_DH1, MG_gamma_DH1, tauN) # Computing the decay lenght for DH1
    MG_Lxy_tot_DH2, MG_Lz_tot_DH2 = cmfp.decaylenght_MG_DH2(MG_px_DH2, MG_py_DH2, MG_pz_DH2, E_DH2, MG_gamma_DH2, tauN) # Computing the decay lenght for DH2

    #HEP data
    data_HEP = cmfp.elem_list(Fichier_HEP[count], Branche_HEP[count]) # Recover public data from ATLAS to compare the results

############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################

###########################################################Computing the efficiencies and ploting the results###############################################################

    if mass_phi >= 400: # Condition if the sample is 'High-ET' or ' Low-ET'
        MG_eff_highETX = cmfp.eff_map_MG_high(MG_pT_DH1, MG_eta_DH1,MG_Lxy_tot_DH1, MG_Lz_tot_DH1, MG_pdg_DH1_1, MG_pT_DH2, MG_eta_DH2, MG_Lxy_tot_DH2, MG_Lz_tot_DH2, MG_pdg_DH2_1, tauN) # Compute the efficiency from MG
        eff_highETX = cmfp.eff_map_High(pT_DH1, eta_DH1, Lxy_tot_DH1, Lz_tot_DH1, abs(pdg_tot_DH1), pT_DH2, eta_DH2, Lxy_tot_DH2, Lz_tot_DH2, abs(pdg_tot_DH2), tauN) # Compute the efficiency from Pythia
        cmfp.plt_eff_high(MG_eff_highETX, eff_highETX, tauN, data_HEP, mass_H, mass_DH ) # Ploting and saving a comparison of all the results of efficiencies
    else:
        MG_eff_lowETX = cmfp.eff_map_MG_low(MG_pT_DH1, MG_eta_DH1,MG_Lxy_tot_DH1, MG_Lz_tot_DH1, MG_pdg_DH1_1, MG_pT_DH2, MG_eta_DH2, MG_Lxy_tot_DH2, MG_Lz_tot_DH2, MG_pdg_DH2_1, tauN)
        eff_lowETX = cmfp.eff_map_Low(pT_DH1, eta_DH1, Lxy_tot_DH1, Lz_tot_DH1, abs(pdg_tot_DH1), pT_DH2, eta_DH2, Lxy_tot_DH2, Lz_tot_DH2, abs(pdg_tot_DH2), tauN)
        cmfp.plt_eff_low(MG_eff_lowETX, eff_lowETX, tauN, data_HEP, mass_H, mass_DH )
    count+=1
