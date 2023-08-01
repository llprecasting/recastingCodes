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
import Computation_Functions_asymmetric as cmfpa
import random
import math
import os

random.seed(123)
hep.style.use("ATLAS") # Define a style for the plots

File_selection = ["./Script_Asymetric_mZp300_mS200/Events/run_01/tag_1_pythia8_events.hepmc.gz", "./Script_Asymetric_mZp150_mS100/Events/run_01/tag_1_pythia8_events.hepmc.gz","./Script_Asymetric_mZp110_mS80/Events/run_01/tag_1_pythia8_events.hepmc.gz", "./Script_Asymetric_mZp80_mS50/Events/run_01/tag_1_pythia8_events.hepmc.gz"] # file that Pythia8 outputs

MG_File_selection = ['./Script_Asymetric_mZp300_mS200/Events/run_01/unweighted_events.lhe.gz', './Script_Asymetric_mZp150_mS100/Events/run_01/unweighted_events.lhe.gz', './Script_Asymetric_mZp110_mS80/Events/run_01/unweighted_events.lhe.gz', './Script_Asymetric_mZp80_mS50/Events/run_01/unweighted_events.lhe.gz'] # file that MadGraph outputs

#Constant
c = 3e8# Light velocity in m/s

Mass_Zp = [300,150,110, 80] # Heavy Neutral Boson mass
Mass_s = [200,100,80, 50] # LLP mass (DH = Dark Higgs)
Nevent = [10000,10000,10000,10000] # Nbr of events
Factor = [1,1,1,1]
os.system("mkdir -p Plots_High")
os.system("mkdir -p Plots_Low")

count=0
tauN=np.geomspace(0.1,1e2,60) # New lifetime range
for file_selection, MG_file_selection,mass_Zp,mass_s, nevent, factor in zip(File_selection, MG_File_selection, Mass_Zp,Mass_s, Nevent, Factor):

    #Pythia
    events = hepmc.HEPMC_EventFile(file_selection) # Open HEPMC file
    px_TOT_DH, py_TOT_DH, pz_TOT_DH, E_TOT_DH, mass_TOT_DH,pdg_TOT_DH = cmfpa.parsing_hepmc_DH(events) # Parsing the HEPMC file
    events = hepmc.HEPMC_EventFile(file_selection) # Open HEPMC file
    px_TOT_Zp, py_TOT_Zp, pz_TOT_Zp, E_TOT_Zp, mass_TOT_Zp, pdg_TOT_Zp = cmfpa.parsing_hepmc_Zp(events) # Parsing the HEPMC file
    P8_px_Zp, P8_py_Zp, P8_pz_Zp, P8_E_Zp, P8_mass_Zp, P8_pdg_Zp = cmfpa.conversion_one_list_Zp(px_TOT_Zp, py_TOT_Zp, pz_TOT_Zp, E_TOT_Zp, mass_TOT_Zp, pdg_TOT_Zp) # Obtaining data in one list
    P8_px_DH, P8_py_DH, P8_pz_DH, P8_E_DH, P8_mass_DH, P8_pdg_DH = cmfpa.conversion_one_list_DH(px_TOT_DH, py_TOT_DH, pz_TOT_DH, E_TOT_DH, mass_TOT_DH, pdg_TOT_DH)
    beta_DH, gamma_DH, pT_DH, eta_DH = cmfpa.kinematics_DH(P8_px_DH, P8_py_DH, P8_pz_DH, P8_E_DH) # Computing kinematics for DH1
    beta_Zp, gamma_Zp, pT_Zp, eta_Zp = cmfpa.kinematics_Zp(P8_px_Zp, P8_py_Zp, P8_pz_Zp, P8_E_Zp) # Computing kinematics for DH2
    Lxy_tot_DH, Lz_tot_DH = cmfpa.decaylenghtDH(P8_px_DH, P8_py_DH, P8_pz_DH, P8_E_DH, gamma_DH, tauN) # Computing the decay lenght for DH1
    Lxy_tot_Zp, Lz_tot_Zp = cmfpa.decaylenghtZp(P8_px_Zp, P8_py_Zp, P8_pz_Zp, P8_E_Zp, gamma_Zp, tauN) # Computing the decay lenght for DH2

    #MG
    MG_events = lhe.EventFile(MG_file_selection) # Open LHE file
    MG_px_DH, MG_py_DH, MG_pz_DH, MG_pdg_DH, MG_E_DH, MG_MASS_DH = cmfpa.parsing_LHE_DH(MG_events) #Parsing the LHE file
    MG_events = lhe.EventFile(MG_file_selection) # Open LHE file
    MG_px_Zp, MG_py_Zp, MG_pz_Zp, MG_pdg_Zp, MG_E_Zp, MG_MASS_Zp = cmfpa.parsing_LHE_Zp(MG_events) #Parsing the LHE file
    MG_pT_DH,MG_eta_DH, MG_gamma_DH = cmfpa.kinematics_MG_DH(MG_px_DH,MG_py_DH,MG_pz_DH,MG_E_DH) # Computing kinematics for DH1
    MG_pT_Zp,MG_eta_Zp, MG_gamma_Zp = cmfpa.kinematics_MG_Zp(MG_px_Zp,MG_py_Zp,MG_pz_Zp,MG_E_Zp) # Computing kinematics for DH2
    MG_Lxy_tot_DH, MG_Lz_tot_DH = cmfpa.decaylenght_MG_DH(MG_px_DH, MG_py_DH, MG_pz_DH, MG_E_DH, MG_gamma_DH, tauN) # Computing decay lenght for DH1
    MG_Lxy_tot_Zp, MG_Lz_tot_Zp = cmfpa.decaylenght_MG_Zp(MG_px_Zp, MG_py_Zp, MG_pz_Zp, MG_E_Zp, MG_gamma_Zp, tauN) # Computing decay lenght for DH2

############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################

###########################################################Computing the efficiencies and ploting the results###########################################################

    if mass_Zp >= 150: # Condition if the sample is 'High-ET' or ' Low-ET'
        MG_eff_highETX = cmfpa.eff_map_MG_high(MG_pT_DH, MG_eta_DH,MG_Lxy_tot_DH, MG_Lz_tot_DH, MG_pdg_DH, MG_pT_Zp, MG_eta_Zp, MG_Lxy_tot_Zp, MG_Lz_tot_Zp, MG_pdg_Zp, tauN, nevent, mass_Zp, mass_s) # Compute the efficiency from MG
        eff_highETX = cmfpa.eff_map_High(pT_DH, eta_DH, Lxy_tot_DH, Lz_tot_DH, P8_pdg_DH, pT_Zp, eta_Zp, Lxy_tot_Zp, Lz_tot_Zp, P8_pdg_Zp, tauN, nevent, mass_Zp, mass_s) # Compute the efficiency from Pythia
        cmfpa.plt_eff_high(MG_eff_highETX, eff_highETX,tauN,  mass_Zp , mass_s) # Ploting and saving a comparison of all the results of efficiencies
        cmfpa.plt_cross_High(eff_highETX, tauN, mass_Zp, mass_s, factor)# Ploting and saving a comparison of the limits obtained with the map and by ATLAS.
    else:
        MG_eff_lowETX = cmfpa.eff_map_MG_low(MG_pT_DH, MG_eta_DH,MG_Lxy_tot_DH, MG_Lz_tot_DH, MG_pdg_DH, MG_pT_Zp, MG_eta_Zp, MG_Lxy_tot_Zp, MG_Lz_tot_Zp, MG_pdg_Zp, tauN, nevent, mass_Zp, mass_s)
        eff_lowETX = cmfpa.eff_map_Low(pT_DH, eta_DH, Lxy_tot_DH, Lz_tot_DH, P8_pdg_DH, pT_Zp, eta_Zp, Lxy_tot_Zp, Lz_tot_Zp, P8_pdg_Zp, tauN,nevent, mass_Zp, mass_s)
        cmfpa.plt_eff_low(MG_eff_lowETX, eff_lowETX, tauN, mass_Zp, mass_s )
        cmfpa.plt_cross_Low(eff_lowETX, tauN, mass_Zp, mass_s, factor)
    count+=1
