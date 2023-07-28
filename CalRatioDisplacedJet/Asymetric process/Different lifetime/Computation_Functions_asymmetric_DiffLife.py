 
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
import uproot
import random
import os
random.seed(123)
hep.style.use("ATLAS")

# set cst
c = 3e8 # Light velocity in m/s

#########################################################################################
#Parsing the hepmc file from the hadronization of the MG outputs to recover the data from the process.
#########################################################################################

def parsing_hepmc_DH(events):

    px_TOT_DH = []
    py_TOT_DH = []
    pz_TOT_DH = []
    E_TOT_DH = []
    mass_TOT_DH = []
    pdg_TOT_DH = []

    for ie , event in enumerate(events):
        count=0
        for id, vertex in event.vertex.items():
            if [p.pdg for p in vertex.incoming] == [35] and [p.pdg for p in vertex.outcoming] == [1, -1] or [p.pdg for p in vertex.incoming] == [35] and [p.pdg for p in vertex.outcoming] == [2, -2] or [p.pdg for p in vertex.incoming] == [35] and [p.pdg for p in vertex.outcoming] == [3, -3] or [p.pdg for p in vertex.incoming] == [35] and [p.pdg for p in vertex.outcoming] == [4, -4] or [p.pdg for p in vertex.incoming] == [35] and [p.pdg for p in vertex.outcoming] == [5, -5] or [p.pdg for p in vertex.incoming] == [35] and [p.pdg for p in vertex.outcoming] == [6, -6] or [p.pdg for p in vertex.incoming] == [35] and [p.pdg for p in vertex.outcoming] == [11, -11] or [p.pdg for p in vertex.incoming] == [35] and [p.pdg for p in vertex.outcoming] == [13, -13] or [p.pdg for p in vertex.incoming] == [35] and [p.pdg for p in vertex.outcoming] == [15, -15]:
                px_TOT_DH.append(list(p.px for p in vertex.incoming))
                py_TOT_DH.append(list(p.py for p in vertex.incoming))
                pz_TOT_DH.append(list(p.pz for p in vertex.incoming))
                E_TOT_DH.append(list(p.E for p in vertex.incoming))
                mass_TOT_DH.append(list(p.mass for p in vertex.incoming))
                pdg_TOT_DH.append(list(p.pdg for p in vertex.outcoming))


    return px_TOT_DH, py_TOT_DH, pz_TOT_DH, E_TOT_DH, mass_TOT_DH, pdg_TOT_DH

def parsing_hepmc_Zp(events):
    px_TOT_Zp = []
    py_TOT_Zp = []
    pz_TOT_Zp = []
    E_TOT_Zp = []
    mass_TOT_Zp = []
    pdg_TOT_Zp = []

    for ie , event in enumerate(events):
        count=0
        for id, vertex in event.vertex.items():
            if [p.pdg for p in vertex.incoming] == [1023] and [p.pdg for p in vertex.outcoming] == [1, -1] or [p.pdg for p in vertex.incoming] == [1023] and [p.pdg for p in vertex.outcoming] == [2, -2] or [p.pdg for p in vertex.incoming] == [1023] and [p.pdg for p in vertex.outcoming] == [3, -3] or [p.pdg for p in vertex.incoming] == [1023] and [p.pdg for p in vertex.outcoming] == [4, -4] or [p.pdg for p in vertex.incoming] == [1023] and [p.pdg for p in vertex.outcoming] == [5, -5] or [p.pdg for p in vertex.incoming] == [1023] and [p.pdg for p in vertex.outcoming] == [6, -6] or [p.pdg for p in vertex.incoming] == [1023] and [p.pdg for p in vertex.outcoming] == [11, -11] or [p.pdg for p in vertex.incoming] == [1023] and [p.pdg for p in vertex.outcoming] == [13, -13] or [p.pdg for p in vertex.incoming] == [1023] and [p.pdg for p in vertex.outcoming] == [15, -15]:
                px_TOT_Zp.append(list(p.px for p in vertex.incoming))
                py_TOT_Zp.append(list(p.py for p in vertex.incoming))
                pz_TOT_Zp.append(list(p.pz for p in vertex.incoming))
                E_TOT_Zp.append(list(p.E for p in vertex.incoming))
                mass_TOT_Zp.append(list(p.mass for p in vertex.incoming))
                pdg_TOT_Zp.append(list(p.pdg for p in vertex.outcoming))

    return px_TOT_Zp, py_TOT_Zp, pz_TOT_Zp, E_TOT_Zp, mass_TOT_Zp, pdg_TOT_Zp

#########################################################################################
#The data recovered are list of list, we need to convert them into one list to be able to separate the contribution of each LLP.
#########################################################################################

def conversion_one_list_Zp(px_TOT_Zp, py_TOT_Zp, pz_TOT_Zp, E_TOT_Zp, mass_TOT_Zp, pdg_TOT_Zp):

    P8_px_Zp = []
    for i in range(len(px_TOT_Zp)):
        for y in range(len(px_TOT_Zp[i])):
            P8_px_Zp.append(px_TOT_Zp[i][y])

    P8_py_Zp = []
    for i in range(len(py_TOT_Zp)):
        for y in range(len(py_TOT_Zp[i])):
            P8_py_Zp.append(py_TOT_Zp[i][y])

    P8_pz_Zp = []
    for i in range(len(pz_TOT_Zp)):
        for y in range(len(pz_TOT_Zp[i])):
            P8_pz_Zp.append(pz_TOT_Zp[i][y])

    P8_E_Zp = []
    for i in range(len(E_TOT_Zp)):
        for y in range(len(E_TOT_Zp[i])):
            P8_E_Zp.append(E_TOT_Zp[i][y])

    P8_mass_Zp = []
    for i in range(len(mass_TOT_Zp)):
        for y in range(len(mass_TOT_Zp[i])):
            P8_mass_Zp.append(mass_TOT_Zp[i][y])

    P8_pdg_Zp = []
    for i in range(len(pdg_TOT_Zp)):
        for y in range(len(pdg_TOT_Zp[i])):
            P8_pdg_Zp.append(pdg_TOT_Zp[i][y])

    P8_px_Zp = np.array(P8_px_Zp)/c
    P8_py_Zp = np.array(P8_py_Zp)/c
    P8_pz_Zp = np.array(P8_pz_Zp)/c
    P8_E_Zp = np.array(P8_E_Zp)
    P8_mass_Zp = np.array(P8_mass_Zp)
    P8_pdg_Zp = np.array(P8_pdg_Zp)

    return P8_px_Zp, P8_py_Zp, P8_pz_Zp, P8_E_Zp, P8_mass_Zp, P8_pdg_Zp

def conversion_one_list_DH(px_TOT_DH, py_TOT_DH, pz_TOT_DH, E_TOT_DH, mass_TOT_DH, pdg_TOT_DH):

    P8_px_DH = []
    for i in range(len(px_TOT_DH)):
        for y in range(len(px_TOT_DH[i])):
            P8_px_DH.append(px_TOT_DH[i][y])

    P8_py_DH = []
    for i in range(len(py_TOT_DH)):
        for y in range(len(py_TOT_DH[i])):
            P8_py_DH.append(py_TOT_DH[i][y])

    P8_pz_DH = []
    for i in range(len(pz_TOT_DH)):
        for y in range(len(pz_TOT_DH[i])):
            P8_pz_DH.append(pz_TOT_DH[i][y])

    P8_E_DH = []
    for i in range(len(E_TOT_DH)):
        for y in range(len(E_TOT_DH[i])):
            P8_E_DH.append(E_TOT_DH[i][y])

    P8_mass_DH = []
    for i in range(len(mass_TOT_DH)):
        for y in range(len(mass_TOT_DH[i])):
            P8_mass_DH.append(mass_TOT_DH[i][y])

    P8_pdg_DH = []
    for i in range(len(pdg_TOT_DH)):
        for y in range(len(pdg_TOT_DH[i])):
            P8_pdg_DH.append(pdg_TOT_DH[i][y])

    P8_px_DH = np.array(P8_px_DH)/c
    P8_py_DH = np.array(P8_py_DH)/c
    P8_pz_DH = np.array(P8_pz_DH)/c
    P8_E_DH = np.array(P8_E_DH)
    P8_mass_DH = np.array(P8_mass_DH)
    P8_pdg_DH = np.array(P8_pdg_DH)

    return P8_px_DH, P8_py_DH, P8_pz_DH, P8_E_DH, P8_mass_DH, P8_pdg_DH

#########################################################################################
# Computation of the kinematics variable for DH and Zp (velocities, beta, gamma, pT the transverse momenta, eta the pseudo-rapidity).
#########################################################################################

def kinematics_DH(P8_px_DH, P8_py_DH, P8_pz_DH, P8_E_DH):

    #computing the velocities to obtain beta, gamma and boost DH
    vx_DH = (P8_px_DH*c**2)/P8_E_DH #compute the velocities in each direction in m/s
    vy_DH = (P8_py_DH*c**2)/P8_E_DH
    vz_DH = (P8_pz_DH*c**2)/P8_E_DH
    beta_DH = np.sqrt(vx_DH**2 + vy_DH**2 + vz_DH**2)/c # compute beta
    gamma_DH = 1/(np.sqrt(1-beta_DH**2)) # compute gamma

    pT_DH = np.sqrt(P8_px_DH**2 + P8_py_DH**2)*c # compute the transverse momenta
    eta_DH = np.arctanh(P8_pz_DH/(np.sqrt(P8_px_DH**2 + P8_py_DH**2 + P8_pz_DH**2)))# compute the pseudorapidity

    return beta_DH, gamma_DH, pT_DH, eta_DH

def kinematics_Zp(P8_px_Zp, P8_py_Zp, P8_pz_Zp, P8_E_Zp):

    vx_Zp = (P8_px_Zp*c**2)/P8_E_Zp
    vy_Zp = (P8_py_Zp*c**2)/P8_E_Zp
    vz_Zp = (P8_pz_Zp*c**2)/P8_E_Zp
    beta_Zp = np.sqrt(vx_Zp**2 + vy_Zp**2 + vz_Zp**2)/c
    gamma_Zp = 1/(np.sqrt(1-beta_Zp**2))

    pT_Zp = np.sqrt(P8_px_Zp**2 + P8_py_Zp**2)*c
    eta_Zp = np.arctanh(P8_pz_Zp/(np.sqrt(P8_px_Zp**2 + P8_py_Zp**2 + P8_pz_Zp**2)))

    return beta_Zp, gamma_Zp, pT_Zp, eta_Zp

#########################################################################################
# lifetime function.
#########################################################################################

def lifetime(avgtau = 4.3):
    import math
    avgtau = avgtau / c
    t = random.random()
    return -1.0 * avgtau * math.log(t)

#########################################################################################
# Decay lenght computation for LLP1.
#########################################################################################

def decaylenghtDH(P8_px_DH, P8_py_DH, P8_pz_DH, P8_E_DH, gamma_DH, tauN_DH):

    Lx_tot_DH = []
    Ly_tot_DH = []
    Lz_tot_DH = []
    Lxy_tot_DH = []

    for ctau in range(len(tauN_DH)):

        Lx_DH = []
        Ly_DH = []
        Lz_DH = []
        Lxy_DH = []

        for i in range(len(gamma_DH)):
            lt = lifetime(tauN_DH[ctau])
            Lx_DH.append((P8_px_DH[i]/P8_E_DH[i])*c**2 * lt * gamma_DH[i])
            Ly_DH.append((P8_py_DH[i]/P8_E_DH[i])*c**2 * lt * gamma_DH[i])
            Lz_DH.append((abs(P8_pz_DH[i])/P8_E_DH[i])*c**2 * lt  * gamma_DH[i])
            Lxy_DH.append(np.sqrt((Lx_DH[i])**2 + (Ly_DH[i])**2))

        Lx_tot_DH.append(Lx_DH)
        Ly_tot_DH.append(Ly_DH)
        Lz_tot_DH.append(Lz_DH)
        Lxy_tot_DH.append(Lxy_DH)

    return Lxy_tot_DH, Lz_tot_DH

#########################################################################################
# Decay lenght computation for LLP2.
#########################################################################################

def decaylenghtZp(P8_px_Zp, P8_py_Zp, P8_pz_Zp, P8_E_Zp, gamma_Zp, tauN_Zp):
    Lx_tot_Zp = []
    Ly_tot_Zp = []
    Lz_tot_Zp = []
    Lxy_tot_Zp = []

    for ctau in range(len(tauN_Zp)):
        Lx_Zp = []
        Ly_Zp = []
        Lz_Zp = []
        Lxy_Zp = []

        for i in range(len(gamma_Zp)):
                lt = lifetime(tauN_Zp[ctau])
                Lx_Zp.append((P8_px_Zp[i]/P8_E_Zp[i])*c**2 * lt * gamma_Zp[i])
                Ly_Zp.append((P8_py_Zp[i]/P8_E_Zp[i])*c**2 * lt * gamma_Zp[i])
                Lz_Zp.append((abs(P8_pz_Zp[i])/P8_E_Zp[i])*c**2 * lt* gamma_Zp[i])
                Lxy_Zp.append(np.sqrt((Lx_Zp[i])**2 + (Ly_Zp[i])**2))

        Lx_tot_Zp.append(Lx_Zp)
        Ly_tot_Zp.append(Ly_Zp)
        Lz_tot_Zp.append(Lz_Zp)
        Lxy_tot_Zp.append(Lxy_Zp)
    return Lxy_tot_Zp, Lz_tot_Zp

#########################################################################################
# Computation of the efficiency with the map from the data obtained with MG+Pythia8 for the high-ET samples (mZp >= 400GeV).
#########################################################################################

def eff_map_High(pT_DH, eta_DH, Lxy_tot_DH, Lz_tot_DH, P8_pdg_DH, pT_Zp, eta_Zp, Lxy_tot_Zp, Lz_tot_Zp, P8_pdg_Zp, tauN_DH, tauN_Zp, nevent, mass_Zp, mass_s):

    eff_highETX = []
    for index_Zp in range(len(tauN_Zp)):
        queryMapResult = []
        eff_lifetime=[]
        for index_DH in tqdm.tqdm(range(len(tauN_DH))):
            for iEvent in range(len(pT_DH)):
                queryMapResult.append(rmN.queryMapFromKinematics(pT_DH[iEvent],
                                                                eta_DH[iEvent],
                                                                Lxy_tot_DH[index_DH][iEvent],
                                                                Lz_tot_DH[index_DH][iEvent],
                                                                abs(P8_pdg_DH[iEvent]),
                                                                pT_Zp[iEvent],
                                                                eta_Zp[iEvent],
                                                                Lxy_tot_Zp[index_Zp][iEvent],
                                                                Lz_tot_Zp[index_Zp][iEvent],
                                                                abs(P8_pdg_Zp[iEvent]),
                                                                selection = "high-ET"))
            eff_lifetime.append(sum(queryMapResult))
        eff_highETX.append(eff_lifetime)
    eff_highETX = np.array(eff_highETX) #convertion into array
    eff_highETX = eff_highETX/nevent #efficiency/(nbr of event)

    Data_Eff_High = np.column_stack(eff_highETX)
    np.savetxt(f'./Plots_High/Efficiencies_Text_{mass_Zp}_{mass_s}.txt', Data_Eff_High)

    return eff_highETX

#########################################################################################
# Computation of the efficiency with the map from the data obtained with MG+Pythia8 for the low-ET samples (mZp <= 400GeV).
#########################################################################################

def eff_map_Low(pT_DH, eta_DH, Lxy_tot_DH, Lz_tot_DH, P8_pdg_DH, pT_Zp, eta_Zp, Lxy_tot_Zp, Lz_tot_Zp, P8_pdg_Zp, tauN_DH, tauN_Zp ,nevent, mass_Zp, mass_s):

    eff_lowETX = []
    for index_Zp in range(len(tauN_Zp)):
        queryMapResult = []
        eff_lifetime=[]
        for index_DH in tqdm.tqdm(range(len(tauN_DH))):
            for iEvent in range(len(pT_DH)):
                queryMapResult.append(rmN.queryMapFromKinematics(pT_DH[iEvent],
                                                                eta_DH[iEvent],
                                                                Lxy_tot_DH[index_DH][iEvent],
                                                                Lz_tot_DH[index_DH][iEvent],
                                                                abs(P8_pdg_DH[iEvent]),
                                                                pT_Zp[iEvent],
                                                                eta_Zp[iEvent],
                                                                Lxy_tot_Zp[index_Zp][iEvent],
                                                                Lz_tot_Zp[index_Zp][iEvent],
                                                                abs(P8_pdg_Zp[iEvent]),
                                                                selection = "low-ET"))
            eff_lifetime.append(sum(queryMapResult))
        eff_lowETX.append(eff_lifetime)
    eff_lowETX = np.array(eff_lowETX) #convertion into array
    eff_lowETX = eff_lowETX/nevent #efficiency/(nbr of event)

    Data_Eff_Low = np.column_stack(eff_lowETX)
    np.savetxt(f'./Plots_Low/Efficiencies_Text_{mass_Zp}_{mass_s}.txt', Data_Eff_Low)

    return eff_lowETX


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################

#########################################################################################
#Parsing the lhe file from the MG output to recover the data from the process.
#########################################################################################


def parsing_LHE_DH(MG_events):

    MG_px_DH = []
    MG_py_DH = []
    MG_pz_DH = []
    MG_pdg_DH = []
    MG_E_DH = []
    MG_MASS_DH = []
    for event in MG_events:
        for particle in event:
            if particle.pdg == 35:
                MG_px_DH.append(particle.px)
                MG_py_DH.append(particle.py)
                MG_pz_DH.append(particle.pz)
                MG_pdg_DH.append(particle.pdg)
                MG_E_DH.append(particle.E)
                MG_MASS_DH.append(particle.mass)

    MG_px_DH = np.array(MG_px_DH)/c # GeV/c = kg m s⁻¹
    MG_py_DH = np.array(MG_py_DH)/c #est bon sinon beta >1
    MG_pz_DH = np.array(MG_pz_DH)/c

    return MG_px_DH, MG_py_DH, MG_pz_DH, MG_pdg_DH, MG_E_DH, MG_MASS_DH


def parsing_LHE_Zp(MG_events):

    MG_px_Zp = []
    MG_py_Zp = []
    MG_pz_Zp = []
    MG_pdg_Zp = []
    MG_E_Zp = []
    MG_MASS_Zp = []
    for event in MG_events:
        for particle in event:
            if particle.pdg == 1023:
                MG_px_Zp.append(particle.px)
                MG_py_Zp.append(particle.py)
                MG_pz_Zp.append(particle.pz)
                MG_pdg_Zp.append(particle.pdg)
                MG_E_Zp.append(particle.E)
                MG_MASS_Zp.append(particle.mass)

    MG_px_Zp = np.array(MG_px_Zp)/c # GeV/c = kg m s⁻¹
    MG_py_Zp = np.array(MG_py_Zp)/c #est bon sinon beta >1
    MG_pz_Zp = np.array(MG_pz_Zp)/c

    return MG_px_Zp, MG_py_Zp, MG_pz_Zp, MG_pdg_Zp, MG_E_Zp, MG_MASS_Zp

#########################################################################################
# Computation of the kinematics variable for LLP1 (velocities, beta, gamma, pT the transverse momenta, eta the pseudo-rapidity).
#########################################################################################

def kinematics_MG_DH(MG_px_DH,MG_py_DH,MG_pz_DH,MG_E_DH):

    MG_vx_DH = (MG_px_DH*c**2)/MG_E_DH # #compute the velocities in each direction in m/s
    MG_vy_DH = (MG_py_DH*c**2)/MG_E_DH
    MG_vz_DH = (MG_pz_DH*c**2)/MG_E_DH
    MG_beta_DH = np.sqrt(MG_vx_DH**2 + MG_vy_DH**2 + MG_vz_DH**2)/c # compute beta
    MG_gamma_DH = 1/(np.sqrt(1-MG_beta_DH**2)) # compute gamma
    MG_pT_DH = np.sqrt(MG_px_DH**2 + MG_py_DH**2)*c # compute the transverse momenta
    MG_eta_DH = np.arctanh(MG_pz_DH/(np.sqrt(MG_px_DH**2 + MG_py_DH**2 + MG_pz_DH**2))) # compute the pseudorapidity

    return MG_pT_DH,MG_eta_DH, MG_gamma_DH

def kinematics_MG_Zp(MG_px_Zp,MG_py_Zp,MG_pz_Zp,MG_E_Zp):

    MG_vx_Zp = (MG_px_Zp*c**2)/MG_E_Zp
    MG_vy_Zp = (MG_py_Zp*c**2)/MG_E_Zp
    MG_vz_Zp = (MG_pz_Zp*c**2)/MG_E_Zp
    MG_beta_Zp = np.sqrt(MG_vx_Zp**2 + MG_vy_Zp**2 + MG_vz_Zp**2)/c
    MG_gamma_Zp = 1/(np.sqrt(1-MG_beta_Zp**2))
    MG_pT_Zp = np.sqrt(MG_px_Zp**2 + MG_py_Zp**2)*c
    MG_eta_Zp = np.arctanh(MG_pz_Zp/(np.sqrt(MG_px_Zp**2 + MG_py_Zp**2 + MG_pz_Zp**2)))

    return MG_pT_Zp,MG_eta_Zp, MG_gamma_Zp

#########################################################################################
# Decay lenght computation for LLP1.
#########################################################################################

def decaylenght_MG_DH(MG_px_DH, MG_py_DH, MG_pz_DH, MG_E_DH, MG_gamma_DH, tauN_DH):

    MG_Lx_tot_DH = []
    MG_Ly_tot_DH = []
    MG_Lz_tot_DH = []
    MG_Lxy_tot_DH = []

    for ctau in range(len(tauN_DH)):

        MG_Lx_DH = []
        MG_Ly_DH = []
        MG_Lz_DH = []
        MG_Lxy_DH = []

        for i in range(len(MG_gamma_DH)):
            MG_lt = lifetime(tauN_DH[ctau])
            MG_Lx_DH.append((MG_px_DH[i]/MG_E_DH[i])*c**2 * MG_lt * MG_gamma_DH[i])
            MG_Ly_DH.append((MG_py_DH[i]/MG_E_DH[i])*c**2 * MG_lt * MG_gamma_DH[i])
            MG_Lz_DH.append((abs(MG_pz_DH[i])/MG_E_DH[i])*c**2 * MG_lt  * MG_gamma_DH[i] )
            MG_Lxy_DH.append(np.sqrt((MG_Lx_DH[i])**2 + (MG_Ly_DH[i])**2))

        MG_Lx_tot_DH.append(MG_Lx_DH)
        MG_Ly_tot_DH.append(MG_Ly_DH)
        MG_Lz_tot_DH.append(MG_Lz_DH)
        MG_Lxy_tot_DH.append(MG_Lxy_DH)


    return MG_Lxy_tot_DH, MG_Lz_tot_DH

#########################################################################################
# Decay lenght computation for LLP2.
#########################################################################################

def decaylenght_MG_Zp(MG_px_Zp, MG_py_Zp, MG_pz_Zp, MG_E_Zp, MG_gamma_Zp, tauN_Zp):

    MG_Lx_tot_Zp = []
    MG_Ly_tot_Zp = []
    MG_Lz_tot_Zp = []
    MG_Lxy_tot_Zp = []

    for ctau in range(len(tauN_Zp)):

        MG_Lx_Zp = []
        MG_Ly_Zp = []
        MG_Lz_Zp = []
        MG_Lxy_Zp = []

        for i in range(len(MG_gamma_Zp)):
            MG_lt = lifetime(tauN_Zp[ctau])
            MG_Lx_Zp.append((MG_px_Zp[i]/MG_E_Zp[i])*c**2 * MG_lt * MG_gamma_Zp[i])
            MG_Ly_Zp.append((MG_py_Zp[i]/MG_E_Zp[i])*c**2 * MG_lt * MG_gamma_Zp[i])
            MG_Lz_Zp.append((abs(MG_pz_Zp[i])/MG_E_Zp[i])*c**2 * MG_lt  * MG_gamma_Zp[i] )
            MG_Lxy_Zp.append(np.sqrt((MG_Lx_Zp[i])**2 + (MG_Ly_Zp[i])**2))

        MG_Lx_tot_Zp.append(MG_Lx_Zp)
        MG_Ly_tot_Zp.append(MG_Ly_Zp)
        MG_Lz_tot_Zp.append(MG_Lz_Zp)
        MG_Lxy_tot_Zp.append(MG_Lxy_Zp)

    return MG_Lxy_tot_Zp, MG_Lz_tot_Zp

#########################################################################################
# Computation of the efficiency with the map from the data obtained with MG for the high-ET samples.
#########################################################################################

def eff_map_MG_high(MG_pT_DH, MG_eta_DH,MG_Lxy_tot_DH, MG_Lz_tot_DH, MG_pdg_DH, MG_pT_Zp, MG_eta_Zp, MG_Lxy_tot_Zp, MG_Lz_tot_Zp, MG_pdg_Zp, tauN_DH, tauN_Zp, nevent, mass_Zp, mass_s):

    MG_eff_highETX = []
    for index_Zp in range(len(tauN_Zp)):
        MG_queryMapResult = []
        MG_eff_lifetime=[]
        for index_DH in tqdm.tqdm(range(len(tauN_DH))):
            for iEvent in range(len(MG_pT_DH)):
                MG_queryMapResult.append(rmN.queryMapFromKinematics(MG_pT_DH[iEvent],
                                                                MG_eta_DH[iEvent],
                                                                MG_Lxy_tot_DH[index_DH][iEvent],
                                                                MG_Lz_tot_DH[index_DH][iEvent],
                                                                abs(MG_pdg_DH[iEvent]),
                                                                MG_pT_Zp[iEvent],
                                                                MG_eta_Zp[iEvent],
                                                                MG_Lxy_tot_Zp[index_Zp][iEvent],
                                                                MG_Lz_tot_Zp[index_Zp][iEvent],
                                                                abs(MG_pdg_Zp[iEvent]),
                                                                selection = "high-ET"))


            MG_eff_lifetime.append(sum(MG_queryMapResult))
        MG_eff_highETX.append(MG_eff_lifetime)
    MG_eff_highETX = np.array(MG_eff_highETX) #convertion into array
    MG_eff_highETX = MG_eff_highETX//nevent #eff/nbrevent

    MG_Data_Eff_High = np.column_stack(MG_eff_highETX)
    np.savetxt(f'./Plots_High/Efficiencies_Text_{mass_Zp}_{mass_s}.txt', MG_Data_Eff_High)

    return MG_eff_highETX

#########################################################################################
# Computation of the efficiency with the map from the data obtained with MG for the low-ET samples.
#########################################################################################

def eff_map_MG_low(MG_pT_DH, MG_eta_DH,MG_Lxy_tot_DH, MG_Lz_tot_DH, MG_pdg_DH, MG_pT_Zp, MG_eta_Zp, MG_Lxy_tot_Zp, MG_Lz_tot_Zp, MG_pdg_Zp, tauN_DH, tauN_Zp, nevent, mass_Zp, mass_s):

    MG_eff_lowETX = []
    for index_Zp in range(len(tauN_Zp)):
        MG_queryMapResult = []
        MG_eff_lifetime=[]
        for index_DH in tqdm.tqdm(range(len(tauN_DH))):
            for iEvent in range(len(MG_pT_DH)):
                MG_queryMapResult.append(rmN.queryMapFromKinematics(MG_pT_DH[iEvent],
                                                                MG_eta_DH[iEvent],
                                                                MG_Lxy_tot_DH[index_DH][iEvent],
                                                                MG_Lz_tot_DH[index_DH][iEvent],
                                                                abs(MG_pdg_DH[iEvent]),
                                                                MG_pT_Zp[iEvent],
                                                                MG_eta_Zp[iEvent],
                                                                MG_Lxy_tot_Zp[index_Zp][iEvent],
                                                                MG_Lz_tot_Zp[index_Zp][iEvent],
                                                                abs(MG_pdg_Zp[iEvent]),
                                                                selection = "low-ET"))


            MG_eff_lifetime.append(sum(MG_queryMapResult))
        MG_eff_lowETX.append(MG_eff_lifetime)
    MG_eff_lowETX = np.array(MG_eff_lowETX) #convertion into array
    MG_eff_lowETX = MG_eff_lowETX/nevent #eff/nbrevent

    MG_Data_Eff_Low = np.column_stack(MG_eff_lowETX)
    np.savetxt(f'./Plots_Low/Efficiencies_Text_{mass_Zp}_{mass_s}.txt', MG_Data_Eff_Low)

    return MG_eff_lowETX


#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################

#########################################################################################
# Plots to compare the results of efficiency obtained with MG, MG+Pythia8 (High-ET).
#########################################################################################

def plt_eff_high(MG_eff_highETX, eff_highETX,tauN,  mass_Zp , mass_s):

    ################## PLOT EFFICIENCY ##################
    fig, ax = plt.subplots()

    ################## Plot efficiency from MG ##################
    plt.plot(tauN,MG_eff_highETX, 'k--', linewidth=2, label = 'MG')

    ################## Plot efficiency from MG+Pythia8 ##################
    plt.plot(tauN,eff_highETX, 'r', linewidth=2, label = 'MG + Pythia')

    ################## Uncertainties from Map ##################
    plt.fill_between(tauN, np.array(eff_highETX) + 0.25* np.array(eff_highETX), np.array(eff_highETX) - 0.25 * np.array(eff_highETX), label='MG+Pythia8, with error bands ', alpha=.7)

    ################## Limits of validity ##################
    ax.hlines(y=(0.25*(max(eff_highETX))), xmin=0, xmax=1e2, linewidth=2, color='g', label = 'Limits of validity' )

    # place a text box in upper left in axes coords
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax.text(0.05, 0.95, f" $ m_Zp $ = {mass_Zp} GeV, $m_S$ = {mass_s} GeV" , transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    x = np.linspace(0,100)
    ax.fill_between(x, 0.25*(max(eff_highETX)), color='black', alpha=.2, hatch="/", edgecolor="black", linewidth=1.0) # adding hatch
    plt.ylim(0) # start at 0

    plt.xscale('log')
    plt.xlabel(r'c$\tau$ [m]', fontsize=20)
    plt.ylabel('Efficiency', fontsize=20 )
    plt.legend(fontsize = 11, loc=1) # set the legend in the upper right corner
    plt.savefig(f"./Plots_High/Efficiency_comparison_mZp{mass_Zp}_mS{mass_s}.png")
    plt.close()


#########################################################################################
# Plots to compared the reasults of efficiency obtained with MG, MG+Pythia8 (Low-ET).
#########################################################################################

def plt_eff_low(MG_eff_lowETX, eff_lowETX,tauN,  mass_Zp , mass_s):

    ################## PLOT EFFICIENCY ##################
    fig, ax = plt.subplots()

    ################## Plot efficiency from MG ##################
    plt.plot(tauN,MG_eff_lowETX, 'k--',linewidth=2, label = 'MG')

    ################## Plot efficiency from MG+Pythia8 ##################
    plt.plot(tauN,eff_lowETX, 'r', linewidth=2 ,label = 'MG + Pythia')

    ################## Uncertainties from Map ##################
    plt.fill_between(tauN, np.array(eff_lowETX) + 0.25* np.array(eff_lowETX), np.array(eff_lowETX) - 0.25*np.array(eff_lowETX), label='MG+Pythia8, with error bands ',alpha=.7)

    ################## Limits of validity ##################
    ax.hlines(y=(0.33*(max(eff_lowETX))), xmin=0, xmax=1e2, linewidth=2, color='g', label = 'Limits of validity' )

    # place a text box in upper left in axes coords
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax.text(0.05, 0.95, f" $ m_Zp $ = {mass_Zp} GeV, $m_S$ = {mass_s} GeV" , transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    x = np.linspace(0,100)
    ax.fill_between(x, 0.33*(max(eff_lowETX)), color='black', alpha=.2, hatch="/", edgecolor="black", linewidth=1.0) # adding hatch
    plt.ylim(0) # start at 0

    plt.xscale('log')
    plt.xlabel(r'c$\tau$ [m]', fontsize=20)
    plt.ylabel('Efficiency', fontsize=20 )
    plt.legend( fontsize = 10, loc=1) # set the legend in the upper right corner
    plt.savefig(f"./Plots_Low/Efficiency_comparison_mZp{mass_Zp}_mS{mass_s}.png")
    plt.close()

#########################################################################################
# Plot limits obtained with the map, to compare with those obtain by ATLAS (High-ET).
#########################################################################################

def plt_cross_High(eff_highETX, tauN, mass_Zp, mass_s, factor):

    fig, ax = plt.subplots()

    Nsobs = 0.5630 * 26 * factor # nbr of observed events = 26 ( factor )

    Crr_Sec_obs = (Nsobs)/((np.array(eff_highETX)) * 139e3 ) # Luminosity = 139e3 fb**(-1)

    plt.plot(tauN, Crr_Sec_obs, 'r', label ='Map results', linewidth = 2)

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'c$\tau$ [m]')
    plt.ylabel(r'95% CL limit on $\sigma \times B$ [pb]')

    # place a text box in upper left in axes coords
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax.text(0.05, 0.95, f" $ m_Zp $ = {mass_Zp} GeV, $m_S$ = {mass_s} GeV" , transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    plt.legend( fontsize = 10, loc=3)
    plt.savefig(f"./Plots_High/Cross_section_mZp{mass_Zp}_mS{mass_s}.png") #create a new fodlder ' Plots ' and save the fig in it
    plt.close()

#########################################################################################
# Plot limits obtained with the map, to those obtain by ATLAS (Low-ET).
#########################################################################################

def plt_cross_Low(eff_lowETX , tauN, mass_Zp, mass_s, factor):

    fig, ax = plt.subplots()

    Nsobs = 0.6592 * 26 * factor # nbr of observed events = 26

    Crr_Sec_obs = (Nsobs)/((np.array(eff_lowETX)) * 139e3 )

    plt.plot(tauN, Crr_Sec_obs, 'r', label ='Map results', linewidth = 2)

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'c$\tau$ [m]')
    plt.ylabel(r'95% CL limit on $\sigma \times B$ [pb]')

    # place a text box in upper left in axes coords
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax.text(0.05, 0.95, f" $ m_Zp $ = {mass_Zp} GeV, $m_S$ = {mass_s} GeV" , transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    plt.legend( fontsize = 10, loc=3)
    plt.savefig(f"./Plots_Low/Cross_section_mZp{mass_Zp}_mS{mass_s}.png") #create a new fodlder ' Plots ' and save the fig in it
    plt.close()

#########################################################################################
# Plot limits obtained with the map, to those obtain by ATLAS (Low-ET).
#########################################################################################

def plt_contour_high(tauN_DH, tauN_Zp, eff_highETX):

    fig, ax = plt.subplots()


    plt.contour(tauN_DH, tauN_Zp, eff_highETX)

    # place a text box in upper left in axes coords
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax.text(0.05, 0.95, f" $ m_Zp $ = {mass_Zp} GeV, $m_S$ = {mass_s} GeV" , transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    plt.colorbar()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('cTau_DH [m]',fontsize = 20)
    plt.ylabel('cTau_Zp [m]',fontsize = 20)
    plt.savefig(f"./Plots_High/Contour_mZp{mass_Zp}_mS{mass_s}.png") #create a new fodlder ' Plots ' and save the fig in it


def plt_contour_low(tauN_DH, tauN_Zp, eff_lowETX):

    fig, ax = plt.subplots()


    # place a text box in upper left in axes coords
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax.text(0.05, 0.95, f" $ m_Zp $ = {mass_Zp} GeV, $m_S$ = {mass_s} GeV" , transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    plt.contour(tauN_DH, tauN_Zp, eff_lowETX)
    plt.colorbar()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('cTau_DH [m]',fontsize = 20)
    plt.ylabel('cTau_Zp [m]',fontsize = 20)
    plt.savefig(f"./Plots_Low/Contour_mZp{mass_Zp}_mS{mass_s}.png") #create a new fodlder ' Plots ' and save the fig in it

