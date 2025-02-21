import sys
import pandas as pd
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import scipy
import time
import tqdm
import mplhep as hep
import lhe_parser as lhe
import hepmc_parser as hepmc
import uproot
import random
import copy
import os
import pickle
random.seed(123)
hep.style.use("ATLAS")

# set cst
c = 3e8# Light velocity in m/s

#All plots are made with 10 000 events, if you want to try with other numbers of events, you will have to change the number of events in lines 285-314-551-580 for the calculation of the efficiency.

#########################################################################################
#Parsing the hepmc file from the hadronization of the MG outputs to recover the data from the process.
#########################################################################################



def parsing_hepmc_generic(events, mediator=25, vectorBoson=-999, LLP=35, LLP2=None, verbose=0):

    variables = ['px', 'py', 'pz', 'E', 'mass', 'decay'] 
    result={
    'Phi':pd.DataFrame(columns=variables),
    'V':pd.DataFrame(columns=variables),
    'llp1':pd.DataFrame(columns=variables),
    'llp2':pd.DataFrame(columns=variables),
    }
    print (f"using mediator={mediator}, LLP={LLP}, LLP2={LLP2}, V={vectorBoson}")
    for ie , event in enumerate(events):
        #if ie > 20 : break
        if verbose: print("============NEW EVENT================")
        llpProdCounter=1
        llpDecayCounter=1
        vectorBosonProd, vectorBosonDecay, mediatorProd, mediatorDecay, LLP1Prod, LLP1Decay, LLP2Prod, LLP2Decay =  None, None, None, None, None, None, None, None
        for id, vertex in event.vertex.items():
            particlesOut = sorted([p for p in vertex.outcoming], key = lambda k : abs(k.pdg))
            particlesIn = sorted([p for p in vertex.incoming], key = lambda k : abs(k.pdg))
            pdgOut = [abs(p.pdg) for p in particlesOut] 
            pdgIn = [abs(p.pdg) for p in particlesIn] 
            pxIn = [p.px for p in particlesIn] 
            pxOut = [p.px for p in particlesOut] 
            # skip the cases where we have a particle both in the ingoing and outgoing particles
            # this is going to be simply a radiation of a gluon or photon...
            if LLP in pdgOut and LLP in pdgIn: continue
            if LLP2 in pdgOut and LLP2 in pdgIn: continue
            if vectorBoson in pdgOut and vectorBoson in pdgIn: continue
            if mediator in pdgOut and mediator in pdgIn: continue
            # also skip vertices where none of the particles we are interested in are involved
            LLP_inNeither = LLP not in pdgOut and LLP not in pdgIn
            LLP2_inNeither = LLP2 not in pdgOut and LLP2 not in pdgIn
            V_inNeither = vectorBoson not in pdgOut and vectorBoson not in pdgIn
            med_inNeither = mediator not in pdgOut and mediator not in pdgIn
            if LLP2 is None:
              if (LLP_inNeither and V_inNeither and med_inNeither): continue
            else: 
              if (LLP2_inNeither and LLP_inNeither and V_inNeither and med_inNeither): continue
            if verbose: print (id, pdgIn, pxIn, "--->", pdgOut, pxOut)


            #if  vectorBoson in pdgOut: # production of vector boson
            #   if verbose: print("--> it's the vector boson production vertex!")
            #   vectorBosonProd  = particlesOut[pdgOut.index(vectorBoson)]
            if  vectorBoson in pdgIn: # decay of vector boson
               if verbose: print("--> it's the vector boson decay vertex!")
               vectorBosonDecay  = particlesOut[0] # always take the first of the two decay particles, should always get the lepton if l,v
               vectorBosonProd =  particlesIn[0]
            
            #if  mediator in pdgOut: # production of mediator
            #   if verbose: print("--> it's the mediator production vertex!")
            #   mediatorProd  = particlesOut[pdgOut.index(mediator)]
            if  mediator in pdgIn: # decay of mediator
               if verbose: print("--> it's the mediator decay vertex!")
               mediatorDecay  = particlesOut[0]
               mediatorProd =  particlesIn[0]
            
            #if  LLP in pdgOut: # production of LLP
            #   if verbose: print("--> it's the LLP production vertex!")
            #   if pdgOut == [LLP, LLP] : # pair production of LLP
            #     LLP1Prod  = particlesOut[0]
            #     LLP2Prod  = particlesOut[1]
            #   else:
            #     if llpProdCounter==1: LLP1Prod  = particlesOut[pdgOut.index(LLP)]
            #     if llpProdCounter==2: LLP2Prod  = particlesOut[pdgOut.index(LLP)]
            #     llpProdCounter+=1
            
            if  LLP in pdgIn: # decay of LLP
               if verbose: print(f"--> it's the {llpDecayCounter}^th LLP decay vertex!",  particlesOut[0].pdg)
               if llpDecayCounter==1: LLP1Decay  = particlesOut[0]
               if llpDecayCounter==2: LLP2Decay  = particlesOut[0]
               if llpDecayCounter==1: LLP1Prod  = particlesIn[0]
               if llpDecayCounter==2: LLP2Prod  = particlesIn[0]
               llpDecayCounter+=1
            if  LLP2 in pdgIn: # decay of LLP
               if verbose: print(f"--> it's the {llpDecayCounter}^th LLP decay vertex!",  particlesOut[0].pdg)
               LLP2Decay  = particlesOut[0]
               LLP2Prod  = particlesIn[0]

        #if LLP2 is None:
        #  if len(result['llp2'])>0 and np.random.random() > 0.5:
        #    tmp = LLP2Decay
        #    LLP2Decay = LLP1Decay
        #    LLP1Decay = tmp
        #    tmp = LLP2Prod
        #    LLP2Prod = LLP1Prod
        #    LLP1Prod = tmp
        #if LLP2 is None:
        #  if len(result['llp2'])>0 and np.random.random() > 0.5:
        #print("LC DEBUG I'm in the EVNET", "LLP1 pT=", LLP1Prod.mass, LLP2Prod.mass)
        if np.random.random() > 0.5 and LLP2 is not None:
            tmp = LLP2Decay
            LLP2Decay = LLP1Decay
            LLP1Decay = tmp
            tmp = LLP2Prod
            LLP2Prod = LLP1Prod
            LLP1Prod = tmp
            #print("LC DEBUG I'm SWAPPING LIKE A CRAZY LOCO", LLP1Prod.mass, LLP2Prod.mass)

        if vectorBosonProd is not None and vectorBosonDecay is not None:
           result['V'].loc[len(result['V'])] = [vectorBosonProd.px, vectorBosonProd.py, vectorBosonProd.pz, vectorBosonProd.E, vectorBosonProd.mass, vectorBosonDecay.pdg]
        if mediatorProd is not None and mediatorDecay is not None :
           result['Phi'].loc[len(result['Phi'])] = [mediatorProd.px, mediatorProd.py, mediatorProd.pz, mediatorProd.E, mediatorProd.mass, mediatorDecay.pdg]
        
        if LLP1Prod is not None and LLP1Decay is not None:
           result['llp1'].loc[len(result['llp1'])] = [LLP1Prod.px, LLP1Prod.py, LLP1Prod.pz, LLP1Prod.E, LLP1Prod.mass, LLP1Decay.pdg]
        
        if LLP2Prod is not None and LLP2Decay is not None:
           result['llp2'].loc[len(result['llp2'])] = [LLP2Prod.px, LLP2Prod.py, LLP2Prod.pz, LLP2Prod.E, LLP2Prod.mass, LLP2Decay.pdg]
        
        
        
    return result

def parsing_hepmc(events):

    px_TOT = []
    py_TOT = []
    pz_TOT = []
    E_TOT = []
    mass_TOT = []
    pdg_TOT = []

    for ie , event in enumerate(events):
        count=0
        for id, vertex in event.vertex.items():
            if [p.pdg for p in vertex.incoming] == [25] and [p.pdg for p in vertex.outcoming] == [35, 35]: # PDGID 25 = Higgs, PDGID 35 Dark Higgs
                px_TOT.append(list(p.px for p in vertex.outcoming)) # recover the x momenta in GeV
                py_TOT.append(list(p.py for p in vertex.outcoming)) # recover the y momenta in GeV
                pz_TOT.append(list(p.pz for p in vertex.outcoming)) # recover the z momenta in GeV
                E_TOT.append(list(p.E for p in vertex.outcoming)) # recover the Energy in GeV
                mass_TOT.append(list(p.mass for p in vertex.outcoming)) # recover the mass in GeV

            if [p.pdg for p in vertex.incoming] == [35]: # PDGID 35 Dark Higgs
                pdg_TOT.append((list(p.pdg for p in vertex.outcoming))) # recover the PDG ID of the particle produced

                count = count+1
                if count==2: ##
                    break
                    pass

    return px_TOT, py_TOT, pz_TOT, E_TOT, mass_TOT, pdg_TOT

def parsing_hepmc_ALP(events):

    px_TOT = []
    py_TOT = []
    pz_TOT = []
    E_TOT = []
    mass_TOT = []
    pdg_TOT = []

    for ie , event in enumerate(events):
        count=0
        for id, vertex in event.vertex.items():
            particles = sorted([p for p in vertex.outcoming], key = lambda k : k.pdg)
            particlesIn = sorted([p for p in vertex.incoming], key = lambda k : k.pdg)
            if 9000005 in [p.pdg for p in particles]  and len(particles)==2 : # PDGID 25 = Higgs, PDGID 35 Dark Higgs
                px_TOT.append(list(p.px for p in particles)) # recover the x momenta in GeV
                py_TOT.append(list(p.py for p in particles)) # recover the y momenta in GeV
                pz_TOT.append(list(p.pz for p in particles)) # recover the z momenta in GeV
                E_TOT.append(list(p.E for p in particles)) # recover the Energy in GeV
                mass_TOT.append(list(p.mass for p in particles)) # recover the mass in GeV

            if [p.pdg for p in vertex.incoming] == [9000005] and len(particles)==2: # PDGID 35 Dark Higgs
                pdg_TOT.append((list(p.pdg for p in vertex.outcoming))) # recover the PDG ID of the particle produced

                count = count+1
                if count==2: ##
                    break
                    pass

    return px_TOT, py_TOT, pz_TOT, E_TOT, mass_TOT, pdg_TOT


#########################################################################################
# Computation of the kinematics variable for LLP1 (velocities, beta, gamma, pT the transverse momenta, eta the pseudo-rapidity).
#########################################################################################

def kinematics(df):
    px = df.px
    py = df.py
    pz = df.pz
    E = df.E
    vx = (px*c)/E #compute the velocities in each direction
    vy = (py*c)/E
    vz = (pz*c)/E
    beta = np.sqrt(vx**2 + vy**2 + vz**2)/c # compute beta
    gamma = 1/(np.sqrt(1-beta**2)) # compute gamma

    pT = np.sqrt(px**2 + py**2) # compute the transverse momenta
    eta = np.arctanh(pz/(np.sqrt(px**2 + py**2 + pz**2))) # compute the pseudorapidity
    df['beta'] = beta
    df['gamma'] = gamma
    df['pT'] = pT
    df['eta'] = eta
    return df


#########################################################################################
# Computation of the kinematics variable for LLP2 (velocities, beta, gamma, pT the transverse momenta, eta the pseudo-rapidity).
#########################################################################################

def kinematics_DH2(px_DH2, py_DH2, pz_DH2, E_DH2):

    vx_DH2 = (px_DH2*c**2)/E_DH2 #compute the velocities in each direction
    vy_DH2 = (py_DH2*c**2)/E_DH2
    vz_DH2 = (pz_DH2*c**2)/E_DH2
    beta_DH2 = np.sqrt(vx_DH2**2 + vy_DH2**2 + vz_DH2**2)/c # compute beta
    gamma_DH2 = 1/(np.sqrt(1-beta_DH2**2)) # compute gamma

    pT_DH2 = np.sqrt(px_DH2**2 + py_DH2**2)*c # compute the transverse momenta
    eta_DH2 = np.arctanh(pz_DH2/(np.sqrt(px_DH2**2 + py_DH2**2 + pz_DH2**2))) # compute the pseudorapidity

    return beta_DH2, gamma_DH2, pT_DH2, eta_DH2

#########################################################################################
# lifetime function.
#########################################################################################

def lifetime(avgtau = 4.3):
    import math
    avgtau = avgtau #/ c
    t = random.random()
    return -1.0 * avgtau * math.log(t)

#########################################################################################
# Decay lenght computation for LLP1.
#########################################################################################

def decayLength(df, tauN):
   
    px = df.px
    py = df.py
    pz = df.pz
    E = df.E
    gamma = df.gamma

    ctLx = []
    ctLy = []
    ctLz = []
    ctLxy = []

    for ctau in range(len(tauN)):

        Lx = []
        Ly = []
        Lz = []
        Lxy = []

        for i in range(len(gamma)):
            lt = lifetime(tauN[ctau]) # set mean lifetime
            Lx.append((px[i]/E[i]) * lt * gamma[i]) # compute the decay lenght in x,y,z
            Ly.append((py[i]/E[i]) * lt * gamma[i])
            Lz.append((abs(pz[i])/E[i]) * lt  * gamma[i] )
            Lxy.append(np.sqrt((Lx[i])**2 + (Ly[i])**2)) # compte the transverse decay lenght

        ctLx.append(Lx)
        ctLy.append(Ly)
        ctLz.append(Lz)
        ctLxy.append(Lxy)
    return [np.array(ctLxy), np.array(ctLz)]

#########################################################################################
# Decay lenght computation for LLP2.
#########################################################################################

def decaylenghtDH2(px_DH2, py_DH2, pz_DH2, E_DH2, gamma_DH2, tauN):

    Lx_tot_DH2 = []
    Ly_tot_DH2 = []
    Lz_tot_DH2 = []
    Lxy_tot_DH2 = []

    for ctau in range(len(tauN)):
        Lx_DH2 = []
        Ly_DH2 = []
        Lz_DH2 = []
        Lxy_DH2 = []

        for i in range(len(gamma_DH2)):
                lt = lifetime(tauN[ctau]) # set mean lifetime
                Lx_DH2.append((px_DH2[i]/E_DH2[i])*c**2 * lt * gamma_DH2[i]) # compute the decay lenght in x,y,z
                Ly_DH2.append((py_DH2[i]/E_DH2[i])*c**2 * lt * gamma_DH2[i])
                Lz_DH2.append((abs(pz_DH2[i])/E_DH2[i])*c**2 * lt* gamma_DH2[i])
                Lxy_DH2.append(np.sqrt((Lx_DH2[i])**2 + (Ly_DH2[i])**2)) # compte the transverse decay lenght

        Lx_tot_DH2.append(Lx_DH2)
        Ly_tot_DH2.append(Ly_DH2)
        Lz_tot_DH2.append(Lz_DH2)
        Lxy_tot_DH2.append(Lxy_DH2)
    return Lxy_tot_DH2, Lz_tot_DH2

#########################################################################################
# Computation of the efficiency with the map from the data obtained with MG+Pythia8 for the high-ET samples (mH >= 400GeV).
#########################################################################################

def eff_map_High(pT_DH1, eta_DH1, Lxy_tot_DH1, Lz_tot_DH1, pdg_tot_DH1, pT_DH2, eta_DH2, Lxy_tot_DH2, Lz_tot_DH2, pdg_tot_DH2, tauN, nevent, mass_phi, mass_s):

    eff_highETX = []

    for index in tqdm.tqdm(range(len(tauN))):
        queryMapResult = []
        for iEvent in range(len(pT_DH1)):
            queryMapResult.append(rmN.queryMapFromKinematics(pT_DH1[iEvent],
                                                            eta_DH1[iEvent],
                                                            Lxy_tot_DH1[index][iEvent],
                                                            Lz_tot_DH1[index][iEvent],
                                                            abs(pdg_tot_DH1[iEvent]),
                                                            pT_DH2[iEvent],
                                                            eta_DH2[iEvent],
                                                            Lxy_tot_DH2[index][iEvent],
                                                            Lz_tot_DH2[index][iEvent],
                                                            abs(pdg_tot_DH2[iEvent]),
                                                            selection = "high-ET"))
        eff_highETX.append(sum(queryMapResult))
    queryMapResult = np.array(queryMapResult) #convertion into array
    eff_highETX = np.array(eff_highETX) #convertion into array
    eff_highETX = eff_highETX/nevent #efficiency/(nbr of event)

    Data_Eff_High = np.column_stack(eff_highETX)
    np.savetxt(f'./Plots_High/Efficiencies_Text_{mass_phi}_{mass_s}.txt', Data_Eff_High)

    return eff_highETX

#########################################################################################
# Computation of the efficiency with the map from the data obtained with MG+Pythia8 for the low-ET samples (mH <= 400GeV).
#########################################################################################

def eff_map_Low(pT_DH1, eta_DH1, Lxy_tot_DH1, Lz_tot_DH1, pdg_tot_DH1, pT_DH2, eta_DH2, Lxy_tot_DH2, Lz_tot_DH2, pdg_tot_DH2, tauN,nevent, mass_phi, mass_s):

    eff_lowETX = []

    for index in tqdm.tqdm(range(len(tauN))):
        queryMapResult = []
        for iEvent in range(len(pT_DH1)):
            queryMapResult.append(rmN.queryMapFromKinematics(pT_DH1[iEvent],
                                                            eta_DH1[iEvent],
                                                            Lxy_tot_DH1[index][iEvent],
                                                            Lz_tot_DH1[index][iEvent],
                                                            abs(pdg_tot_DH1[iEvent]),
                                                            pT_DH2[iEvent],
                                                            eta_DH2[iEvent],
                                                            Lxy_tot_DH2[index][iEvent],
                                                            Lz_tot_DH2[index][iEvent],
                                                            abs(pdg_tot_DH2[iEvent]),
                                                            selection = "low-ET"))
        eff_lowETX.append(sum(queryMapResult))
    queryMapResult = np.array(queryMapResult) #convertion into array
    eff_lowETX = np.array(eff_lowETX) #convertion into array
    eff_lowETX = eff_lowETX/nevent #efficiency/(nbr of event)

    Data_Eff_Low = np.column_stack(eff_lowETX)
    np.savetxt(f'./Plots_Low/Efficiencies_Text_{mass_phi}_{mass_s}.txt', Data_Eff_Low)

    return eff_lowETX

    
def eff_bdt_tauN(values, tauN,  sel, do2D=False):
   eff = []
   for index in tqdm.tqdm(range(len(tauN))):
     for index2 in tqdm.tqdm(range(len(tauN))):
        if not do2D: 
          if index2 > 0 : continue
          index2=index 

        thesevalues =copy.deepcopy(values)
        thesevalues['llp1_Lxy']=values['llp1_Lxy'][index]
        thesevalues['llp2_Lxy']=values['llp2_Lxy'][index2]
        thesevalues['llp1_Lz']=values['llp1_Lz'][index]
        thesevalues['llp2_Lz']=values['llp2_Lz'][index2]
        keys = thesevalues.keys()
        factor=1
        for k in keys:
          thesevalues[k] = np.array(thesevalues[k])
        #if 0: thisEffA =-1
        factor=1
        #if    np.append(thesevalues['llp1_Lxy'][thesevalues['llp1_eta'] > 1.5], thesevalues['llp2_Lxy'][thesevalues['llp2_eta'] > 1.5]).mean()*factor < 0.25/factor : thisEffA= -1
        #elif  np.append(thesevalues['llp1_Lxy'][thesevalues['llp1_eta'] > 1.5], thesevalues['llp2_Lxy'][thesevalues['llp2_eta'] > 1.5]).mean()*factor > 30*factor : thisEffA= -1
        #elif  np.append(thesevalues['llp1_Lz'][thesevalues['llp1_eta'] < 1.5], thesevalues['llp2_Lz'][thesevalues['llp2_eta'] < 1.5]).mean()*factor  < 0.75/factor : thisEffA= -1
        #elif  np.append(thesevalues['llp1_Lz'][thesevalues['llp1_eta'] < 1.5], thesevalues['llp2_Lz'][thesevalues['llp2_eta'] < 1.5]).mean()*factor  > 60*factor : thisEffA= -1
        if    thesevalues['llp1_Lxy'][thesevalues['llp1_eta'] > 1.5].mean()*factor < 0.25/factor : thisEffA= -1
        elif  thesevalues['llp1_Lxy'][thesevalues['llp1_eta'] > 1.5].mean()*factor > 16*factor : thisEffA= -1
        elif  thesevalues['llp1_Lz'][thesevalues['llp1_eta'] < 1.5].mean()*factor  < 0.55/factor : thisEffA= -1
        elif  thesevalues['llp1_Lz'][thesevalues['llp1_eta'] < 1.5].mean()*factor  > 28*factor : thisEffA= -1
        elif  thesevalues['llp2_Lxy'][thesevalues['llp2_eta'] > 1.5].mean()*factor < 0.25/factor : thisEffA= -1
        elif  thesevalues['llp2_Lxy'][thesevalues['llp2_eta'] > 1.5].mean()*factor > 16*factor : thisEffA= -1
        elif  thesevalues['llp2_Lz'][thesevalues['llp2_eta'] < 1.5].mean()*factor  < 0.75/factor : thisEffA= -1
        elif  thesevalues['llp2_Lz'][thesevalues['llp2_eta'] < 1.5].mean()*factor  > 28*factor : thisEffA= -1
        else:
          thesevalues['llp1_Lxy'] = thesevalues['llp1_Lxy']
          thesevalues['llp2_Lxy'] = thesevalues['llp2_Lxy']
          thesevalues['llp1_Lz'] = thesevalues['llp1_Lz']
          thesevalues['llp2_Lz'] = thesevalues['llp2_Lz']
          thesevalues['llp1_ET'] = thesevalues['llp1_ET']
          thesevalues['llp1_pT'] = thesevalues['llp1_pT']
          thesevalues['llp2_ET'] = thesevalues['llp2_ET']
          thesevalues['llp2_pT'] = thesevalues['llp2_pT']
          thesevalues['V_pt'] = thesevalues['V_pt']
          mean, std, var, clf =  load_model(sel)
          thisEffA = bdt_eval(pd.DataFrame(thesevalues), sel=sel, mean=mean, std=std, var=var, clf=clf) 
        if thisEffA==-1 : thisEffA=1e-6
        eff.append(thisEffA)
   if do2D: eff = np.array(eff).reshape(len(tauN),len(tauN))
   return eff

def load_model(sel):
   modelDir = "models"
   scaler_mean =  np.load(f"{modelDir}/{sel}_scaler_mean.npy")
   scaler_std =  np.load(f"{modelDir}/{sel}_scaler_std.npy")
   with open(f"{modelDir}/{sel}_features.txt") as g: var=eval(g.read())
   f = open(f'{modelDir}/{sel}_model.pkl', 'rb')
   clf = pickle.load(f)
   newVar =[]
   for varName in var:
     if "W" in varName: newVar+=[varName.replace("W","V")]
     elif "Z" in varName: newVar+=[varName.replace("Z","V")]
     else: newVar+=[varName]
   return scaler_mean, scaler_std, newVar, clf

def bdt_eval(values, clf=None, mean=None, std=None, var =None, sel=None):
  if clf is None:
    mean, std, var, clf =  load_model(sel)
  X = np.array(values[var])
  X = (X - mean) /std
  pred_proba = clf.predict_proba(X).T
  den = len(pred_proba[1])
  return sum(pred_proba[1])/den #, sum(pred_proba[1])/den,  sum(pred_proba[2])/den ,  sum(pred_proba[3])/den, sum(pred_proba[4])/den


   
def bdt_eval_list(values, clf=None, mean=None, std=None, var =None, sel=None):
  if clf is None:
    mean, std, var, clf =  load_model(sel)
  X = np.array(values[var])
  X = (X - mean) /std
  pred_proba = clf.predict_proba(X).T
  den = len(pred_proba[1])
  return pred_proba[1]

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################

#########################################################################################
#Parsing the lhe file from the MG output to recover the data from the process.
#########################################################################################


def parsing_LHE(MG_events, mediator=25, vectorBoson=-999, LLP=35):
    variables = ['px', 'py', 'pz', 'E', 'mass', 'decay'] 
    result={
    'Phi':pd.DataFrame(columns=variables),
    'V':pd.DataFrame(columns=variables),
    'llp1':pd.DataFrame(columns=variables),
    'llp2':pd.DataFrame(columns=variables),
    }
    for event in MG_events:
        counter=1
        for particle in event:
            if particle.pdg == mediator: # PDGID 35 Dark Higgs
                result['Phi'].loc[len(result['Phi'])] = [particle.px, particle.py, particle.pz, particle.E, particle.mass, particle.pdg]
            if particle.pdg == vectorBoson: # PDGID 35 Dark Higgs
                result['V'].loc[len(result['V'])] = [particle.px, particle.py, particle.pz, particle.E, particle.mass, particle.pdg]
            if particle.pdg == LLP: # PDGID 35 Dark Higgs
                result[f'llp{counter}'].loc[len(result[f'llp{counter}'])] = [particle.px, particle.py, particle.pz, particle.E, particle.mass, 5]
                counter+=1

    return result 


#########################################################################################
# Recovering the data from MG (LHE) (PDG ID, px,py,pz,E,mass).
#########################################################################################

def recover_MG_DH1(px, py, pz, E, MASS, pdg):

    MG_pdg_DH1_1 = []
    for i in range(5,len(pdg),9):
        MG_pdg_DH1_1.append(pdg[i]) #List with the PDG ID of the particle produced by the decay of the LLP1

    MG_E_DH1 = []
    for i in range(0,len(px),2):
        MG_E_DH1.append(E[i]) #List with the energy of the LLP1

    MG_px_DH1 = []
    for i in range(0,len(px),2):
        MG_px_DH1.append(px[i]) #List with x momenta from LLP1

    MG_py_DH1 = []
    for i in range(0,len(px),2):
        MG_py_DH1.append(py[i]) #List with y momenta from LLP1

    MG_pz_DH1 = []
    for i in range(0,len(px),2):
        MG_pz_DH1.append(pz[i]) #List with z momenta from LLP1

    MG_mass_DH1 = []
    for i in range(0,len(px),2):
        MG_mass_DH1.append(MASS[i]) #List with the mass from LLP1

    MG_px_DH1 = np.array(MG_px_DH1) # convertion into arrays
    MG_py_DH1 = np.array(MG_py_DH1)
    MG_pz_DH1 = np.array(MG_pz_DH1)
    MG_E_DH1 = np.array(MG_E_DH1)
    MG_mass_DH1 = np.array(MG_mass_DH1)

    return MG_px_DH1, MG_py_DH1,MG_pz_DH1,MG_E_DH1,MG_mass_DH1,MG_pdg_DH1_1

#########################################################################################
# Computation of the kinematics variable for LLP1 (velocities, beta, gamma, pT the transverse momenta, eta the pseudo-rapidity).
#########################################################################################

def kinematics_MG_DH1(MG_px_DH1,MG_py_DH1,MG_pz_DH1,MG_E_DH1 ):

    MG_vx_DH1 = (MG_px_DH1*c**2)/MG_E_DH1 #compute the velocities in each direction
    MG_vy_DH1 = (MG_py_DH1*c**2)/MG_E_DH1
    MG_vz_DH1 = (MG_pz_DH1*c**2)/MG_E_DH1
    MG_beta_DH1 = np.sqrt(MG_vx_DH1**2 + MG_vy_DH1**2 + MG_vz_DH1**2)/c # compute beta
    MG_gamma_DH1 = 1/(np.sqrt(1-MG_beta_DH1**2)) # compute gamma
    MG_pT_DH1 = np.sqrt(MG_px_DH1**2 + MG_py_DH1**2)*c # compute the transverse momenta
    MG_eta_DH1 = np.arctanh(MG_pz_DH1/(np.sqrt(MG_px_DH1**2 + MG_py_DH1**2 + MG_pz_DH1**2))) # compute the pseudorapidity

    return MG_pT_DH1,MG_eta_DH1, MG_gamma_DH1

#########################################################################################
# Recovering the data from LLP2 (PDG ID, px,py,pz,E,mass).
#########################################################################################

def recover_MG_DH2(px, py, pz, E, MASS, pdg):

    MG_pdg_DH2_1 = []
    for i in range(7,len(pdg),9):
        MG_pdg_DH2_1.append(pdg[i]) #List with the PDG ID of the particle produced by the decay of the LLP2

    MG_E_DH2 = []
    for i in range(1,len(px),2):
        MG_E_DH2.append(E[i]) #List with the energy of the LLP2

    MG_px_DH2 = []
    for i in range(1,len(px),2):
        MG_px_DH2.append(px[i]) #List with x momenta from LLP2

    MG_py_DH2 = []
    for i in range(1,len(px),2):
        MG_py_DH2.append(py[i]) #List with y momenta from LLP2

    MG_pz_DH2 = []
    for i in range(1,len(px),2):
        MG_pz_DH2.append(pz[i]) #List with z momenta from LLP2

    MG_mass_DH2 = []
    for i in range(1,len(px),2):
        MG_mass_DH2.append(MASS[i]) #List with the mass from LLP2

    MG_px_DH2 = np.array(MG_px_DH2) # convertion into arrays
    MG_py_DH2 = np.array(MG_py_DH2)
    MG_pz_DH2 = np.array(MG_pz_DH2)
    MG_E_DH2 = np.array(MG_E_DH2)
    MG_mass_DH2 = np.array(MG_mass_DH2)

    return MG_px_DH2, MG_py_DH2,MG_pz_DH2,MG_E_DH2,MG_mass_DH2,MG_pdg_DH2_1

#########################################################################################
# Computation of the kinematics variable for LLP2 (velocities, beta, gamma, pT the transverse momenta, eta the pseudo-rapidity).
#########################################################################################

def kinemamtics_MG_DH2(MG_px_DH2,MG_py_DH2,MG_pz_DH2,MG_E_DH2):

    MG_vx_DH2 = (MG_px_DH2*c**2)/MG_E_DH2 #compute the velocities in each direction
    MG_vy_DH2 = (MG_py_DH2*c**2)/MG_E_DH2
    MG_vz_DH2 = (MG_pz_DH2*c**2)/MG_E_DH2
    MG_beta_DH2 = np.sqrt(MG_vx_DH2**2 + MG_vy_DH2**2 + MG_vz_DH2**2)/c # compute beta
    MG_gamma_DH2 = 1/(np.sqrt(1-MG_beta_DH2**2)) # compute gamma

    MG_pT_DH2 = np.sqrt(MG_px_DH2**2 + MG_py_DH2**2)*c # compute the transverse momenta
    MG_eta_DH2 = np.arctanh(MG_pz_DH2/(np.sqrt(MG_px_DH2**2 + MG_py_DH2**2 + MG_pz_DH2**2))) # compute the pseudorapidity

    return MG_pT_DH2,MG_eta_DH2, MG_gamma_DH2

#########################################################################################
# Decay lenght computation for LLP1.
#########################################################################################

def decaylenght_MG_DH1(MG_px_DH1, MG_py_DH1, MG_pz_DH1, E_DH1, MG_gamma_DH1, tauN):

    MG_Lx_tot_DH1 = []
    MG_Ly_tot_DH1 = []
    MG_Lz_tot_DH1 = []
    MG_Lxy_tot_DH1 = []

    for ctau in range(len(tauN)):

        MG_Lx_DH1 = []
        MG_Ly_DH1 = []
        MG_Lz_DH1 = []
        MG_Lxy_DH1 = []

        for i in range(len(MG_gamma_DH1)):
            MG_lt = lifetime(tauN[ctau]) # set the mean lifetime
            MG_Lx_DH1.append((MG_px_DH1[i]/E_DH1[i])*c**2 * MG_lt * MG_gamma_DH1[i]) # compute the decay lenght in x,y,z
            MG_Ly_DH1.append((MG_py_DH1[i]/E_DH1[i])*c**2 * MG_lt * MG_gamma_DH1[i])
            MG_Lz_DH1.append((abs(MG_pz_DH1[i])/E_DH1[i])*c**2 * MG_lt  * MG_gamma_DH1[i] )
            MG_Lxy_DH1.append(np.sqrt((MG_Lx_DH1[i])**2 + (MG_Ly_DH1[i])**2)) # compute the transverse decay lenght

        MG_Lx_tot_DH1.append(MG_Lx_DH1) # convertion into arrays
        MG_Ly_tot_DH1.append(MG_Ly_DH1)
        MG_Lz_tot_DH1.append(MG_Lz_DH1)
        MG_Lxy_tot_DH1.append(MG_Lxy_DH1)

    return MG_Lxy_tot_DH1, MG_Lz_tot_DH1

#########################################################################################
# Decay lenght computation for LLP2.
#########################################################################################

def decaylenght_MG_DH2(MG_px_DH2, MG_py_DH2, MG_pz_DH2, E_DH2, MG_gamma_DH2, tauN):

    MG_Lx_tot_DH2 = []
    MG_Ly_tot_DH2 = []
    MG_Lz_tot_DH2 = []
    MG_Lxy_tot_DH2 = []

    for ctau in range(len(tauN)):

        MG_Lx_DH2 = []
        MG_Ly_DH2 = []
        MG_Lz_DH2 = []
        MG_Lxy_DH2 = []

        for i in range(len(MG_gamma_DH2)):
            MG_lt = lifetime(tauN[ctau]) # set the mean lifetime
            MG_Lx_DH2.append((MG_px_DH2[i]/E_DH2[i])*c**2 * MG_lt * MG_gamma_DH2[i]) # compute the decay lenght in x,y,z
            MG_Ly_DH2.append((MG_py_DH2[i]/E_DH2[i])*c**2 * MG_lt * MG_gamma_DH2[i])
            MG_Lz_DH2.append((abs(MG_pz_DH2[i])/E_DH2[i])*c**2 * MG_lt  * MG_gamma_DH2[i] )
            MG_Lxy_DH2.append(np.sqrt((MG_Lx_DH2[i])**2 + (MG_Ly_DH2[i])**2)) # compute the transverse decay lenght

        MG_Lx_tot_DH2.append(MG_Lx_DH2)
        MG_Ly_tot_DH2.append(MG_Ly_DH2)
        MG_Lz_tot_DH2.append(MG_Lz_DH2)
        MG_Lxy_tot_DH2.append(MG_Lxy_DH2)

    return MG_Lxy_tot_DH2, MG_Lz_tot_DH2

#########################################################################################
# Computation of the efficiency with the map from the data obtained with MG for the high-ET samples (mH <= 400GeV).
#########################################################################################

def eff_map_MG_high(MG_pT_DH1, MG_eta_DH1,MG_Lxy_tot_DH1, MG_Lz_tot_DH1, MG_pdg_DH1_1, MG_pT_DH2, MG_eta_DH2, MG_Lxy_tot_DH2, MG_Lz_tot_DH2, MG_pdg_DH2_1, tauN, nevent, mass_phi, mass_s):

    MG_eff_highETX = []

    for index in tqdm.tqdm(range(len(tauN))):
        MG_queryMapResult = []
        for iEvent in range(len(MG_pT_DH1)):
            MG_queryMapResult.append(rmN.queryMapFromKinematics(MG_pT_DH1[iEvent],
                                                            MG_eta_DH1[iEvent],
                                                            MG_Lxy_tot_DH1[index][iEvent],
                                                            MG_Lz_tot_DH1[index][iEvent],
                                                            abs(MG_pdg_DH1_1[iEvent]),
                                                            MG_pT_DH2[iEvent],
                                                            MG_eta_DH2[iEvent],
                                                            MG_Lxy_tot_DH2[index][iEvent],
                                                            MG_Lz_tot_DH2[index][iEvent],
                                                            abs(MG_pdg_DH2_1[iEvent]),
                                                            selection = "high-ET"))
        MG_eff_highETX.append(sum(MG_queryMapResult))
    MG_queryMapResult = np.array(MG_queryMapResult) # convertion into arrays
    MG_eff_highETX = np.array(MG_eff_highETX) # convertion into arrays
    MG_eff_highETX = MG_eff_highETX/nevent #eff/nbrevent

    MG_Data_Eff_High = np.column_stack(MG_eff_highETX)
    np.savetxt(f'./Plots_High/Efficiencies_Text_{mass_phi}_{mass_s}.txt', MG_Data_Eff_High)

    return MG_eff_highETX

#########################################################################################
# Computation of the efficiency with the map from the data obtained with MG for the low-ET samples (mH <= 400GeV).
#########################################################################################

def eff_map_MG_low(MG_pT_DH1, MG_eta_DH1,MG_Lxy_tot_DH1, MG_Lz_tot_DH1, MG_pdg_DH1_1, MG_pT_DH2, MG_eta_DH2, MG_Lxy_tot_DH2, MG_Lz_tot_DH2, MG_pdg_DH2_1, tauN, nevent, mass_phi, mass_s):

    MG_eff_lowETX = []

    for index in tqdm.tqdm(range(len(tauN))):
        MG_queryMapResult = []
        for iEvent in range(len(MG_pT_DH1)):
            MG_queryMapResult.append(rmN.queryMapFromKinematics(MG_pT_DH1[iEvent],
                                                            MG_eta_DH1[iEvent],
                                                            MG_Lxy_tot_DH1[index][iEvent],
                                                            MG_Lz_tot_DH1[index][iEvent],
                                                            abs(MG_pdg_DH1_1[iEvent]),
                                                            MG_pT_DH2[iEvent],
                                                            MG_eta_DH2[iEvent],
                                                            MG_Lxy_tot_DH2[index][iEvent],
                                                            MG_Lz_tot_DH2[index][iEvent],
                                                            abs(MG_pdg_DH2_1[iEvent]),
                                                            selection = "low-ET"))
        MG_eff_lowETX.append(sum(MG_queryMapResult))
    MG_queryMapResult = np.array(MG_queryMapResult) # convertion into arrays
    MG_eff_lowETX = np.array(MG_eff_lowETX) # convertion into arrays
    MG_eff_lowETX = MG_eff_lowETX/nevent #eff/nbrevent

    MG_Data_Eff_Low = np.column_stack(MG_eff_lowETX)
    np.savetxt(f'./Plots_Low/Efficiencies_Text_{mass_phi}_{mass_s}.txt', MG_Data_Eff_Low)

    return MG_eff_lowETX


#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################

#########################################################################################
# Computing HEP data
#########################################################################################

def elem_list(HEP, File_HEP_limit) :

    file_HEP = uproot.open(HEP) # open the file from HEP data for the efficiency
    data_HEP = file_HEP[file_HEP.keys()[1]] # open the branch

    file_HEP_limit = uproot.open(File_HEP_limit) # open the file from HEP data for the limits
    branch_HEP_limit = file_HEP_limit[file_HEP_limit.keys()[2]] # open the branch

    return data_HEP, branch_HEP_limit

#########################################################################################
# Plots to compare the results of efficiency obtained with MG, MG+Pythia8 (High-ET).
#########################################################################################

def plt_eff2D(eff_highETX, tauN,  mass_phi , mass_s, mass_s2=None, model="HSS", sel="", data_HEP=None) :
    
    slug = f"{sel}_{model}_mH{mass_phi}_mS{mass_s}"
    if model=="HS1S2": slug = f"{sel}_{model}_mH{mass_phi}_mS1{mass_s}_mS2{mass_s2}"
    if model=="ALP": slug = f"{sel}_{model}_mALP{mass_s}"
    np.save(f"./Plots/values2D_{slug}_ctau.npy", tauN)
    np.save(f"./Plots/values2D_{slug}_eff.npy", eff_highETX)
    ctau1 =  tauN
    eff1 = eff_highETX
    mask = np.max(eff1, axis=0) >0 
    # Generate some data
    x1, y1 = np.meshgrid(ctau1[mask], ctau1[mask])
    z1  =  eff1[mask].T[mask]
    z_min, z_max = 0, z1.max()
    
    # Create the figure and subplots
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    c = ax.pcolormesh(x1, y1, z1, cmap='Blues', vmin=z_min, vmax=z_max)
    ax.set_title('pcolormesh')
    # set the limits of the plot to the limits of the data
    ax.axis([x1.min(), x1.max(), y1.min(), y1.max()])
    fig.colorbar(c, ax=ax)
    
    # Upper subplot with two lines and uncertainty bands
    #ax1.plot(x, y1, label=f'{label1}')
    s = ax.contour(x1, y1, z1, [0.01, 0.05])
    #ax.contour(cs, colors='k')
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_ylabel("$c\\tau$ [m]")
    ax.set_xlabel("$c\\tau$ [m]")
    
    
    # Show the plot
    plt.tight_layout()
    plt.savefig(f"Plots/eff2D_{slug}.png")
    plt.savefig(f"Plots/eff2D_{slug}.pdf")
  


def plt_eff(eff_highETX, tauN,  mass_phi , mass_s, mass_s2=None, model="HSS", sel="", data_HEP=None) :
    
    eff_highETX,tauN = np.array(eff_highETX), np.array(tauN)
    # only plot within valid range
    mask =  eff_highETX >= 0
    
    eff_highETX = eff_highETX[mask]
    tauN = tauN[mask]
    ################## PLOT EFFICIENCY ##################
    fig, ax = plt.subplots()

    ################## Plot efficiency from MG+Pythia8 ##################
    plt.plot(tauN,eff_highETX, 'r', linewidth=2, label = 'BDT ')
    
    if data_HEP is not None:
       ################## Plot efficiency from HEP data ##################
       plt.plot(data_HEP.values(axis='both')[0],data_HEP.values(axis='both')[1], 'b')

       ################ Uncertainties from HEP ##################
       plt.fill_between(data_HEP.values(axis='both')[0], data_HEP.values(axis='both')[1] +  data_HEP.errors('high')[1] , data_HEP.values(axis='both')[1] - data_HEP.errors('high')[1] , color = 'blue', label = r'ATLAS, with $\pm$ 1 $\sigma$ error bands',alpha=.7)
 
    ################## Uncertainties from Map ##################
    plt.fill_between(tauN, np.array(eff_highETX) + 0.30* np.array(eff_highETX), np.array(eff_highETX) - 0.30 * np.array(eff_highETX), label='30% error bands ', color='r', alpha=.7)

    ################## Limits of validity ##################
    # place a text box in upper left in axes coords
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    if model=="HS1S2": ax.text(0.05, 0.95, f"{sel} $ m_Φ $ = {mass_phi} GeV, $m_S1$ = {mass_s} GeV, $m_S2$ = {mass_s2} GeV" , transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    if model=="HSS": ax.text(0.05, 0.95, f"{sel} $ m_Φ $ = {mass_phi} GeV, $m_S$ = {mass_s} GeV" , transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    if model=="ALP": ax.text(0.05, 0.95, f"{sel} $m_a$ = {mass_s} GeV" , transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    #x = np.linspace(0.5,10)
    #ax.fill_between(x, max(eff_highETX), color='orange', alpha=.2, hatch="/", edgecolor="black", linewidth=1.0, label = 'validity window') # adding hatch

    plt.xscale('log')
    plt.ylim([0, max(eff_highETX)*2]) # start at 0
    #plt.ylim([1e-4,1]) # start at 0
    #plt.yscale('log')
    plt.xlabel(r'c$\tau$ [m]', fontsize=20)
    plt.ylabel('Efficiency', fontsize=20 )
    plt.legend(fontsize = 11) # set the legend in the upper right corner
    slug = f"{sel}_{model}_mH{mass_phi}_mS{mass_s}"
    if model=="HS1S2": slug = f"{sel}_{model}_mH{mass_phi}_mS1{mass_s}_mS2{mass_s2}"
    if model=="ALP": slug = f"{sel}_{model}_mALP{mass_s}"
    plt.savefig(f"./Plots/Efficiency_comparison_{slug}.png")
    plt.savefig(f"./Plots/Efficiency_comparison_{slug}.pdf")
    print(f"./Plots/Efficiency_comparison_{slug}.png")
    plt.close()
    np.save(f"./Plots/values_{slug}_ctau.npy", tauN)
    np.save(f"./Plots/values_{slug}_eff.npy", eff_highETX)

def plt_multi_eff(eff_ctau_pairs_dict, ref=None,  model="HSS", sel=""):
    

    ################## PLOT EFFICIENCY ##################
    fig, ax = plt.subplots()
    maxVal = -1 
    
    plt.plot(-1, -1, 'k-', label ="$\mathbfit{ATLAS}$ full analysis")
    plt.plot(-1, -1, 'k--', label ="MG5 + Py8 + BDT")
    
    colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
    colors *= 10
    cIndex = -1
    for name, values in eff_ctau_pairs_dict.items():
      cIndex +=1
     
      for refName in ref.keys():
        #print (f"{name} in {refName.replace('0p','0.')}?")
        if name in refName or name in refName.replace("GeV",".0GeV") or name in refName.replace("0p","0."):
          #print("->yes!")
          plt.plot(ref[refName][1],ref[refName][0], '-', color=colors[cIndex], linewidth=2)

      eff, ctau = values
      if maxVal < max(eff): maxVal=max(eff)
      name = name.replace("mALP", "$m_\mathrm{ALP}=$ ")
      name = name.replace(".0",  "")
      name = name +" GeV"
      if "*" in model: model =  model.split("_")[0]
      if name[0] != "m": name="mH"+name
      plt.plot(ctau,eff, '--', linewidth=2, color=colors[cIndex], label = name)
      plt.fill_between(ctau, np.array(eff) * 1.30, np.array(eff) * 0.7, color=colors[cIndex], alpha=.7)
    

    #ax.text(0.95, 0.95, "$\mathbfit{ATLAS}$ $\mathit{Simulation}$ $\mathit{Internal}$" , transform=ax.transAxes, fontsize=17, verticalalignment='top', horizontalalignment='right')
    ax.text(0.05, 0.95, f"Model: {model}" , transform=ax.transAxes, fontsize=14, verticalalignment='top', horizontalalignment='left')
    ax.text(0.05, 0.9, f"Selection: {sel}" , transform=ax.transAxes, fontsize=14, verticalalignment='top', horizontalalignment='left')

    plt.xscale('log')
    plt.ylim([0, maxVal*2]) # start at 0
    plt.xlabel(r'c$\tau$ [m]', fontsize=20)
    plt.ylabel('Efficiency', fontsize=20 )
    plt.legend(fontsize = 14, loc='upper right') # set the legend in the upper right corner
    plt.savefig(f"./PlotsSummary/Efficiency_summary_{sel}_{model}.png")
    plt.savefig(f"./PlotsSummary/Efficiency_summary_{sel}_{model}.pdf")
    print(f"./PlotsSummary/Efficiency_summary_{sel}_{model}.png")


#########################################################################################
# Plot limits obtained with the map, to compare with those obtain by ATLAS (High-ET).
#########################################################################################

def plt_cross(eff_highETX, tauN, mass_phi, mass_s, branch_HEP_limit, factor, hepdata_eff=None, label=""):

    fig, ax = plt.subplots()

    Nsobs = 0.5630 * 26 * factor # nbr of observed events = 26 ( factor )

    Crr_Sec_obs = (Nsobs)/((np.array(eff_highETX)) * 139e3 ) # Luminosity = 139e3 fb**(-1)

    plt.plot(tauN, Crr_Sec_obs, 'r', label ='Map results', linewidth = 2)
    plt.fill_between(tauN,  1.25* np.array(Crr_Sec_obs), 0.75 * np.array(Crr_Sec_obs), label='25\% error bands ', color='r', alpha=.7)
    if hepdata_eff is not None: 
       Crr_Sec_obs_hepdata_eff = (Nsobs)/((np.array(hepdata_eff.values(axis='both')[1])) * 139e3 ) # Luminosity = 139e3 fb**(-1)
       plt.plot(hepdata_eff.values(axis='both')[0], Crr_Sec_obs_hepdata_eff, 'g', label ='Observed', linewidth = 2)
    if branch_HEP_limit is not None:
      plt.plot(np.array(branch_HEP_limit.values(axis='both')[0]), np.array(branch_HEP_limit.values(axis='both')[1]), 'b', label ='Observed', linewidth = 2)

    x = np.linspace(0.5,10)
    ax.fill_between(x, max(Crr_Sec_obs)*1.1, color='orange', alpha=.2, hatch="/", edgecolor="black", linewidth=1.0, label = 'validity window') # adding hatch
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([1e-3, 1e3])
    plt.xlabel(r'c$\tau$ [m]')
    plt.ylabel(r'95% CL limit on $\sigma \times B$ [pb]')

    # place a text box in upper left in axes coords
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax.text(0.05, 0.95, f" $ m_Φ $ = {mass_phi} GeV, $m_S$ = {mass_s} GeV" , transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    plt.legend( fontsize = 10, loc=3)
    plt.savefig(f"./Plots/Cross_section_{label}_mH{mass_phi}_mS{mass_s}.png") #create a new fodlder ' Plots ' and save the fig in it
    plt.savefig(f"./Plots/Cross_section_{label}_mH{mass_phi}_mS{mass_s}.pdf") #create a new fodlder ' Plots ' and save the fig in it
    print(f"./Plots/Cross_section_{label}_mH{mass_phi}_mS{mass_s}.png") #create a new fodlder ' Plots ' and save the fig in it
    plt.close()

