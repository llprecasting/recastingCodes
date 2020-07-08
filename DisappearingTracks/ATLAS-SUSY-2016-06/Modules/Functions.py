#!/usr/bin/env python

## Last version of Function.py

import os, copy
import numpy as np
import random
import math
from ROOT import TLorentzVector
from __future__ import print_function

def PrintEvent(particle,jet,electron,muon,met):

    print('\n----- Parton level -----')
    for p in particle:
    ## jet at parton level
       if p.Status == 23 and p.PID < 23 :
	  print( " PID = {:.0f}        Stat={:.0f}  Px={:.2f}  Py={:.2f}  PT={:.2f}".format(p.PID, p.Status, p.P4().Px(), p.P4().Py(), (p.P4().Px()**2 + p.P4().Py()**2)**0.5))
    ## neutralinos and charginos at parton level
       if( (p.Status == 23 and p.PID == 1000022) or (p.Status == 22 and abs(p.PID) == 1000024) ):
  	  print( " PID = {:.0f}  Stat={:.0f}  Px={:.2f}  Py={:.2f}  PT={:.2f}".format(p.PID, p.Status, p.P4().Px(), p.P4().Py(), (p.P4().Px()**2 + p.P4().Py()**2)**0.5))
    print( '----- reconstruction level -----')
    for i,j in enumerate(jet):
	print( '{:.0f} Jets: flavor = {:.0f} Px = {:.2f}  Py = {:.2f}  PT = {:.2f} '.format(i+1,j.Flavor,j.P4().Px(), j.P4().Py(), j.PT))
    if ( electron.GetEntries()==0 ): print( 'no electrons')
    for i,e in enumerate(electron):
	print( '{:.0f} recons Elect: Px = {:.2f}  Py = {:.2f}  PT = {:.2f} '.format(i+1,e.P4().Px(), e.P4().Py(), e.PT))
    if ( muon.GetEntries()==0 ): print( 'no muon')
    for i,m in enumerate(muon):
	print( '{:.0f} recons Muons: Px = {:.2f}  Py = {:.2f}  PT = {:.2f} '.format(i+1,m.P4().Px(), m.P4().Py(), m.PT))
    for i,n in enumerate(particle):
      if( n.Status==1 and (abs(p.PID) == 12 or abs(p.PID) == 14 or abs(p.PID) == 16) ):
	print( '{:.0f} recons nuetrinos: Px = {:.2f}  Py = {:.2f}  PT = {:.2f} '.format(i+1,m.P4().Px(), m.P4().Py(), m.PT))
    for me in met:
       print('Delphes MET: Px = {:.2f}  Py = {:.2f}  MET = {:.2f}'.format(me.P4().Px(), me.P4().Py(), me.MET))

def PrintParticle(objects):
    for i,obj in enumerate(objects):
       print( "{}  PID ={:.0f}  Stat={}  Px={:.2f}  Py={:.2f}  PT={:.2f}".format(i,obj.PID, obj.Status, obj.Px, obj.Py, math.sqrt(obj.Px**2 + obj.Py**2)))

def PrintJet(objects):
    Pxjet = 0
    Pyjet = 0
    for obj in objects:
       print( "PT ={} Eta={}  Phi={}".format(obj.PT,obj.Eta,obj.Phi))
       Pxjet = Pxjet + obj.P4().Px()
       Pyjet = Pyjet + obj.P4().Py()
    print( 'Jet PT = ', (Pxjet**2 + Pyjet**2)**0.5)

def selectChargino(objects, Status, PID=None):
    """ Select chargino from Particle branch and store the information into a new branch """
    if(PID is not None):
        PID_char = PID
    else:
        PID_char = 1000024
    new_objects = copy.deepcopy(objects)
    for obj in new_objects:
	if( abs(obj.PID) != PID_char or obj.Status != Status ):
		new_objects.RecursiveRemove(obj)
    return new_objects

def selectNeutralino(objects, Status, PID=None):
    """ Select Neutralino from Particle branch and store the information into a new branch"""
    if(PID is not None):
        PID_neutr = PID
    else:
        PID_neutr = 1000022
    new_objects = copy.deepcopy(objects)
    for obj in new_objects:
	if( abs(obj.PID) != PID_neutr or obj.Status != Status ):
		new_objects.RecursiveRemove(obj)
    return new_objects

def selectJet(objects, Status):
    """ Select Jet at parton level from Particle branch and store the information into a new branch"""
    new_objects = copy.deepcopy(objects)
    for obj in new_objects:
	if( abs(obj.PID) > 30 or obj.Status != Status ):
		new_objects.RecursiveRemove(obj)
    return new_objects

def selectNeutrino(objects, status):
    """ Select neutrino from Particle branch and store the information into a new branch """
    new_objects = copy.deepcopy(objects)
    for obj in new_objects:
	if( abs(obj.PID) != 12 and abs(obj.PID) != 14 and abs(obj.PID) != 16 or obj.Status != status ):
	        new_objects.RecursiveRemove(obj)
    return new_objects

def filterPhaseSpace(objects,pt,eta):
    """ filter for pt and eta """
    for obj in objects:
        if( obj.PT < pt or abs(obj.Eta) > eta  ):
	    objects.RecursiveRemove(obj)
	else:
	    continue #print(obj.PT, obj.Eta)
    return objects

def TriggerMET(objects,MET_max):
    """ filter for met """
    for obj in objects:
        if( obj.MET < MET_max ):
	    objects.RecursiveRemove(obj)
	else:
	    continue #print(obj.PT, obj.Eta)
    return objects

def overlapRemoval(candidates,neighbours,dR):
    if neighbours.GetEntries()==0: return candidates
    for candidate in candidates:
        for neighbour in neighbours:
            if candidate.P4().DeltaR(neighbour.P4()) < dR:
                candidates.RecursiveRemove(candidate)
    return candidates

def overlapRemovalMS(candidates,neighbours,dR):
    if neighbours.GetEntries()==0: return candidates
    for candidate in candidates:
        for neighbour in neighbours:
	   if(abs(neighbour.PID)==13):
             if candidate.P4().DeltaR(neighbour.P4()) < dR:
                candidates.RecursiveRemove(candidate)
    return candidates

def overlapRemovalUpTo(candidates,neighbours,dR,n):
    if neighbours.GetEntries()==0: return candidates
    for candidate in candidates:
        for j,neighbour in enumerate(neighbours):
 	    if j<n:
	       if candidate.P4().DeltaR(neighbour.P4()) < dR:
                  candidates.RecursiveRemove(candidate)
    return candidates

def overlapRemovalInterval(candidates,neighbours,dR_min,dR_max):
    if neighbours.GetEntries()==0: return candidates
    for candidate in candidates:
        for neighbour in neighbours:
            DR = candidate.P4().DeltaR(neighbour.P4())
	    if (DR > dR_min and DR < dR_max):
	       candidates.RecursiveRemove(candidate)
    return candidates

def overlapRemovalTracks(tracks,d0,zsinT,Rmin):
    for tr in tracks:
      if(tr.D0 > d0):
        if( abs(tr.DZ)*math.sin(tr.P4().Theta()) > 10 ):
	  R = tr.L*math.sin(tr.P4().Theta())
          if( R < Rmin ):
            tracks.RecursiveRemove(tr)
    return tracks

def isolated(candidate,neighbours,dr):
    if neighbours.GetEntries()==0: return True
    nPT = 0
    for n in neighbours:
	if( candidate.P4().DeltaR(n.P4()) < dr ):
	  if (n.PT > 1 and abs(n.DZ*np.sin(n.P4().Theta()))<3 ):
	     nPT = nPT + n.PT
    CONE = nPT/candidate.PT
    if (CONE < 0.04): return True
    else: return False

def separated(candidate,neighbours,dr):
    if neighbours.GetEntries()==0: return True
    counter = 0
    for n in neighbours:
	if( candidate.P4().DeltaR(n.P4()) < dr ):
	  counter = counter + 1
    if (counter == 0): return True
    else: return False

def separatedMS(candidate,neighbours,dr):
    if neighbours.GetEntries()==0: return True
    counter = 0
    for n in neighbours:
     if(abs(n.PID)==13):
	if( candidate.P4().DeltaR(n.P4()) < dr ):
	  counter = counter + 1
    if (counter == 0): return True
    else: return False

def separatedUpTo(candidate,neighbours,dr,num):
    if neighbours.GetEntries()==0: return True
    counter = 0
    for i,n in enumerate(neighbours):
      if(i<num):
	if( candidate.P4().DeltaR(n.P4()) < dr ):
	  counter = counter + 1
    if (counter == 0): return True
    else: return False

def dPhiMin(jets, MET, n):
    dPhimin = 10
    for m in MET:
      for i,jet in enumerate(jets):
	if i<n:
	   dPhimin = min( abs( m.P4().DeltaPhi(jet.P4()) ) , dPhimin )
    return dPhimin

def METGL(neutralinos, charginos, MET):
	met_n = TLorentzVector() #Define TLorentzVector for neutralinos and charginos
 	met_c = TLorentzVector()

	MET_nX = 0
	MET_nY = 0
	for n in neutralinos:
		MET_nX = MET_nX + n.P4().Px()
		MET_nY = MET_nY + n.P4().Py()
	met_n.SetPxPyPzE(MET_nX,MET_nY,0,0)

	MET_cX = 0
	MET_cY = 0
	for c in charginos:
		MET_cX = MET_cX + c.P4().Px()
		MET_cY = MET_cY + c.P4().Py()
	met_c.SetPxPyPzE(MET_cX,MET_cY,0,0)

    ## Sum of the vector components
	met_gl = met_n + met_c

        for m in MET:
	    m.MET = abs(met_gl.Mag())
	    m.Phi = met_gl.Phi()
	return MET

def METGL2(neutralinos, charginos, neutrinos, MET): # neutrinos,
    #Define TLorentzVector for neutralinos and charginos
	met_n = TLorentzVector()
 	met_c = TLorentzVector()
 	met_nt = TLorentzVector()

	MET_nX = 0
	MET_nY = 0
	for n in neutralinos:
		MET_nX = MET_nX + n.P4().Px()
		MET_nY = MET_nY + n.P4().Py()
	met_n.SetPxPyPzE(MET_nX,MET_nY,0,0)

	MET_cX = 0
	MET_cY = 0
	for c in charginos:
		MET_cX = MET_cX + c.P4().Px()
		MET_cY = MET_cY + c.P4().Py()
	met_c.SetPxPyPzE(MET_cX,MET_cY,0,0)

	MET_ntX = 0
	MET_ntY = 0
	for nt in neutrinos:
		MET_ntX = MET_ntX + nt.P4().Px()
		MET_ntY = MET_ntY + nt.P4().Py()
	met_nt.SetPxPyPzE(MET_ntX,MET_ntY,0,0)

    ## Sum of the vector components
	met_gl = met_n + met_c + met_nt

        for m in MET:
	    m.MET = abs(met_gl.Mag())
	    m.Phi = met_gl.Phi()
	return MET

def OffMissingET(jets, electrons, muons, MET):
    """ Calculate the MissingET from jets, electrons and muons """
    #Define TLorentzVector for jets, electrons and muons
    met_j = TLorentzVector()
    met_e = TLorentzVector()
    met_mu = TLorentzVector()

    ## The x and y components must be with (-1)
    MET_jX = 0
    MET_jY = 0
    for j in jets:
        MET_jX = j.P4().Px() + MET_jX
        MET_jY = j.P4().Py() + MET_jY
    met_j.SetPxPyPzE(-MET_jX,-MET_jY,0,0)

    MET_eX = 0
    MET_eY = 0
    for e in electrons:
        MET_eX = e.P4().Px() + MET_eX
        MET_eY = e.P4().Py() + MET_eY
    met_e.SetPxPyPzE(-MET_eX,-MET_eY,0,0)

    MET_mX = 0
    MET_mY = 0
    for m in muons:
        MET_mX = m.P4().Px() + MET_mX
        MET_mY = m.P4().Py() + MET_mY
    met_mu.SetPxPyPzE(-MET_mX,-MET_mY,0,0)

    ##sum of the vector components
    offmet = met_j + met_e + met_mu

    for m in MET:
        m.MET = abs(offmet.Mag())
        m.Phi = offmet.Phi()
    return MET

def smear(candidates,delta):
    """Smear the 4-momentum of a particle according to Sasha's function"""
    PI = math.pi
    for c in candidates:
	EI = c.E
	U1 = random.random()
	U2 = random.random()
	RGEN1 = math.sqrt(-2.*math.log(U1))
	RGEN2 = 2*PI*U2
	RN = RGEN1*math.cos(RGEN2)
	EIP = 1. + delta/math.sqrt(EI)*RN

	c.Px = c.Px*EIP
	c.Py = c.Py*EIP
	c.Pz = c.Pz*EIP
	c.E  = c.E*EIP
	c.PT = c.PT*EIP
    return candidates

def smearingPT(candidates):
    """Smear the 4-momentum of a particle according to ATLAS function"""
    alpha = 1.67
    beta  = -1.72
    sigma = 13.2
    for c in candidates:
      DeltacPT = c.Charge/c.PT
      z = (DeltacPT - beta)/sigma

      if (z<-alpha):
         funcz = math.exp(alpha*(z+alpha/2))
      elif (z>-alpha and z < alpha):
         funcz = math.exp(-z**2/2)
      elif (z>alpha):
         funcz = math.exp(-alpha*(z-alpha/2))

      c.PT = c.PT * funcz
    return candidates

def CDF(tauo, P, M):  # (lifetime in particle frame, module of 3-momentum, mass)
    """ Cumulative Distribution Function C.D.F. """
    c = 2.998e+11
    U = random.random()
    l = -(c*tauo*P/M)*math.log(1-U)
    return l

def charginoIndex(objects):
	## select up to 1 chargino in an event ##
	chargino_Index = -2
	char_number = objects.GetEntries()
	if (char_number == 0): 	chargino_Index = -1
 	if (char_number == 1): 	chargino_Index = 0
	if (char_number == 2):
	   # select randomly (assuming up to two charginos)
	   U = random.uniform(0,1)*10
	   if (U<5.0): chargino_Index = 0
	   else: chargino_Index = 1
	return chargino_Index
