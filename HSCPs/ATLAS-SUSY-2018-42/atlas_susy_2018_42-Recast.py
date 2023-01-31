#!/usr/bin/env python3

import os,glob
import numpy as np
import pandas as pd
import glob
import pyslha
import time
import progressbar as P
from ATLAS_data.effFunctions import (getMuonRecoEff,getTriggerEff,getTrackEff,
                                     getSelectionEff,getTargetMass,getMassSelEff,
                                     massLong,massShort)

delphesDir = os.path.abspath("./DelphesHSCP")
os.environ['ROOT_INCLUDE_PATH'] = os.path.join(delphesDir,"external")

import ROOT
import xml.etree.ElementTree as ET


ROOT.gSystem.Load(os.path.join(delphesDir,"libDelphes.so"))

ROOT.gInterpreter.Declare('#include "classes/SortableObject.h"')
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')


def getHSCPCandidates(hscps,daughters):

    candidates = []
    for ip in range(hscps.GetEntries()):
        p = hscps.At(ip)
        if p.Status == 1:
            p.r_decay = np.inf
            p.z_decay = np.inf
            candidates.append(p)
            continue
        
        # If HSCP is unstable compute its decay position from daughter
        d1 = p.D1
        daughter = daughters.At(d1)
        # Get the HSCP decay radius from the first daughter production vertex:
        p.r_decay = np.sqrt(daughter.X**2+daughter.Y**2)      
        p.z_decay = daughter.Z        
        candidates.append(p)

    for p in candidates:
        # Add beta and gamma*beta
        trimom = np.sqrt(p.Px**2 +p.Py**2 + p.Pz**2)
        p.beta = trimom/p.E
        p.gbeta = trimom/p.Mass

    return candidates

def applyHSCPSelection(hscpList,pT=50.,eta=2.4,r=500.):

    selHSCPs = []
    for hscp in hscpList:
        if hscp.PT < pT: continue 
        if abs(hscp.Eta) > eta: continue
        if hscp.r_decay < r: continue
        selHSCPs.append(hscp)
    
    return selHSCPs

def applyIsolation(hscpList,pTmax=5.0):

    isoHSCPs = []
    # Apply isolation requirement for HSCP tracks
    for hscp in hscpList:
        sumPT = hscp.SumPtCharged
        if sumPT > pTmax: continue
        isoHSCPs.append(hscp)
    return isoHSCPs


def applyMuonTagging(hscpList):

    """
    Computes the probability of reconstructing the hscp as a muon.
    :param hscpList: List of GenParticle objects
    """

    muonsLLP = []
    for hscp in hscpList:
        if hscp.r_decay < 3.9e3 and hscp.z_decay < 6e3: # Skip decays before MS
            continue
        
        beta = hscp.beta
        eta = abs(hscp.Eta)
        eff = getMuonRecoEff(beta,eta,hscp.PID)
        # Randomly reconstrunct the HSCP as a muon
        if np.random.uniform() < eff:
            continue
        muonsLLP.append(hscp)
    
    return muonsLLP

def removeFromMET(particles,METobj):
    """
    Removes the contribution from the particles in the list
    to the total MET.
    """

    metx = METobj.MET*np.cos(METobj.Phi)
    mety = METobj.MET*np.sin(METobj.Phi)

    if particles:
        # Remove particles from MET:            
        pxTot = sum([p.Px for p in particles])
        pyTot = sum([p.Py for p in particles])        
        metx = (metx-pxTot)
        mety = (mety-pyTot)
    
    return [metx,mety]

def applyMTCut(hscps,METvector):
    """
    Remove tracks which have mT < 130 GeV
    """

    selHSCPs = []
    met = np.sqrt(METvector[0]**2 + METvector[1]**2)    
    for hscp in hscps:
        pThscp = [hscp.Px,hscp.Py]
        cosdphi = np.dot(pThscp,METvector)/(hscp.PT*met)
        mT = np.sqrt(2*hscp.PT*met*(1-cosdphi))
        if mT < 130.: continue
        selHSCPs.append(hscp)
    
    return selHSCPs

def getModelDict(inputFiles):

    modelInfoDict = {}
    f = inputFiles[0]
    if not os.path.isfile(f):
        print('File %s not found' %f)
        raise OSError()
    parsDict = {}    
    for banner in glob.glob(os.path.join(os.path.dirname(f),'*banner*txt')):
        with open(banner,'r') as ff:
            slhaData = ff.read().split('<slha>')[1].split('</slha>')[0]
            slhaData = pyslha.readSLHA(slhaData)
    parsDict = {}
    parsDict['mLLP'] = slhaData.blocks['MASS'][1000024]
    parsDict['mLSP'] = slhaData.blocks['MASS'][1000022]
    parsDict['width'] = slhaData.decays[1000024].totalwidth
    if parsDict['width']:
        parsDict['tau_ns'] = (6.582e-25/parsDict['width'])*1e9
    else:
        parsDict['tau_ns'] = np.inf    

    modelInfoDict.update(parsDict)
    print('mLLP = ',parsDict['mLLP'])
    print('width (GeV) = ',parsDict['width'])
    print('tau (ns) = ',parsDict['tau_ns'])

    return modelInfoDict

# ### Define dictionary to store data
def getRecastData(inputFiles,pTcut=60.,normalize=False):

    if len(inputFiles) > 1:
        print('Combining files:')
        for f in inputFiles:
            print(f)

    modelDict = getModelDict(inputFiles)
    if not modelDict:
        modelDict = {}

    # Select mass windows accoding to lifetime
    # (use long lifetime by default)
    if 'width' not in modelDict or (6.582e-25/modelDict['width']) < 1e-9:
        massWindows = massShort
        lifetimeRegime = 'Short Lifetime'
    else:
        massWindows = massLong
        lifetimeRegime = 'Long Lifetime'

    massWindows = massWindows[['Mass_window_Low','Mass_window_High','Target_Mass_GeV']]

    modelDict['Total MC Events'] = 0

    nevtsDict = {}
    # Get total number of events:
    for inputFile in inputFiles:
        f = ROOT.TFile(inputFile,'read')
        tree = f.Get("Delphes")
        nevts = tree.GetEntries()
        modelDict['Total MC Events'] += nevts        
        nevtsDict[inputFile] = nevts
        f.Close()


    lumi = 139.0
    yields = {}
    targetMasses = []
    totalweightPB = 0.0
    # Keep track of yields for each dataset
    cutFlow = { "Total" : 0.0,
                "Trigger" : 0.0,
                "ETmiss > 170 GeV" : 0.0,
                "pT > 50 GeV" : 0.0,
                "Track isolation" : 0.0,
                "pT > 120 GeV" : 0.0,
                "eta < 1.8" : 0.0,
                "mT(track,pTmiss) > 130 GeV" : 0.0,
                "(Acceptance)" : 0.0,
                "(SR-Low - no mass Window)" : 0.0,
                "(SR-High - no mass Window)" : 0.0          
                }


    progressbar = P.ProgressBar(widgets=["Reading %i Events: " %modelDict['Total MC Events'], 
                                P.Percentage(),P.Bar(marker=P.RotatingMarker()), P.ETA()])
    progressbar.maxval = modelDict['Total MC Events']
    progressbar.start()

    ntotal = 0
    for inputFile in inputFiles:
        f = ROOT.TFile(inputFile,'read')
        tree = f.Get("Delphes")
        nevts = tree.GetEntries()
        if normalize:
            norm =nevtsDict[inputFile]/modelDict['Total MC Events']
        else:
            norm = 1.0

        for ievt in range(nevts):    
            
            ntotal += 1
            progressbar.update(ntotal)
            tree.GetEntry(ievt)        

            metCalo = tree.MissingETCalo.At(0).MET
            hscpCandidates = getHSCPCandidates(tree.isoHSCPs,tree.hscpDaughters)
            dmParticles = [tree.dmParticles.At(idm) for idm in range(tree.dmParticles.GetEntries())]

            try:
                weightPB = tree.Weight.At(1).Weight
            except:
                weightPB = tree.Weight.At(0).Weight
            weightPB = weightPB*norm
            ns = weightPB*1e3*lumi # number of signal events
            totalweightPB += weightPB


            # For chargino -> LSP + pion the LSP+LSP pT is ≃ the gaugino system (C1N1 or C1C1) pT
            gauginoPT = sum([p.Px for p in dmParticles])**2
            gauginoPT += sum([p.Py for p in dmParticles])**2
            gauginoPT = np.sqrt(gauginoPT)
            # Apply pT cut on HSCP system to reproduce ATLAS selection
            if gauginoPT < pTcut:
                continue


            hscpsFilter = applyHSCPSelection(hscpCandidates,pT=120.,eta=1.8,r=500.0)
            if hscpsFilter:
                cutFlow['(Acceptance)'] += ns

            muonsLLP = applyMuonTagging(hscpCandidates)
            hscps = [hscp for hscp in hscpCandidates if hscp not in muonsLLP]
            newMETv = removeFromMET(muonsLLP,tree.MissingET.At(0))
            newMET = np.sqrt(newMETv[0]**2+newMETv[1]**2)
            triggerEff = getTriggerEff(metCalo)
            
            cutFlow['Total'] += ns
            ns = ns*triggerEff
            if not ns: continue
            cutFlow['Trigger'] += ns
            # Apply event selection efficiency
            eventEff = getSelectionEff(newMET)
            ns = ns*eventEff
            if not ns: continue
            cutFlow["ETmiss > 170 GeV"] += ns

            # Apply selection to HSCP candidates (following ATLAS snippet)
            hscps = applyHSCPSelection(hscps,pT=50.,eta=3.0,r=500.0)
            if not hscps: continue
            cutFlow['pT > 50 GeV'] += ns        
            # hscps = applyIsolation(hscps,tree.Track) # Already included in trackEff
            if not hscps: continue
            cutFlow['Track isolation'] += ns     
            hscps = applyHSCPSelection(hscps,pT=120.)
            if not hscps: continue
            cutFlow['pT > 120 GeV'] += ns  
            if not hscps: continue
            hscps = applyHSCPSelection(hscps,pT=120.,eta=1.8)
            if not hscps: continue
            cutFlow['eta < 1.8'] += ns
            # hscps = applyMTCut(hscps,newMETv) # Already included in trackEff
            if not hscps: continue
            cutFlow['mT(track,pTmiss) > 130 GeV'] += ns

            gbetas = [h.gbeta for h in hscps]
            trackEffHigh = getTrackEff(gbetas,sr='High')
            trackEffLow =  getTrackEff(gbetas,sr='Low')
            cutFlow['(SR-High - no mass Window)'] += ns*(1-np.prod(1.0-trackEffHigh))
            cutFlow['(SR-Low - no mass Window)'] += ns*(1-np.prod(1.0-trackEffLow))

            # Assume hscp masses can be approximated
            # by the closest target mass for computing the
            # mass window efficiency:
            masses = [getTargetMass(hscp.Mass) for hscp in hscps]
            if all(m is None for m in masses):
                continue
            # masses = [hscp.Mass for hscp in hscps] # Use target mass or real HSCP mass?
            # Not sure if we can use the same mass window efficiencies if the masses
            # are different!
            wmassSRHigh = getMassSelEff(masses,sr='High')
            wmassSRLow = getMassSelEff(masses,sr='Low')
               
            yieldHigh = ns*(1-np.prod(1.0-trackEffHigh*wmassSRHigh))            
            yieldLow = ns*(1-np.prod(1.0-trackEffLow*wmassSRLow))
        
            # In case there are distinct masses, use the largest one (correct?)
            # (if there are no good masses, store in 0.)
            mTarget = max([m for m in masses if m])
            if not mTarget in yields:
                yields[mTarget] = {'SR-Inclusive_Low': [], 'SR-Inclusive_High' : []}
                
            # Store event for a given target window
            yields[mTarget]['SR-Inclusive_Low'].append(yieldLow)
            yields[mTarget]['SR-Inclusive_High'].append(yieldHigh)
            
        f.Close()
    progressbar.finish()

    modelDict['Total xsec-pTcut (pb)'] = cutFlow['Total']/(1e3*lumi)
    # Store total (combined xsec)
    modelDict['Total xsec (pb)'] = totalweightPB
    print('\nCross-section (pb) = %1.3e\n' %totalweightPB)

    # Compute normalized cutflow
    for key,val in cutFlow.items():
        if key == 'Total':
            continue
        valRound = float('%1.3e' %val)
        valNorm = float('%1.3e' %(val/cutFlow['Total']))
        cutFlow[key] = (valRound,valNorm)
    cutFlow['Total'] = (float('%1.3e' %cutFlow['Total']),1.0)

    
    # Create a dictionary for storing data
    dataDict = {'Luminosity (1/fb)' : lumi, 'Regime' : lifetimeRegime}
    # Signal regions
    if not yields:
        dataDict['SR'] = ['SR-Inclusive_Low', 'SR-Inclusive_High']
        dataDict['Target Mass [GeV]'] = [0.0,0.0]
        dataDict['$N_s$'] = [0.0,0.0]
        dataDict['$\sigma_{Ns}$'] = [0.0,0.0]
    else:
        dataDict['SR'] = []
        dataDict['Target Mass [GeV]'] = []
        dataDict['$N_s$'] = []
        dataDict['$\sigma_{Ns}$'] = []

        # Total signal regions = mass targets * (Low, High)
        for targetMass in yields.keys():
            for sr in yields[targetMass]:
                ns = np.array(yields[targetMass][sr])
                nsTot = sum(ns)
                nsError = np.sqrt(sum(ns**2))
                if not nsTot:
                    continue # Skip empty bins
                dataDict['SR'].append(sr)
                dataDict['Target Mass [GeV]'].append(targetMass)
                dataDict['$N_s$'].append(nsTot)
                dataDict['$\sigma_{Ns}$'].append(nsError)    

    # Expand cutflow and modelDict to match number of rows in dataDict:
    for key,val in cutFlow.items():
        cutFlow[key] = [val]*len(dataDict['SR'])
    for key,val in modelDict.items():
        modelDict[key] = [val]*len(dataDict['SR'])

    # Create a dictionary for storing data
    dataDict.update(modelDict)
    dataDict.update(cutFlow)
   

    return dataDict


if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser( description=
            "Run the recasting for ATLAS-SUSY-2018-42 using one or multiple Delphes ROOT files as input. "
            + "If multiple files are given as argument, add them (the samples weights will be normalized if -n is given)."
            + " Store the cutflow and SR bins in a pickle (Pandas DataFrame) file." )
    ap.add_argument('-f', '--inputFile', required=True,nargs='+',
            help='path to the ROOT event file(s) generated by Delphes.', default =[])
    ap.add_argument('-o', '--outputFile', required=False,
            help='path to output file storing the DataFrame with the recasting data.'
                 + 'If not defined, will use the name of the first input file', 
            default = None)
    ap.add_argument('-pt', '--pTcut', required=False,default=60.0,type=float,
            help='Gen level cut on gaugino pT for computing partial cross-sections.')
    ap.add_argument('-n', '--normalize', required=False,action='store_true',
            help='If set, the input files will be considered to refer to multiple samples of the same process and their weights will be normalized.')
    ap.add_argument('-v', '--verbose', default='info',
            help='verbose level (debug, info, warning or error). Default is info')


    # First make sure the correct env variables have been set:
    import subprocess
    import sys
    LDPATH = subprocess.check_output('echo $LD_LIBRARY_PATH',shell=True,text=True)
    ROOTINC = subprocess.check_output('echo $ROOT_INCLUDE_PATH',shell=True,text=True)
    pythiaDir = os.path.abspath('./MG5/HEPTools/pythia8/lib')
    delphesDir = os.path.abspath('./DelphesHSCP/external')
    if pythiaDir not in LDPATH or delphesDir not in ROOTINC:
        print('Enviroment variables not properly set. Run source setenv.sh first.')
        sys.exit()


    t0 = time.time()

    # # Set output file
    args = ap.parse_args()
    inputFiles = args.inputFile
    outputFile = args.outputFile
    if outputFile is None:
        outputFile = inputFiles[0].replace('delphes_events.root','atlas_2018_42.pcl')

    if os.path.splitext(outputFile)[1] != '.pcl':
        outputFile = os.path.splitext(outputFile)[0] + '.pcl'

    modelDict = getRecastData(inputFiles,args.pTcut,args.normalize)

    # #### Create pandas DataFrame
    df = pd.DataFrame.from_dict(modelDict)

    # ### Save DataFrame to pickle file
    print('Saving to',outputFile)
    df.to_pickle(outputFile)

    print("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
