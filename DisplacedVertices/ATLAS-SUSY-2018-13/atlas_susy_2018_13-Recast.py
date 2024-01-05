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

delphesDir = os.path.abspath("./DelphesLLP")
os.environ['ROOT_INCLUDE_PATH'] = os.path.join(delphesDir,"external")

import ROOT
import xml.etree.ElementTree as ET


ROOT.gSystem.Load(os.path.join(delphesDir,"libDelphes.so"))

ROOT.gInterpreter.Declare('#include "classes/SortableObject.h"')
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')


def getLLPCandidates(llps,daughters):

    candidates = []
    for ip in range(llps.GetEntries()):
        p = llps.At(ip)
        if p.Status == 1:
            p.r_decay = np.inf
            p.z_decay = np.inf
            candidates.append(p)
            continue
        
        # If LLP is unstable compute its decay position from daughter
        d1 = p.D1
        daughter = daughters.At(d1)
        # Get the LLP decay radius from the first daughter production vertex:
        p.r_decay = np.sqrt(daughter.X**2+daughter.Y**2)      
        p.z_decay = daughter.Z        
        candidates.append(p)

    for p in candidates:
        # Add beta and gamma*beta
        trimom = np.sqrt(p.Px**2 +p.Py**2 + p.Pz**2)
        p.beta = trimom/p.E
        p.gbeta = trimom/p.Mass

    return candidates

def applyLLPSelection(llpList,pT=50.,eta=2.4,r=500.):

    selLLPs = []
    for llp in llpList:
        if llp.PT < pT: continue 
        if abs(llp.Eta) > eta: continue
        if llp.r_decay < r: continue
        selLLPs.append(llp)
    
    return selLLPs

def applyIsolation(llpList,pTmax=5.0):

    isoLLPs = []
    # Apply isolation requirement for LLP tracks
    for llp in llpList:
        sumPT = llp.SumPtCharged
        if sumPT > pTmax: continue
        isoLLPs.append(llp)
    return isoLLPs


def applyMuonTagging(llpList):

    """
    Computes the probability of reconstructing the llp as a muon.
    :param llpList: List of GenParticle objects
    """

    muonsLLP = []
    for llp in llpList:
        if llp.r_decay < 3.9e3 and llp.z_decay < 6e3: # Skip decays before MS
            continue
        
        beta = llp.beta
        eta = abs(llp.Eta)
        eff = getMuonRecoEff(beta,eta,llp.PID)
        # Randomly reconstrunct the LLP as a muon
        if np.random.uniform() < eff:
            continue
        muonsLLP.append(llp)
    
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


def applyMTCut(llps,METvector):
    """
    Remove tracks which have mT < 130 GeV
    """

    selLLPs = []
    met = np.sqrt(METvector[0]**2 + METvector[1]**2)    
    for llp in llps:
        pTllp = [llp.Px,llp.Py]
        cosdphi = np.dot(pTllp,METvector)/(llp.PT*met)
        mT = np.sqrt(2*llp.PT*met*(1-cosdphi))
        if mT < 130.: continue
        selLLPs.append(llp)
    
    return selLLPs

def getModelDict(inputFiles,model):

    if model == 'wino':
        LLP = 1000024
        LSP = 1000022
    elif model == 'stau':
        LLP = 1000015
        LSP = 1000039
    elif model == 'gluino':
        LLP = 1000021
        LSP = 1000022

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
    parsDict['mLLP'] = slhaData.blocks['MASS'][LLP]
    parsDict['mLSP'] = slhaData.blocks['MASS'][LSP]
    parsDict['width'] = slhaData.decays[LLP].totalwidth
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
def getRecastData(inputFiles,pTcut=60.,normalize=False,model='wino'):

    if len(inputFiles) > 1:
        print('Combining files:')
        for f in inputFiles:
            print(f)

    modelDict = getModelDict(inputFiles,model)
    if not modelDict:
        modelDict = {}

    # Select mass windows accoding to lifetime
    # (use long lifetime by default)
    if 'width' not in modelDict or (modelDict['width'] > 0 and (6.582e-25/modelDict['width']) < 1e-9):
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
    totalweightPB = 0.0
    # Keep track of yields for each dataset
    cutFlow = { "Total" : 0.0,
                "Trigger" : 0.0,
                "$E_{T}^{miss}>170$ GeV" : 0.0,
                "$p_{T} > 50$ GeV" : 0.0,
                "Track isolation" : 0.0,
                "$p_{T} > 120$ GeV" : 0.0,
                "$|\eta|<1.8$" : 0.0,
                "$m_{T}({track},{p}_{{T}}^{{ miss}}) > 130$ GeV" : 0.0,
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
            llpCandidates = getLLPCandidates(tree.llps,tree.llpDaughters)
            dmParticles = [tree.dmParticles.At(idm) for idm in range(tree.dmParticles.GetEntries())]

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
            help='Gen level cut on gaugino pT for computing partial cross-sections [default = 60].')
    ap.add_argument('-n', '--normalize', required=False,action='store_true',
            help='If set, the input files will be considered to refer to multiple samples of the same process and their weights will be normalized.')
    ap.add_argument('-m', '--model', required=False,type=str,default='wino',
            help='Defines which model should be considered for extracting model parameters (stau,wino,gluino).')

    ap.add_argument('-v', '--verbose', default='info',
            help='verbose level (debug, info, warning or error). Default is info')


    # First make sure the correct env variables have been set:
    import subprocess
    import sys
    LDPATH = subprocess.check_output('echo $LD_LIBRARY_PATH',shell=True,text=True)
    ROOTINC = subprocess.check_output('echo $ROOT_INCLUDE_PATH',shell=True,text=True)
    pythiaDir = os.path.abspath('./MG5/HEPTools/pythia8/lib')
    delphesDir = os.path.abspath('./DelphesLLP/external')
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

    modelDict = getRecastData(inputFiles,args.pTcut,args.normalize,args.model)

    # #### Create pandas DataFrame
    df = pd.DataFrame.from_dict(modelDict)

    # ### Save DataFrame to pickle file
    print('Saving to',outputFile)
    df.to_pickle(outputFile)

    print("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
