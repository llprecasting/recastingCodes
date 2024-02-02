#!/usr/bin/env python3

import os,glob
from typing import Any
import numpy as np
import pandas as pd
import glob
import pyslha
import time
import progressbar as P
from helper import LLP
from ATLAS_data.effFunctions import eventEff,vertexEff

delphesDir = os.path.abspath("../DelphesLLP")
os.environ['ROOT_INCLUDE_PATH'] = os.path.join(delphesDir,"external")

import ROOT
import xml.etree.ElementTree as ET


ROOT.gSystem.Load(os.path.join(delphesDir,"libDelphes.so"))

ROOT.gInterpreter.Declare('#include "classes/SortableObject.h"')
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')

def getLLPs(llpList,directDaughters,finalDaughters,maxMomViolation=1e-2):

    llps = []
    for ip in range(llpList.GetEntries()):
        p = llpList.At(ip)        
        # Get direct daughters
        llp_ddaughters = []        
        for d in directDaughters:
            if d.M1 == ip:
                llp_ddaughters.append(d)
        # Get final daughters
        llp_fdaughters = []
        for d in finalDaughters:
            if d.M1 == ip:
                llp_fdaughters.append(d)
        
        # Get final daughters
        llps.append(LLP(p,llp_ddaughters,llp_fdaughters,maxMomViolation))
        
    return llps

def getJets(jets,pTmin=50.0,etaMax=5.0):
    """
    Select jets with pT > pTmin and |eta| < etaMax
    """

    jetsSel = []
    for ijet in range(jets.GetEntries()):
        jet = jets.At(ijet)
        if jet.PT < pTmin:
            continue
        if abs(jet.Eta) > etaMax:
            continue
        jetsSel.append(jet)
    
    return jetsSel

def getDisplacedJets(jets,llps,skipPIDs=[1000022]):
    """
    Select from the list of all jets, the displaced jets associated
    with a LLP decay.
    """

    displacedJets = []
    for jet in jets:
        deltaRmin = 0.3
        llpMatch = None
        for llp in llps: 
            for daughter in llp.directDaughters:
                if abs(daughter.PID) in skipPIDs:
                    continue
                deltaR = np.sqrt((jet.Eta-daughter.Eta)**2 + (jet.Phi-daughter.Phi)**2)
                if deltaR < deltaRmin:
                    deltaRmin = deltaR
                    llpMatch = llp # Store LLP parent
        
        jet.llp = llpMatch
        if llpMatch is not None:
            R = np.sqrt(llpMatch.Xd**2 + llpMatch.Yd**2 + llpMatch.Zd**2)
            if R > 3870:
                continue
            displacedJets.append(jet)
    
    return displacedJets

def eventAcc(jets,met,metCut=200.0,
             maxJetChargedPT=np.inf,
             minJetPt1=0.,minJetPt2=0.,
             minPVdistance=0.0):
    
    if met < metCut:
        return 0.0

    # Split analysis is two bunchs: 75% and 25%
    lumCut = np.random.random()
    if lumCut > 0.75:
        return 1.0
    
    passAcc = 0.0
    # Apply jet cuts
    good_jets = []
    for jet in jets:
        pTCharged = 0.0
        for ip in range(jet.Constituents.GetEntries()):
            particle = jet.Constituents.At(ip)
            if particle.Charge == 0:
                continue
            r_prod = np.sqrt(particle.X**2 + particle.Y**2)
            if r_prod  > minPVdistance:
                continue
            pTCharged += particle.PT
        if pTCharged > maxJetChargedPT:
            continue
        good_jets.append(jet)

    
    good_jets = sorted(good_jets[:], key = lambda j: j.PT, reverse=True)

    if len(good_jets) > 0 and good_jets[0].PT > minJetPt1:
        passAcc = 1.0
    elif len(good_jets) > 1 and  good_jets[1].PT > minJetPt2:
        passAcc = 1.0
    
    return passAcc

def vertexAcc(llp,Rmax=np.inf,zmax=np.inf,Rmin=0.0,d0min=0.0,nmin=0,mDVmin=0):
    
    passAcc = 1.0
    
    if np.sqrt(llp.Xd**2 + llp.Yd**2) > Rmax or abs(llp.Zd) > zmax:
        passAcc = 0.0
        
    if np.sqrt(llp.Xd**2 + llp.Yd**2) < Rmin:
        passAcc = 0.0
    
    maxD0 = 0.0
    for d in llp.finalDaughters:
        if d.Charge == 0: # Skip neutral
            continue
        d0 = abs((d.Y*d.Px - d.X*d.Px)/d.PT)
        maxD0 = max(maxD0,d0)
    if maxD0 < d0min:
        passAcc = 0.0
            
    if llp.nTracks < nmin:
        passAcc = 0.0
        
    if llp.mDV < mDVmin:
        passAcc = 0.0
        
    return passAcc
    
def getModelDict(inputFiles,model):

    if model == 'ewk':
        LLP = 1000022
        LSP = 1000024
    elif model == 'strong':
        LLP = 1000022
        LSP = 1000021
    elif model == 'gluino':
        LLP = 1000021
        LSP = 1000022
    elif model == 'sbottom':
        LLP = 1000005
        LSP = 1000022        
    else:
        raise ValueError("Unreconized model %s" %model)

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
def getRecastData(inputFiles,normalize=False,model='strong'):

    if len(inputFiles) > 1:
        print('Combining files:')
        for f in inputFiles:
            print(f)

    modelDict = getModelDict(inputFiles,model)
    if not modelDict:
        modelDict = {}

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


    lumi = 32.8
    totalweightPB = 0.0
    # Keep track of yields for each dataset
    cutFlow = { "Total" : 0.0,
                "Jet+MET selection" : 0.0,
                # "$R_{xy},z <$ 300 mm" : 0.0,
                # "$R_{DV} > 4$ mm" : 0.0,
                # "$nTracks >= 5$" : 0.0,
                # "mDV > 10 GeV" : 0.0
                "DV selection" : 0.0
                }
    

    progressbar = P.ProgressBar(widgets=["Reading %i Events: " %modelDict['Total MC Events'], 
                                P.Percentage(),P.Bar(marker=P.RotatingMarker()), P.ETA()])
    progressbar.maxval = modelDict['Total MC Events']
    progressbar.start()

    ntotal = 0
    totalweightPB = 0.0
    for inputFile in inputFiles:
        f = ROOT.TFile(inputFile,'read')
        tree = f.Get("Delphes")
        nevts = tree.GetEntries()
        # if normalize:
        #     norm =nevtsDict[inputFile]/modelDict['Total MC Events']
        # else:
        #     norm = 1.0/modelDict['Total MC Events']
        norm = 1.0

        for ievt in range(nevts):    
            
            ntotal += 1
            progressbar.update(ntotal)
            tree.GetEntry(ievt)   
            weightPB = tree.Weight.At(1).Weight     
            weightPB = weightPB*norm
            totalweightPB += weightPB
            ns = weightPB*1e3*lumi # number of signal events

            jets = getJets(tree.GenJet,pTmin=25.,etaMax=5.0)
            met = tree.GenMissingET.At(0).MET

            # Event acceptance
            evt_acc = eventAcc(jets,met,metCut=200.0,
                               maxJetChargedPT=5.0,minJetPt1=70.,
                               minJetPt2=25.,minPVdistance=4.0)

            cutFlow["Total"] += ns
            if (not evt_acc):
                continue

            cutFlow["Jet+MET selection"] += ns*evt_acc

            llps = getLLPs(tree.bsm,tree.bsmDirectDaughters,tree.bsmFinalDaughters)
            # Vertex acceptances:
            v_acc = np.array([vertexAcc(llp,Rmax=300.0,zmax=300.0,Rmin=4.0,
                                        d0min=2.0,nmin=5,mDVmin=10.0)  for llp in llps])
            good_llps = np.array(llps)[v_acc > 0.0]
            if len(good_llps) == 0:
                continue
            # Event efficiency
            evt_eff = eventEff(met,good_llps)
            
            # Vertex efficiencies:
            v_eff = np.array([vertexEff(llp) for llp in llps])
            
            wvertex = 1.0-np.prod(1.0-v_acc*v_eff)
            
            # Add to the total weight in each SR:
            cutFlow["DV selection"] += ns*evt_acc*evt_eff*wvertex

        f.Close()
    progressbar.finish()

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
    dataDict = {}
    dataDict['Luminosity (1/fb)'] = [lumi]
    dataDict['$N_s$'] = [cutFlow["DV selection"][0]]
    dataDict['AccEff'] = [cutFlow["DV selection"][1]]
    for cut,val in cutFlow.items():
        dataDict.setdefault(cut,[val])

    # Create a dictionary for storing data
    dataDict.update(modelDict)
   

    return dataDict


if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser( description=
            "Run the recasting for ATLAS-SUSY-2018-13 using one or multiple Delphes ROOT files as input. "
            + "If multiple files are given as argument, add them (the samples weights will be normalized if -n is given)."
            + " Store the cutflow and SR bins in a pickle (Pandas DataFrame) file." )
    ap.add_argument('-f', '--inputFile', required=True,nargs='+',
            help='path to the ROOT event file(s) generated by Delphes.', default =[])
    ap.add_argument('-o', '--outputFile', required=False,
            help='path to output file storing the DataFrame with the recasting data. '
                 + 'If not defined, will use the name of the first input file', 
            default = None)
    ap.add_argument('-n', '--normalize', required=False,action='store_true',
            help='If set, the input files will be considered to refer to multiple samples of the same process and their weights will be normalized.')
    ap.add_argument('-m', '--model', required=False,type=str,default='sbottom',
            help='Defines which model should be considered for extracting model parameters (strong,ewk,gluino,sbottom).')

    ap.add_argument('-v', '--verbose', default='info',
            help='verbose level (debug, info, warning or error). Default is info')


    # First make sure the correct env variables have been set:
    import subprocess
    import sys
    LDPATH = subprocess.check_output('echo $LD_LIBRARY_PATH',shell=True,text=True)
    ROOTINC = subprocess.check_output('echo $ROOT_INCLUDE_PATH',shell=True,text=True)
    pythiaDir = os.path.abspath('../MG5/HEPTools/pythia8/lib')
    delphesDir = os.path.abspath('../DelphesLLP/external')
    if pythiaDir not in LDPATH or delphesDir not in ROOTINC:
        print('Enviroment variables not properly set. Run source setenv.sh first.')
        sys.exit()


    t0 = time.time()

    # # Set output file
    args = ap.parse_args()
    inputFiles = args.inputFile
    outputFile = args.outputFile
    if outputFile is None:
        outputFile = inputFiles[0].replace('delphes_events.root','atlas_2016_08.pcl')

    if os.path.splitext(outputFile)[1] != '.pcl':
        outputFile = os.path.splitext(outputFile)[0] + '.pcl'

    dataDict = getRecastData(inputFiles,args.normalize,args.model)
    if args.verbose in ['info','debug']:
        for k,v in dataDict.items():
            print(k,v)

    # #### Create pandas DataFrame
    df = pd.DataFrame.from_dict(dataDict)

    # ### Save DataFrame to pickle file
    print('Saving to',outputFile)
    df.to_pickle(outputFile)

    print("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
