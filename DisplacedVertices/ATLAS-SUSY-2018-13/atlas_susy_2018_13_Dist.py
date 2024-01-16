#!/usr/bin/env python3

import os,glob
from typing import Any
import numpy as np
import pandas as pd
import glob
import pyslha
import time
import progressbar as P
from helper import LLP,BinnedData
from ATLAS_data.effFunctions import eventEff,vertexEff
from ATLAS_data.atlasBins import atlas_bins
from atlas_susy_2018_13_Recast import getLLPs,getJets,getDisplacedJets,eventAcc,vertexAcc,getModelDict


delphesDir = os.path.abspath("./DelphesLLP")
os.environ['ROOT_INCLUDE_PATH'] = os.path.join(delphesDir,"external")

import ROOT
import xml.etree.ElementTree as ET


ROOT.gSystem.Load(os.path.join(delphesDir,"libDelphes.so"))

ROOT.gInterpreter.Declare('#include "classes/SortableObject.h"')
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')


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

    
    binnedData = {}
    for sr in ['HighPT', 'Trackless']:
        xbins = atlas_bins['HighPT']['nTracks']
        ybins = atlas_bins['HighPT']['mDV']
        binnedData[sr] = BinnedData(xbins,ybins)

    lumi = 139.0
    totalweightPB = 0.0
    # Keep track of yields for each dataset
    cutFlowHighPT = { "Total" : 0.0,
                "Jet selection" : 0.0,
                # "$R_{xy},z <$ 300 mm" : 0.0,
                # "$R_{DV} > 4$ mm" : 0.0,
                # "$nTracks >= 5$" : 0.0,
                # "mDV > 10 GeV" : 0.0
                "DV selection" : 0.0
                }
    
    cutFlowTrackless = {k : v for k,v in cutFlowHighPT.items()}
    cutFlows = {'HighPT' : cutFlowHighPT, 'Trackless' : cutFlowTrackless}


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
        if normalize:
            norm =nevtsDict[inputFile]/modelDict['Total MC Events']
        else:
            norm = 1.0/modelDict['Total MC Events']

        for ievt in range(nevts):    
            
            ntotal += 1
            progressbar.update(ntotal)
            tree.GetEntry(ievt)   
            weightPB = tree.Weight.At(0).Weight     
            weightPB = weightPB*norm
            totalweightPB += weightPB
            ns = weightPB*1e3*lumi # number of signal events

            llps = getLLPs(tree.llps,tree.llpDaughters)
            jets = getJets(tree.GenJet,pTmin=25.,etaMax=5.0)
            jetsDisp = getDisplacedJets(jets,llps)
            
            # Event acceptance:
            ev_acc = {}
            ev_acc['HighPT'] = eventAcc(jets,jetsDisp,sr='HighPT')
            ev_acc['Trackless'] = eventAcc(jets,jetsDisp,sr='Trackless')                        


            cutFlowHighPT["Total"] += ns
            cutFlowTrackless["Total"] += ns
            if sum(ev_acc.values()) == 0.::
                continue

            cutFlowHighPT["Jet selection"] += ns*ev_acc['HighPT']
            cutFlowTrackless["Jet selection"] += ns*ev_acc['Trackless']
            
            # Event efficiency
            ev_eff = {}
            ev_eff['HighPT'] = eventEff(jets,llps,sr='HighPT')
            ev_eff['Trackless'] = eventEff(jets,llps,sr='Trackless')

            # Vertex efficiencies:
            v_eff = np.array([vertexEff(llp) for llp in llps])            

            # Vertex acceptances:
            nTrack_mDV_pairs = np.array[(llp.nTracks,llp.mDV) for llp in llps]
            for sr in ['HighPT','Trackless']:
                binData = binnedData[sr]
                for nmin in binData.xbins:
                    if max(nTrack_mDV_pairs[:,0]) < nmin:
                        continue
                    for mDVmin in binData.ybinx:
                        if max(nTrack_mDV_pairs[:,1]) < mDVmin:
                            continue
                        
                        v_acc = np.array([vertexAcc(llp,Rmax=300.0,zmax=300.0,Rmin=4.0,
                                        d0min=2.0,nmin=nmin,mDVmin=mDVmin)  for llp in llps])
                
                        wvertex = 1.0-np.prod(1.0-v_acc*v_eff)
                    binData.fill(nmin,mDVmin,ns*ev_acc[sr]*ev_eff[sr]*wvertex)

            
                    if nmin == 5 and mDVmin == 10.0:
                        # Add to the total weight in each SR:
                        cutFlow[sr]["DV selection"] += ns*ev_acc[sr]*ev_eff[sr]*wvertex

        f.Close()
    progressbar.finish()

    modelDict['Total xsec (pb)'] = totalweightPB
    print('\nCross-section (pb) = %1.3e\n' %totalweightPB)

    # Compute normalized cutflow
    for cutFlow in [cutFlowHighPT,cutFlowTrackless]:
        for key,val in cutFlow.items():
            if key == 'Total':
                continue
            valRound = float('%1.3e' %val)
            valNorm = float('%1.3e' %(val/cutFlow['Total']))
            cutFlow[key] = (valRound,valNorm)
        cutFlow['Total'] = (float('%1.3e' %cutFlow['Total']),1.0)

    
    # Create a dictionary for storing data
    dataDict = {}
    dataDict['Luminosity (1/fb)'] = []
    dataDict['SR'] = []
    dataDict['$N_s$'] = []
    # Signal regions
    cutFlows = {'HighPT' : cutFlowHighPT, 'Trackless' : cutFlowTrackless}
    for sr,cutFlow in cutFlows.items():
        dataDict['Luminosity (1/fb)'].append(lumi)
        dataDict['SR'].append(sr)
        dataDict['$N_s$'].append(cutFlow["DV selection"][0])
        for cut in cutFlow:
            if cut not in dataDict:
                dataDict[cut] = []
            dataDict[cut].append(cutFlow[cut])

    # Expand modelDict to match number of rows in dataDict:
    for key,val in modelDict.items():
        modelDict[key] = [val]*len(dataDict['SR'])

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
    ap.add_argument('-m', '--model', required=False,type=str,default='strong',
            help='Defines which model should be considered for extracting model parameters (strong,ewk,gluino).')

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
        outputFile = inputFiles[0].replace('delphes_events.root','atlas_2018_13.pcl')

    if os.path.splitext(outputFile)[1] != '.pcl':
        outputFile = os.path.splitext(outputFile)[0] + '.pcl'

    dataDict = getRecastData(inputFiles,args.normalize,args.model)
    if args.verbose == 'debug':
        for k,v in dataDict.items():
            print(k,v)

    # #### Create pandas DataFrame
    df = pd.DataFrame.from_dict(dataDict)

    # ### Save DataFrame to pickle file
    print('Saving to',outputFile)
    df.to_pickle(outputFile)

    print("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
