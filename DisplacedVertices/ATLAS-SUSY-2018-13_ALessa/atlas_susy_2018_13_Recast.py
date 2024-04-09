#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import time
import progressbar as P
import sys
sys.path.append('../')
from helper import getLLPs,getJets,getDisplacedJets,getModelDict,splitModels
from ATLAS_data.effFunctions import eventEff,vertexEff

delphesDir = os.path.abspath("./DelphesLLP")
os.environ['ROOT_INCLUDE_PATH'] = os.path.join(delphesDir,"external")

import ROOT

ROOT.gSystem.Load(os.path.join(delphesDir,"libDelphes.so"))

ROOT.gInterpreter.Declare('#include "classes/SortableObject.h"')
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')



def eventAcc(jets,jetsDisp,sr):
    passAcc = 0.0
    if sr == 'HighPT':
        # Apply HighPT jet selection    
        njet250 = len([j for j in jets if j.PT > 250.0])
        njet195 = len([j for j in jets if j.PT > 195.0])
        njet116 = len([j for j in jets if j.PT > 116.0])
        njet90 = len([j for j in jets if j.PT > 90.0])
        if (njet250 >= 4) or (njet195 >= 5) or (njet116 >= 6) or (njet90 >= 7):
            passAcc = 1.0
    elif sr == 'Trackless':    
        # Apply Trackless jet selection (only if HighPT has failed)    
        njet137 = len([j for j in jets if j.PT > 137.0])
        njet101 = len([j for j in jets if j.PT > 101.0])
        njet83 = len([j for j in jets if j.PT > 83.0])
        njet55 = len([j for j in jets if j.PT > 55.0])
        njetDisp70 = len([j for j in jetsDisp if j.PT > 70.0])
        njetDisp50 = len([j for j in jetsDisp if j.PT > 50.0])
        if (njet137 >= 4) or (njet101 >= 5) or (njet83 >= 6) or (njet55 >= 7):
            if (njetDisp70 >=1) or (njetDisp50 >= 2):
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
    
def getRecastData(inputFiles,model='strong',modelDict=None,addweights=False):

    if len(inputFiles) > 1:
        print('Combining files:')
        for f in inputFiles:
            print(f)

    if modelDict is None:
        modelDict = getModelDict(inputFiles[0],model)
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


    lumi = 139.0
    totalweightPB = 0.0
    # Keep track of yields for each dataset
    keys = ["Total", "Jet selection", "DV selection"]
    cutFlowHighPT = { k : np.zeros(2) for k in keys}    
    cutFlowTrackless = {k : np.zeros(2) for k in keys}

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
        # If addweights = Fakse: 
        # assume multiple files correspond to equivalent samplings
        # of the same distributions
        # If addweights = True: directly add events
        if not addweights:
            norm =nevtsDict[inputFile]/modelDict['Total MC Events']
        else:
            norm = 1.0

        for ievt in range(nevts):    
            
            ntotal += 1
            progressbar.update(ntotal)
            tree.GetEntry(ievt)   
            weightPB = tree.Event.At(0).Weight/nevts
            weightPB = weightPB*norm
            totalweightPB += weightPB
            ns = weightPB*1e3*lumi # number of signal events

            llps = getLLPs(tree.bsm,tree.bsmDirectDaughters,tree.bsmFinalDaughters,maxMomViolation=6e-2)
            jets = getJets(tree.GenJet,pTmin=25.,etaMax=5.0)
            jetsDisp = getDisplacedJets(jets,llps)
            
            # Event acceptance:
            highPT_acc = eventAcc(jets,jetsDisp,sr='HighPT')
            trackless_acc = eventAcc(jets,jetsDisp,sr='Trackless')                        


            cutFlowHighPT["Total"] += (ns,ns**2)
            cutFlowTrackless["Total"] += (ns,ns**2)

            if (not highPT_acc) and (not trackless_acc):
                continue

            cutFlowHighPT["Jet selection"] += (ns*highPT_acc,(ns*highPT_acc)**2)
            cutFlowTrackless["Jet selection"] += (ns*trackless_acc,(ns*trackless_acc)**2)

            if len(llps) == 0:
                continue
            
            # Event efficiency
            highPT_eff = eventEff(jets,llps,sr='HighPT')
            trackless_eff = eventEff(jets,llps,sr='Trackless')

            # Vertex acceptances:
            v_acc = np.array([vertexAcc(llp,Rmax=300.0,zmax=300.0,Rmin=4.0,
                                        d0min=2.0,nmin=5,mDVmin=10.0)  for llp in llps])

            
            # Vertex efficiencies:
            v_eff = np.array([vertexEff(llp) for llp in llps])

            
            wvertex = 1.0-np.prod(1.0-v_acc*v_eff)
            
            # Add to the total weight in each SR:
            cutFlowHighPT["DV selection"] += (ns*highPT_acc*highPT_eff*wvertex,(ns*highPT_acc*highPT_eff*wvertex)**2)
            cutFlowTrackless["DV selection"] += (ns*trackless_acc*trackless_eff*wvertex,(ns*trackless_acc*trackless_eff*wvertex)**2)

        f.Close()
    progressbar.finish()

    modelDict['Total xsec (pb)'] = totalweightPB
    print('\nCross-section (pb) = %1.3e\n' %totalweightPB)


    cutFlowDicts = {'HighPT' : cutFlowHighPT, 'Trackless' : cutFlowTrackless}
    cutFlowErrDict = {}
    for sr, cutFlow in cutFlowDicts.items():
        cutFlowErrDict[sr] = {k : np.sqrt(v[1]) for k,v in cutFlow.items()}
        cutFlowDicts[sr] = {k : v[0] for k,v in cutFlow.items()}

    # Compute normalized cutflow
    for sr,cutFlow in cutFlowDicts.items():
        cutFlowErr = cutFlowErrDict[sr]
        for key,val in cutFlow.items():
            if key == 'Total':
                continue
            valRound = float('%1.3e' %val)
            valNorm = float('%1.3e' %(val/cutFlow['Total']))
            cutFlow[key] = (valRound,valNorm)
            errRound = float('%1.3e' %cutFlowErr[key])
            errNorm = float('%1.3e' %(cutFlowErr[key]/cutFlow['Total']))
            cutFlowErr[key] = (errRound,errNorm)
        cutFlow['Total'] = (float('%1.3e' %cutFlow['Total']),1.0)
        cutFlowErr['Total'] = (float('%1.3e' %cutFlowErr['Total']),0.0)

    
    # Create a dictionary for storing data
    dataDict = {}
    dataDict['Luminosity (1/fb)'] = []
    dataDict['SR'] = []
    dataDict['$N_s$'] = []
    dataDict['$N_s$ Err'] = []
    dataDict['AccEff'] = []
    dataDict['AccEffErr'] = []
    for sr,cutFlow in cutFlowDicts.items():
        cutFlowErr = cutFlowErrDict[sr]
        dataDict['Luminosity (1/fb)'].append(lumi)
        dataDict['SR'].append(sr)
        dataDict['$N_s$'].append(cutFlow["DV selection"][0])
        dataDict['$N_s$ Err'].append(cutFlowErr["DV selection"][0])
        dataDict['AccEff'].append(cutFlow["DV selection"][1])
        dataDict['AccEffErr'].append(cutFlowErr["DV selection"][1])
        for cut in cutFlow:
            if cut not in dataDict:
                dataDict[cut] = []
                dataDict[cut+' Error'] = []
            dataDict[cut].append(cutFlow[cut])
            dataDict[cut+' Error'].append(cutFlowErr[cut])

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
    ap.add_argument('-A', '--add', required=False,action='store_true',default=False,
            help='If set, the input files will be considered to refer to samples of the orthogonal processes and their weights will be added.')    
    ap.add_argument('-m', '--model', required=False,type=str,default='sbottom',
            help='Defines which model should be considered for extracting model parameters (strong,ewk,gluino,sbottom).')
    ap.add_argument('-U', '--update', required=False,action='store_true',
            help='If the flag is set only the model points containing data newer than the dataframe will be read.')
    
    ap.add_argument('-v', '--verbose', default='info',
            help='verbose level (debug, info, warning or error). Default is info')


    # First make sure the correct env variables have been set:
    import subprocess
    import sys
    from datetime import datetime as dt
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
    
    # Split input files by distinct models and get recast data for
    # the set of files from the same model:
    for fileList,mDict in splitModels(inputFiles,args.model):

        if outputFile is None:
            outFile = fileList[0].replace('delphes_events.root','atlas_2018_13.pcl')
        else:
            outFile = outputFile[:]

        if os.path.splitext(outFile)[1] != '.pcl':
            outFile = os.path.splitext(outFile)[0] + '.pcl'

        skipModel = False
        if args.update and os.path.isfile(outFile):
            outFile_date = dt.fromtimestamp(os.path.getctime(outFile))
            inputFiles_date = max([dt.fromtimestamp(os.path.getctime(f)) for f in fileList])
            if inputFiles_date <= outFile_date:
                skipModel = True
        if skipModel:
            print('\nSkipping',mDict,'\n')
            # print('files=',fileList)
            # sys.exit()
            continue

        print('----------------------------------')
        print('\t Model: %s (%i files)' %(mDict,len(fileList)))

        dataDict = getRecastData(fileList,args.model,mDict,addweights=args.add)
        if args.verbose == 'debug':
            for k,v in dataDict.items():
                print(k,v)

        

        # #### Create pandas DataFrame
        df = pd.DataFrame.from_dict(dataDict)
        # ### Save DataFrame to pickle file
        print('Saving to',outFile)
        df.to_pickle(outFile)
        print('\n')

    print("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
