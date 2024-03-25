#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import time
import progressbar as P
import sys
sys.path.append('../')
from helper import getLLPs,getJets,getModelDict,splitModels
from ATLAS_data.effFunctions import eventEff,vertexEff

delphesDir = os.path.abspath("./DelphesLLP")
os.environ['ROOT_INCLUDE_PATH'] = os.path.join(delphesDir,"external")

import ROOT


ROOT.gSystem.Load(os.path.join(delphesDir,"libDelphes.so"))

ROOT.gInterpreter.Declare('#include "classes/SortableObject.h"')
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')



def eventAcc(jets,met,metCut=200.0,
             maxJetChargedPT=np.inf,
             minJetPt1=0.,minJetPt2=0):
    
    if met < metCut:
        return 0.0

    # Split analysis is two bunchs: 75% and 25%
    lumCut = np.random.random()
    if lumCut > 0.75:
        return 1.0
    
    # Apply jet cuts
    good_jets = []
    jets_sorted = sorted(jets, key = lambda j: j.PT, reverse=True)
    good_jets = [j for j in jets_sorted if j.ChargedPTPV < maxJetChargedPT]
    
    if len(good_jets) > 0 and good_jets[0].PT > minJetPt1:
        return 1.0
    elif len(good_jets) > 1 and  good_jets[1].PT > minJetPt2:
        return 1.0
    else:
        return 0.0

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
    

def getRecastData(inputFiles,model='sbottom',modelDict=None,effStrategy='official',mDVcut=10.0,addweights=False):

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


    lumi = 32.8
    totalweightPB = 0.0
    # Keep track of yields for each dataset
    keys = ["Total","Jet+MET selection","DV selection"]
    cutFlow = {k  : np.zeros(2) for k in keys}    

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
            weightPB = tree.Weight.At(1).Weight     
            weightPB = weightPB*norm
            totalweightPB += weightPB
            ns = weightPB*1e3*lumi # number of signal events

            jets = getJets(tree.GenJet,pTmin=25.,etaMax=5.0)
            met = tree.GenMissingET.At(0).MET

            # Event acceptance
            evt_acc = eventAcc(jets,met,metCut=200.0,
                               maxJetChargedPT=5.0,minJetPt1=70.,
                               minJetPt2=25.)

            cutFlow["Total"] += (ns,ns**2)
            if (not evt_acc):
                continue

            ns = ns*evt_acc
            cutFlow["Jet+MET selection"] += (ns,ns**2)

            llps = getLLPs(tree.bsm,tree.bsmDirectDaughters,tree.bsmFinalDaughters)
            # Vertex acceptances:
            v_acc = np.array([vertexAcc(llp,Rmax=300.0,zmax=300.0,Rmin=4.0,
                                        d0min=2.0,nmin=5,mDVmin=mDVcut)  for llp in llps])
            good_llps = np.array(llps)[v_acc > 0.0]
            if len(good_llps) == 0:
                continue
            # Event efficiency
            evt_eff = eventEff(met,good_llps)

            ns = ns*evt_eff
            
            # Vertex efficiencies:
            v_eff = np.array([vertexEff(llp,strategy=effStrategy) for llp in llps])
            
            wvertex = 1.0-np.prod(1.0-v_acc*v_eff)

            ns = ns*wvertex
            
            # Add to the total weight in each SR:
            cutFlow["DV selection"] += (ns,ns**2)

        f.Close()
    progressbar.finish()

    modelDict['Total xsec (pb)'] = totalweightPB
    print('\nCross-section (pb) = %1.3e\n' %totalweightPB)

    cutFlowErr = {k : np.sqrt(v[1]) for k,v in cutFlow.items()}
    cutFlow = {k : v[0]  for k,v in cutFlow.items()}
    
    # Compute normalized cutflow
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
    dataDict['Luminosity (1/fb)'] = [lumi]
    dataDict['$N_s$'] = [cutFlow["DV selection"][0]]
    dataDict['$N_s$ Err'] = [cutFlowErr["DV selection"][0]]
    dataDict['AccEff'] = [cutFlow["DV selection"][1]]
    dataDict['AccEffErr'] = [cutFlowErr["DV selection"][1]]
    dataDict['VertexEff Strategy'] = [effStrategy]
    dataDict['mDV cut'] = [mDVcut]
    for cut,val in cutFlow.items():
        dataDict.setdefault(cut,[val])
        dataDict.setdefault(cut+' Error',[cutFlowErr[cut]])

    # Create a dictionary for storing data
    dataDict.update(modelDict)
   

    return dataDict


    


if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser( description=
            "Run the recasting for ATLAS-SUSY-2016-08 using one or multiple Delphes ROOT files as input. "
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
    ap.add_argument('-S', '--effstrategy', required=False,type=str,default='official',
            help='Defines which strategy to use for the applying the vertex efficiencies (official, nearest, average).')
    ap.add_argument('-mDV', '--mDVcut', required=False,type=float,default=10.0,
            help='Value for the mDV cut.')
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

    # Set random seed
    # np.random.seed(22)
    np.random.seed(15)

    t0 = time.time()

    # # Set output file
    args = ap.parse_args()
    inputFiles = args.inputFile
    outputFile = args.outputFile

    if args.effstrategy not in  ['official','average','nearest']:
        print("Select a valid vertex efficiency strategy")
        sys.exit()

    if args.update:
        print('\n\n======= Updating files with new event data ==========\n')

    # Split input files by distinct models and get recast data for
    # the set of files from the same model:
    for fileList,mDict in splitModels(inputFiles,args.model):

        # Set output file
        if outputFile is None:
            if args.effstrategy == 'official':
                outFile = fileList[0].replace('delphes_events.root','atlas_2016_08.pcl')
            elif args.effstrategy == 'average':
                outFile = fileList[0].replace('delphes_events.root','atlas_2016_08_average.pcl')
            elif args.effstrategy == 'nearest':
                outFile = fileList[0].replace('delphes_events.root','atlas_2016_08_nearest.pcl')
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

        dataDict = getRecastData(fileList,args.model,mDict,
                                 effStrategy=args.effstrategy,mDVcut=args.mDVcut,addweights=args.add)
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
