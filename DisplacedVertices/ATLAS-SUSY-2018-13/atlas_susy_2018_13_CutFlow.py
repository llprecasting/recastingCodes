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
from atlas_susy_2018_13_Recast import eventAcc, vertexAcc

delphesDir = os.path.abspath("../DelphesLLP")
os.environ['ROOT_INCLUDE_PATH'] = os.path.join(delphesDir,"external")

import ROOT

ROOT.gSystem.Load(os.path.join(delphesDir,"libDelphes.so"))

ROOT.gInterpreter.Declare('#include "classes/SortableObject.h"')
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')


# ### Define dictionary to store data
def getCutFlow(inputFiles,model='sbottom',sr='HighPT',nevtsMax=-1,modelDict=None,addweights=False):

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
        if nevtsMax > 0:
            nevts = min(nevtsMax,nevts)
        modelDict['Total MC Events'] += nevts        
        nevtsDict[inputFile] = nevts
        f.Close()


    lumi = 139.0
    totalweightPB = 0.0
    # Keep track of yields for each dataset
    cutFlow = {}
    if model == 'strong':
        cutFlow['$m_{\\tilde g} (GeV)$'] = modelDict['mLSP']
        cutFlow['$m_{\\tilde \chi_1^0} (GeV)$'] = modelDict['mLLP']
        cutFlow['$\\tau(\\tilde \chi_1^0) (ns)$'] = modelDict['tau_ns']
    elif model == 'ewk':
        cutFlow['$m_{\\tilde \chi_1^0} (GeV)$'] = modelDict['mLLP']
        cutFlow['$\\tau(\\tilde \chi_1^0) (ns)$'] = modelDict['tau_ns']
    elif model == 'bb':
        cutFlow['$m_{\\tilde b_1} (GeV)$'] = modelDict['mLLP']
        cutFlow['$\\tau(\\tilde b_1) (ns)$'] = modelDict['tau_ns']

    keys = ["Total", "Jet selection", "$R_{xy},z <$ 300 mm", 
            "$R_{DV} > 4$ mm", "$d_0 > 2$ mm", 
            "$nTracks >= 5$", "$mDV > 10$ GeV", "final Acc*Eff"]
    for k in keys:
        cutFlow[k] = np.zeros(2)
    
    progressbar = P.ProgressBar(widgets=["Reading %i Events: " %modelDict['Total MC Events'], 
                                P.Percentage(),P.Bar(marker=P.RotatingMarker()), P.ETA()])
    progressbar.maxval = modelDict['Total MC Events']
    progressbar.start()

    ntotal = 0
    totalweightPB = 0.0
    for inputFile in inputFiles:
        f = ROOT.TFile(inputFile,'read')
        tree = f.Get("Delphes")
        nevts = nevtsDict[inputFile]
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

            llps = getLLPs(tree.bsm,tree.bsmDirectDaughters,tree.bsmFinalDaughters)
            jets = getJets(tree.GenJet,pTmin=25.,etaMax=5.0)
            jetsDisp = getDisplacedJets(jets,llps)
            
            # Event acceptance:
            jet_acc = eventAcc(jets,jetsDisp,sr=sr)

            cutFlow["Total"] += (ns,ns**2)
            if (not jet_acc): continue
            ns = ns*jet_acc

            cutFlow["Jet selection"] += (ns,ns**2)
            
            llpsSel = [llp for llp in llps if vertexAcc(llp,Rmax=300.0,zmax=300.0)]
            if not llpsSel: continue
            cutFlow["$R_{xy},z <$ 300 mm"] += (ns,ns**2)

            llpsSel = [llp for llp in llps if vertexAcc(llp,Rmax=300.0,zmax=300.0,Rmin=4.0)]
            if not llpsSel: continue
            cutFlow["$R_{DV} > 4$ mm"] += (ns,ns**2)

            llpsSel = [llp for llp in llps if vertexAcc(llp,Rmax=300.0,zmax=300.0,Rmin=4.0,d0min=2.0)]
            if not llpsSel: continue
            cutFlow["$d_0 > 2$ mm"] += (ns,ns**2)

            llpsSel = [llp for llp in llps if vertexAcc(llp,Rmax=300.0,zmax=300.0,Rmin=4.0,d0min=2.0,nmin=5)]
            if not llpsSel: continue
            cutFlow["$nTracks >= 5$"] += (ns,ns**2)

            llpsSel = [llp for llp in llps if vertexAcc(llp,Rmax=300.0,zmax=300.0,Rmin=4.0,d0min=2.0,nmin=5,mDVmin=10.0)]
            if not llpsSel: continue
            cutFlow["$mDV > 10$ GeV"] += (ns,ns**2)

            # Event efficiency
            ev_eff = eventEff(jets,llps,sr=sr)
            # Vertex acceptances:
            v_acc = np.array([vertexAcc(llp,Rmax=300.0,zmax=300.0,Rmin=4.0,d0min=2.0,nmin=5,mDVmin=10.0) for llp in llps])
            # Vertex efficiencies:
            v_eff = np.array([vertexEff(llp) for llp in llps])

            wvertex = 1.0-np.prod(1.0-v_acc*v_eff)

            ns = ns*ev_eff*wvertex
            cutFlow['final Acc*Eff'] += (ns, ns**2)

        f.Close()
    progressbar.finish()

    modelDict['Total xsec (pb)'] = totalweightPB
    print('\nCross-section (pb) = %1.3e\n' %totalweightPB)

    cutFlowErr = {k : np.sqrt(v[1]) for k,v in cutFlow.items()}
    cutFlow = {k : v[0]  for k,v in cutFlow.items()}


    # Compute normalized cutflow
    for key,val in cutFlow.items():
        if key == 'Total' or '(GeV)' in key or '(ns)' in key:
            continue
        valNorm = float('%1.3e' %(val/cutFlow['Total']))
        errNorm = float('%1.3e' %(cutFlowErr[key]/cutFlow['Total']))
        cutFlow[key] = valNorm
        cutFlowErr[key] = errNorm
    cutFlow['Total'] = 1.0
    cutFlowErr['Total'] = 0.0

    print('Acceptance for %s:' %sr)
    for k,v in cutFlow.items():
        if v != 0.0:
            print('%s : %1.3e +- %1.1f%%' %(k,v,1e2*cutFlowErr[k]/v))
        else:
            print('%s : %1.3e +- ??' %(k,v))

    return modelDict,cutFlow,cutFlowErr


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
    ap.add_argument('-s', '--SR', required=False,type=str,default='HighPT',
            help='Defines which signal region should be considered for the cutflow (HighPT or Trackless).')
    ap.add_argument('-N', '--nevts', required=False,type=int,default=-1,
            help='Maximum number of events to use.')


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

    # Split input files by distinct models and get recast data for
    # the set of files from the same model:
    for fileList,mDict in splitModels(inputFiles,args.model):
        modelDict,cutFlow,cutFlowErr = getCutFlow(fileList,args.model,modelDict=mDict,sr=args.SR,addweights=args.add)
        dataDict = {key : [val] for key,val in modelDict.items()}
        for key,val in cutFlow.items():
            dataDict[key] = [(val,cutFlowErr[key])]

        if outputFile is None:
            outFile = fileList[0].replace('delphes_events.root','atlas_2018_13_cutflow.pcl')
        else:
            outFile = outputFile[:]

        if os.path.splitext(outFile)[1] != '.pcl':
            outFile = os.path.splitext(outFile)[0] + '.pcl'

        # #### Create pandas DataFrame
        df = pd.DataFrame.from_dict(dataDict)

        # ### Save DataFrame to pickle file
        print('Saving to',outFile)
        df.to_pickle(outFile)
        print('\n\n')

    print("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
