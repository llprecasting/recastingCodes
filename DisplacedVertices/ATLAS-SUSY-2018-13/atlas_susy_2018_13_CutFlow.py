#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import time
import progressbar as P
from ATLAS_data.effFunctions import eventEff,vertexEff
from atlas_susy_2018_13_Recast import (getLLPs, getJets, getDisplacedJets, eventAcc, 
                                       vertexAcc, getModelDict)

delphesDir = os.path.abspath("./DelphesLLP")
os.environ['ROOT_INCLUDE_PATH'] = os.path.join(delphesDir,"external")

import ROOT

ROOT.gSystem.Load(os.path.join(delphesDir,"libDelphes.so"))

ROOT.gInterpreter.Declare('#include "classes/SortableObject.h"')
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')


# ### Define dictionary to store data
def getCutFlow(inputFiles,normalize=False,model='strong',sr='HighPT',nevtsMax=-1):

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
        if nevtsMax > 0:
            nevts = min(nevtsMax,nevts)
        modelDict['Total MC Events'] += nevts        
        nevtsDict[inputFile] = nevts
        f.Close()


    lumi = 139.0
    totalweightPB = 0.0
    # Keep track of yields for each dataset
    cutFlowAcceptance = {}
    if model == 'strong':
        cutFlowAcceptance['$m_{\\tilde g} (GeV)$'] = modelDict['mLSP']
        cutFlowAcceptance['$m_{\\tilde \chi_1^0} (GeV)$'] = modelDict['mLLP']
        cutFlowAcceptance['$\\tau(\\tilde \chi_1^0) (ns)$'] = modelDict['tau_ns']
    elif model == 'ewk':
        cutFlowAcceptance['$m_{\\tilde \chi_1^0} (GeV)$'] = modelDict['mLLP']
        cutFlowAcceptance['$\\tau(\\tilde \chi_1^0) (ns)$'] = modelDict['tau_ns']
    cutFlowAcceptance.update({"Total" : 0.0,
                "Jet selection" : 0.0,
                "$R_{xy},z <$ 300 mm" : 0.0,
                "$R_{DV} > 4$ mm" : 0.0,
                "$d_0 > 2$ mm" : 0.0,
                "$nTracks >= 5$" : 0.0,
                "mDV > 10 GeV" : 0.0,
                "final Acc*Eff" : 0.0
                })
    
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

            llps = getLLPs(tree.bsm,tree.bsmDirectDaughters,tree.bsmFinalDaughters)
            jets = getJets(tree.GenJet,pTmin=25.,etaMax=5.0)
            jetsDisp = getDisplacedJets(jets,llps)
            
            # Event acceptance:
            jet_acc = eventAcc(jets,jetsDisp,sr=sr)

            cutFlowAcceptance["Total"] += ns
            if (not jet_acc): continue
            ns = ns*jet_acc

            cutFlowAcceptance["Jet selection"] += ns
            
            llpsSel = [llp for llp in llps if vertexAcc(llp,Rmax=300.0,zmax=300.0)]
            if not llpsSel: continue
            cutFlowAcceptance["$R_{xy},z <$ 300 mm"] += ns

            llpsSel = [llp for llp in llps if vertexAcc(llp,Rmax=300.0,zmax=300.0,Rmin=4.0)]
            if not llpsSel: continue
            cutFlowAcceptance["$R_{DV} > 4$ mm"] += ns

            llpsSel = [llp for llp in llps if vertexAcc(llp,Rmax=300.0,zmax=300.0,Rmin=4.0,d0min=2.0)]
            if not llpsSel: continue
            cutFlowAcceptance["$d_0 > 2$ mm"] += ns

            llpsSel = [llp for llp in llps if vertexAcc(llp,Rmax=300.0,zmax=300.0,Rmin=4.0,d0min=2.0,nmin=5)]
            if not llpsSel: continue
            cutFlowAcceptance["$nTracks >= 5$"] += ns

            llpsSel = [llp for llp in llps if vertexAcc(llp,Rmax=300.0,zmax=300.0,Rmin=4.0,d0min=2.0,nmin=5,mDVmin=10.0)]
            if not llpsSel: continue
            cutFlowAcceptance["mDV > 10 GeV"] += ns

            # Event efficiency
            ev_eff = eventEff(jets,llps,sr=sr)
            # Vertex acceptances:
            v_acc = np.array([vertexAcc(llp,Rmax=300.0,zmax=300.0,Rmin=4.0,d0min=2.0,nmin=5,mDVmin=10.0) for llp in llps])
            # Vertex efficiencies:
            v_eff = np.array([vertexEff(llp) for llp in llps])

            wvertex = 1.0-np.prod(1.0-v_acc*v_eff)

            cutFlowAcceptance['final Acc*Eff'] += ns*ev_eff*wvertex

        f.Close()
    progressbar.finish()

    modelDict['Total xsec (pb)'] = totalweightPB
    print('\nCross-section (pb) = %1.3e\n' %totalweightPB)

    # Compute normalized cutflow
    for key,val in cutFlowAcceptance.items():
        if key == 'Total' or '(GeV)' in key or '(ns)' in key:
            continue
        valNorm = float('%1.3e' %(val/cutFlowAcceptance['Total']))
        cutFlowAcceptance[key] = valNorm
    cutFlowAcceptance['Total'] = 1.0

    print('Acceptance for %s:' %sr)
    for k,v in cutFlowAcceptance.items():
        print('%s : %1.3f%%' %(k,v*1e2))
    
    return cutFlowAcceptance


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
        outputFile = inputFiles[0].replace('delphes_events.root','atlas_2018_42_cutflow.pcl')

    if os.path.splitext(outputFile)[1] != '.pcl':
        outputFile = os.path.splitext(outputFile)[0] + '.pcl'

    cutFlowAcceptance = getCutFlow(inputFiles,args.normalize,args.model,args.SR,args.nevts)

    # ### Save DataFrame to pickle file
    # #### Create pandas DataFrame
    df = pd.DataFrame.from_dict(cutFlowAcceptance, orient='index')
    print('Saving to',outputFile)
    df.to_pickle(outputFile)

    print("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
