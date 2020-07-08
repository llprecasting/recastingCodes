import sys, os, copy
import ROOT
#import Functions as Ft
## Define Delphes libraries path ##
Delphes_libraries = '/Users/jzurita/Physics/MG5_aMC_v2_7_3/Delphes/libDelphes'
ROOT.gSystem.Load(Delphes_libraries)
from ROOT import TTree, TFile, gROOT, gDirectory, TLorentzVector
import numpy as np
import random
import math

exec(open(Modules_dir+'Functions.py').read())

######## root input file ########
rootFilePath=current_dir+'/'+rootFile
######## output file ######
OutputFilePath = current_dir+'/'+outputFile

######## kfactor option #####
if(kfactor is not None): kfac = kfactor
else: kfac = 1.0

fAccEff = open(OutputFilePath,'w')
fAccEff.write(' Char_Mass(GeV)   Neutr_Mass(GeV)   tau(ns)      EAxEE          TAxTE         EAxEExTAxTE    xs100(fb)    Ntr        MCev  \n')

## Event and Tracklet Selection file
exec(open(Modules_dir+'Event_selection.py').read())
## Calculate Tracklet Acceptance x Efficiency and write the final result
exec(open(Modules_dir+'Acc_Eff_Write.py').read())
