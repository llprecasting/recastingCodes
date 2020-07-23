#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import time,os,sys
# import random
# random.seed(10) #Set the random seeed to generate reproduciable results

start_time = time.time()

current_dir = os.getcwd()
Modules_dir = os.path.abspath(os.path.join(current_dir,'Modules'))

sys.path.append(Modules_dir)
from eventSelector import EventSelector

### Define variables ##
# User can provide different tau values in order to study other regions of the parameter space. The range of tau values used for the analysis is (0.01 - 10) nano-seconds
tau_array = np.array([0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])

#Create an event selector:
evtSelector = EventSelector('input_param.dat',tau_array)

#Set Delphes path and load required libraries:
evtSelector.loadDelphesLib('/home/felipe/Documents/Programas/MG5_aMC_v2_6_5/Delphes')

#Load ROOT/Delphes events file:
evtSelector.loadRootFile()

#Select events and compute weights for the lifetimes in tau_array:
evtSelector.selectEvents()

#Save output to file:
evtSelector.writeResults()


final_time = time.time() - start_time
if (final_time < 60.0): print("time final =",final_time , "seconds")
if (final_time >= 60.0 and final_time < 3600): print("time final =",final_time/60.0 , "minutes")
if (final_time >= 3600): print("time final =",final_time/3600.0 , "hous")
