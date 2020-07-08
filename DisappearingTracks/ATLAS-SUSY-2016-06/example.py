#!/usr/bin/env python

import numpy as np
import time
from __future__ import print_function

start_time = time.clock()

current_dir = os.getcwd()
Modules_dir = current_dir+'/Modules/'

### Define .root input file ###
rootFile = 'example1.root'
### Define name of output file ###
outputFile = 'example_output.dat'

### Define variables ##
# User can provide different tau values in order to study other regions of the parameter space. The range of tau values used for the analysis is (0.01 - 10) nano-seconds
tau_array = np.array([0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
### Define Luminosity ###
Lum = 36.1
### Define chargino/neutralino PID ### (if PID_chargino/PID_neutralino = None, default value for SUSY models is taken Chargino_PID=1000024 and Neutralino_PID=1000022)
PID_chargino = None
PID_neutralino = None
### Define k-factor of sample ### (if kfactor = None, default value is 1)
kfactor = 1.227

exec(open(Modules_dir+'LLP_DM.py').read())

final_time = time.clock() - start_time
if (final_time < 60.0): print("time final =",final_time , "seconds")
if (final_time >= 60.0 and final_time < 3600): print("time final =",final_time/60.0 , "minutes")
if (final_time >= 3600): print("time final =",final_time/3600.0 , "hous")
