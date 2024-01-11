#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import interp1d,LinearNDInterpolator
import os,glob
atlasDir = os.path.dirname(os.path.abspath(__file__))


functions_event_eff = {}
for f in glob.glob(os.path.join(atlasDir,'./HEPData-ins2628398-v1-csv/event_efficiency_*.csv')):
    fname = os.path.basename(f)
    sr = fname.split('_')[2].lower()
    if not sr in functions_event_eff:
        functions_event_eff[sr] = {}
    Rmin = int(fname.split('_')[4])
    if len(fname.split('_')) == 7:
        Rmax = int(fname.split('_')[5])
    elif Rmin == 3870:
        Rmax = np.inf
    elif Rmin == 1150:
        Rmin = 0
        Rmax = 1150
        
    pts = np.genfromtxt(f,names=True,skip_header=10,delimiter=',')
    eff_F = interp1d(pts['Sumpt_GeV'],pts['Efficiency'],fill_value=(0.0,pts['Efficiency'][-1]),bounds_error=False)
    functions_event_eff[sr][(Rmin,Rmax)] = eff_F


def eventEff(jets,llps,sr):
    
    jetPT = 0.0
    for j in jets:
        jetPT += j.PT
    r_max = max([np.sqrt(llp.Xd**2 + llp.Yd**2) 
                 for llp in llps])
    
    # Select eff function
    eff_F = [f for Rint,f in functions_event_eff[sr.lower()].items() 
             if Rint[0] < r_max <= Rint[1]]
    if not eff_F:
        return 0.0
    else:
        eff_F = eff_F[0]
        return eff_F(jetPT/1e3) # jetPT in TeV!


functions_vertex_eff = {}
for f in glob.glob(os.path.join(atlasDir,'./HEPData-ins2628398-v1-csv/vertex_efficiency_*.csv')):
    fname = os.path.basename(f)
    Rmin = int(fname.split('_')[3])    
    if len(fname.split('_')) == 6:
        Rmax = int(fname.split('_')[4])
    elif Rmin == 22:
        Rmin = 4
        Rmax = 22
    pts = np.genfromtxt(f,names=True,skip_header=10,delimiter=',')
    eff_F = LinearNDInterpolator((pts['m_DV_GeV'],pts['n_tracks']),pts['Efficiency'],fill_value=0.0)
    functions_vertex_eff[(Rmin,Rmax)] = eff_F


def vertexEff(llp):
    
    r = np.sqrt(llp.Xd**2 + llp.Yd**2)
    mDV = llp.mDV
    n = llp.nTracks
    
    # Select eff function
    eff_F = [f for Rint,f in functions_vertex_eff.items() if Rint[0] < r <= Rint[1]]
    if not eff_F:
        return 0.0
    else:
        eff_F = eff_F[0]
        return eff_F(mDV,n)
