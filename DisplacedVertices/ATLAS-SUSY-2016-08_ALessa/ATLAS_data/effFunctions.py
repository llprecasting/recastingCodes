#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import interp1d,LinearNDInterpolator,NearestNDInterpolator
import os
atlasDir = os.path.dirname(os.path.abspath(__file__))


functions_event_eff = {}
files = {(0.,1150.) : "Table22.csv", (1150.,3870.) :  "Table23.csv", (3870.,np.inf) : "Table24.csv"}
for (Rmin,Rmax),f in files.items():
    fname = os.path.join(atlasDir,'HEPData-ins1630632-v2-csv',f)        
    pts = np.genfromtxt(fname,names=True,
                        comments='#',delimiter=',',skip_header=10)
    eff_F = interp1d(pts['Truth_MET_GeV'],pts['Event_selection_efficiency'],
                     fill_value=(0.0,pts['Event_selection_efficiency'][-1]),
                     bounds_error=False,
                      kind='linear')
                    #  kind='nearest')
    functions_event_eff[(Rmin,Rmax)] = eff_F


def eventEff(met,llps):

    r_max = max([np.sqrt(llp.Xd**2 + llp.Yd**2) 
                 for llp in llps])
    
    # Select eff function
    eff_F = [f for Rint,f in functions_event_eff.items() 
             if Rint[0] < r_max <= Rint[1]]
    if not eff_F:
        return 0.0
    else:
        eff_F = eff_F[0]
        return eff_F(met)

functions_vertex_eff = {}
files = {(4.,22.) : "Table25.csv", 
         (22.,25.) : "Table26.csv", 
         (25.,29.) : "Table27.csv", 
         (29.,38.) : "Table28.csv", 
         (38.,46.) : "Table29.csv", 
         (46.,73.) : "Table30.csv", 
         (73.,84.) : "Table31.csv", 
         (84.,111.) : "Table32.csv", 
         (111.,120.) : "Table33.csv", 
         (120.,145.) : "Table34.csv", 
         (145.,180.) : "Table35.csv", 
         (180.,300.) : "Table36.csv"}

for (Rmin,Rmax),f in files.items():
    fname = os.path.join(atlasDir,'HEPData-ins1630632-v2-csv',f)        
    pts = np.genfromtxt(fname,names=True,
                        comments='#',delimiter=',',skip_header=10)

    eff_F = LinearNDInterpolator((pts['m_Truth_vertex_GeV'],pts['Number_of_tracks_Truth_vertex']),pts['Vertex_selection_efficiency'],fill_value=0.0)
    # eff_F = NearestNDInterpolator((pts['m_Truth_vertex_GeV'],pts['Number_of_tracks_Truth_vertex']),pts['Vertex_selection_efficiency'])
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
