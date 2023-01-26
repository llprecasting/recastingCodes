#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import interp1d,griddata,LinearNDInterpolator
import os
atlasDir = os.path.dirname(os.path.abspath(__file__))

inputFile = os.path.join(atlasDir,'./HEPData-ins2080541-v1-csv/MuonReconstructionEfficiencydistribution.csv')
gridPtsReco = np.genfromtxt(inputFile,names=True,skip_header=7,delimiter=',')
muonEff_F = LinearNDInterpolator((gridPtsReco['beta'],gridPtsReco['eta']),gridPtsReco['Efficiency'],fill_value=0.0)

def getMuonRecoEff(beta,eta,pid=None):
    """
    Return the interpolated muon reconstruction efficiency
    as a function of beta and eta. If the pid is 
    
    :param beta: HSCP beta
    :param eta: HSCP eta
    
    :return: Efficiency
    """
    
    eta = abs(eta)
    
    eff = muonEff_F(beta,eta)
    
    return float(eff)


inputFile = os.path.join(atlasDir,'./HEPData-ins2080541-v1-csv/TriggerEfficiencydistribution.csv')
gridPtsTrig = np.genfromtxt(inputFile,names=True,skip_header=7,delimiter=',')
triggerEff_F = interp1d(gridPtsTrig['Truth_Level_ETmiss_calo_GeV'],gridPtsTrig['Efficiency'],
                      fill_value=0.0,bounds_error=False)

def getTriggerEff(metCalo):
    
    if metCalo > 1000.:
        return 1.0
    eff = triggerEff_F(metCalo)
    return float(eff)        


inputFile = os.path.join(atlasDir,'./HEPData-ins2080541-v1-csv/EventSelectionEfficiencydistribution.csv')
gridPtsSel = np.genfromtxt(inputFile,names=True,skip_header=7,delimiter=',')
selectionEff_F = interp1d(gridPtsSel['Truth_Level_ETmiss_Total_GeV'],gridPtsSel['Efficiency'],
                      fill_value=0.0,bounds_error=False)

def getSelectionEff(met):
    
    if met > 1000.:
        return 1.0
    eff = selectionEff_F(met)
    return float(eff)

inputFile = os.path.join(atlasDir,'./HEPData-ins2080541-v1-csv/TrackSelectionEfficiencydistribution.csv')
gridPtsHigh = np.genfromtxt(inputFile,names=True,skip_header=34,delimiter=',')
# Central value
trackEffHigh_F = interp1d(gridPtsHigh['beta__gamma'],gridPtsHigh['Efficiency'],
                      fill_value=0.0,bounds_error=False)
# Upper value                    
trackEffHigh_Fp = interp1d(gridPtsHigh['beta__gamma'],gridPtsHigh['error_'],
                      fill_value=0.0,bounds_error=False)
# Lower value                      
trackEffHigh_Fm = interp1d(gridPtsHigh['beta__gamma'],gridPtsHigh['error__1'],
                      fill_value=0.0,bounds_error=False)


gridPtsLow = np.genfromtxt(inputFile,names=True,skip_header=8,skip_footer=32-8,delimiter=',')
# Central value
trackEffLow_F = interp1d(gridPtsLow['beta__gamma'],gridPtsLow['Efficiency'],
                      fill_value=0.0,bounds_error=False)
# Upper value
trackEffLow_Fp = interp1d(gridPtsLow['beta__gamma'],gridPtsLow['error_'],
                      fill_value=0.0,bounds_error=False)
# Lower value
trackEffLow_Fm = interp1d(gridPtsLow['beta__gamma'],gridPtsLow['error__1'],
                      fill_value=0.0,bounds_error=False)
@np.vectorize
def getTrackEff(gbeta,sr,use='central'):
    
    eff = 0.0
    if sr.lower() == 'high':
        if use == 'central':
            eff = trackEffHigh_F(gbeta)
        elif use == 'lower':
            eff = trackEffHigh_F(gbeta)+trackEffHigh_Fm(gbeta)
        elif use == 'higher':
            eff = trackEffHigh_F(gbeta)+trackEffHigh_Fp(gbeta)
        elif use == 'smear':
            mu = trackEffHigh_F(gbeta)
            sigma = (trackEffHigh_Fp(gbeta)-trackEffHigh_Fm(gbeta))/2.0
            eff = np.random.normal(mu,sigma,1)[0]
        else:
            raise ValueError()
    elif sr.lower() == 'low':
        if use == 'central':
            eff = trackEffLow_F(gbeta)
        elif use == 'lower':
            eff = trackEffLow_F(gbeta)+trackEffLow_Fm(gbeta)
        elif use == 'higher':
            eff = trackEffLow_F(gbeta)+trackEffLow_Fp(gbeta)
        elif use == 'smear':
            mu = trackEffLow_F(gbeta)
            sigma = (trackEffLow_Fp(gbeta)-trackEffLow_Fm(gbeta))/2.0
            eff = np.random.normal(mu,sigma,1)[0]
        else:
            raise ValueError()
    else:
        raise ValueError()

    return float(eff)


# Load mass bins for long lifetimes
inputFile = os.path.join(atlasDir,'./HEPData-ins2080541-v1-csv/p0-valuesandmodel-independentlimits,longregime.csv')
massLongLT = np.genfromtxt(inputFile,names=True,skip_header=9,skip_footer=397-73,delimiter=',',dtype='f8,U15')
# Create new array with lower bin and higher bin values
new_dt = np.dtype(massLongLT.dtype.descr + [('Mass_window_Low', 'f8'),('Mass_window_High', 'f8')])
massLong = np.empty(massLongLT.shape, dtype=new_dt)
for c in massLongLT.dtype.names:
    massLong[c] = massLongLT[c]
massLong['Mass_window_Low'] = [eval(x.strip().replace(" ",","))[0] for x in massLongLT['Mass_window_GeV']]
massLong['Mass_window_High'] = [eval(x.strip().replace(" ",","))[1] for x in massLongLT['Mass_window_GeV']]

# Load mass bins for short lifetimes
inputFile = os.path.join(atlasDir,'./HEPData-ins2080541-v1-csv/p0-valuesandmodel-independentlimits,shortregime.csv')
massShortLT = np.genfromtxt(inputFile,names=True,skip_header=9,skip_footer=280-64,delimiter=',',dtype='f8,U15')
# Create new array with lower bin and higher bin values
new_dt = np.dtype(massShortLT.dtype.descr + [('Mass_window_Low', 'f8'),('Mass_window_High', 'f8')])
massShort = np.empty(massShortLT.shape, dtype=new_dt)
for c in massShortLT.dtype.names:
    massShort[c] = massShortLT[c]
massShort['Mass_window_Low'] = [eval(x.strip().replace(" ",","))[0] for x in massShortLT['Mass_window_GeV']]
massShort['Mass_window_High'] = [eval(x.strip().replace(" ",","))[1] for x in massShortLT['Mass_window_GeV']]

def getTargetMass(trueMass,trueWidth=0.0):
    """
    Returns the target mass (which falls within a given mass window)
    as a function of the true LLP mass and its decay regime (long-lived or short-lived).
    According to the ATLAS guide, we should assume the long-lived regime (default)
    """
    
    # The long regime is intended for particles with tau > 1 ns
    # while the short regime windows are intended for particles with tau <= 1 ns.
    if trueWidth == 0.0 or (6.582e-25/trueWidth) > 1e-9:
        massWindow = massLong
    else:
        massWindow = massShort
    
    targetMass = None
    for pt in massWindow:
        if pt['Mass_window_Low'] <= trueMass < pt['Mass_window_High']:
            targetMass = pt['Target_Mass_GeV']
            
            
    return targetMass

@np.vectorize
def getMassSelEff(targetMass,sr):            

    if sr.lower() == 'high':
        wMass = max(0.,0.74 - 0.052*(targetMass/1000.))
    else:
        wMass = 0.6

    return wMass                    