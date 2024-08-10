#!/usr/bin/env python3

from typing import Any,Union
import numpy as np
import os
import glob
import pyslha

class LLP(object):
    """
    Convenience class for holding the LLP, its daughters and several
    other useful methods
    """

    def __init__(self, candidate, direcdaughters, finaldaughters, mothers=[],
                 maxMomViolation=5e-2,trackEff = 1.0) -> None:
        self._candidate = candidate
        self.directDaughters = direcdaughters[:]
        self.finalDaughters = finaldaughters[:]
        self.mothers = mothers[:]
        self._selectedDecays = []
        self.nTracks = None
        self.mDV = None
        self.Xd = np.inf
        self.Yd = np.inf
        self.Zd = np.inf
        self.r_decay = np.inf
        self.z_decay = np.inf
        trimom = np.sqrt(self.Px**2 +self.Py**2 + self.Pz**2)
        self.beta = trimom/self.E
        self.gbeta = trimom/self.Mass

        if candidate.Status == 102: # LLP has not decayed (use final R-Hadron as stable LLP)
            # Get final R-hadron:
            finalRHadrons = [d for d in finaldaughters if d.Mass >= candidate.Mass]
            if len(finalRHadrons) != 1:
                print("Error getting stable R-Hadron (%i)" %len(finalRHadrons))                
            # Use final daughter as stable LLP
            self._candidate = finalRHadrons[0]
            self.directDaughters = []
            self.finalDaughters = []
            self.mothers = [candidate]

        for d in self.finalDaughters:
            if d.Status != 1:
                continue
            if d.Charge == 0:
                continue
            if d.PT < 1.0:
                continue
            pTratio = abs(d.PT/d.Charge)
            if pTratio < 1.0:
                continue
            
            if trackEff < 1.0 and np.random.random() > trackEff: # Apply random efficiency for track reco?
                continue
            self._selectedDecays.append(d)
            

        # Force the computation of mDV and nTracks
        pTot = np.zeros(4)
        for d in self._selectedDecays:
            p = np.array([0.,d.Px,d.Py,d.Pz])
            mpion = 0.140
            p[0] = np.sqrt(mpion**2 + np.dot(p[1:],p[1:]))
            pTot = pTot + p
        mDV = np.sqrt(pTot[0]**2 - np.dot(pTot[1:],pTot[1:]))
        self.mDV = mDV
        self.nTracks = len(self._selectedDecays)

        # Consistency checks:
        if self.directDaughters:
            pTot = np.array([self.E,self.Px,self.Py,self.Pz])
            pNorm = np.linalg.norm(pTot)
            for d in self.directDaughters:                
                pTot -= np.array([d.E,d.Px,d.Py,d.Pz])
            if self.directDaughters and np.linalg.norm(pTot)/pNorm > maxMomViolation/1e3: # Be more strict about direct daughters
                # raise ValueError("Error getting direct daughters, momentum conservation violated by %1.1e!" %(np.linalg.norm(pTot)/pNorm))
                print("Error getting direct daughters, momentum conservation violated by %1.1e!" %(np.linalg.norm(pTot)/pNorm))

            rList = [np.sqrt(d.X**2 + d.Y**2 + d.Z**2) for d in self.directDaughters]

            if max(rList) != min(rList) and abs(max(rList)-min(rList))/(max(rList)+min(rList)) > 0.01:
                raise ValueError("Direct daughters do not have the same production vertex!")
            
            daughter = self.directDaughters[0]
            self.Xd = daughter.X
            self.Yd = daughter.Y
            self.Zd = daughter.Z
            self.r_decay = np.sqrt(self.Xd**2+self.Yd**2)
            self.z_decay = self.Zd

        if self.finalDaughters:      
            pTot = np.array([self.E,self.Px,self.Py,self.Pz])
            pNorm = np.linalg.norm(pTot)
            for d in self.finalDaughters:                
                pTot -= np.array([d.E,d.Px,d.Py,d.Pz])
            if self.finalDaughters and np.linalg.norm(pTot)/pNorm > maxMomViolation:
                # raise ValueError("Error getting final daughters, momentum conservation violated by %1.1e! (%s)" %(np.linalg.norm(pTot)/pNorm,str(pTot)))
                print("Error getting final daughters, momentum conservation violated by %1.1e! (%s)" %(np.linalg.norm(pTot)/pNorm,str(pTot)))
            

    def __getattr__(self, attr: str) -> Any:
        try:
            return self.__getattribute__(attr)
        except AttributeError:
            try:
                return self._candidate.__getattribute__(attr)
            except:
                raise AttributeError("Could not get attribute %s" %attr)
            
    def __str__(self):
        return "%s (%s)" %(self.PID,self.Status)
    
    def __repr__(self) -> str:
        return str(self)   

def getLLPs(llpList,directDaughters,finalDaughters,mothers=[],maxMomViolation=5e-2,trackEff=1.0):

    llps = []
    for ip,p in enumerate(llpList):
        # Get direct daughters
        llp_ddaughters = []        
        for d in directDaughters:
            if d.M1 == ip:
                llp_ddaughters.append(d)
        # Get final daughters
        llp_fdaughters = []
        for d in finalDaughters:
            if d.M1 == ip:
                llp_fdaughters.append(d)
        # Get mothers
        llp_mothers = list(mothers)[p.M1:p.M2+1]        
        # Get final daughters
        llps.append(LLP(p,llp_ddaughters,llp_fdaughters,llp_mothers,maxMomViolation,trackEff))
        
    return llps

def getJets(jets,pTmin=50.0,etaMax=5.0):
    """
    Select jets with pT > pTmin and |eta| < etaMax
    """

    jetsSel = []
    for ijet in range(jets.GetEntries()):
        jet = jets.At(ijet)
        if jet.PT < pTmin:
            continue
        if abs(jet.Eta) > etaMax:
            continue
        jetsSel.append(jet)
    
    return jetsSel

def getDisplacedJets(jets,llps,skipPIDs=[1000022]):
    """
    Select from the list of all jets, the displaced jets associated
    with a LLP decay.
    """

    displacedJets = []
    for jet in jets:
        deltaRmin = 0.3
        llpMatch = None
        for llp in llps: 
            for daughter in llp.directDaughters:
                if abs(daughter.PID) in skipPIDs:
                    continue
                deltaR = np.sqrt((jet.Eta-daughter.Eta)**2 + (jet.Phi-daughter.Phi)**2)
                if deltaR < deltaRmin:
                    deltaRmin = deltaR
                    llpMatch = llp # Store LLP parent
        
        jet.llp = llpMatch
        if llpMatch is not None:
            R = np.sqrt(llpMatch.Xd**2 + llpMatch.Yd**2 + llpMatch.Zd**2)
            if R > 3870:
                continue
            displacedJets.append(jet)
    
    return displacedJets

def getModelDict(inputFile,model,verbose=True):

    if model == 'ewk':
        LLP = 1000022
        LSP = 1000024
    elif model == 'strong':
        LLP = 1000022
        LSP = 1000021
    elif model == 'gluino':
        LLP = 1000021
        LSP = 1000022
    elif model == 'sbottom':
        LLP = 1000005
        LSP = 1000022
    elif model == 'hs':
        LLP = 35
        LSP = 35     
    else:
        raise ValueError("Unreconized model %s" %model)

    modelInfoDict = {}
    f = inputFile
    if not os.path.isfile(f):
        print('File %s not found' %f)
        raise OSError()
    parsDict = {}    
    for banner in glob.glob(os.path.join(os.path.dirname(f),'*banner*txt')):
        with open(banner,'r') as ff:
            slhaData = ff.read().split('<slha>')[1].split('</slha>')[0]
            slhaData = pyslha.readSLHA(slhaData)
    parsDict = {}
    parsDict['mLLP'] = slhaData.blocks['MASS'][LLP]
    parsDict['mLSP'] = slhaData.blocks['MASS'][LSP]
    parsDict['width'] = slhaData.decays[LLP].totalwidth
    if parsDict['width']:
        parsDict['tau_ns'] = (6.582e-25/parsDict['width'])*1e9
    else:
        parsDict['tau_ns'] = np.inf    

    modelInfoDict.update(parsDict)
    if verbose:
        print('mLLP = ',parsDict['mLLP'])
        print('width (GeV) = ',parsDict['width'])
        print('tau (ns) = ',parsDict['tau_ns'])

    return modelInfoDict

def splitModels(inputFiles,model='ewk'):

    # First extract the model for each file
    fileModels = {}
    for f in inputFiles:
        modelDict = getModelDict(f,model=model,verbose=False)
        fileModels[f] = modelDict

    uniqueModels = []
    for mDict in fileModels.values():
        if mDict not in uniqueModels:
            uniqueModels.append(mDict)
    uniqueModels = sorted(uniqueModels, key= lambda d: list(d.values()))

    print('Splitting %i files into %i models\n\n\n' %(len(inputFiles),len(uniqueModels)))
    # Now group together files with the same model
    # Get list of unique dictionaries
    for mDict in uniqueModels:        
        fileList = [f for f in fileModels if fileModels[f] == mDict]
        # Sort files by modified time, so older comes first
        fileList.sort(key=os.path.getmtime)
        yield (fileList,mDict)
