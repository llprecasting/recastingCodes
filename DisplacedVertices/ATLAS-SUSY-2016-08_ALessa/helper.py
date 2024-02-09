#!/usr/bin/env python3

import os,glob
from typing import Any
import numpy as np
import glob
import itertools


class LLP(object):
    """
    Convenience class for holding the LLP, its daughters and several
    other useful methods
    """

    def __init__(self, candidate, direcdaughters, finaldaughters,
                 maxMomViolation=5e-2) -> None:
        self._candidate = candidate
        self.directDaughters = direcdaughters[:]
        self.finalDaughters = finaldaughters[:]
        self._selectedDecays = []
        self.nTracks = None
        self.mDV = None

        for d in self.finalDaughters:
            if d.Charge == 0:
                continue
            if d.Status != 1:
                continue
            if d.PT < 1.0:
                continue
            pTratio = abs(d.PT/d.Charge)
            if pTratio < 1.0:
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
        pTot = np.array([self.E,self.Px,self.Py,self.Pz])
        pNorm = np.linalg.norm(pTot)
        for d in self.directDaughters:                
            pTot -= np.array([d.E,d.Px,d.Py,d.Pz])
        if np.linalg.norm(pTot)/pNorm > maxMomViolation/1e3: # Be more strict about direct daughters
            raise ValueError("Error getting direct daughters, momentum conservation violated!")
      
        pTot = np.array([self.E,self.Px,self.Py,self.Pz])
        pNorm = np.linalg.norm(pTot)
        for d in self.finalDaughters:                
            pTot -= np.array([d.E,d.Px,d.Py,d.Pz])
        if np.linalg.norm(pTot)/pNorm > maxMomViolation:
            raise ValueError("Error getting final daughters, momentum conservation violated! (%s)" %(str(pTot)))        
        
        rList = [np.sqrt(d.X**2 + d.Y**2 + d.Z**2) for d in self.directDaughters]

        if max(rList) != min(rList) and abs(max(rList)-min(rList))/(max(rList)+min(rList)) > 0.01:
            raise ValueError("Direct daughters do not have the same production vertex!")
        
        daughter = self.directDaughters[0]
        self.Xd = daughter.X
        self.Yd = daughter.Y
        self.Zd = daughter.Z
        self.r_decay = np.sqrt(self.Xd**2+self.Yd**2)
        self.z_decay = self.Zd

        trimom = np.sqrt(self.Px**2 +self.Py**2 + self.Pz**2)
        self.beta = trimom/self.E
        self.gbeta = trimom/self.Mass


    def __getattr__(self, attr: str) -> Any:
        try:
            return self.__getattribute__(attr)
        except AttributeError:
            try:
                return self._candidate.__getattribute__(attr)
            except:
                raise AttributeError("Could not get attribute %s" %attr)


