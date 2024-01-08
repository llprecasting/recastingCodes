#!/usr/bin/env python3

import os,glob
from typing import Any
import numpy as np
import pandas as pd
import glob
import pyslha
import time



class LLP(object):
    """
    Convenience class for holding the LLP, its daughters and several
    other useful methods
    """

    def __init__(self, candidate, daughters=[]) -> None:
        self._candidate = candidate
        self.directDaughters = []
        self.finalDaughters = []
        self._selectedDecays = []
        self._nTracks = None
        self._mDV = None

        if daughters:
            for d in daughters:
                if d.M1 < 0 and d.Status == 1:
                    self.finalDaughters.append(d)
                    if d.Charge == 0:
                        continue
                    if d.Statue != 1:
                        continue
                    pTratio = abs(d.PT/d.Charge)
                    if pTratio < 1.0:
                        continue
                    self._selectedDecays.append(d)
                else:
                    self.directDaughters.append(d)

            # Consistency checks:
            pTot = np.array([self.E,self.Px,self.Py,self.Pz])
            for d in self.directDaughters:
                pTot -= np.array([d.E,d.Px,d.Py,d.Pz])
            if np.linalg.norm(pTot) < 1e-4:
                raise ValueError("Error getting direct daughters, momentum conservation violated!")

            pTot = np.array([self.E,self.Px,self.Py,self.Pz])
            for d in self.finalDaughters:
                pTot -= np.array([d.E,d.Px,d.Py,d.Pz])
            if np.linalg.norm(pTot) < 1e-4:
                raise ValueError("Error getting final daughters, momentum conservation violated!")
            
            rList = [np.sqrt(d.X**2 + d.Y**2 + d.Z**2) for d in self.directDaughters]
            if max(rList)/min(rList) > 1.001:
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



    def __getattribute__(self, attr: str) -> Any:
        if hasattr(self,attr):
            return getattr(self,attr)
        else:
            return getattr(self._candidate,attr)
        
    @property
    def mDV(self):
        if self._mDV is not None:
            return self._mDV
        else:
            pTot = np.zeros(4)
            for d in self._selectedDecays:
                p = np.array([0.,d.Px,d.Py,d.Pz])
                mpion = 0.140
                p[0] = np.sqrt(mpion**2 + np.dot(p[1:],p[1:]))
                pTot += p
            mDV = np.sqrt(pTot[0]**2 - np.dot(pTot[1:],pTot[1:]))
            self._mDV = mDV
            return mDV
        
    @property
    def nTracks(self):
        return len(self._selectedDecays)
