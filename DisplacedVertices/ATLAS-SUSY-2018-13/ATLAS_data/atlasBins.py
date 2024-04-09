#!/usr/bin/env python3

import numpy as np
import os

atlasDir = os.path.dirname(os.path.abspath(__file__))

atlas_bins = {}

# High PT:
data = np.genfromtxt(os.path.join(atlasDir,'HEPData-ins2628398-v1-csv/yields_highpt_SR_observed.csv'),comments='#',delimiter=',',skip_header=10)
nTrack_bins = np.append(data[:,1],data[-1,2])
nTrack_bins = np.sort(np.unique(nTrack_bins))
mDV_bins = np.append(data[:,4],data[-1,5])
mDV_bins = np.sort(np.unique(mDV_bins))
n_obs = data[:,-1]

atlas_bins['HighPT'] = {'mDV' : mDV_bins, 'nTracks' : nTrack_bins, 'nobs' : n_obs}



# Trackless:
data = np.genfromtxt(os.path.join(atlasDir,'./HEPData-ins2628398-v1-csv/yields_trackless_SR_observed.csv'),comments='#',delimiter=',',skip_header=10)
nTrack_bins = np.append(data[:,1],data[-1,2])
nTrack_bins = np.sort(np.unique(nTrack_bins))
mDV_bins = np.append(data[:,4],data[-1,5])
mDV_bins = np.sort(np.unique(mDV_bins))
n_obs = data[:,-1]

atlas_bins['Trackless'] = {'mDV' : mDV_bins, 'nTracks' : nTrack_bins, 'nobs' : n_obs}