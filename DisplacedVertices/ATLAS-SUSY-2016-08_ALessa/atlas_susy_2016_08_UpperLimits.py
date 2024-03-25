#!/usr/bin/env python3

import pandas as pd
import time,os,sys
import pandas as pd
import numpy as np

def computeULs(inputFile,outputFile):

    # ### Load Recast Data
    recastData = pd.read_pickle(inputFile)

    # First compute ULs for each individual SR
    # ATLAS limits:
    atlasUL = {'S95_obs' : 3.0, 'S95_exp' : 3.0} # Limits are the same for mDV > 5 GeV, nTracks >=5 (see getUL.ipynb)    
    # Set r-value
    robs = []
    rexp = []
    # Set r-value
    for _,row in recastData.iterrows():
        S95obs = atlasUL['S95_obs']
        S95exp = atlasUL['S95_exp']
        robs.append((row['$N_s$']/S95obs,row['$N_s$ Err']/S95obs))
        rexp.append((row['$N_s$']/S95exp,row['$N_s$ Err']/S95exp))

    robs = np.array(robs)
    rexp = np.array(rexp)
    recastData['robs'] = robs[:,0]
    recastData['rexp'] = rexp[:,0]
    recastData['robsErr'] = robs[:,1]
    recastData['rexpErr'] = rexp[:,1]



    # Store all data to the output file
    recastData.to_pickle(outputFile)

   

if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser(description=
            "Compute the r-values for ATLAS-SUSY-2016-08 for the recast data stored in the input file. \
            The points without any signal will be removed from the output.")
    ap.add_argument('-f', '--inputFile', required=True,
            help='path to the pickle file containing the Pandas DataFrame with the recasting results for the models')
    ap.add_argument('-o', '--outputFile', required=False,
            help='path to output file. If not defined the upper limits will be stored in the input file.',
            default = None)

    t0 = time.time()

    # # Set output file
    args = ap.parse_args()
    inputFile = args.inputFile
    if not os.path.isfile(inputFile):
        print("File %s not found" %inputFile)
        sys.exit()
    outputFile = args.outputFile
    if outputFile is None:
        outputFile = inputFile

    computeULs(inputFile,outputFile)
    
    
    print("\n\nDone in %3.2f min" %((time.time()-t0)/60.))





