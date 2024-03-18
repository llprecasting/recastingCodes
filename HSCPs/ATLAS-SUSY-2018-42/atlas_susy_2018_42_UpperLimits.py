#!/usr/bin/env python3

import pandas as pd
import time,os,sys
import pandas as pd


def computeULs(inputFile,outputFile):


    # ### Load ATLAS data
    atlasFile = './ATLAS_data/modelIndependentLimits.pcl'
    atlasDF = pd.read_pickle(atlasFile)

    # ### Load Recast Data
    recastDF = pd.read_pickle(inputFile)

    # Merge data frames based on Target Mass and SR
    allDF = pd.merge(recastDF,atlasDF,on=['Target Mass [GeV]','SR', 'Regime'],how="inner",
                     suffixes=(None,"del_me"))

    allDF.drop(columns=[c for c in allDF.columns if "del_me" in c],inplace=True)
    allDF['robs'] = allDF['$N_s$']/allDF['S95_obs']
    allDF['rexp'] = allDF['$N_s$']/allDF['S95_exp']

    # Filter data, so only the highest robs is kept
    mLLP = 'mLLP'
    mn1= 'mLSP'
    tau = 'tau_ns'
    mCols = [mLLP,mn1,tau,'robs']
    allDF.sort_values(mCols,inplace=True,ascending=True)
    recastDF2 = allDF.drop_duplicates(subset=[mLLP,mn1,tau],keep="last")

    # Store all data to the output file
    recastDF2.to_pickle(outputFile)



if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser(description=
            "Compute the upper limit on mu for ATLAS-SUSY-2018-42 for the recast data stored in the input file. \
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





