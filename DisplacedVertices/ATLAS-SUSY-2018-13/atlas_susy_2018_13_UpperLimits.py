#!/usr/bin/env python3

import pandas as pd
import time,os,sys
import pandas as pd


def computeULs(inputFile,outputFile):


    # ATLAS limits:
    atlasUL = {'HighPT' : {'S95_obs' : 3.8, 'S95_exp' : 3.1},
               'Trackless' : {'S95_obs' : 3.0, 'S95_exp' : 3.4}}

    # ### Load Recast Data
    recastDF = pd.read_pickle(inputFile)

    # Set r-value
    robs = []
    rexp = []
    # Set r-value
    for _,row in recastDF.iterrows():
        S95obs = atlasUL[row['SR']]['S95_obs']
        S95exp = atlasUL[row['SR']]['S95_exp']
        robs.append(row['$N_s$']/S95obs)
        rexp.append(row['$N_s$']/S95exp)

    recastDF['robs'] = robs
    recastDF['rexp'] = rexp

    # Store all data to the output file
    recastDF.to_pickle(outputFile)



if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser(description=
            "Compute the r-values for ATLAS-SUSY-2018-13 for the recast data stored in the input file. \
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





