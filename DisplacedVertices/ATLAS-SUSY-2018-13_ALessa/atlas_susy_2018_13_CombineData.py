#!/usr/bin/env python3


import pandas as pd


def combineRecastData(files,outputFile,axis):

    # Get files
    if not files:
        print('No valid files found')
    else:
        print('Combining %i files' %len(files))

    allData = pd.read_pickle(files[0])
    for f in files[1:]:
        recastData = pd.read_pickle(f)
        allData = pd.concat((allData,recastData),ignore_index=True,
                            axis=axis)

       

    allData.to_pickle(outputFile)



if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser(description=
            "Merge individual DataFrames for model points to a single DataFrame (pickle file).")
    ap.add_argument('-f', '--inputFiles', required=True,nargs='+',
            help='list of pickle files to be merged', default =[])
    ap.add_argument('-a', '--axis', required=False,default=0,
            help='Which axis to merge the files (set to 1 when combining cutflows)')
    
    ap.add_argument('-o', '--outputFile', required=True, help='output file.')


    args = ap.parse_args()
    inputFiles = args.inputFiles
    outputFile = args.outputFile
    
    combineRecastData(inputFiles,outputFile,args.axis)







