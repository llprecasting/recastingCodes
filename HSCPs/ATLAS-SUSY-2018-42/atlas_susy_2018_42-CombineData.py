#!/usr/bin/env python3


import pandas as pd


def combineRecastData(files,outputFile):

    # Get files
    if not files:
        print('No valid files found')
    else:
        print('Combining %i files' %len(files))

    allData = pd.read_pickle(files[0])
    for f in files[1:]:
        recastData = pd.read_pickle(f)
        allData = pd.concat((allData,recastData),ignore_index=True)

        
    allColumns = allData.columns.tolist()
    orderColumns = ['mLLP','width','tau_ns','Regime','Target Mass [GeV]', 'SR'] 
    orderColumns += [c for c in allColumns if not c in orderColumns]
    allData = allData[orderColumns]
    allData.sort_values(['mLLP','width','tau_ns','Regime','Target Mass [GeV]', 'SR'],inplace=True,
                ascending=[True,False,True,True,True,True],ignore_index=True)        

    allData.to_pickle(outputFile)



if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser(description=
            "Merge individual DataFrames for model points to a single DataFrame (pickle file).")
    ap.add_argument('-f', '--inputFiles', required=True,nargs='+',
            help='list of pickle files to be merged', default =[])
    ap.add_argument('-o', '--outputFile', required=True, help='output file.')


    args = ap.parse_args()
    inputFiles = args.inputFiles
    outputFile = args.outputFile
    
    combineRecastData(inputFiles,outputFile)







