#!/usr/bin/env python3

import json
import argparse,os,glob

def combineFiles(effFiles,outputFile):

    combinedData = list()
    for effFile in effFiles:
        with open(effFile,'r') as f:
            combinedData.append(json.load(f))

    with open(outputFile,'w') as f:
        json.dump(combinedData,f,indent=4)

if __name__ == "__main__":

    ap = argparse.ArgumentParser(description=
            "Merge individual json files generated with computeEfficiencies.py.")
    ap.add_argument('-i','--input', help='Path to the folder containing the json files.')

    ap.add_argument('-o', '--outputFile', required=False, help='output file [atlas_susy_2018_14_effs.json].',
                    default='atlas_susy_2018_14_effs.json')


    args = ap.parse_args()
    inputFiles = []
    if os.path.isdir(args.input):
        # Find root files:
        pattern = os.path.join(args.input, "**", f"*.json")
        inputFiles = list(glob.glob(pattern, recursive=True))
        if not inputFiles:
            print(f"No .json files found in {args.input}!")
            raise ValueError()
    else:
        print(f"Folder {args.input} not found!")
        raise ValueError()

    combineFiles(inputFiles,args.outputFile)







