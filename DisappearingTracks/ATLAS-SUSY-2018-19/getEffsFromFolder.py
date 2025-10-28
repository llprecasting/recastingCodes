#!/usr/bin/env python3
import numpy as np
import os,sys,glob
from getEffsFromROOT import logger,main
import subprocess
import logging
import pyslha
import multiprocessing


def getModelInfo(bannerFile,llpPDG):

    with open(bannerFile,'r') as ff:
        slhaData = ff.read()
        if '<slha>' in slhaData:
            slhaData = slhaData.split('<slha>')[1].split('</slha>')[0]
            
        slhaData = pyslha.readSLHA(slhaData)
    
    if slhaData is None:
       logger.error(f'Error reading banner file {bannerFile}')
       raise ValueError()
    
    modelInfoDict = {}
    modelInfoDict['mLLP'] = slhaData.blocks['MASS'][llpPDG]
    modelInfoDict['tau0_ns'] = (6.58212e-16/slhaData.decays[llpPDG].totalwidth)
    
    return modelInfoDict

if __name__ == "__main__":
      
  import argparse
  parser = argparse.ArgumentParser(description='Analyse delphesLLP output to produce efficiencies for ATLAS-SUSY-2019-18 DT search')
  parser.add_argument('-F','--folder', help='Path to the folder containing Delphes ROOT files.')
  parser.add_argument('-l','--llpPDG',help='LLP PDG [1000024]',type=int, required=False, default=1000024)
  parser.add_argument('-tauF','--tau_file',metavar='tau_file', help='CSV file containing the lifetime values (in ns) used for reweighting',
                      type=str, required=False, default=None)
  parser.add_argument('-n', '--ncpus',type=int,default=1,help='number of parallel jobs to run.')
  parser.add_argument('-v', '--verbose', default='info',
                      help='verbose level (debug, info, warning or error). Default is info')

  args = parser.parse_args()


  # First make sure the correct env variables have been set:
  LDPATH = subprocess.check_output('echo $LD_LIBRARY_PATH',shell=True,text=True)
  ROOTINC = subprocess.check_output('echo $ROOT_INCLUDE_PATH',shell=True,text=True)
  pythiaDir = os.path.abspath('./MG5/HEPTools/pythia8/lib')
  delphesDir = os.path.abspath('./DelphesLLP/external')
  if pythiaDir not in LDPATH or delphesDir not in ROOTINC:
      print('Enviroment variables not properly set. Run source setenv.sh first.')
      sys.exit()

  level = args.verbose
  levels = { "debug": logging.DEBUG, "info": logging.INFO,
              "warn": logging.WARNING,
              "warning": logging.WARNING, "error": logging.ERROR }
  if level in levels:       
      logger.setLevel(level = levels[level])

  if not os.path.isdir(args.folder):
     logger.error(f"Folder {args.folder} not found!")
     raise ValueError()
  
  # Find root files:
  pattern = os.path.join(args.folder, "**", f"*.root")
  found_files = list(glob.glob(pattern, recursive=True))
  if not found_files:
     logger.error(f"No .root files found in {args.folder}!")
     raise ValueError()
  else:
     logger.info(f"Found {len(found_files)} files")

  banner_files = []
  for f in found_files:
    d = os.path.dirname(f)
    b_files = list(glob.glob(os.path.join(d,'*_banner.txt')))
    if not b_files:
      logger.error(f"No banner files found in {d}!")
      raise ValueError()
    else:
      banner_files.append(b_files[0])
  
  ncpus = min(len(found_files),args.ncpus)
  pool = multiprocessing.Pool(processes=ncpus)
  children = []
  ijob = 0
  for rootFile,bannerFile in zip(found_files,banner_files):
    tau0 = getModelInfo(bannerFile,args.llpPDG)['tau0_ns']
    p = pool.apply_async(main, args=(rootFile,tau0,args.tau_file,ijob,))
    ijob += 1
    children.append(p)

  logger.info(f'Running {ijob} jobs in {ncpus} instances')
  for p in children: 
    p.get()

  

