#!/usr/bin/env python3

import os, time, sys
import numpy as np
import readMapNew as rmN
import pyhepmc
import random
import logging
import configparser
import multiprocessing
import progressbar as P

c = 3e8

FORMAT = '%(levelname)s: %(message)s at %(asctime)s'
logging.basicConfig(format=FORMAT,datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger()   


random.seed(123)

def getLLPdecays(event,llps=[35]):

    vertexDict = {}
    for ivertex,vertex in enumerate(event.vertices):
        if len(vertex.particles_in) != 1:
            continue
        if len(vertex.particles_out) == 1:
            continue
        part_in = vertex.particles_in[0]
        # Check if incoming particle appear in pdgList
        if abs(part_in.pid) not in llps:
            continue

        part_out_list = vertex.particles_out
        # Skip vertices which could correspond to FSR (initial state appears in final state)
        if any(p.pid == part_in.pid for p in part_out_list):
            continue

        vertexDict[ivertex] = vertex
    
    return vertexDict

def getDataFrom(event,llps=[35],invisibles=[12,14,16]):

    vertexDict = getLLPdecays(event,llps=llps)
    eventDict = {}
    if len(vertexDict) != 2:
        erroMsg = f"{len(vertexDict)} LLP decay vertices found (can only handle 2 decay vertices)"
        logger.error(erroMsg)
        raise ValueError(erroMsg)
    for ivertex,vertex in vertexDict.items():
        visible_particles = [p for p in vertex.particles_out 
                                if abs(p.pid) not in invisibles]
        invisible_particles = [p for p in vertex.particles_out 
                                if abs(p.pid) in invisibles]
        # In case the LLP decays to a single visible particle (e.g. higgs), replace it
        # by its visible decays (e.g. higgs -> b bar)
        if len(visible_particles) == 1:
            p = visible_particles[0]
            if len(p.children) != 2:
                erroMsg = f"Single visible particle from LLP decay with {len(p.children)} daughters (can only handle 2)"
                logger.error(erroMsg)
                raise ValueError(erroMsg)
            else:
                visible_particles = [d for d in p.children
                                     if abs(d.pid) not in invisibles]
        if len(visible_particles) != 2:
            errorMsg = f"{len(visible_particles)} visible particles found in LLP decay (can only handle 2 particles)"
            logger.error(errorMsg)
            raise ValueError(errorMsg)
        p_visible = pyhepmc.FourVector(0.,0.,0.,0.)
        p_invisible = pyhepmc.FourVector(0.,0.,0.,0.)
        for p in visible_particles:
            p_visible = p_visible + p.momentum
        for p in invisible_particles:
            p_invisible = p_invisible + p.momentum

        visParticle = pyhepmc.GenParticle(momentum = p_visible, pid = ivertex)
        invisParticle = pyhepmc.GenParticle(momentum = p_invisible, pid = ivertex)
        eventDict[ivertex] = {'visible' : visParticle, 
                              'invisible' : invisParticle,
                              'parent' : vertex.particles_in[0],
                              'visiblePDGs' : [p.pid for p in visible_particles],
                              'invisiblePDGs' : [p.pid for p in invisible_particles],
                              }

    return eventDict

def lifetime(avgtau = 4.3):
    avgtau = avgtau / c
    t = random.random()
    return -1.0 * avgtau * np.log(t)

def getDecayLength(genParticle,tau):
    
    lt = lifetime(tau) # set mean lifetime
    gamma = genParticle.momentum.e/genParticle.momentum.m()
    L_vec = (genParticle.momentum/genParticle.momentum.e)*c*lt*gamma

    return L_vec.pt(),abs(L_vec.pz)

def passHTmissCut(eventDict):

    # Compute HT as the scalar sum of the visible pT out 
    # of each LLP decay
    HT = 0.0
    for vertex in eventDict.values():
        HT += vertex['visible'].momentum.pt()

    # Compute HTmiss as the total invisible pT from both decays
    pMiss = pyhepmc.FourVector(0.,0.,0.,0.)
    for vertex in eventDict.values():
        pMiss = pMiss + vertex['invisible'].momentum
    HTmiss = pMiss.pt()

    if (not HT) or HTmiss/HT > 0.6:
        return False
    
    return True

def getEffFor(event,tauList,llps,invisibles,applyHTcut):


    evt_effs = {"low-ET" : np.zeros(len(tauList)),
                "high-ET"  : np.zeros(len(tauList))}
    # Extract necessary data from event
    eventDict = getDataFrom(event,llps=llps,invisibles=invisibles)
    if len(eventDict) != 2:
        errorMsg = f"{len(eventDict)} LLP decay vertices found (can only handle 2 decay vertices)"
        logger.error(errorMsg)
        # raise ValueError(errorMsg)
        return evt_effs # return zero effs

    if applyHTcut and (not passHTmissCut(eventDict)):
        return evt_effs # return zero effs
    
    llpList = [d['parent'] for d in eventDict.values()]
    visList = [d['visible'] for d in eventDict.values()]
    visPDGList = [d['visiblePDGs'] for d in eventDict.values()]
    
    # Get total momentum of LLP pair (equals momentum of parent)
    pTot = llpList[0].momentum + llpList[1].momentum
    if pTot.m() >= 400:
        sr = "high-ET"
    else:
        sr = "low-ET"

    
    p1 = visList[0].momentum        
    p1_pt = p1.pt()
    p1_eta = p1.eta()
    p1_pdgs = set([abs(pdg) for pdg in visPDGList[0]])
    if len(p1_pdgs) != 1:
        errorMsg = f'Can not handle decays to distinct pair of fermions (e.g. {p1_pdgs})'
        logger.error(errorMsg)
        # raise ValueError(errorMsg)
        return evt_effs # return zero effs
    else:
        p1_pdgs = list(p1_pdgs)[0]
    
    p2 = visList[1].momentum
    p2_pt = p2.pt()
    p2_eta = p2.eta()
    p2_pdgs = set([abs(pdg) for pdg in visPDGList[1]])
    if len(p2_pdgs) != 1:
        errorMsg = f'Can not handle decays to distinct pair of fermions (e.g. {p2_pdgs})'
        # raise ValueError(errorMsg)
        return evt_effs # return zero effs
    else:
        p2_pdgs = list(p2_pdgs)[0]
    
    
    
    for i,tau in enumerate(tauList):
        L1xy,L1z = getDecayLength(llpList[0],tau)
        L2xy,L2z = getDecayLength(llpList[1],tau)

        eff = rmN.queryMapFromKinematics(p1_pt,p1_eta,L1xy,L1z,p1_pdgs,
                                            p2_pt,p2_eta,L2xy,L2z,p2_pdgs,
                                            selection = sr)
        evt_effs[sr][i] = eff
        
    return evt_effs


def getEfficiencies(inputFile,tauList,
                    llps=[35],invisibles=[12,14,16],applyHTcut=True):

    # Load hepmc File
    if not os.path.isfile(inputFile):
        logger.error(f"File {inputFile} not found")
        return None
    
    try:
        f = pyhepmc.open(inputFile) # Open HEPMC file
    except:
        logger.error(f"Error reading {inputFile}")
        return None
    
    # Total efficiencies (for each tau value)
    effsDict = {"high-ET" : np.zeros(len(tauList)),
                "low-ET" : np.zeros(len(tauList)),
                "Nevents" : 0}
    
    # Get next event
    event = f.read()
    while event is not None:
        effsDict['Nevents'] += 1
        evt_effs = getEffFor(event,tauList,llps,invisibles,applyHTcut)
        # Add event efficiency to total efficiencies 
        # #for the given selection (for each tau value)
        for sr in evt_effs:
            effsDict[sr] += np.array(evt_effs[sr])
        
        event  = f.read()
        
    f.close()
    # Divide the total efficiency by the number of events:
    effsDict["high-ET"] = effsDict["high-ET"]/effsDict['Nevents']
    effsDict["low-ET"] = effsDict["low-ET"]/effsDict['Nevents']
    # Store ctau list
    effsDict['ctau'] = tauList[:]
    # Store input file name
    effsDict['inputFile'] = inputFile

    return effsDict

def saveOutput(effsDict,outputFile):
        
    tauList = effsDict['ctau']
    data = np.array(list(zip(tauList,effsDict['low-ET'],
                                effsDict['high-ET'])))

    
    np.savetxt(outputFile, data, 
                header=f'Input file: {effsDict['inputFile']}\nNumber of events: {effsDict['Nevents']}\nctau(m),eff(low-ET),eff(high-ET)',
                delimiter=',',fmt='%1.3e')
    


if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser( description=
            "Compute the efficiencies for a given input HepMC file" )
    ap.add_argument('-f', '--inputfile',nargs='+',
            help='path to the HepMC file(s).')
    ap.add_argument('-p', '--parfile',default='parameters_getEff.ini',
            help='parameters file name, where additional options are defined [parameters_getEff.ini].')
    ap.add_argument('-v', '--verbose', default='info',
            help='verbose level (debug, info, warning or error). Default is info')
    

    args = ap.parse_args()

    level = args.verbose
    levels = { "debug": logging.DEBUG, "info": logging.INFO,
               "warn": logging.WARNING,
               "warning": logging.WARNING, "error": logging.ERROR }
    if level in levels:       
        logger.setLevel(level = levels[level])



    parser = configparser.ConfigParser(inline_comment_prefixes="#")   
    ret = parser.read(args.parfile)
    if ret == []:
        logger.error(f"No such file or directory: {args.parfile}")
        sys.exit()

    
    try:
        tauList = parser.get("options","ctau_values")
        tauList = eval(tauList)
    except:
        logger.error("Error getting list of ctau values")
        tauList = np.geomspace(args.tmin,args.tmax,args.ntau)

    try:
        llpPDGs = parser.get("model","llp_pdgs").split(',')
        llpPDGs = [int(pdg) for pdg in llpPDGs]
    except:
        logger.error("Error getting list of LLP PDG values")
        llpPDGs = [35]

    try:
        invisiblePDGs = parser.get("model","invisible_pdgs").split(',')
        invisiblePDGs = [int(pdg) for pdg in invisiblePDGs]
    except:
        logger.error("Error getting list of LLP PDG values")
        invisiblePDGs = [12,14,16]

    t0 = time.time()

    # Check whether to apply the HTmiss cut
    if parser.has_option("options","applyHTcut"):
        applyHTcut = bool(parser.get("options","applyHTcut"))
    else:
        applyHTcut = True

    # Check for output file suffix
    output_suffix = ""
    if parser.has_option("options","output_suffix"):
        output_suffix = parser.get("options","output_suffix")
    
    if not output_suffix:
        if applyHTcut:
            output_suffix = "_HTcut"
        else:
            output_suffix = "_noHTcut"

    # Split input files by distinct models and get recast data for
    # the set of files from the same model:
    ncpus = int(parser.get("options","ncpus"))
    inputFileList = args.inputfile
    ncpus = min(ncpus,len(inputFileList))
    pool = multiprocessing.Pool(processes=ncpus)
    children = []
    for inputfile in inputFileList:
        p = pool.apply_async(getEfficiencies, args=(inputfile,tauList,
                           llpPDGs,invisiblePDGs,applyHTcut,))
        children.append(p)

    nfiles = len(inputFileList)
    progressbar = P.ProgressBar(widgets=[f"Reading {nfiles} files", 
                                P.Percentage(),P.Bar(marker=P.RotatingMarker()), P.ETA()])
    progressbar.maxval = nfiles
    progressbar.start()
    ndone = 0
    for p in children: 
        effsDict = p.get()
        outFile = effsDict['inputFile'].split('.hepmc')[0]
        outFile = outFile + output_suffix +'_effs.csv'
        saveOutput(effsDict,outFile)
        ndone += 1
        progressbar.update(ndone)
        
    logger.info("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
