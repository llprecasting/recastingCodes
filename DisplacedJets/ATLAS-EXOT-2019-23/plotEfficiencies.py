#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import uproot
import logging


FORMAT = '%(levelname)s: %(message)s at %(asctime)s'
logging.basicConfig(format=FORMAT,datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger()   


try:
    import mplhep as hep
    hep.style.use("ATLAS")
except:
    logger.info("mplhep not found. Falling back to default plot style")


def getATLASdata(mPhi,mS):

    massPair = (int(mPhi),int(mS))
    effFileDict = {(1000,275) : "./ATLAS_data/HEPData-ins2043503-v3-Figure_2e_of_Aux._Mat._1000_275.root",
                    (1000,475) : "./ATLAS_data/HEPData-ins2043503-v3-Figure_2e_of_Aux._Mat._1000_475.root",
                    } #HEP data files

    effBranchDict = {(1000,275) : "Figure 2e of Aux. Mat. 1000_275/Graph1D_y1;1",
                    (1000,475) : "Figure 2e of Aux. Mat. 1000_475/Graph1D_y1;1",
                    } # Efficiencies from HEP data

    limitFileDict = {(1000,275) : "./ATLAS_data/HEP_Limits/HEPData-ins2043503-v3-Figure_6f_of_Aux._Mat._1000_275.root",
                    (1000,475) : "./ATLAS_data/HEP_Limits/HEPData-ins2043503-v3-Figure_6e_of_Aux._Mat._1000_475.root",
                    }

    limitBranchExpDict = {(1000,275) : "Figure 6f of Aux. Mat./Graph1D_y1;1",
                        (1000,475) : "Figure 6e of Aux. Mat./Graph1D_y1;1",
                    }

    limitBranchObsDict = {(1000,275) : "Figure 6f of Aux. Mat./Graph1D_y2;1",
                        (1000,475) : "Figure 6e of Aux. Mat./Graph1D_y2;1",
                    }
    
    if massPair not in effFileDict:
        logger.debug(f"mPhi = {mPhi} and mS = {mS} not found in HEP data")
        return None,None,None

    HEP = effFileDict[massPair]
    file_HEP = uproot.open(HEP) # open the file from HEP data for the efficiency
    Branch_HEP = effBranchDict[massPair]
    data_HEP = file_HEP[Branch_HEP] # open the branch

    File_HEP_limit = limitFileDict[massPair]
    file_HEP_limit = uproot.open(File_HEP_limit) # open the file from HEP data for the limits
    Branch_HEP_limit = limitBranchExpDict[massPair]
    branch_HEP_limit_exp = file_HEP_limit[Branch_HEP_limit] # open the branch
    Branch_HEP_limit = limitBranchObsDict[massPair]
    branch_HEP_limit_obs = file_HEP_limit[Branch_HEP_limit]

    return data_HEP, branch_HEP_limit_exp, branch_HEP_limit_obs



def plotEffs(effFile,mPhi,mS,outFile,sr='high-ET'):

    data_HEP, _, _ = getATLASdata(mPhi,mS)
    
    data_recast = np.genfromtxt(effFile,skip_header=2,delimiter=',',names=True)
    tauN = data_recast['ctaum']
    eff = data_recast['eff'+sr.replace('-','')]

    _, ax = plt.subplots()

    ################## Plot efficiency from MG+Pythia8 ##################
    plt.plot(tauN,eff, 'r', linewidth=2, label = 'MG + Pythia')

    ################## Plot efficiency from HEP data ##################
    if data_HEP is not None:
        plt.plot(data_HEP.values(axis='both')[0],data_HEP.values(axis='both')[1], 'b')

        ################ Uncertainties from HEP ##################
        plt.fill_between(data_HEP.values(axis='both')[0], data_HEP.values(axis='both')[1] +  data_HEP.errors('high')[1] , data_HEP.values(axis='both')[1] - data_HEP.errors('high')[1] , color = 'blue', label = r'ATLAS, with $\pm$ 1 $\sigma$ error bands',alpha=.7)

    ################## Uncertainties from Map ##################
    plt.fill_between(tauN, np.array(eff) + 0.25* np.array(eff), np.array(eff) - 0.25 * np.array(eff), label='MG+Pythia8, with error bands ', alpha=.7)

    ################## Limits of validity ##################
    ax.hlines(y=(0.25*(max(eff))), xmin=0, xmax=1e2, linewidth=2, color='g', label = 'Limits of validity' )

    # place a text box in upper left in axes coords
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax.text(0.05, 0.95, f" $ m_Φ $ = {mPhi} GeV, $m_S$ = {mS} GeV" , transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    x = np.linspace(0,100)
    ax.fill_between(x, 0.25*(max(eff)), color='black', alpha=.2, hatch="/", edgecolor="black", linewidth=1.0) # adding hatch
    plt.ylim(0) # start at 0

    plt.xscale('log')
    plt.xlabel(r'c$\tau$ [m]', fontsize=20)
    plt.ylabel('Efficiency', fontsize=20 )
    plt.legend(fontsize = 11, loc=1) # set the legend in the upper right corner
    plt.savefig(outFile)


def plotXsecLimit(effFile,mPhi,mS,outFile,sr='high-ET',factor=1):

    _, branch_HEP_limit_exp, branch_HEP_limit_obs = getATLASdata(mPhi,mS)
        
    data_recast = np.genfromtxt(effFile,skip_header=2,delimiter=',',names=True)
    tauN = data_recast['ctaum']
    eff = data_recast['eff'+sr.replace('-','')]
    
    
    _, ax = plt.subplots()

    # nsUL_obs = 0.5630 * 26 * factor (outdated value)
    # Upper limits computed using abcd_pyhf and the numbers from Table 4 in 2203.01009
    if sr == 'high-ET':
        nsUL_obs = 27.44 * factor
        nsUL_exp = 17.31 * factor
    else:
        nsUL_obs = 32.98 * factor
        nsUL_exp = 21.51 * factor

    Crr_Sec_obs = (nsUL_obs)/((np.array(eff)) * 139e3 ) # Luminosity = 139e3 fb**(-1)
    Crr_Sec_exp = (nsUL_exp)/((np.array(eff)) * 139e3 ) # Luminosity = 139e3 fb**(-1)

    plt.plot(tauN, Crr_Sec_obs, 'r', label ='Observed (Recast)', linewidth = 2)
    # plt.fill_between(tauN, Crr_Sec_obs/1.1,Crr_Sec_obs/0.9, label=r'Observed (Recast) $\pm 10$%', 
                    #  alpha=0.4, color='r')
    plt.plot(tauN, Crr_Sec_exp, 'r', label ='Expected (Recast)', linewidth = 2, linestyle='dashed')

    if branch_HEP_limit_exp is not None:
        plt.plot(branch_HEP_limit_exp.values(axis='both')[0], branch_HEP_limit_exp.values(axis='both')[1], 'b', label ='Expected (HEPData)', linewidth = 2, linestyle='dashed')
        plt.plot(branch_HEP_limit_obs.values(axis='both')[0], branch_HEP_limit_obs.values(axis='both')[1], 'b', label ='Observed (HEPData)', linewidth = 2)
    

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'c$\tau$ [m]')
    plt.ylabel(r'95% CL limit on $\sigma \times B$ [pb]')
    plt.xlim(0.1,100)
    plt.ylim(1e-4,100)

    # place a text box in upper left in axes coords
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax.text(0.05, 0.95, f" $ m_Φ $ = {mPhi} GeV, $m_S$ = {mS} GeV" , transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    plt.legend( fontsize = 10, loc=3)
    plt.grid()
    plt.savefig(outFile)




if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser(description=
            "Plot the efficiencies for a given input efficiency file generated with getEfficiencies.py" )
    ap.add_argument('-f', '--inputfile',
            help='path to the efficiencies file.')
    ap.add_argument('-sr', '--signalregion', default = 'high-ET',
            help='signal selection (high-ET or low-ET) [high-ET].')
    ap.add_argument('-M', '--mPhi', default = 1000.0,
            help='parent (s-channel) mass [1000.0].')
    ap.add_argument('-m', '--mS', default = 275.0,
            help='LLP mass [275.0].')
    ap.add_argument('-e', '--effPlot', default = None,
            help='File name for the efficiency plot. If not defined, it will not be plotted.')
    ap.add_argument('-x', '--xsecPlot', default = None,
            help='File name for the xsection limit plot. If not defined, it will not be plotted.')
    ap.add_argument('-v', '--verbose', default='info',
            help='verbose level (debug, info, warning or error). Default is info')
    

    args = ap.parse_args()

    level = args.verbose
    levels = { "debug": logging.DEBUG, "info": logging.INFO,
               "warn": logging.WARNING,
               "warning": logging.WARNING, "error": logging.ERROR }
    if level in levels:       
        logger.setLevel(level = levels[level])


    if args.effPlot is not None:
        plotEffs(args.inputfile,args.mPhi,args.mS,outFile=args.effPlot,sr=args.signalregion)
        logger.info(f"Efficiency plot saved to {args.effPlot}")
    if args.xsecPlot is not None:
        plotXsecLimit(args.inputfile,args.mPhi,args.mS,outFile=args.xsecPlot,sr=args.signalregion)
        logger.info(f"Cross-section limit plot saved to {args.xsecPlot}")
