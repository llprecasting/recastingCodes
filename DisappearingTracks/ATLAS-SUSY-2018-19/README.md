# Disappearing Track Recast ([ATLAS-SUSY-2018-19](https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-19/))


## Authors: ##
[Lucas Magno](mailto:lucas.magno.ramos@usp.br) and [Andre Lessa](mailto:andre.lessa@ufabc.edu.br)

The recast code and results are based on the [recast note](https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-19/hepdata_info.pdf) and the [SimpleAnalysis code](https://gitlab.cern.ch/atlas-sa/simple-analysis/-/blob/master/SimpleAnalysisCodes/src/ANA-SUSY-2018-19_TrackletAcc.cxx).


Validation of the results can be found in the [validation folder](./validation).

**Important Note**: *the recasting code implements only the model-independent signal regions, which includes a pT > 60 GeV cut on the LLPs. As a result the efficiencies are slightly smaller than the model dependent approach used to constrain the wino scenario.*

## Pre-Requisites and Installation ##

The following pre-requisites must be installed before running the main code:

  * [DelphesLLP](https://github.com/llprecasting/recastingCodes/tree/main/Delphes_LLP)
  * [pyROOT](https://root.cern.ch/)

It might also be necessary to set environment variables, so the required libraries can be found (see [setenv.sh](./setenv.sh)).

An installation script ([installer.sh](./installer.sh)) is provided, which will try to install [MadGraph5](https://launchpad.net/mg5amcnlo) and [DelphesLLP](https://github.com/llprecasting/recastingCodes/tree/main/Delphes_LLP).
A ROOT installation must already be present in the system.

**Note**: *It is not a given that the installation script will work for all systems. Also, the MadGraph version may have to be updated. Nonetheless, it may serve as a guide for a manual installation.*

## Running the recast code


### Generating Input Files

For running the recast code or reproducing the results in the [validation folder](./validation/) one must first:

 1. Generate HepMC events with MadGraph5 plus Pythia8.
 2. Run Delphes on the HepMC events using a [modified version of Delphes](https://github.com/llprecasting/recastingCodes/tree/main/Delphes_LLP), which stores LLPs and their decays in the output ROOT file.

Examples of cards for generating events for the wino model can be found in the [validation/Cards](./validation/Cards/) folder.

#### Running a scan over model space

A scan script [runScanMG5.py](./runScanMG5.py) is provided for running the MG5+Pythia8+DelphesLLP pipeline over a set of model parameters. The output is a set of Delphes/ROOT files containing the minimal information required for
computing the efficiencies.
The usage of the scan code is:


```
./runScanMG5.py -p <parameters_file>
```

where the parameter defines all the necessary input and options. Examples of parameter files can be found [here](./validation/scan_parameters_wino_chgchg.ini) for generating events for chargino pair production and [here](./validation/scan_parameters_wino_chgn1.ini) for chargino-neutralino associated production.


### Computing efficiencies

After generating ROOT file(s) using the modified Delphes code (see above),  the efficiencies can be computed running:

```
./computeEfficiencies.py -i <input> -l <LLP PDG> -tauF <lifetime reweighting file> -n <number of parallel jobs> -v <verbose>
```

The options include:

*  `-h` or `--help` : show the help message and exit
*  `-i` or `--input` : Path to Delphes ROOT file or to a folder containing Delphes ROOT files with the event samples to be analysed.
*  `-l` or `--llpPDG` : LLP PDG [default = 1000024]
*  `-tauF` or `--tau_file` : CSV file containing the lifetime values (in ns) used for reweighting. If empty or file not found, it will not apply reweighting [default=tau_list.csv].
*  `-n` or `--ncpus` : number of parallel jobs to run when running over multiple files [default=1].
*  `-v` or `--verbose` : verbose level (debug, info, warning or error). If debug, it will also print the cutflows [default=warning].

An example of the lifetime CSV file which can be used for lifetime reweighting can be found [here](./tau_list.csv).
*One must be careful when applying the lifetime reweighting to lifetime values too far away from the one assumed for event generation,
since it may result in large uncertainties.*

If running over multiple files, the code allows to run the efficiency calculation in parallel (the number of parallel runs is set by the `ncpus` option).

For each ROOT file, the output is stored in a JSON file in the same folder containing the input ROOT file, with the suffix `_effs.json`. This file contains basic information about the input file as
well as the efficiencies for each signal region and each lifetime.
An example is shown below:

```json
{
    "totalweight": 1179718.542,
    "Nevents": 54986,
    "inputFile": "pp2chgchg/Events/run_08/wino_100GeV_10.000ns_delphes_events.root",
    "tau0_ns": 10.0,
    "llpPDG": 1000024,
    "Number of Events": 125000.0,
    "Integrated weight (pb)": 21.45489,
    "Matched Integrated weight (pb)": 9.43774971,
    "mLLP": 100.0,
    "Cross-Section (pb)": 9.4377,
    "Efficiencies": [
        {
            "tau_ns": 0.01,
            "EWK SR": 8.740e-09,
            "EWK SR Error": 7.02719e-09,
            "Strong SR": 6.0124e-15,
            "Strong SR Error": 6.01235e-15
        },
        {
            "tau_ns": 0.012,
            "EWK SR": 7.771599e-08,
            "EWK SR Error": 5.84265e-08,
            "Strong SR": 5.00076e-13,
            "Strong SR Error": 5.000444e-13
        }
    ]
}

```

These efficiencies can then be used to constrain the model.
For convenience a script ([collectData.py](./collectData.py)) is provided to collect all the results from a set of json files into a single output file, which can then be used to analyse the parameter space. The usage is:

```
./collectData.py -i <folder containing json files> -o output file [default=atlas_susy_2018_19_effs.json]
```

The output file will then contain a list of dictionaries contanining the information of each `*_effs.json` file found in the input folder (or its subfolders).


## Plotting the results

Validation plots similar to the ones in [arXiv:2412.13976](https://arxiv.org/pdf/2412.13976) can be generated running:

```
./plotEfficiencies.py -f INPUTFILE -M MPHI -m MS -e EFFPLOT -x XSECPLOT
```

where INPUTFILE is the one of the csv files containing the signal efficiencies, (MPHI,MS) are
the corresponding ( $m_{\Phi}$ , $m_{S}$ ) mass values assumed by ATLAS (for comparison against the official ATLAS curves) and EFFPLOT and XSECPLOT are the names of the efficiency plot and cross-section upper limit plot, respectively.


## Validation

More information about validation can be found in the [validation](./validation/) folder.
