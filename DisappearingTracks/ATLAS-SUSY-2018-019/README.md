# CalRatio Recast ([ATLAS-SUSY-2018-19](https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-19/))


## Authors: ##
[Lucas Magno](mailto:lucas.magno.ramos@usp.br) and [Andre Lessa](mailto:andre.lessa@ufabc.edu.br)

The recast code and results are based on the [recast note](https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-19/hepdata_info.pdf) and the [SimpleAnalysis code](https://gitlab.cern.ch/atlas-sa/simple-analysis/-/blob/master/SimpleAnalysisCodes/src/ANA-SUSY-2018-19_TrackletAcc.cxx).


Validation of the results can be found in the [validation folder](./validation).

**Important Note**: *the recasting code requires the pair production of LLPs, where each LLP decays to a pair of fermions (or a pair of fermions plus invisible particles) or to a Higgs (or a Higgs plus invisible particles). For the latter, the Higgs is assumed to decay to a pair of b quarks.*

## Pre-Requisites and Installation ##

The following pre-requisites must be installed before running the main code:

  * [DelphesLLP](https://github.com/llprecasting/recastingCodes/tree/main/Delphes_LLP)
  * [pyROOT](https://root.cern.ch/)


## Running the recast code

For running the recast code or reproducing the results in the [validation folder](./validation/) one must:

 1. Generate HepMC events with MadGraph5 plus Pythia8.
 2. Run Delphes on the HepMC events using a [modified version of Delphes](https://github.com/llprecasting/recastingCodes/tree/main/Delphes_LLP), which stores LLPs and their decays in the output ROOT file.

Examples of cards for generating events for the wino model can be found in the [validation/Cards](./validation/Cards/) folder.


After generating a ROOT file using the modified Delphes code (see above),  the efficiencies can be computed running:

```
./getEffsFromROOT.py -p <parameters_file> -f <path-to-root-file(s)>
```
An example of the parameters file can be found [here](./parameters_getEff.ini).
This file defines the PDGs for the LLP and what should be considered as invisible in the LLP decays.
It also provides a list of $c\tau$ values for which the efficiencies will be computed and other options.
The output will consist of a csv file (stored in the same folder as the input file) for each ROOT file containing the efficiencies for
the provided $c\tau$ values.

If running over multiple files, the code allows to run the efficiency calculation in parallel (the number of parallel runs is set by the `ncpus` option).

## Plotting the results

Validation plots similar to the ones in [arXiv:2412.13976](https://arxiv.org/pdf/2412.13976) can be generated running:

```
./plotEfficiencies.py -f INPUTFILE -M MPHI -m MS -e EFFPLOT -x XSECPLOT
```

where INPUTFILE is the one of the csv files containing the signal efficiencies, (MPHI,MS) are
the corresponding ( $m_{\Phi}$ , $m_{S}$ ) mass values assumed by ATLAS (for comparison against the official ATLAS curves) and EFFPLOT and XSECPLOT are the names of the efficiency plot and cross-section upper limit plot, respectively.


## Validation

More information about validation can be found in the [validation](./validation/) folder and in the [recast note](https://arxiv.org/pdf/2412.13976).
