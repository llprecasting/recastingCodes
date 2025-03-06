# CalRatio Recast ([ATLAS-EXOT-2019-23](https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/EXOT-2019-23/))


## Authors: ##
[Louie Corpe](mailto:l.corpe@cern.ch)

[Andre Lessa](mailto:andre.lessa@ufabc.edu.br)

The recast code and results are based on [arXiv:2412.13976](https://arxiv.org/pdf/2412.13976), the auxiliary material
provided in [HepDATA](https://www.hepdata.net/record/ins2043503) and the [code developed by Louie Corpe](https://github.com/llprecasting/recastingCodes/tree/main/DisplacedJets/ATLAS-EXOT-2019-23).


Validation of the results can be found in the [validation folder](./validation)

**Important Note**: *the recasting code requires the pair production of LLPs, where each LLP decays to a pair of fermions (or a pair of fermions plus invisible particles) or to a Higgs (or a Higgs plus invisible particles). For the latter, the Higgs is assumed to decay to a pair of b quarks.*

## Pre-Requisites and Installation ##

The following pre-requisites must be installed before compiling the main code:

  * [PyYAML](https://pypi.org/project/PyYAML/)
  * [pyhepmc](https://pypi.org/project/pyhepmc/) (if using HepMC files as input)
  * [pyROOT](https://root.cern.ch/) (if using ROOT files as input)


## Running the recast code

There are two paths for running the recast code or reproducing the results in the [validation folder](./validation/):

 1. Generate HepMC events with MadGraph5 plus Pythia8, which is taken as input by [getEffsfromHepMC.py](./getEffsFromHepMC.py) to compute the efficiencies or
 2. Run MadGraph5 plus DelphesPythia8, which stores the necessary information in a ROOT file, which is taken as input by [getEffsFromROOT.py](./getEffsFromROOT.py) to compute the efficiencies.

Method 1 requires requires running only MadGraph5 and Pythia8, while Method 2 requires running in addition a [modified version of Delphes](https://github.com/llprecasting/recastingCodes/tree/main/Delphes_LLP), which stores LLPs and their decays in the output ROOT file.
However, Method 2 can be about 5x faster.

Examples of cards for generating events for the HAHM model using either of the methods can be found in the [validation/Cards](./validation/Cards/) folder.

### Method 1

After generating a HepMC file containing the events the efficiencies can be computed running:

```
./getEffsFromHepMC.py -p <parameters_file> -f <path-to-hepmc-file(s)>
```
where one example of the parameters file can be found [here](./parameters_getEff.ini).
This file defines the PDGs for the LLP and what should be considered as invisible in the LLP decays.
It also provides a list of $c\tau$ values for which the efficiencies will be computed and other options.
The output will consist of a csv file (stored in the same folder as the input file) for each HepMC file containing the efficiencies for
the provided $c\tau$ values.

If running over multiple files, the code allows to run the efficiency calculation in parallel (the number of parallel runs is set by the `ncpus` option).


### Method 2

After generating a ROOT file using the modified Delphes code (see above),  the efficiencies can be computed running:

```
./getEffsFromROOT.py -p <parameters_file> -f <path-to-root-file(s)>
```
where one example of the parameters file can be found [here](./parameters_getEff.ini).
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
the corresponding ( $m_{\Phi}$ , $m_{S}$ ) mass values assumed by ATLAS (for comparing against the official ATLAS curves) and EFFPLOT and XSECPLOT are the names of the efficiency plot and cross-section upper limit plot, respectively.


## Validation

More information about validation can be found in the [validation](./validation/) folder and in the [recast note](https://arxiv.org/pdf/2412.13976).