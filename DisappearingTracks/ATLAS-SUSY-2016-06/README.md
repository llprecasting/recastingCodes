# Disappearing Tracks Recasting #

## Authors: ##
[Felipe Rojas](mailto:astrofis.rojas@gmail.com)

[Jose Zurita](mailto:jose.zurita@kit.edu)


This folder holds the main code for recasting the 13 TeV disappearing-track signature in pp collisions with the ATLAS detector ([ATLAS-SUSY-2016-06](http://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-06/)).
This recast follows the "EW-production" analysis, hence it is valid for models with EW charged multiplets.

## Pre-Requisites ##

The following pre-requisites must be installed before compiling the main code:

  * Python (>= v2.7)
  * [ROOT](https://root.cern/) (with Python bindings)
  * [Delphes](https://cp3.irmp.ucl.ac.be/projects/delphes)

## Installation ##

Once ROOT (compatible with Python 2) and Delphes are available, make sure to set the enviroment variable:

```
ROOT_INCLUDE_PATH=$DELPHESPATH/external
```

where $DELPHESPATH is the Delphes directory.

## Running ##

For running the code, the user must provide:
 1. the input parameters which must to be written inside the file input_param.dat. The minimun amount of parameters the used must to provide are:
   * DelphesPath: the path to the Delphes folder
   * rootFile: the path to the .root file with Delphes processed events
   * outputFile: the path to the output result of the analysis
   * Luminosity: the integrated luminosity in inverse femtobarns
   * Weight: the way events are weighted. It can be either 'same' weight for all events or 'root' for take the weight from root file itself. If Weight = same, the user must also provided the total cross section of the events with the parameter InitXs.
   * tau_values: values for the LLP lifetime for which the user wants the efficiencies to be computed.
   * Aditional options are available:
     * the user can define PID number to identify the chargino/neutraino particle for a given model with the flag: PIDchargino, PIDneutralino.
     * NLO-QCD corrections can be added with parameter kfactor.
     * In order to increase the tracklets statistic the user can use the parameter nloop to an integer larger than 1. The parameter nloop controls how many times the decay distribution is sampled for a given MC event.

Once the above parameters are set, the efficiencies can be computed running:

```
python ./example.py
```

## Output ##

The output of the program is a .dat file with 9 columns defined as:

 * Char_Mass(GeV): "Chargino" mass of the sample
 * Neut_Mass(GeV): "Neutralino" mass of the sample
 * tau(ns): List of lifetime values defined in tau_array flag inside example.py file.
 * EAxEE: Event Acceptance, number of events pasing the Event-Selection level over the total number of events.
 * TAxTE: Tracklet Acceptance, number of tracklets passing the reconstruction Level over number of tracklets passing the Event-Selection level
 * EAxEExTAxTE: Total Acceptance x Efficiency.
 * total_eff: Total efficiency after the analysis, where Init_xs * total_eff = xs100.
 * xs100(fb): Cross section after the event selection in fb (the 100 refers to the cut on "chargino" pT > 100 GeV).
 * xslim(fb): is the upper limit LO production cross section cross section of the point (in fb).
 * Ntr: Expected number of tracklets based on the value of the cross section after event selection ( xs100 x Luminosity )
 * r: Fraction between the cross section after the analysis over the upper limit on the model-independent visible cross-section at 95% CL in fb. from Table 4, article 1712.02118
 * MCev: Number of MonteCarlo events


## Notes ##


 * For a succesful run (namely avoiding the "ImportError: No module named ROOT" message we had to set the following environmental variables:
    1. ROOT_INCLUDE_PATH=$DELPHESPATH/external (where $DELPHESPATH is the Delphes directory)
    2. PYTHONPATH=$ROOTSYS/lib:$ROOT_INCLUDE_PATH
 * The LLP decays in the event file are not used and the charged LLP is always assumed to decay to invisible and soft particles.

## Validation ##

The validation of the signal efficiencies (efficiency times acceptance) and upper limits are documented in [arXiv:2008.08581](https://arxiv.org/abs/2008.08581).
