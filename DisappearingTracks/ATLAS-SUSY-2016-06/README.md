# Disappearing Tracks Recasting #

## Authors: ##
[Felipe Rojas](mailto:astrofis.rojas@gmail.com)

[Jose Zurita](mailto:jose.zurita@kit.edu)


This folder holds the main code for recasting the 13 TeV disappearing-track signature in pp collisions with the ATLAS detector ([ATLAS-SUSY-2016-06](http://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-06/)).
This recast follows the "EW-production" analysis, hence it is valid for models with EW charged multiplets.

## Pre-Requisites ##

The following pre-requisites must be installed before compiling the main code:

  * Python2
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
 1. the .root file with Delphes processed events.  The path to this file is set by 'rootFile' inside the example.py program.  By default the "example1.root" file in the main directory is used.
 2. the path to the Delphes folder
 3. additonal (optional input), such as the lifetimes for which to compute the efficiencies

After this, write in a terminal:

```
python ./example.py
```

## Output ##

The output of the program is a .dat file with 9 columns defined as:

 * Char_Mass(GeV): "Chargino" mass of the sample
 * Neut_Mass(GeV): "Neutralino" mass of the sample
 * tau(ns): List of lifetime values defined in tau_array flag inside Variable.py file.
 * EAxEE: Event Acceptance, number of events pasing the Event-Selection level over the total number of events.
 * TAxTE: Tracklet Acceptance, number of tracklets pasing the reconstruction Level over number of tracklets pasing the Event-Selection level
 * EAxEExTAxTE: Total Acceptance x Efficiency.
 * total_eff: Total efficiency after the analysis, where Init_xs * total_eff = xs100.
 * xs100(fb): Cross section after the event selection in fb (the 100 refers to the cut on "chargino" pT > 100 GeV).
 * xslim(fb): Cross section limit is the cross section that should have the sample to be excluded (in fb).
 * Ntr: Expected number of tracklets based on the value of the cross section (xs100) x Luminosity
 * r: Fraction between the number of tracklets after the analysis over the observed number of tracklets.
 * MCev: Number of MonteCarlo events


## Notes ##


 * For a succesful run (namely avoiding the "ImportError: No module named ROOT" message we had to set the following environmental variables:
    1. ROOT_INCLUDE_PATH=$DELPHESPATH/external (where $DELPHESPATH is the Delphes directory)
    2. PYTHONPATH=$ROOTSYS/lib:$ROOT_INCLUDE_PATH
 * The charged and neutral particle are generically called "Chargino" and "Neutralino", a MSSM inspired notation. To employ another model the CharginoPID (1000024) and NeutralinoPID (1000022) variables must be changed in ./modules/Variables.py<sup>*</sup>
 * For a given "trial" value of the lifetime we increase the statistics by sampling its value Ntrial times from an exponential distribution. The user can adjust the Ntrial value in ./modules/Variables.py<sup>*</sup>


<sup>*</sup> *This feature has not been implemented yet.*

## Validation ##
