# ATLAS Displaced Vertex 13 TeV Recast #

## Authors: ##
[Giovanna Cottin](mailto:gfcottin@gmail.com)

This repository holds the main code for recasting the 13 TeV ATLAS search for displaced vertices
plus missing transverse momenta ([ATLAS-SUSY-2016-08](https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-08/))
using the parametrized efficiencies for event and displaced vertex reconstruction provided [here](https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-08/hepdata_info.pdf).

This code was used in:

* https://arxiv.org/abs/1803.10379

Please cite if you use this code :-)

## Pre-Requisites ##

The following pre-requisites must be installed before compiling the main code:

  * python
  * [FastJet](http://fastjet.fr/)
  * [Pythia8](http://home.thep.lu.se/~torbjorn/pythia8) (compiled with FastJet)

## Installation/Compiling ##

In order to compile the recasting code, run:

```
make 
```

inside the pythiaCode folder, where a Makefile example can be found

## Running ##

You can run with:

```
./displacedTruthVertex_paramEff_RhadID
```
For the [LesHouches2017](https://arxiv.org/abs/1803.10379) reinterpretation. 
The basic required input is a (parton level) LHE or SLHA file.
The output are efficiency files and cutflows (save to the "truthDV_data" folder inside Plots)

## pythiaCode ##

In this folder you will find:

* Makefile -- paths should be edited accordingly      
* gluino_tau_0.0.slha -- the Split SUSY model spectra          
* pythia8.cfg -- pythia configuration parameters     
* parametrized_truthEff.cc -- the ATLAS digitized efficiencies
* parametrized_truthEff.h  
* displacedTruthVertex_paramEff_RhadID.cc -- the recasting code for each mass, varying the lifetime
* displacedTruthVertex_varyMass.cc -- the recasting code for for a fixed lifetime, varying the gluino mass 
* ToyDetector-ATLAS-tracklessjet.cc -- the custom made detector simulation for this recast
* ToyDetector-ATLAS-tracklessjet.h 

## Plots ##

Comparison plots with the recasted limits and the ATLAS observed limits.

## ATLASDV_MET note ##

ATLASDV_MET_recast.pdf is a note detailing the recasting procedure. An older recasted version (without the parametrized efficiencies) is also mentioned. 