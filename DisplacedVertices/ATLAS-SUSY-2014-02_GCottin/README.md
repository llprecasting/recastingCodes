# ATLAS Displaced Vertex plus jets 8 TeV Recast #

## Authors: ##
[Giovanna Cottin](mailto:gfcottin@gmail.com)

This repository holds the main code for recasting the 8 TeV ATLAS search for displaced vertices
in association with jets ([arXiv:1504.05162](https://arxiv.org/abs/1504.05162)). See also [ATLAS-SUSY-2014-02](https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2014-02/). It includes
a simple displaced vertex reconstruction algorithm based on a functional form for tracking efficiency.

This code was validated in:

* https://arxiv.org/abs/1606.03099

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
./displacedRecoVertex
```
An example Makefile can be found in DisplacedVertices/ATLAS-SUSY-2016-08_GCottin/pythiaCode/Makefile
The basic required input is a (parton level) LHE or SLHA file.
The output are efficiency files and cutflows

## pythiaCode ##

In this folder you will find:

* GGM2_softsusy_sdecay.slha  -- the ATLAS model benchmark we generated
* GGM_softsusy_sdecay.slha -- the ATLAS model benchmark we generated
* RPV_softsusy_sdecay.slha -- the ATLAS model benchmark we generated
* displacedBM.cmnd -- pythia configuration parameters 
* displacedRecoVertex.cc -- the main displaced vertex plus jets code
* ToyDetector-ATLAS.cc -- the custom made detector simulation for this recast
* ToyDetector-ATLAS.h

## Plots ##

In this folder you will find:

* dgsPhenoPlots.py -- the main plotting script for making the validation plots shown in [arXiv:1606.03099](https://arxiv.org/abs/1606.03099)
* Validation plots showing efficiency against proper lifetime for the three ATLAS benchmarks
