# Displaced Vertex Recasting #

## Authors: ##
[Andre Lessa](mailto:andre.lessa@ufabc.edu.br)

This repository holds the main code for recasting the 13 TeV ATLAS search for displaced vertices
plus missing energy ([ATLAS-SUSY-2016-08](https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-08/))
using the efficiency grids  for event and vertex reconstruction selection provided [here](https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-08/hepdata_info.pdf).

## Pre-Requisites ##

The following pre-requisites must be installed before compiling the main code:

  * [C++ boost library](https://www.boost.org/)
  * [FastJet](http://fastjet.fr/)
  * [Pythia8](http://home.thep.lu.se/~torbjorn/pythia8) (compiled with FastJet and GZIP support)

## Installation/Compiling ##

In order to compile the recasting code, run:

```
make main_atlas_susy_2016_08.exe pythia8path=<PATH-TO-PYTHIA8FOLDER> fastjetpath=<PATH-TO-FASTJETFOLDER>
```

where <PATH-TO-PYTHIA8FOLDER> (<PATH-TO-FASTJETFOLDER>) should point to the Pythia8 (FastJet) folder containing the lib and include folders.

## Running ##

Instructions for running the main code can be obtained running:

```
./main_atlas_susy_2016_08.exe -h
```

The basic required input is a SLHA file.
For instance, running:

```
./main_atlas_susy_2016_08.exe -f example.slha -n 100 -o example.eff
```

should generate 100 events and display the efficiencies for the given input file.

## Validation ##

The validation of the signal efficiencies (efficiency times acceptance)
for direct production of staus can be found in the [validation folder](validation).
The output (.eff files) was generated using the SLHA files and pythia8.cfg and parameters.ini files stored in the folder.
The following validation plot can be generated running this [ipython notebook](validation/validation.ipynb):


![Alt text](validation/validationPlot.png?raw=true "Validation Plot")

