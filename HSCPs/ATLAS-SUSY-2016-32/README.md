# HSCP Recasting #

## Authors: ##
[Andre Lessa](mailto:andre.lessa@ufabc.edu.br)

This repository holds the main code for recasting the 13 TeV ATLAS heavy stable charged particle
search ([ATLAS-SUSY-2016-32](http://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-32/))
based on the recasting details and code provided [here](http://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-32/hepdata_info.pdf) .

## Pre-Requisites ##

The following pre-requisites must be installed before compiling the main code:

  * [Pythia8](http://home.thep.lu.se/~torbjorn/pythia8/) (with ROOT bindings)

## Installation/Compiling ##

In order to compile the recasting code for the R-hadron searches, run:

```
make main_Rhadron.exe pythia8path=<PATH-TO-PYTHIA8FOLDER>
```

where PATH-TO-PYTHIA8FOLDER should point to the Pythia8 folder containing the Pythia lib and include folders.
The recasting code for the color neutral (HSCP) searches can be compiled through:

```
make main_HSCP.exe pythia8path=<PATH-TO-PYTHIA8FOLDER>
```


## Running ##

Instructions for running the main code can be obtained running:

```
./main_HSCP.exe -h  or ./main_Rhadron.exe -h
```

The basic required input is a (parton level) LHE or SLHA file.
For instance, running:

```
./main_HSCP.exe -f example.slha -n 100
```

should generate 100 events and display the efficiencies for the given input file.
Notice that for R-hadron searches the Pythia configuration file (pythia8.cfg) may have to be adapted
in order to properly specify the hadronization properties of the new long-lived colored BSM particles.


## Validation ##

The validation of the signal efficiencies (efficiency times acceptance) can be found in the [validation folder](validation).
The efficiencies were computed using the SLHA files and the pythia8_xx.cfg file stored in the folder.
The following validation plot can be generated running this [ipython notebook](validation/validationGluino.ipynb):


![Alt text](validation/gluino_eff.png?raw=true "Validation Plot for Efficiencies")

![Alt text](validation/gluino_UL.png?raw=true "Validation Plot for Upper Limits")

* *The values provided by ATLAS for small lifetime values (< 50 ns) for the gluino R-hadron benchmark could not be validated.
  Therefore we recommend to use the results with care. The comparison between the ATLAS efficiencies (Table 6) and 
  the ones obtained through the recasting code for a few benchmarks is shown below:

  | gluino Mass | gluino lifetime | ATLAS eff. | Recasting eff. |
  | ----------- | --------------- | ---------- | -------------- |
  |   1 TeV     |      10 ns      |   0.065    |   0.015        |
  |   1 TeV     |      30 ns      |   0.121    |   0.076        |
  |   1 TeV     |      50 ns      |   0.125    |   0.101        |
  |   2 TeV     |      10 ns      |   0.060    |   0.023        |
  |   2 TeV     |      30 ns      |   0.132    |   0.116        |
  |   2 TeV     |      50 ns      |   0.146    |   0.160        |
  | ----------- | --------------- | ---------- | -------------- |







