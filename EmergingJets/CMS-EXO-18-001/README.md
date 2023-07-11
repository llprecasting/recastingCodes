# Emerging Jets Recasting #

## Authors: ##
[Juliana Mara Carrasco](jmcarrasco@ific.uv.es)

[Jose Zurita](mailto:jose.zurita@kit.edu)


This folder holds the main code for recasting the 13 TeV emerging jet signature in pp collisions with the CMS detector ([CMS-EXO-18-001](https://cms-results.web.cern.ch/cms-results/public-results/publications/EXO-18-001/)).

## Pre-Requisites ##

The following pre-requisites must be installed before compiling the main code:

  * Pythia8.230 (or higher)

## Installation ##

To compile the code, run:

```
make maindEJ pythia8path=<PATH-TO-PYTHIA8FOLDER>
```


## Running ##

The code requires as input a pythia configuration card, a .cmnd file with the model parameters. We provide two examples, one used by the CMS Collaboration (cardEJcms.cmd) based on the model from Schwaller et.al. arXiv:1502.05409 and another one for SM Higgs decaying into two dark quarks (cardHiggs.cmnd) making use of the python package from Knapen et.al. arxiv:2103.01238.
 
You can run, for example, with:

```
./maindEJ cardEJcms.cmnd
```

## Output ##

As output, the program prints, for each signal region, the signal acceptances and the respective ratio of the calculated number of signal events over the excluded one by CMS.


## Validation ##

The validation of the signal acceptances and exclusion limits is shown in section 3 of arxiv:xxxx.xxxx.
