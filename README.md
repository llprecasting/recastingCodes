# LLP Recasting Repository

This repository holds example codes for recasting long-lived particle (LLP) searches.

## Adding your recasting code ##

This is an open repository and if you have developed a code for recasting a LLP analysis, 
we encourage you to include it here.
Please contact <llp-recasting@googlegroups.com> and we will provide you with the necessary information
for including your code.

## Repository Structure ##

The repository folder structure is organized according to the type of LLP signature and the
corresponding analysis and authors:

  * [Displaced Vertices](DisplacedVertices)  
    * [13 TeV ATLAS Displaced Vertex plus MET by ALessa](DisplacedVertices/ATLAS-SUSY-2016-08)
    * [13 TeV ATLAS Displaced Vertex plus MET by GCottin](DisplacedVertices/ATLAS-SUSY-2016-08_GCottin)
    * [8 TeV ATLAS Displaced Vertex plus jets by GCottin](DisplacedVertices/ATLAS-SUSY-2014-02_GCottin)
  * [Heavy Stable Charged Particles](HSCPs)
    * [8 TeV CMS HSCP](HSCPs/CMS-EXO-12-026)
    * 13 TeV CMS HSCP (CMS-PAS-EXO-16-036)
  * [Displaced Leptons](DisplacedLeptons)
  * [Displaced Jets](DisplacedJet)

A README file can be found inside each folder with the required dependencies
and basic instructions on how to run the recasting codes.

## Running the Recasting Code ##

A Makefile is provided inside each analysis folder which
compiles the main executable once the code pre-requisites
have been installed.
For instance, the recasting of the 8 TeV CMS search
requires [Pythia8](http://home.thep.lu.se/~torbjorn/pythia8/).
After downloading and compiling Pythia 8, the main recasting code
can be compiled with the following steps:

   1. Go to the [HSCPs/CMS-EXO-12-026](HSCPs/CMS-EXO-12-026) folder
   2. Make sure Pythia 8 is available and run  `make main_hscp.exe -pythia8path=<path-to-pythia>`

Finally the compiled code can be run and its options displayed running:

```
./main_hscp.exe --help
```

## Citation ##

If you use any of the recasting code stored in this repository please cite ...

## Contact ##

If you have any questions, comments or want to contact the repository maintainers,
please send an e-mail to <lp-recasting@googlegroups.com>
