# LLP Recasting Repository

This repository holds example codes for recasting long-lived particle (LLP) searches.
The code authors and repository maintainers are not responsible for how the code is used and the user should
use discretion when applying it to new models.

## Adding your recasting code ##

This is an open repository and if you have developed a code for recasting a LLP analysis,
we encourage you to include it here.
Please contact <llp-recasting@googlegroups.com> and we will provide you with the necessary information
for including your code.

## Repository Structure ##

The repository folder structure is organized according to the type of LLP signature and the
corresponding analysis and authors:


  * [Displaced Vertices](DisplacedVertices)
    * [13 TeV ATLAS Displaced Jets](DisplacedVertices/ATLAS-SUSY-2018-13)
    * [13 TeV ATLAS Displaced Vertex plus MET by ALessa](DisplacedVertices/ATLAS-SUSY-2016-08_ALessa)
    * [13 TeV ATLAS Displaced Vertex plus MET by GCottin](DisplacedVertices/ATLAS-SUSY-2016-08_GCottin)
    * [8 TeV ATLAS Displaced Vertex plus jets by GCottin](DisplacedVertices/ATLAS-SUSY-2014-02_GCottin)
  * [Displaced Jets](DisplacedJets) (CalRatio)
    * [13 TeV ATLAS Displaced Jets in the calorimeter](DisplacedJets/ATLAS-EXOT-2019-23)
    * [13 TeV ATLAS Displaced Jets plus X in the calorimeter](DisplacedJets/ATLAS-EXOT-2022-04)
  * [Emerging Jets](EmergingJets/CMS-EXO-18-001)    
  * [Heavy Stable Charged Particles](HSCPs)
    * [13 TeV ATLAS HSCP - 139/fb](HSCPs/ATLAS-SUSY-2018-42)
    * [13 TeV ATLAS HSCP - 31.6/fb](HSCPs/ATLAS-SUSY-2016-32)
    * [8 TeV CMS HSCP](HSCPs/CMS-EXO-12-026)    
  * [Disappearing Tracks](DisappearingTracks/ATLAS-SUSY-2016-06)

A README file can be found inside each folder with the required dependencies
and basic instructions on how to run the recasting codes.

## Running the Recasting Code ##

Instructions on how to compile and run the main executable are provided in every recasting folder.
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


## Contact ##

If you have any questions, comments or want to contact the repository maintainers,
please send an e-mail to <lp-recasting@googlegroups.com>
