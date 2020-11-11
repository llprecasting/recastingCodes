# LLP Recasting Repository - LLP Workshop

This repository holds codes for recasting long-lived particle (LLP) searches.

**Not all the code versions available in this branch are final/stable. Use the master branch
for the stable releases.**

The code authors and repository maintainers are not responsible for how the code is used and the user should
use discretion when applying it to new models.


## Repository Structure ##

The repository folder structure is organized according to the type of LLP signature and the
corresponding analysis and authors:


  * [Displaced Vertices](DisplacedVertices)
    * [13 TeV ATLAS Displaced Vertex plus MET by ALessa](DisplacedVertices/ATLAS-SUSY-2016-08_ALessa)
    * [13 TeV ATLAS Displaced Vertex plus MET by GCottin](DisplacedVertices/ATLAS-SUSY-2016-08_GCottin)
    * [8 TeV ATLAS Displaced Vertex plus jets by GCottin](DisplacedVertices/ATLAS-SUSY-2014-02_GCottin)
  * [Heavy Stable Charged Particles](HSCPs)
    * [8 TeV CMS HSCP](HSCPs/CMS-EXO-12-026)
    * [13 TeV ATLAS HSCP](HSCPs/ATLAS-SUSY-2016-32)
  * [Disappearing Tracks](DisappearingTracks)
    * [13 TeV ATLAS DT](DisappearingTracks/ATLAS-SUSY-2016-06)


A README file can be found inside each folder with the required dependencies
and basic instructions on how to run the recasting codes.

## Requirements ##

Most of the recasting codes present here requires a local installation of [Pythia8](http://home.thep.lu.se/Pythia/) and [Delphes](https://cp3.irmp.ucl.ac.be/projects/delphes).
For some specific analyses [FastJet](http://fastjet.fr/) and [ROOT](https://root.cern/) are also requirements.
It is recommended that Pythia, Delphes and FastJet are installed in the top folder.


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
