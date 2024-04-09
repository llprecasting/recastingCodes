# Code for recasting 2301.13866

**By Kingman Cheung, Fei-Tung Chung, Giovanna Cottin, and Zeren Simon Wang**


[![arXiv](http://img.shields.io/badge/arXiv-2404.XXXXX-B31B1B.svg)](https://arxiv.org/abs/2404.XXXXX)
[![license](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/llprecasting/recastingCodes/XXX/LICENSE.md)


## Introduction

We develop a code for recasting the LLP search reported in 2301.13866.


## Paper

This code is used for recasting the LLP search reported in 2301.13866. Its validation and application in one BSM theoretical scenario are presented in [![arXiv](http://img.shields.io/badge/arXiv-2404.XXXXX-B31B1B.svg)](https://arxiv.org/abs/2404.XXXXX).


## Compilation

The code is written in C++. There is a Makefile. In the first lines of Makefile, put in the paths to your own Pythia8 and FastJet installation. Pythia8 should have been compiled with a link to FastJet. Run "make -jN" in the terminal to compile, where N is the number of CPU cores you want to use.



## Input cards

We provide the MadGraph5 input cards and validation results.


## Run the code

Run the code with ./recast_2301_13866 and the path to the input pythia card. In the cmnd card, specify the path to the input LHE files.


## Questions and bug reports

If you have any questions, or would like to report any bug, you are welcome to write to any of us:

- Kingman Cheung: [cheung@phys.nthu.edu.tw](mailto:cheung@phys.nthu.edu.tw)
- Fei-Tung Chung: [feitung.chung@gapp.nthu.edu.tw](mailto:feitung.chung@gapp.nthu.edu.tw)
- Giovanna Cottin: [gfcottin@uc.cl](mailto:gfcottin@uc.cl)
- Zeren Simon Wang: [wzs@mx.nthu.edu.tw](mailto:wzs@mx.nthu.edu.tw)
