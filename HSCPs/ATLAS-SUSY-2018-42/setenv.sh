#!/bin/sh

currentDIR="$( pwd )"
delphesDIR=$currentDIR/DelphesHSCP
pythiaDIR=$currentDIR/MG5/HEPTools/pythia8
#Make sure pythia can be found by Delphes
export LD_LIBRARY_PATH=$pythiaDIR/lib:$LD_LIBRARY_PATH
#Make sure Delphes can be found by ROOT
export ROOT_INCLUDE_PATH=$delphesDIR/external

