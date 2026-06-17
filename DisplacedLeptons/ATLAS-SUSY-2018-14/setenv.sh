#!/bin/sh

currentDIR="$( pwd )"
pythiaDIR=$currentDIR/MG5/HEPTools/pythia8
lhapdfDIR=$currentDIR/MG5/HEPTools/lhapdf6_py3
export LD_LIBRARY_PATH=$pythiaDIR/lib:$lhapdfDIR/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$currentDIR/MG5/HEPTools/lhapdf6_py3/local/lib/python3.12/dist-packages:$PYTHONPATH
export ROOT_INCLUDE_PATH=$currentDIR/DelphesLLP/external:$ROOT_INCLUDE_PATH
