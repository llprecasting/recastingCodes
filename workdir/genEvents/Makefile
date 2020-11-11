# Set the shell.
SHELL=/usr/bin/env bash

# Include the configuration.
# -include Makefile.inc

homeDIR = $(shell pwd)

DVPATH = ../../DisplacedVertices/ATLAS-SUSY-2016-08_ALessa/

PYTHIA8 := $(homeDIR)/../pythia8
HEPMC := $(homeDIR)/../HepMC
FASTJET := $(homeDIR)/../fastjet-3.3.4

CXX      := g++
XMLDOC   := $(PYTHIA8)/share/Pythia8/xmldoc
CXXFLAGS := -O3 -std=c++11 -I$(PYTHIA8)/include -I$(PYTHIA8)/include/Pythia8/ -I$(HEPMC)/include -I$(ROOTSYS)/include  -I$(FASTJET)/include
LDFLAGS  := -L$(PYTHIA8)/lib/ -L$(PYTHIA8)/lib -Wl,-rpath,$(PYTHIA8)/lib -L$(FASTJET)/lib -Wl,-rpath,$(FASTJET)/lib -L$(HEPMC)/lib -Wl,-rpath,$(HEPMC)/lib

all: main_pythia2HepMC.exe main_atlas_susy_2016_08_cutFlow

clean:
	rm main_pythia2HepMC.exe main_atlas_susy_2016_08_cutFlow


main_pythia2HepMC.exe: main_pythia2HepMC.cc 
	echo $(XMLDOC) > xml.doc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@  main_pythia2HepMC.cc -lpythia8 -ldl -DGZIPSUPPORT -lz -L$(HEPMC)/lib -Wl,-rpath,$(HEPMC)/lib -lHepMC -L$(ROOTSYS)/lib -Wl,-rpath,$(ROOTSYS)/lib `$(ROOTSYS)/bin/root-config --glibs`
	
	
main_atlas_susy_2016_08_cutFlow.exe: $(DVPATH)/recastCode/main_atlas_susy_2016_08_cutFlow.cc $(DVPATH)/recastCode/helperFunctions.h $(DVPATH)/recastCode/dataReader.h
	echo $(XMLDOC) > xml.doc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(DVPATH)/recastCode/main_atlas_susy_2016_08_cutFlow.cc -lpythia8 -ldl -DGZIPSUPPORT -lz -lfastjet -lHepMC

