homeDIR = $(shell pwd)

#######################################
pythia8path=
fastjetpath=
#######################################

CXX      := g++
CXXFLAGS := -O3 -std=c++11 -DGZSTREAM_NAMESPACE=GZ -I$(pythia8path)/include -I$(pythia8path)/include/Pythia8/ -I$(ROOTSYS)/include 
LDFLAGS  := -L$(pythia8path)/lib/ -L$(pythia8path)/lib -Wl,-rpath,$(pythia8path)/lib 

GZSTREAM_NAMESPACE=GZ


all: create.exe


create.exe: 	create.cc
	$(CXX)   -c -o gzstream.o src/gzstream.C $(CXXFLAGS) -iquote./include
	$(CXX)   -c -o create.o create.cc $(CXXFLAGS) -iquote./include 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ create.o gzstream.o -iquote. -iquote./include -lpythia8 -ldl -DGZIPSUPPORT -lz
