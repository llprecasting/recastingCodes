## Makefile, static linking. Pythia8 is fetched and built
## automatically.



homeDIR = $(shell pwd)
pythia8path = $(homeDIR)/../pythia8 ## Set the path to your local pythia folder

CXX      := g++

CXXFLAGS := -O3 -I$(pythia8path)/include -I$(pythia8path)/include/Pythia8/
LDFLAGS  := -L$(pythia8path)/lib/ -L$(pythia8path)/lib -Wl,-rpath,$(pythia8path)/lib
XMLDOC   := $(pythia8path)/share/Pythia8/xmldoc

all: maindEJ


clean:
	rm maindEJ.cc

maindEJ: maindEJ.cc
	echo $(XMLDOC) > xml.doc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ maindEJ.cc -lpythia8 -ldl -DGZIPSUPPORT -lz

