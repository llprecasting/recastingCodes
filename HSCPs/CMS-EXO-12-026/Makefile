
CXX      := g++
CXXFLAGS := -O3 -std=c++11 -I$(pythia8path)/include -I$(pythia8path)/include/Pythia8/ -I$(ROOTSYS)/include
LDFLAGS  := -L$(pythia8path)/lib/ -L$(pythia8path)/lib -Wl,-rpath,$(pythia8path)/lib
XMLDOC   := $(pythia8path)/share/Pythia8/xmldoc

all: main_hscp.exe

clean:
	rm main_hscp.exe

main_hscp.exe: recastCode/main_hscp.cc recastCode/helperFunctionsHSCP.h 
	echo $(XMLDOC) > xml.doc
    
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@  recastCode/main_hscp.cc -lpythia8 -ldl -DGZIPSUPPORT -lz -L$(ROOTSYS)/lib -Wl,-rpath,$(ROOTSYS)/lib `$(ROOTSYS)/bin/root-config --glibs`

