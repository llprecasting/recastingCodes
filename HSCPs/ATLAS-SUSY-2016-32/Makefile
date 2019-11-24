
CXX      := g++
CXXFLAGS := -O3 -std=c++11 -I$(pythia8path)/include -I$(pythia8path)/include/Pythia8/ -I$(ROOTSYS)/include
LDFLAGS  := -L$(pythia8path)/lib/ -L$(pythia8path)/lib -Wl,-rpath,$(pythia8path)/lib
XMLDOC   := $(pythia8path)/share/Pythia8/xmldoc

all: main_Rhadron.exe main_HSCP.exe

clean:
	rm main_Rhadron.exe main_HSCP.exe

main_Rhadron.exe: recastCode/main_Rhadron.cc recastCode/helperFunctions.h 
	echo $(XMLDOC) > xml.doc
    
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@  recastCode/main_Rhadron.cc -lpythia8 -ldl -L$(ROOTSYS)/lib -Wl,-rpath,$(ROOTSYS)/lib `$(ROOTSYS)/bin/root-config --glibs`


main_HSCP.exe: recastCode/main_HSCP.cc recastCode/helperFunctions.h 
	echo $(XMLDOC) > xml.doc
    
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@  recastCode/main_HSCP.cc -lpythia8 -ldl -L$(ROOTSYS)/lib -Wl,-rpath,$(ROOTSYS)/lib `$(ROOTSYS)/bin/root-config --glibs`

