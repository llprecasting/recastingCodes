###### The following paths need to be specified:

hepmcpath=
pythia8path=
fastjetpath=

YODApath=
#### Here a subtlety, need e.g.:
#### YAMLpath=<where your YODA code is>/YODA-1.8.5/src/yamlcpp
YAMLpath=



#################################################################
# These are standard defaults
# You can change the compiler but changing the lib extension requires
# more effort.

CXX      := g++

LIBEXTENSION=so


#################################################################

homeDIR = $(shell pwd)
HEPDATAINCLUDEDIR=$(homeDIR)/DATA/




CXXFLAGS := -O3 -g  -std=c++11 -DGZSTREAM_NAMESPACE=GZ -DHEPDATAINCLUDEDIR=\"$(HEPDATAINCLUDEDIR)\" -I$(pythia8path)/include -I$(pythia8path)/include/Pythia8 -I$(fastjetpath)/include -I$(hepmcpath)/include -I$(YODApath)/include -I$(YAMLpath)


PLDFLAGS  := -L$(pythia8path)/lib -L$(fastjetpath)/lib -Wl,-rpath,$(pythia8path)/lib -Wl,-rpath,$(fastjetpath)/lib -L$(YODApath)/lib -Wl,-rpath,$(YODApath)/lib
HLDFLAGS  := -L$(fastjetpath)/lib -Wl,-rpath,$(fastjetpath)/lib -L$(hepmcpath)/lib -Wl,-rpath,$(hepmcpath)/lib -L$(YODApath)/lib -Wl,-rpath,$(YODApath)/lib


HEADERS=$(wildcard *.h) $(wildcard Analyses/*.h) $(wildcard include/*.h)  $(wildcard include/ExHEPUtils/*.h)
CSOURCES =src/gzstream.C
CCSOURCES=$(wildcard src/*.cc) $(wildcard Analyses/*.cc)

PY_SOURCE=analysePYTHIA.cc
PYLHE_SOURCE=analysePYTHIA_LHE.cc
HEPMC_SOURCE=analyseHEPMC.cc


OBJS=$(CCSOURCES:.cc=.o) $(CSOURCES:.C=.o)


PY_OBJ=$(PY_SOURCE:.cc=.o)
PYLHE_OBJ=$(PYLHE_SOURCE:.cc=.o)
HEPMC_OBJ=$(HEPMC_SOURCE:.cc=.o)
GZSTREAM_NAMESPACE=GZ

MY_LIB=HackAnalysis
MY_LIB_FULL=lib$(MY_LIB).$(LIBEXTENSION)



all: analysePYTHIA.exe analysePYTHIA_LHE.exe analyseHEPMC.exe

analysePYTHIA.exe: $(MY_LIB_FULL) $(A_LIB_FULL) $(PY_OBJ)
	$(CXX) -o $@ $(PY_OBJ) -L$(homeDIR)/lib -Wl,-rpath,$(homeDIR)/lib $(PLDFLAGS) -fopenmp -DGZSTREAM_NAMESPACE=GZ -DGZIPSUPPORT -lpythia8 -ldl -lz -lfastjet -lfastjettools -l$(MY_LIB) -lYODA 

analysePYTHIA_LHE.exe:	$(MY_LIB_FULL) $(A_LIB_FULL) $(PYLHE_OBJ)
	$(CXX) -o $@ $(PYLHE_OBJ) -L$(homeDIR)/lib -Wl,-rpath,$(homeDIR)/lib $(PLDFLAGS) -fopenmp -DGZSTREAM_NAMESPACE=GZ -DGZIPSUPPORT -lpythia8 -ldl -lz -lfastjet -lfastjettools -l$(MY_LIB) -lYODA

analyseHEPMC.exe:	$(MY_LIB_FULL) $(A_LIB_FULL) $(HEPMC_OBJ) 
	$(CXX) -o $@ $(HEPMC_OBJ) -L$(homeDIR)/lib -Wl,-rpath,$(homeDIR)/lib $(HLDFLAGS) -lHepMC  -DGZSTREAM_NAMESPACE=GZ -DGZIPSUPPORT -ldl -lz -lfastjet -lfastjettools -l$(MY_LIB) -lYODA


%.o: %.C $(HEADERS)
	$(CXX)   -c -o $@ $< $(CXXFLAGS) -fPIC -fopenmp -DGZSTREAM_NAMESPACE=GZ -iquote./include

%.o: %.cc $(HEADERS) 
	$(CXX) -c -o $@ $< $(CXXFLAGS) -fPIC -fopenmp -DGZSTREAM_NAMESPACE=GZ -iquote.

$(MY_LIB_FULL): $(OBJS)
	$(CXX) -shared -fopenmp $(CXXFLAGS) -L. -Wl,-rpath,. -o $(homeDIR)/lib/$@ $^

clean:
	rm -f *.o *.so *.exe include/*.o src/*.o Analyses/*.o lib/*.so lib/*.a

