# Makefile for displaced searches for MAC desktop
PYTHIA8=/PATH/TO/pythia8308
FASTJET3=/PATH/TO/fastjet-3.4.0-install
PYTHIALIB=$(PYTHIA8)/lib
PYTHIAINC=$(PYTHIA8)/include/
FJINC=$(FASTJET3)/include
FJLIB=$(FASTJET3)/lib
OBJECTS=Simpler_ToyDetector.o recast_2301_13866.o efficiency.o
INCLUDE=-I$(PYTHIAINC) -I./ -I$(FJINC)


#Need to do in MAC
# export DYLD_LIBRARY_PATH=/PATH/TO/pythia8308/lib:...../PATH/TO/fastjet-3.4.0-install/lib

recast_2301_13866: $(OBJECTS)
	$(CXX) -o recast_2301_13866 $(OBJECTS) -L$(PYTHIALIB) -lpythia8 \
	-L$(FJLIB) -lfastjet -lz -Wl,-rpath $(PYTHIALIB) -Wl,-rpath $(FJLIB)

%.o: %.cc
	$(CXX) -c $^ $(INCLUDE) -std=c++11 


