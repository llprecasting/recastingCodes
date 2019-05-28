PYTHIA8=/Users/ChubiCottin/HEPSoftware/pythia8215
FASTJET3=/Users/ChubiCottin/HEPSoftware/fastjet-3.1.3/fastjet-install
PYTHIALIB=$(PYTHIA8)/lib
PYTHIAINC=$(PYTHIA8)/include/
FJINC=$(FASTJET3)/include
FJLIB=$(FASTJET3)/lib
INCLUDE=-I$(PYTHIAINC) -I./ -I$(FJINC)
OBJECTS=displacedTruthVertex_paramEff_RhadID.o ToyDetector-ATLAS-tracklessjet.o parametrized_truthEff.o
# OBJECTS=displacedTruthVertex_varyMass.o ToyDetector-ATLAS-tracklessjet.o parametrized_truthEff.o

displacedTruthVertex_paramEff_RhadID: $(OBJECTS)
	$(CXX) $(CXX_COMMON) -o displacedTruthVertex_paramEff_RhadID $(OBJECTS) -L$(PYTHIALIB) -lpythia8 \
	-L$(FJLIB) -lfastjet -lz
%.o: %.cc
	$(CXX) $(CXX_COMMON) -c $^ $(INCLUDE) 

# displacedTruthVertex_varyMass: $(OBJECTS)
# 	$(CXX) $(CXX_COMMON) -o displacedTruthVertex_varyMass $(OBJECTS) -L$(PYTHIALIB) -lpythia8 \
# 	-L$(FJLIB) -lfastjet -lz
# %.o: %.cc
# 	$(CXX) $(CXX_COMMON) -c $^ $(INCLUDE) 
