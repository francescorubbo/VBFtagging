# --------------------------------------------- #
# Makefile for Recluster code                        #
# Pascal Nef, March 6th 2014                    #
#                                               #
# Note: source setup.sh before make             #
# --------------------------------------------- #

CXXFLAGS = -O2 -Wall -Wextra -std=c++0x -g

.PHONY: clean debug all

all: setup VBFTagging

setup:
	mkdir -p lib

VBFTagging:  lib/VBFTagging.so lib/VBFTaggingAnalysis.so lib/Configuration.so
	$(CXX) lib/VBFTagging.so lib/VBFTaggingAnalysis.so lib/Configuration.so -o $@ \
	$(CXXFLAGS) -Wno-shadow  \
	`root-config --glibs` -lEG -lEGPythia8 \
	-I./include -L./lib \
	-L$(FASTJETLOCATION)/lib `$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins ` -lSubjetJVF  -lVertexJets \
	-L$(PYTHIA8LOCATION)/lib -lpythia8 -llhapdfdummy \
	-L$(BOOSTLIBLOCATION) -lboost_program_options 

lib/VBFTagging.so: src/VBFTagging.C lib/VBFTaggingAnalysis.so
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` -lSubjetJVF -lVertexJets \
	-I./include -L./lib \
	-I$(PYTHIA8LOCATION)/include \
	-I $(BOOSTINCDIR) \
	`root-config --cflags` 

lib/VBFTaggingAnalysis.so : src/VBFTaggingAnalysis.cc include/VBFTaggingAnalysis.h 
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` -lSubjetJVF -lVertexJets \
	-I./include \
	-I$(PYTHIA8LOCATION)/include \
	`root-config --cflags --libs` 

lib/Configuration.so : src/Configuration.cc include/Configuration.h
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` -lSubjetJVF -lVertexJets \
	-I./include \
	-I$(PYTHIA8LOCATION)/include \
	`root-config --cflags --libs`

clean:
	rm -rf VBFTagging
	rm -rf lib
	rm -f *~

install:
	install VBFTagging -t ${HOME}/local/bin
	install setup.sh -t ${HOME}/local/bin
	install scripts/Timing.sh -t ${HOME}/local/bin
