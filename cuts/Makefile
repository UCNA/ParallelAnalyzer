# Makefile from Jianglai

ObjSuf        = o
SrcSuf        = C
DllSuf        = so
OutPutOpt     = -o  
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
ROOTCFLAGS    = $(shell root-config --cflags)

VPATH = ./:../include/

# Work with Linux with egcs	
CXX           = g++ 
CXXFLAGS      = -O2 -Wall -fPIC -std=c++11
LD            = g++
SOFLAGS       = -shared
LIBS          = $(ROOTLIBS) $(ROOTGLIBS)
CXXFLAGS     += $(ROOTCFLAGS) -I../include
LIBS         += -lSpectrum -lMinuit

objects1 = DataTree.o
objects2 = MWPCPositionResponse.o 

.PHONY: all
all: beamDropCuts.exe 

beamDropCuts.exe: beamDropCuts.o $(objects1)
	$(CXX) $(CXXFLAGS) -o beamDropCuts.exe beamDropCuts.o $(objects1) $(LIBS)
	@echo "Done"

cathodeTriggerCuts.exe: cathodeTriggerCuts.o $(objects2)
	$(CXX) $(CXXFLAGS) -o cathodeTriggerCuts.exe cathodeTriggerCuts.o $(objects2) $(LIBS)
	@echo "Done"

DataTree.o: DataTree.hh

MWPCPositionResponse.o: MWPCPositionResponse.hh

# -------------------------------------------------------------------------------
#  Generic compilation and linking step to make an executable from
#  a single *.cc file
#
#%: %.$(SrcSuf)
#	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)
#	@echo "$@ done"

clean:
		@rm -f *.o *~  core
