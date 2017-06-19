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
CXXFLAGS      = -O2 -Wall -fPIC --std=c++11
LD            = g++
SOFLAGS       = -shared
LIBS          = $(ROOTLIBS) $(ROOTGLIBS)
CXXFLAGS     += $(ROOTCFLAGS) -I../include
LIBS         += -lSpectrum -lMinuit

objects = DataTree.o MBUtils.o BetaDecayTools.o

endpointGain.exe: endpointGain.o $(objects)
	$(CXX) $(CXXFLAGS) -o endpointGain.exe endpointGain.o $(objects) $(LIBS)
	@echo "Done"

endpointTracker.exe: endpointTracker.o $(objects)
	$(CXX) $(CXXFLAGS) -o endpointTracker.exe endpointTracker.o $(objects) $(LIBS)
	@echo "Done"

DataTree.o: DataTree.hh

MBUtils.o: MBUtils.hh

BetaDecayTools.o: BetaDecayTools.hh

# -------------------------------------------------------------------------------
#  Generic compilation and linking step to make an executable from
#  a single *.cc file
#
#%: %.$(SrcSuf)
#	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)
#	@echo "$@ done"

clean:
		@rm -f *.o *~  core
