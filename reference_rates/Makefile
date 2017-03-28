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
CXXFLAGS      = -O2 -Wall -fPIC -std=c++11 -g
LD            = g++
SOFLAGS       = -shared
LIBS          = $(ROOTLIBS) $(ROOTGLIBS)
CXXFLAGS     += $(ROOTCFLAGS) -I../include
LIBS         += -lSpectrum -lMinuit

.PHONY: all

objects = MBUtils.o DataTree.o calibrationTools.o

all: backgroundSpectra.exe foregroundSpectra.exe SIM_foregroundSpectra.exe

backgroundSpectra.exe: backgroundSpectra.o $(objects)
	$(CXX) $(CXXFLAGS) -o backgroundSpectra.exe backgroundSpectra.o $(objects) $(LIBS)
	@echo "Done"

foregroundSpectra.exe: foregroundSpectra.o $(objects)
	$(CXX) $(CXXFLAGS) -o foregroundSpectra.exe foregroundSpectra.o $(objects) $(LIBS)
	@echo "Done"

SIM_foregroundSpectra.exe: SIM_foregroundSpectra.o MBUtils.o calibrationTools.o 
	$(CXX) $(CXXFLAGS) -o SIM_foregroundSpectra.exe SIM_foregroundSpectra.o MBUtils.o calibrationTools.o $(LIBS)
	@echo "Done"

MBUtils.o: MBUtils.hh

DataTree.o: DataTree.hh

calibrationTools.o: calibrationTools.hh



# -------------------------------------------------------------------------------
#  Generic compilation and linking step to make an executable from
#  a single *.cc file
#
#%: %.$(SrcSuf)
#	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)
#	@echo "$@ done"

clean:
		@rm -f *.o *~  core
