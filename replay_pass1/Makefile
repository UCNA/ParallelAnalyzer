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

objects = WireChamberResponse.o DataTree.o MWPCPositionResponse.o
.PHONY: all
all: replay_pass1.exe

replay_pass1.exe: replay_pass1.o $(objects)
	$(CXX) $(CXXFLAGS) -o replay_pass1.exe replay_pass1.o $(objects) $(LIBS)
	@echo "Done"

WireChamberResponse.o: WireChamberResponse.h

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
