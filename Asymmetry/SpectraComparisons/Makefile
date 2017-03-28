ObjSuf		= o
SrcSuf		= C
DllSuf		= so
OutPutOpt	= -o
ROOTLIBS	= $(shell root-config --libs)
ROOTGLIBS	= $(shell root-config --glibs)
ROOTCFLAGS	= $(shell root-config --cflags)

VPATH = ./:../../include

CXX		= g++
CXXFLAGS	= -O2 -Wall -fPIC -std=c++11 -g
LD		= g++
SOFLAGS		= -shared
LIBS		= $(ROOTLIBS) $(ROOTGLIBS)
CXXFLAGS       += $(ROOTCFLAGS) -I../../include
LIBS	       += -lSpectrum -lMinuit

objects = MBUtils.o


.PHONY: all
all: EastWestComparisons.exe

EastWestComparisons.exe: EastWestComparisons.o $(objects) 
	$(CXX) $(CXXFLAGS) -o EastWestComparisons.exe EastWestComparisons.o $(objects) $(LIBS) 
	@echo "Done"

MBUtils.o: MBUtils.hh

clean:
		@rm -f *.o *~  core
