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

objects = DataTree.o MBUtils.o positionMapHandler.o ../source_peaks/peaks.o

position_map.exe: position_map.o $(objects)
	$(CXX) $(CXXFLAGS) -o position_map.exe position_map.o $(objects) $(LIBS)
	@echo "Done"

plot_position_map.exe: plot_position_map.o $(objects)
	$(CXX) $(CXXFLAGS) -o plot_position_map.exe plot_position_map.o $(objects) $(LIBS)
	@echo "Done"

plot_swank_map.exe: plot_swank_map.o $(objects)
	$(CXX) $(CXXFLAGS) -o plot_swank_map.exe plot_swank_map.o $(objects) $(LIBS)
	@echo "Done"

compare_swank_maps.exe: compare_swank_maps.o $(objects)
	$(CXX) $(CXXFLAGS) -o compare_swank_maps.exe compare_swank_maps.o $(objects) $(LIBS)
	@echo "Done"

DataTree.o: DataTree.hh

peaks.o: peaks.hh

MBUtils.o: MBUtils.hh

positionMapHandler.o: positionMapHandler.hh MBUtils.hh

# -------------------------------------------------------------------------------
#  Generic compilation and linking step to make an executable from
#  a single *.cc file
#
#%: %.$(SrcSuf)
#	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)
#	@echo "$@ done"

clean:
		@rm -f *.o *~  core
