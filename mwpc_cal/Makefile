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

objects = DataTree.o MBUtils.o positionMapHandler.o peaks.o MWPCPositionResponse.o

cathode_gain.exe: cathode_gain.o $(objects)
	$(CXX) $(CXXFLAGS) -o cathode_gain.exe cathode_gain.o $(objects) $(LIBS)
	@echo "Done"

cathode_model.exe: cathode_model.o $(objects)
	$(CXX) $(CXXFLAGS) -o cathode_model.exe cathode_model.o $(objects) $(LIBS)
	@echo "Done"

MWPC_position_map.exe: MWPC_position_map.o $(objects)
	$(CXX) $(CXXFLAGS) -o MWPC_position_map.exe MWPC_position_map.o $(objects) $(LIBS)
	@echo "Done"

plot_MWPC_position_map.exe: plot_MWPC_position_map.o $(objects)
	$(CXX) $(CXXFLAGS) -o plot_MWPC_position_map.exe plot_MWPC_position_map.o $(objects) $(LIBS)
	@echo "Done"

MWPC_calibration.exe: MWPC_calibration.o $(objects)
	$(CXX) $(CXXFLAGS) -o MWPC_calibration.exe MWPC_calibration.o $(objects) $(LIBS)
	@echo "Done"

BackscatterSeparator.exe: BackscatterSeparator.o $(objects)
	$(CXX) $(CXXFLAGS) -o BackscatterSeparator.exe BackscatterSeparator.o $(objects) $(LIBS)
	@echo "Done"

DataTree.o: DataTree.hh

MBUtils.o: MBUtils.hh

positionMapHandler.o: positionMapHandler.hh MBUtils.o

peaks.o: peaks.hh MBUtils.o

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
