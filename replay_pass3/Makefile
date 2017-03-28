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

objects = DataTree.o positionMapHandler.o MBUtils.o calibrationTools.o  MWPCPositionResponse.o

replay_pass3.exe: replay_pass3.o $(objects)
	$(CXX) $(CXXFLAGS) -o replay_pass3.exe replay_pass3.o $(objects) $(LIBS)
	@echo "Done"

DataTree.o: DataTree.hh

positionMapHandler.o: positionMapHandler.hh MBUtils.o

MBUtils.o: MBUtils.hh

calibrationTools.o: calibrationTools.hh

MWPCPositionResponse.o: MWPCPositionResponse.hh


#../include/positionMapHandler.o: ../include/positionMapHandler.hh ../include/MBUtils.o
#	$(CXX) $(CXXFLAGS) -c ../include/positionMapHandler.o ../include/positionMapHandler.cc ../include/MBUtils.o $(LIBS)
#	@echo "Done"
# -------------------------------------------------------------------------------
#  Generic compilation and linking step to make an executable from
#  a single *.cc file
#
#%: %.$(SrcSuf)
#	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)
#	@echo "$@ done"

clean:
		@rm -f *.o *~  core
