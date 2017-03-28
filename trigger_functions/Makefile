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

triggerThreshold.exe: triggerThreshold.o positionMapHandler.o SQLinterface.o calibrationTools.o MBUtils.o TriggerMap.o
	$(CXX) $(CXXFLAGS) -o triggerThreshold.exe triggerThreshold.o MBUtils.o calibrationTools.o TriggerMap.o $(LIBS)
	@echo "Done"

findADCthreshold.exe: findADCthreshold.o
	$(CXX) $(CXXFLAGS) -o findADCthreshold.exe findADCthreshold.o $(LIBS)
	@echo "Done"

findADCthreshold_singleRun.exe: findADCthreshold_singleRun.o DataTree.o
	$(CXX) $(CXXFLAGS) -o findADCthreshold_singleRun.exe findADCthreshold_singleRun.o DataTree.o $(LIBS)
	@echo "Done"

findDiscriminatorThresh.exe: findDiscriminatorThresh.o MBUtils.o 
	$(CXX) $(CXXFLAGS) -o findDiscriminatorThresh.exe findDiscriminatorThresh.o MBUtils.o $(LIBS)
	@echo "Done"

MBUtils.o: MBUtils.hh

calibrationTools.o: calibrationTools.hh

TriggerMap.o: TriggerMap.hh MBUtils.o

DataTree.o: DataTree.hh



# -------------------------------------------------------------------------------
#  Generic compilation and linking step to make an executable from
#  a single *.cc file
#
#%: %.$(SrcSuf)
#	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)
#	@echo "$@ done"

clean:
		@rm -f *.o *~  core
