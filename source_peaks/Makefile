# Makefile from Jianglai

ObjSuf        = o
SrcSuf        = C
DllSuf        = so
OutPutOpt     = -o  
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
ROOTCFLAGS    = $(shell root-config --cflags)

VPATH = ./:../include

# Work with Linux with egcs	
CXX           = g++ 
CXXFLAGS      = -O2 -Wall -fPIC -std=c++11
LD            = g++
SOFLAGS       = -shared
LIBS          = $(ROOTLIBS) $(ROOTGLIBS)
CXXFLAGS     += $(ROOTCFLAGS) -I../include
LIBS         += -lSpectrum -lMinuit

objects = DataTree.o MBUtils.o peaks.o positionMapHandler.o calibrationTools.o MWPCPositionResponse.o

.PHONY: all
all: srcCalReplay.exe sim_source_peak_fitter.exe
#source_peaks.exe source_peaks_EvisPMTbyPMT.exe sim_source_peaks_EvisPMTbyPMT.exe sim_source_peaks_EnergyPeak.exe source_peaks_EnergyPeak.exe 

source_peaks.exe: source_peaks.o $(objects)
	$(CXX) $(CXXFLAGS) -o source_peaks.exe source_peaks.o $(objects) $(LIBS)
	@echo "Done"

srcCalReplay.exe: srcCalReplay.o $(objects)
	$(CXX) $(CXXFLAGS) -o srcCalReplay.exe srcCalReplay.o $(objects) $(LIBS)
	@echo "Done"

sim_source_peak_fitter.exe: sim_source_peak_fitter.o $(objects)
	$(CXX) $(CXXFLAGS) -o sim_source_peak_fitter.exe sim_source_peak_fitter.o $(objects) $(LIBS)
	@echo "Done"


source_peaks_postCal.exe: source_peaks_postCal.o $(objects)
	$(CXX) $(CXXFLAGS) -o source_peaks_postCal.exe source_peaks_postCal.o $(objects) $(LIBS)
	@echo "Done"

sim_source_peaks.exe: sim_source_peaks.o $(objects)
	$(CXX) $(CXXFLAGS) -o sim_source_peaks.exe sim_source_peaks.o $(objects) $(LIBS)
	@echo "Done"

sim_source_peaks_postCal.exe: sim_source_peaks_postCal.o $(objects)
	$(CXX) $(CXXFLAGS) -o sim_source_peaks_postCal.exe sim_source_peaks_postCal.o $(objects) $(LIBS)
	@echo "Done"

source_peaks_EvisPMTbyPMT.exe: source_peaks_EvisPMTbyPMT.o $(objects)
	$(CXX) $(CXXFLAGS) -o source_peaks_EvisPMTbyPMT.exe source_peaks_EvisPMTbyPMT.o $(objects) $(LIBS)
	@echo "Done"

sim_source_peaks_EvisPMTbyPMT.exe: sim_source_peaks_EvisPMTbyPMT.o
	$(CXX) $(CXXFLAGS) -o sim_source_peaks_EvisPMTbyPMT.exe sim_source_peaks_EvisPMTbyPMT.o $(LIBS)
	@echo "Done"

sim_source_peaks_EnergyPeak.exe: sim_source_peaks_EnergyPeak.o
	$(CXX) $(CXXFLAGS) -o sim_source_peaks_EnergyPeak.exe sim_source_peaks_EnergyPeak.o $(LIBS)
	@echo "Done"

source_peaks_EnergyPeak.exe: source_peaks_EnergyPeak.o $(objects)
	$(CXX) $(CXXFLAGS) -o source_peaks_EnergyPeak.exe source_peaks_EnergyPeak.o $(objects) $(LIBS)
	@echo "Done"

peaks.o: peaks.hh

MWPCPositionResponse.o: MWPCPositionResponse.hh

DataTree.o: DataTree.hh

MBUtils.o: MBUtils.hh

positionMapHandler.o: positionMapHandler.hh MBUtils.o

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
