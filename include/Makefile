# Makefile from Jianglai

ObjSuf        = o
SrcSuf        = C
DllSuf        = so
OutPutOpt     = -o  
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
ROOTCFLAGS    = $(shell root-config --cflags)

VPATH = ./:../include/:../Asymmetry/

# Work with Linux with egcs	
CXX           = g++ 
CXXFLAGS      = -O2 -Wall -fPIC
LD            = g++
SOFLAGS       = -shared
LIBS          = $(ROOTLIBS) $(ROOTGLIBS)
CXXFLAGS     += $(ROOTCFLAGS) -I../include
LIBS         += -lSpectrum -lMinuit

positionMapHandler.exe: positionMapHandler.o positionMapHandler.o SQLinterface.o ../include/MBUtils.o
	$(CXX) $(CXXFLAGS) -o positionMapHandler.exe positionMapHandler.o MBUtils.o ../Asymmetry/SQLinterface.o $(LIBS)
	@echo "Done"

#../include/positionMapHandler.o: positionMapHandler.hh MBUtils.o

MBUtils.o: MBUtils.hh

../Asymmetry/SQLinterface.o: ../Asymmetry/SQLinterface.hh

#posMapReader.o: posMapReader.h



# -------------------------------------------------------------------------------
#  Generic compilation and linking step to make an executable from
#  a single *.cc file
#
#%: %.$(SrcSuf)
#	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)
#	@echo "$@ done"

clean:
		@rm -f *.o *~  core
