#include <vector>
#include <TF1.h>

static const bool useRC = false;

class PositionMap {

public:

  PositionMap(Double_t bin_width, Double_t r_max); //Constructor
  ~PositionMap();

  void readPositionMap(Int_t XeRunPeriod); //Read in a trigger map
  //void writePositionMap(Int_t XeRunPeriod); //Write trigger map
  void setPositionMapPoint(Int_t xBin, Int_t yBin, std::vector <Double_t> vals); //Set a trigger map val

  void setRCflag(bool rc) { bUseRC = rc; }

  std::vector <Double_t> getInterpolatedEta(Double_t xE, Double_t yE, Double_t xW, Double_t yW); //return interpolated value from position map

  Int_t getBinNumber(Double_t pos); //returns bin number for position

  Double_t getBinWidth() {return binWidth;}
  Int_t getNbins() {return nBinsTotal;}
  Int_t getNbinsXY() {return nBinsXY;}
  Double_t getBinUpper(Int_t bin) {return (bin >= 0 && bin < nBinsTotal) ? xyBinUpper[bin] : -1000. ;}
  Double_t getBinLower(Int_t bin) {return (bin >= 0 && bin < nBinsTotal) ? xyBinLower[bin] : -1000.;}
  Double_t getBinCenter(Int_t bin) {return (bin >= 0 && bin < nBinsTotal) ? xyBinCenter[bin] : -1000.;}
  Int_t getCurrentXeRunPeriod() {return XePeriod;}

private:

  void setBinValues();

  Double_t binWidth; //x and y width of position bins
  Double_t rmax; //Maximum r where an interpolation can faithfully be constructed
  bool bUseRC; // Flag which defaults to const value at top of file, but can be set to whatever
  
  Int_t nBinsXY, nBinsTotal;
  Int_t XePeriod;
  std::vector <Double_t> xyBinLower;
  std::vector <Double_t> xyBinUpper;
  std::vector <Double_t> xyBinCenter;
  std::vector <TF1*> funcs;
  //std::vector <Double_t> yBinLower;
  //std::vector <Double_t> yBinUpper;
  //std::vector <Double_t> yBinCenter;

  std::vector < std::vector < std::vector < Double_t > > > posMap;
 
};

class MWPCPositionMap {

public:

  MWPCPositionMap(Double_t bin_width, Double_t r_max); //Constructor
  ~MWPCPositionMap();

  void readMWPCPositionMap(Int_t XeRunPeriod, double elow, double ehigh); //Read in a trigger map
  //void writeMWPCPositionMap(Int_t XeRunPeriod); //Write trigger map
  void setMWPCPositionMapPoint(Int_t xBin, Int_t yBin, std::vector <Double_t> vals); //Set a trigger map val
  
  std::vector <Double_t> getInterpolatedEta(Double_t xE, Double_t yE, Double_t xW, Double_t yW); //return interpolated value from position map

  Int_t getBinNumber(Double_t pos); //returns bin number for position

  Double_t getBinWidth() {return binWidth;}
  Int_t getNbins() {return nBinsTotal;}
  Int_t getNbinsXY() {return nBinsXY;}
  Double_t getBinUpper(Int_t bin) {return (bin >= 0 && bin < nBinsTotal) ? xyBinUpper[bin] : -1000. ;}
  Double_t getBinLower(Int_t bin) {return (bin >= 0 && bin < nBinsTotal) ? xyBinLower[bin] : -1000.;}
  Double_t getBinCenter(Int_t bin) {return (bin >= 0 && bin < nBinsTotal) ? xyBinCenter[bin] : -1000.;}
  Int_t getCurrentXeRunPeriod() {return XePeriod;}

private:

  void setBinValues();

  Double_t binWidth; //x and y width of position bins
  Double_t rmax; //Maximum r where an interpolation can faithfully be constructed
  
  Int_t nBinsXY, nBinsTotal;
  Int_t XePeriod;
  std::vector <Double_t> xyBinLower;
  std::vector <Double_t> xyBinUpper;
  std::vector <Double_t> xyBinCenter;
  std::vector <TF1*> funcs;
  //std::vector <Double_t> yBinLower;
  //std::vector <Double_t> yBinUpper;
  //std::vector <Double_t> yBinCenter;

  std::vector < std::vector < std::vector < Double_t > > > posMap;
 
};
