#ifndef TRIGGERMAP_HH
#define TRIGGERMAP_HH

#include <vector>
#include <TF1.h>

class TriggerMap {

public:

  TriggerMap(Double_t bin_width, Int_t nParams); //Constructor
  ~TriggerMap();

  void readTriggerMap(Int_t XeRunPeriod); //Read in a trigger map
  Double_t returnTriggerProbabilityEast(Double_t x, Double_t y, Double_t En);
  Double_t returnTriggerProbabilityWest(Double_t x, Double_t y, Double_t En);
  //std::vector <Double_t> getTriggerFunction(Double_t x, Double_t y); //returns the trigger function for the given event position
  void writeTriggerMap(Int_t XeRunPeriod); //Write trigger map
  void setTriggerMapPointEast(Int_t xBin, Int_t yBin, std::vector <Double_t> params); //Set a trigger map val
  void setTriggerMapPointWest(Int_t xBin, Int_t yBin, std::vector <Double_t> params); //Set a trigger map val

  void setTriggerFunc(TF1* f);
  void setTriggerFunc(TF1 f);

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
  Int_t nBinsXY, nBinsTotal;
  Int_t XePeriod;
  Int_t nParams;
  std::vector <Double_t> xyBinLower;
  std::vector <Double_t> xyBinUpper;
  std::vector <Double_t> xyBinCenter;
  //std::vector <Double_t> yBinLower;
  //std::vector <Double_t> yBinUpper;
  //std::vector <Double_t> yBinCenter;

  TF1 *func;

  std::vector < std::vector < std::vector <Double_t> > > triggMapE;
  std::vector < std::vector < std::vector <Double_t> > > triggMapW;
  
 
};

  
#endif
