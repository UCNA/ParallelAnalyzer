#include "TriggerMap.hh"

#include <string>
#include <fstream>
#include <cstdlib>
#include "MBUtils.hh"

#include <TString.h>

TriggerMap::TriggerMap(Double_t bin_width, Int_t nparams) : binWidth(bin_width), nParams(nparams), func(NULL) {
  
  Int_t nBinsHold = (int)(120./binWidth);
  nBinsXY = ((int)nBinsHold%2)==0 ? (int)nBinsHold+1 : (int)nBinsHold;
  nBinsTotal = nBinsXY*nBinsXY;

  xyBinLower.resize(nBinsXY,0.);
  xyBinUpper.resize(nBinsXY,0.);
  xyBinCenter.resize(nBinsXY,0.);
  

  triggMapE.resize(nBinsXY, std::vector < std::vector < Double_t> > (nBinsXY, std::vector <Double_t> (nParams,0.)));
  triggMapW.resize(nBinsXY, std::vector < std::vector <Double_t> > (nBinsXY, std::vector <Double_t> (nParams,0.)));
  
  setBinValues();

};

TriggerMap::~TriggerMap() {
  if (func) delete func;
};

void TriggerMap::setBinValues() {
  for (Int_t k=0; k<nBinsXY; k++) {
    xyBinLower[k]     = -(double)nBinsXY*binWidth/2. + ((double) k)*binWidth;
    xyBinUpper[k]     = -(double)nBinsXY*binWidth/2. + ((double) k)*binWidth + binWidth;
    xyBinCenter[k]    = (xyBinLower[k] + xyBinUpper[k])/2.;
    
    //cout << xBinLower[k] << " " << intXBinCenter[k] << " " << xBinUpper[k] << endl;
  }
};

void TriggerMap::readTriggerMap(Int_t XeRunPeriod) {
  XePeriod = XeRunPeriod;
  std::string file = std::string(getenv("TRIGGER_FUNC")) + "/trigger_functions_XePeriod_"+itos(XePeriod)+"_"+ftos(binWidth)+"mm_East" + ".dat";
  std::ifstream infile(file.c_str());

  Int_t binValx, binValy;
  std::vector <Double_t> p(8,0.);
  
  while (infile >> binValx >> binValy >> p[0] >> p[1] >> p[2] >> p[3] >> p[4] >> p[5] >> p[6] >> p[7])
    {
      setTriggerMapPointEast(getBinNumber(binValx), getBinNumber(binValy), p);
    }
  infile.close();

  file = std::string(getenv("TRIGGER_FUNC")) + "/trigger_functions_XePeriod_"+itos(XePeriod)+"_"+ftos(binWidth)+"mm_West" + ".dat";
  infile.open(file.c_str());
  
  while (infile >> binValx >> binValy >> p[0] >> p[1] >> p[2] >> p[3] >> p[4] >> p[5] >> p[6] >> p[7])
    {
      setTriggerMapPointWest(getBinNumber(binValx), getBinNumber(binValy), p);
    }
  infile.close();

  std::cout << "Read in Trigger Map from Xe Period " << file << std::endl;
};

void TriggerMap::setTriggerFunc(TF1* f) {
  if (func) delete func;
  func = new TF1(*f); 
  
};
void TriggerMap::setTriggerFunc(TF1 f) {
  if (func) delete func;
  func = new TF1(f);
};

Double_t TriggerMap::returnTriggerProbabilityEast(Double_t x, Double_t y, Double_t En) {
  
  if (!func) { std::cout << "No Function Set..\n"; return 0.; }

  if (En>80.) return 1.;

  else {
    Int_t binX = getBinNumber(x);
    Int_t binY = getBinNumber(y);
    if (binX==-1 || binY==-1) return 0.;
    else return func->EvalPar(&En, &triggMapW[binX][binY][0]);
  }

}

Double_t TriggerMap::returnTriggerProbabilityWest(Double_t x, Double_t y, Double_t En) {
  
  if (!func) { std::cout << "No Function Set..\n"; return 0.; }

  if (En>70.) return 1.;

  else {
    Int_t binX = getBinNumber(x);
    Int_t binY = getBinNumber(y);
    if (binX==-1 || binY==-1) return 0.;
    else return func->EvalPar(&En, &triggMapW[binX][binY][0]);
  }
}

Int_t TriggerMap::getBinNumber(Double_t pos) {
  for (Int_t m=0; m<nBinsXY; m++) {
    if ( (pos >= xyBinLower[m]) && (pos < xyBinUpper[m]) )  return m;
  }
  return -1;
};

void TriggerMap::setTriggerMapPointEast(Int_t xBin, Int_t yBin, std::vector <Double_t> params) {
  triggMapE[xBin][yBin] = params;
};

void TriggerMap::setTriggerMapPointWest(Int_t xBin, Int_t yBin, std::vector <Double_t> params) {
  triggMapW[xBin][yBin] = params;
};


void TriggerMap::writeTriggerMap(Int_t XeRunPeriod) {
  XePeriod = XeRunPeriod;
  std::string file = std::string(getenv("TRIGGER_FUNC")) + "/trigger_functions_XePeriod_"+itos(XePeriod)+"_"+ftos(binWidth)+"mm_East" + ".dat";
  std::cout << file << std::endl;
  std::ofstream ofile(file.c_str());
  
  for (Int_t xbin=0; xbin<nBinsXY; xbin++) {
    for (Int_t ybin=0; ybin<nBinsXY; ybin++) {

      ofile << xyBinCenter[xbin] << " " 
	    << xyBinCenter[ybin] << " ";
      
      for (Int_t n=0; n<nParams; n++) {
	
	ofile << triggMapE[xbin][ybin][n];
	if (n!=nParams-1) ofile << " ";
      }
      
      if (xbin!=nBinsXY-1 && ybin!=nBinsXY-1) ofile << std::endl;
    }	   
  }
  ofile.close();

  file = std::string(getenv("TRIGGER_FUNC")) + "/trigger_functions_XePeriod_"+itos(XePeriod)+"_"+ftos(binWidth)+"mm_West" + ".dat";
  ofile.open(file.c_str());
  
  for (Int_t xbin=0; xbin<nBinsXY; xbin++) {
    for (Int_t ybin=0; ybin<nBinsXY; ybin++) {

      ofile << xyBinCenter[xbin] << " " 
	    << xyBinCenter[ybin] << " ";
      
      for (Int_t n=0; n<nParams; n++) {
	
	ofile << triggMapW[xbin][ybin][n];
	if (n!=nParams-1) ofile << " ";
      }
      
      if (xbin!=nBinsXY-1 && ybin!=nBinsXY-1) ofile << std::endl;
    }	   
  }
  ofile.close();
}
