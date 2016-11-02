#ifndef CALTOOLS_HH
#define CALTOOLS_HH

#include <TF1.h>
#include <TMath.h>
#include <vector>
#include <fstream>
#include <cstdlib>


class LinearityCurve {
 
public:

  LinearityCurve(Int_t srcPeriod, bool useTanhSmear=true);  
  ~LinearityCurve(); 
  void readLinearityCurveParams(Int_t srcPeriod);
  std::vector < std::vector < Double_t > > returnLinCurveParams() { return pmtParams; };
  Double_t applyLinCurve(Int_t pmt, Double_t x); 
  Double_t applyInverseLinCurve(Int_t pmt, Double_t x); 
  Int_t getCurrentSrcCalPeriod() { return sourceCalPeriod; }; 
  
private:                                                                                                                      
  Int_t sourceCalPeriod;  
  bool useTanh;
  std::vector < std::vector < Double_t > > pmtParams;                             
  Int_t nParams;
  TF1 *linCurve;                                                                    
};        


class TriggerFunctions {

public:

  TriggerFunctions(Int_t run);
  ~TriggerFunctions();

  void readTriggerFunctions(Int_t run);
  Int_t returnCurrentRun() { return currentRun; }
  bool decideEastTrigger(std::vector<Double_t> ADC, std::vector<Double_t> rands);
  bool decideWestTrigger(std::vector<Double_t> ADC, std::vector<Double_t> rands);

private:

  Int_t currentRun;
  TF1 *func;
  
  std::vector < std::vector < Double_t > > triggerFuncs;


};
#endif