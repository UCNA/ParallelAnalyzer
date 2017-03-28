#ifndef CALTOOLS_HH
#define CALTOOLS_HH

#include <TF1.h>
#include <TString.h>
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

class WirechamberCal {
 
public:

  WirechamberCal(Int_t run);  
  ~WirechamberCal(); 
  void readWirechamberCalParams(Int_t run);
  std::vector < std::vector < Double_t > > returnWirechamberCalParams() { return _params; };
  Double_t applyCal(Int_t side, Double_t adc); 
  Int_t getCurrentRun() { return _run; }; 
  
private:                                                                                                                      
  Int_t _run;  
  std::vector < std::vector < Double_t > > _params;                             
  Int_t _nParams;
  TF1 *_calFunc;
  TF1 *_extrap; // Function to extrapolate to zero ADC
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

class EreconParameterization {
 
public:
  EreconParameterization(Int_t runNumber);  
  ~EreconParameterization() {}; 
  
  std::vector < std::vector < std::vector < Double_t > > > returnLinCurveParams() { return params; };
  Double_t getErecon(Int_t side, Int_t type, Double_t Evis); 
  TString getCurrentGeometry() { return geometry; }; 
  
private:
  TString geometry;
  std::vector < std::vector < std::vector < Double_t > > >  params;                             
  Int_t nParams;

  void readParams();
};


class BackscatterSeparator {

public:
  BackscatterSeparator();
  ~BackscatterSeparator() {};
  void LoadCutCurve(int run);

  Int_t separate23(Double_t en); 

private:
  Int_t _run;            // Current run
  TString _geometry;   // Geometry being used
  TF1 *_func;          // Function which parametrizes the cut
  std::vector <Double_t> params; // Parameters for the func
  
};
  

#endif
