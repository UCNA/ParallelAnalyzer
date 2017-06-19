// Wrapper class to fill a histogram and iteratively fit the 
// histogram to get the best fit possible

#ifndef PEAKS_HH
#define PEAKS_HH

#include <TH1D.h>
#include <TF1.h>
#include <TString.h>

#include <string>

class SinglePeakHist {

public:

  SinglePeakHist(TH1D* h, Double_t rangeLow, Double_t rangeHigh, bool autoFit=true, Int_t iterations=5, Double_t minScaleFac=1.25, Double_t maxScaleFac=1.75, bool landau=false);
  ~SinglePeakHist();

  void FitHist(Double_t meanGuess=0., Double_t sigGuess=1., Double_t heightGuess=1.);

  Double_t ReturnMean()  { return mean; }
  Double_t ReturnMeanError()  { return meanErr; }
  Double_t ReturnSigma() { return sigma; }
  Double_t ReturnSigmaError() { return sigmaErr; }
  Double_t ReturnScale() { return scale; }
  bool isGoodFit()       { return goodFit_new; }
  
  void SetRangeMin(Double_t m) { min = m; }
  void SetRangeMax(Double_t m) { max = m; }

private:

  TH1D *hist;
  TF1 *func;
  Double_t mean;
  Double_t meanErr;
  Double_t sigma;
  Double_t sigmaErr;
  Double_t scale;
  Double_t min; 
  Double_t max;
  Double_t iters;
  bool goodFit_prev;  // the previous fit is good or bad
  bool goodFit_new;  // the new fit is good or bad
  TString status;
  Double_t minScaleFactor;
  Double_t maxScaleFactor;
  Bool_t _bLandau;
};


class DoublePeakHist {

public:

  DoublePeakHist(TH1D* h, Double_t rangeLow, Double_t rangeHigh, bool autoFit=true, Int_t iterations=6);
  ~DoublePeakHist();

  void FitHist(Double_t meanGuess1=0., Double_t sigGuess1=1., Double_t heightGuess1=1., Double_t meanGuess2=1., Double_t sigGuess2=1., Double_t heightGuess2=1.);

  Double_t ReturnMean1()  { return mean1; } //This is the highest peak (Upper Bi)
  Double_t ReturnMean1Error()  { return mean1Err; }
  Double_t ReturnSigma1() { return sigma1; }
  Double_t ReturnSigma1Error() { return sigma1Err; }
  Double_t ReturnScale1() { return scale1; }
  Double_t ReturnMean2()  { return mean2; }
  Double_t ReturnMean2Error()  { return mean2Err; }
  Double_t ReturnSigma2() { return sigma2; }
  Double_t ReturnSigma2Error() { return sigma2Err; }
  Double_t ReturnScale2() { return scale2; }
  bool isGoodFit()       { return goodFit_new; }

  void SetRangeMin(Double_t m) { min = m; }
  void SetRangeMax(Double_t m) { max = m; }

private:

  TH1D *hist;
  TF1 *func;
  Double_t mean1;
  Double_t mean1Err;
  Double_t sigma1;
  Double_t sigma1Err;
  Double_t scale1;
  Double_t mean2;
  Double_t mean2Err;
  Double_t sigma2;
  Double_t sigma2Err;
  Double_t scale2;
  Double_t min; 
  Double_t max;
  Double_t iters;
  bool goodFit_prev;  // the previous fit is good or bad
  bool goodFit_new;
  TString status;
};

#endif
