
#include <TMath.h>
#include <TRandom3.h>
#include <TMinuit.h>
#include <TString.h>
#include "peaks.hh"


SinglePeakHist::SinglePeakHist(TH1D* h, Double_t rangeLow, Double_t rangeHigh, bool autoFit, Int_t iterations, Double_t minScaleFac, Double_t maxScaleFac, bool landau) :
  hist(h), func(NULL), mean(0.), meanErr(0.), sigma(0.), sigmaErr(0.), scale(0.), min(rangeLow), max(rangeHigh), iters(iterations), goodFit_prev(false), goodFit_new(false), status("NO FITS") {

  minScaleFactor = minScaleFac;
  maxScaleFactor = maxScaleFac;

  _bLandau = landau;

  if (autoFit) {
    hist->GetXaxis()->SetRangeUser(min,max);
    mean = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
    scale = hist->GetBinContent(hist->GetMaximumBin());
    hist->GetXaxis()->SetRange(0,hist->GetNbinsX());
    sigma = 50.;
    Int_t it = 0;
    while (it<iters) {
      FitHist(mean, sigma, scale);
      if (goodFit_new && goodFit_prev && it>1) break;
      goodFit_prev = goodFit_new;
      it++;
    }

    if (!goodFit_new && meanErr==0.) meanErr = 100.;
    
    std::cout << "Finished Fitting after " << it+1 << " iterations with Status = " << status << "\n"; 
  }
}

SinglePeakHist::~SinglePeakHist() {
 
}

void SinglePeakHist::FitHist(Double_t meanGuess, Double_t sigGuess, Double_t heightGuess) {
  
  if ( _bLandau ) func = new TF1("func","[0]*TMath::Landau(x,[1],[2])",min, max);
  else func = new TF1("func","gaus(0)",min, max);

  func->SetParameters(heightGuess, meanGuess, sigGuess);
  
  hist->Fit(func, "LRQ");

  Double_t scaleCheck = func->GetParameter(0);  
  Double_t meanCheck = func->GetParameter(1);
  Double_t sigmaCheck = TMath::Abs(func->GetParameter(2));
  
  status = gMinuit->fCstatu;

  if (status==TString("CONVERGED ")) { 
    //std::cout << gMinuit->fCstatu << std::endl;
    goodFit_new=true;
    scale = scaleCheck;
    mean = meanCheck;
    meanErr = func->GetParError(1);
    sigma = sigmaCheck;
    sigmaErr = func->GetParError(2);
    min = mean-minScaleFactor*sigma;
    max = mean+maxScaleFactor*sigma;
  }

  else if ((TMath::Abs(meanCheck)>50.?TMath::Abs(meanCheck):50.)/sigmaCheck > 0.5 && (TMath::Abs(meanCheck)>50.?TMath::Abs(meanCheck):50.)/sigmaCheck < 500. && meanCheck > min && meanCheck < max) { 

    goodFit_new=false;
    scale = scaleCheck;
    mean = meanCheck;
    meanErr = func->GetParError(1);
    sigma = sigmaCheck;
    sigmaErr = func->GetParError(2);
    min = mean-sigma;
    max = mean+sigma;
    
  }
  else {
    goodFit_new = false;
    TRandom3 rand(0);
    mean = mean + rand.Gaus(mean, 20.); 
  }

  delete func;
}




DoublePeakHist::DoublePeakHist(TH1D* h, Double_t rangeLow, Double_t rangeHigh, bool autoFit, Int_t iterations) :
  hist(h), func(NULL), mean1(0.), mean1Err(0.), sigma1(0.), sigma1Err(0.), scale1(0.), mean2(0.), mean2Err(0.), sigma2(0.), sigma2Err(0.), scale2(0.), min(rangeLow), max(rangeHigh), iters(iterations), goodFit_prev(false), goodFit_new(false), status("NO FITS") {

  if (autoFit) {
    mean1 = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
    scale1 = hist->GetBinContent(hist->GetMaximumBin());
    sigma1 = 50.;

    mean2 = mean1*0.5;
    scale2 = scale1*0.5; //hist->GetBinContent(hist->GetXaxis()->FindBin(mean2));
    sigma2 = 50.;

    Int_t it = 0;
    while (it<iters) {
      FitHist(mean1, sigma1, scale1, mean2, sigma2, scale2);
      if (goodFit_new && goodFit_prev) break;
      goodFit_prev = goodFit_new;
      it++;
    }
    
    if (!goodFit_new && mean1Err==0.) mean1Err = 100.;
    if (!goodFit_new && mean2Err==0.) mean2Err = 100.;

    std::cout << "Finished Fitting after " << it << " iterations with Status = " << status << "\n"; 
  }
}

DoublePeakHist::~DoublePeakHist() {
 
}

void DoublePeakHist::FitHist(Double_t meanGuess1, Double_t sigGuess1, Double_t heightGuess1, Double_t meanGuess2, Double_t sigGuess2, Double_t heightGuess2) {
  
  func = new TF1("func","gaus(0)+gaus(3)",min, max);

  func->SetParameters(heightGuess1, meanGuess1, sigGuess1, heightGuess2, meanGuess2, sigGuess2);
  
  hist->Fit(func, "LRQ");

  Int_t highPeak = func->GetParameter(1) > func->GetParameter(4) ? 0 : 1;

  Double_t scale1Check,mean1Check, sigma1Check, scale2Check, mean2Check, sigma2Check; 

  if (highPeak==0) {
    scale1Check = func->GetParameter(0);  
    mean1Check = func->GetParameter(1);
    sigma1Check = TMath::Abs(func->GetParameter(2));
    scale2Check = func->GetParameter(3);  
    mean2Check = func->GetParameter(4);
    sigma2Check = TMath::Abs(func->GetParameter(5));
  }
  else {
    scale1Check = func->GetParameter(3);  
    mean1Check = func->GetParameter(4);
    sigma1Check = TMath::Abs(func->GetParameter(5));
    scale2Check = func->GetParameter(0);  
    mean2Check = func->GetParameter(1);
    sigma2Check = TMath::Abs(func->GetParameter(2));
  }

  status = gMinuit->fCstatu;

  if (status==TString("CONVERGED ")) { 
    //std::cout << gMinuit->fCstatu << std::endl;
    goodFit_new=true;
    scale1 = scale1Check;
    mean1 = mean1Check;
    mean1Err = highPeak==0 ? func->GetParError(1) : func->GetParError(4);
    sigma1 = sigma1Check;
    sigma1Err = highPeak==0 ? func->GetParError(2) : func->GetParError(5);
    scale2 = scale2Check;
    mean2 = mean2Check;
    mean2Err = highPeak==0 ? func->GetParError(4) : func->GetParError(1);
    sigma2 = sigma2Check;
    sigma2Err = highPeak==0 ? func->GetParError(5) : func->GetParError(2);
    min = mean2-2.25*sigma2;
    max = mean1+2.0*sigma1;    

  }

  else if ((TMath::Abs(mean1Check) > 50. ? (TMath::Abs(mean1Check)) : 50.)/sigma1Check > 0.5 && (TMath::Abs(mean1Check)>50. ? (TMath::Abs(mean1Check)) : 50.)/sigma1Check < 500. && mean1Check > min && mean1Check < max
	   && (TMath::Abs(mean2Check)>50. ? (TMath::Abs(mean2Check)) : 50.)/sigma2Check > 0.5 && (TMath::Abs(mean2Check)>50. ? (TMath::Abs(mean2Check)) : 50.)/sigma2Check < 500. && mean2Check > min && mean2Check < max) { 

    goodFit_new=false;
    scale1 = scale1Check;
    mean1 = mean1Check;
    mean1Err = highPeak==0 ? func->GetParError(1) : func->GetParError(4);
    sigma1 = sigma1Check;
    sigma1Err = highPeak==0 ? func->GetParError(2) : func->GetParError(5);
    scale2 = scale2Check;
    mean2 = mean2Check;
    mean2Err = highPeak==0 ? func->GetParError(4) : func->GetParError(1);
    sigma2 = sigma2Check;
    sigma2Err = highPeak==0 ? func->GetParError(5) : func->GetParError(2);
    min = mean2-2.25*sigma2;
    max = mean1+2.0*sigma1;    

  }

  else {
    goodFit_new = false;
    TRandom3 rand(0);
    mean1 = mean1 + rand.Gaus(mean1, 20.); 
    mean2 = mean2 + rand.Gaus(mean2, 20.); 
  }
  delete func;
}
