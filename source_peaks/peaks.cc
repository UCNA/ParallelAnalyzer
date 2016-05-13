
#include <TMath.h>
#include <TRandom3.h>

#include "peaks.hh"


SinglePeakHist::SinglePeakHist(TH1D* h, Double_t rangeLow, Double_t rangeHigh, bool autoFit, Int_t iterations) :
  hist(h), func(NULL), mean(0.), sigma(0.), scale(0.), min(rangeLow), max(rangeHigh), iters(iterations), goodFit(false) {

  if (autoFit) {
    mean = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
    scale = hist->GetBinContent(hist->GetMaximumBin());
    sigma = 50.;
    for (Int_t it=0; it<iters; it++) {
      FitHist(mean, sigma, scale);
    }
    
    std::cout << "Finished Fitting\n"; 
  }
}

SinglePeakHist::~SinglePeakHist() {
 
}

void SinglePeakHist::FitHist(Double_t meanGuess, Double_t sigGuess, Double_t heightGuess) {
  
  func = new TF1("func","gaus(0)",min, max);

  func->SetParameters(heightGuess, meanGuess, sigGuess);
  
  hist->Fit(func, "R");

  Double_t scaleCheck = func->GetParameter(0);  
  Double_t meanCheck = func->GetParameter(1);
  Double_t sigmaCheck = TMath::Abs(func->GetParameter(2));
  delete func;

  if ((TMath::Abs(meanCheck)>50.?TMath::Abs(meanCheck):50.)/sigmaCheck > 0.5 && (TMath::Abs(meanCheck)>50.?TMath::Abs(meanCheck):50.)/sigmaCheck < 500. && meanCheck > min && meanCheck < max) { 

    goodFit=true;
    scale = scaleCheck;
    mean = meanCheck;
    sigma = sigmaCheck;
    min = mean-sigma;
    max = mean+sigma;
    
  }
  else {
    goodFit = false;
    TRandom3 rand(0);
    mean = mean + rand.Gaus(mean, 20.); 
  }

}




DoublePeakHist::DoublePeakHist(TH1D* h, Double_t rangeLow, Double_t rangeHigh, bool autoFit, Int_t iterations) :
  hist(h), func(NULL), mean1(0.), sigma1(0.), scale1(0.), mean2(0.), sigma2(0.), scale2(0.), min(rangeLow), max(rangeHigh), iters(iterations), goodFit(false) {

  if (autoFit) {
    mean1 = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
    scale1 = hist->GetBinContent(hist->GetMaximumBin());
    sigma1 = 50.;

    mean2 = mean1*0.5;
    scale2 = scale1*0.5; //hist->GetBinContent(hist->GetXaxis()->FindBin(mean2));
    sigma2 = 50.;

    for (Int_t it=0; it<iters; it++) {
      FitHist(mean1, sigma1, scale1, mean2, sigma2, scale2);
    }
    
    std::cout << "Finished Fitting\n"; 
  }
}

DoublePeakHist::~DoublePeakHist() {
 
}

void DoublePeakHist::FitHist(Double_t meanGuess1, Double_t sigGuess1, Double_t heightGuess1, Double_t meanGuess2, Double_t sigGuess2, Double_t heightGuess2) {
  
  func = new TF1("func","gaus(0)+gaus(3)",min, max);

  func->SetParameters(heightGuess1, meanGuess1, sigGuess1, heightGuess2, meanGuess2, sigGuess2);
  
  hist->Fit(func, "R");

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

  delete func;

  if ((TMath::Abs(mean1Check)>50.?TMath::Abs(mean1Check):50.)/sigma1Check > 0.5 && (TMath::Abs(mean1Check)>50.?TMath::Abs(mean1Check):50.)/sigma1Check < 500. && mean1Check > min && mean1Check < max
      && (TMath::Abs(mean2Check)>50.?TMath::Abs(mean2Check):50.)/sigma2Check > 0.5 && (TMath::Abs(mean2Check)>50.?TMath::Abs(mean2Check):50.)/sigma2Check < 500. && mean2Check > min && mean2Check < max) { 

    goodFit=true;
    scale1 = scale1Check;
    mean1 = mean1Check;
    sigma1 = sigma1Check;
    scale2 = scale2Check;
    mean2 = mean2Check;
    sigma2 = sigma2Check;
    min = mean2-sigma2;
    max = mean1+sigma1;    

  }

  else {
    goodFit = false;
    TRandom3 rand(0);
    mean1 = mean1 + rand.Gaus(mean1, 20.); 
    mean2 = mean2 + rand.Gaus(mean2, 20.); 
  }

}