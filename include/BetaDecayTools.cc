#include "BetaDecayTools.hh"

#include <iostream>
#include <fstream>


Double_t calcTotalEn_Wilk(Double_t KE) { return KE>0. ? (m_e + KE)/m_e : 1.; };
Double_t calcMomentum_Wilk(Double_t W) { return TMath::Sqrt( W*W - 1 ); };

KurieFitter::KurieFitter() : _kuriePlot(NULL), _alpha(1.), _W0(0.), _W0err(0.), _actualW0(beta_W0) { };

KurieFitter::~KurieFitter() {

  if ( _kuriePlot ) delete _kuriePlot;

};

void KurieFitter::FitSpectrum(const TH1D* spec, Double_t min, Double_t max, Double_t alpha) {

  if ( _kuriePlot ) delete _kuriePlot; // Start with fresh kurie histogram

  _alpha = alpha;
  
  Int_t nBins = spec->GetNbinsX();

  // For plotting, these are the value and the error
  std::vector < std::vector< Double_t> > W_vals( 2, std::vector<Double_t>(nBins,0.) );
  std::vector < std::vector< Double_t> > K_vals( 2, std::vector<Double_t>(nBins,0.) );
  

  for ( Int_t bin=1; bin<=nBins; ++bin ) {
    
    Double_t W = calcTotalEn_Wilk( alpha * spec->GetXaxis()->GetBinCenter(bin) );
    Double_t p = calcMomentum_Wilk( W );
    
    K_vals[0][bin-1] = ( spec->GetBinContent(bin)>0. && p*W>0 ? 
			 TMath::Sqrt( spec->GetBinContent(bin) / ( p*W ) ) :
			 0. );
    K_vals[1][bin-1] = ( K_vals[0][bin-1]>0. ? 
			 TMath::Abs( 0.5*spec->GetBinError(bin) / TMath::Sqrt(TMath::Abs(spec->GetBinContent(bin))*p*W) )
			 : 0. );
    W_vals[0][bin-1] = W;
    W_vals[1][bin-1] = 0.;

    /*std::cout << "**** Bin " << bin << " ****\n"
	      << "TotalEn W = " << W << "\nMomentum p = " << p << "\n"
	      << "Kurie Val = " << K_vals[0][bin-1] << " +/- " << K_vals[1][bin-1] << "\n" ;*/
    
  }

  _kuriePlot = new TGraphErrors(nBins,&W_vals[0][0],&K_vals[0][0],&W_vals[1][0],&K_vals[1][0]);
  
  TF1 *lin = new TF1("lin","[0]*([1]-x)",calcTotalEn_Wilk(min),calcTotalEn_Wilk(max));
  lin->SetParameter(1,_actualW0);
  lin->SetParameter(0,50.);
  
  _kuriePlot->Fit(lin,"R");
  _W0 = lin->GetParameter(1);
  _W0err = lin->GetParError(1);
  

};

void KurieFitter::IterativeKurie(const TH1D* spec, Double_t min, Double_t max,
				 Double_t alphaStart, Double_t delta) {

  Double_t alpha_new;
  alpha_new = alphaStart;
  Double_t diff = 100.;
  Int_t iters = 0;

  while ( diff > delta && iters<10 ) {

    std::cout << "\n\nNew Alpha: " << alpha_new << "Diff: " << diff << std::endl;
    
    FitSpectrum(spec,min,max,alpha_new);
    
    alpha_new = ( _actualW0 - 1 ) / ( _W0 - 1 ) * _alpha;
    diff = TMath::Abs( 1. - alpha_new/_alpha );    
    iters++;
  }

  std::cout << "Finished after " << iters << " iterations\n";
  
};
