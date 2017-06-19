#ifndef BETADECAYTOOLS_HH
#define BETADECAYTOOLS_HH

#include <TH1D.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TMath.h>
#include <TF1.h>

const double neutronBetaEp = 782.347;			///< neutron beta decay endpoint, keV
const double m_e = 511.00;						///< electron mass, keV/c^2
const double m_p = 938272.046;					///< proton mass, keV/c^2
const double m_n = m_p+m_e+neutronBetaEp;		///< neutron mass, keV/c^2
const double fs_alpha = 1./137.036;				///< fine structure constant
const double lambda = fabs(-1.2723);			///< +/-0.0023, PDG 2016 value, Wilkinson sign convention
//const double A0_PDG = -0.1173;					///< +/-0.0013, PDG 2010 value
const double A0_PDG = -0.1184;                                  ///< +/-0.0010, PDG 2016 value
const double beta_W0 = (neutronBetaEp+m_e)/m_e;	///< beta spectrum endpoint, "natural" units
const double neutron_R0 = 0.0025896*1.2;		///< neutron and proton radius approximation, in "natural" units (1.2fm)/(hbar/m_e*c)
const double proton_M0 = m_p/m_e;				///< proton mass, "natural" units
const double neutron_M0 = m_n/m_e;				///< neutron mass, "natural" units
const double gamma_euler = 0.577215;			///< Euler's constant

extern Double_t calcTotalEn_Wilk(Double_t KE); 
extern Double_t calcMomentum_Wilk(Double_t W); 


class KurieFitter {
public:

  KurieFitter();
  ~KurieFitter();

  //void SetSpectrum(const TH1D* h);
  void FitSpectrum(const TH1D* spec, Double_t min, Double_t max, Double_t alpha=1.);
  void IterativeKurie(const TH1D* spec, Double_t min, Double_t max, Double_t alphaStart=0.9,
		      Double_t delta=1.e-3); // This will iterate the curie fit and calculate
                                             // the appropriate alpha at the delta level

  TGraphErrors returnKuriePlot() { return *_kuriePlot; };
  Double_t returnAlpha() { return _alpha; };
  Double_t returnW0() { return _W0; };
  Double_t returnW0err() { return _W0err; };
  Double_t returnK0() { return m_e*(_W0-1); };
  Double_t returnK0err() { return m_e*_W0err; };

  void setActualW0(Double_t w) { _actualW0 = w; }; // This is the W0 to compare against

private:

  TGraphErrors * _kuriePlot;
  Double_t _alpha;
  Double_t _W0, _W0err;
  Double_t _actualW0; // This defaults to neutronBetaEp

};

#endif
