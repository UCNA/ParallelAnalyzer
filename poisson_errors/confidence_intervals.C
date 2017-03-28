//This code integrates both the poisson and gaussian distributions to find the
// 68% confidence limits on each, and plots them against each other.

//Note that I am using the relationship between poisson ditributions and the
// gamma function to return values for the cdf of continuous counts.

double poiss_cdf(double x, double nu) {
  double a = x + 1.0;
  return ROOT::Math::gamma_cdf_c(nu,a,1.0);
}

double poiss_cdf_c(double x, double nu) {
  double a = x + 1.0;
  return ROOT::Math::gamma_cdf(nu,a,1.0);
}

double lowerBound(TString type, double n, double confLevel = 0.68) {
  double step = 0.1*n;
  double mean = n;
  double check = type==TString("poiss") ? poiss_cdf(n,mean) : ROOT::Math::normal_cdf(n,TMath::Sqrt(mean),mean) ;
  
  while ( check > (0.5 - confLevel/2.) ) {
    //cout << check << endl;
    n -= step;
    check = type==TString("poiss") ? poiss_cdf(n,mean) : ROOT::Math::normal_cdf(n,TMath::Sqrt(mean),mean) ;
  }

  return n;
}

double upperBound(TString type, double n, double confLevel = 0.68) {
  double step = 0.1*n;
  double mean = n;
  double check = type==TString("poiss") ? poiss_cdf_c(n,mean) : ROOT::Math::normal_cdf_c(n,TMath::Sqrt(mean),mean) ;
  cout << check << endl;
  while ( check > (0.5 - confLevel/2.) ) {
    //cout << check << endl;
    n += step;
    check = type==TString("poiss") ? poiss_cdf_c(n,mean) : ROOT::Math::normal_cdf_c(n,TMath::Sqrt(mean),mean) ;
  }

  return n;
}
  

void confidence_intervals()
{

  //make a vector of all x values to calculate confidence intervals on
  std::vector <double> x(26,0.);
  std::vector < std::vector <double> > poiss_conf(26, std::vector<double>(2,0.));
  std::vector < std::vector <double> > gauss_conf(26, std::vector<double>(2,0.));

  //begin by calculating the zeroth and first confidence interval at physical boundary
  x[0] = 0.;
  x[1] = 1.;
  poiss_conf[0][0] = poiss_conf[1][0] = 0.;
  poiss_conf[0][1] = poiss_conf[1][1] = 1.;
  //special case for first 0,1
  double i=0.;
  double ch = 1.;
  while ( ch > (1. - 0.68) ) {
    i += 0.001;
    ch = poiss_cdf_c(i,1.);
    poiss_conf[0][1] = poiss_conf[1][1] = i;
  }
    
  gauss_conf[0][0] = lowerBound("gauss",0.);
  gauss_conf[0][1] = upperBound("gauss",0.);
  gauss_conf[1][0] = lowerBound("gauss",1.);
  gauss_conf[1][1] = upperBound("gauss",1.);
  
  for (int n = 2; n<26; n++) {
    
    double dN = (double) n;
    x[n] = dN;
    poiss_conf[n][0] = lowerBound("poiss",dN);
    poiss_conf[n][1] = upperBound("poiss",dN);
    gauss_conf[n][0] = lowerBound("gauss",dN);
    gauss_conf[n][1] = upperBound("gauss",dN);
    
  }

  
  cout << "n\t\tgaussian\t\tpoisson\n";

  for (int n = 0; n<26; n++) {
    cout << n << "\t\t" << gauss_conf[n][0] << "\t" << gauss_conf[n][1] << "\t" << poiss_conf[n][0] << "\t" << poiss_conf[n][1] << "\n";
  }
  
}

  
