#include "BetaDecayTools.hh"
#include "DataTree.hh"
#include <TString.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TMath.h>

#include <iostream>
#include <fstream>
#include <algorithm>

std::vector <Int_t> badOct {7,9,59,60,61,62,63,64,65,66,67,91,93,101,107,121};

std::vector <int>  readOctetFile(int octet); 

int main(int argc, char *argv[]) {

  if (argc!=4) {
    std::cout << "Usage: ./endpointTracker.exe [octet start] [octet end] [bool simulation]\n";
    //std::cout << "The code will produce comparisons for every octet in the range given,\non an octet-by-octet basis, and as a whole, using the Super-Sum\n";
    exit(0);
  }
  
  //Show that we can do the Kurie fit...

  int octetMin = atoi(argv[1]);
  int octetMax = atoi(argv[2]);

  bool sim = ( std::string(argv[3])=="true" ) ? true : false;
  std::string prefix = sim ? "SIM" : "UK";

  int nBins = 100;
  double minRange = 0.;
  double maxRange = 1000.;

  ///// Loop over all octets in range

  

  for ( int octet=octetMin ; octet<=octetMax ; ++octet ) {

    if ( std::find(badOct.begin(), badOct.end(),octet) != badOct.end() ) continue;  //Checking if octet should be ignored for data quality reasons
  
    std::vector <int> runs = readOctetFile(octet);

  
    // Vectors for creating the individual runs events and errors
    std::vector < std::vector <Double_t> > eastSpectrum;
    std::vector < std::vector <Double_t> > eastSpectrumErr;
    
    eastSpectrum.resize(runs.size(),std::vector<Double_t>(nBins,0.));
    eastSpectrumErr.resize(runs.size(),std::vector<Double_t>(nBins,0.));

    std::vector < std::vector <Double_t> > westSpectrum;
    std::vector < std::vector <Double_t> > westSpectrumErr;
    
    westSpectrum.resize(runs.size(),std::vector<Double_t>(nBins,0.));
    westSpectrumErr.resize(runs.size(),std::vector<Double_t>(nBins,0.));

    // Now loop over each run to determine the individual spectra, then fill their appropriate bins in the vector

    Double_t en, fg, bg, time, bg_time;
    std::string txt_hold;
    Int_t bin = 0;

    int nRun = 0;
  
    for ( auto rn : runs ) {

      std::ifstream rateFile;

      //Read in East rate files and fill vectors for the correct number of bins

      rateFile.open(TString::Format("../Asymmetry/BinByBinComparison/%s_Erun%i_anaChD.dat", prefix.c_str(),rn));
      rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
      rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
    
      bin = 0; 
    
      while ( rateFile >> en >> fg >> bg && bin<nBins )  { 
	eastSpectrum[nRun][bin] = (fg - bg);
	eastSpectrumErr[nRun][bin] = sqrt( eastSpectrum[nRun][bin]/time );
	++bin;
      }
    
      rateFile.close();
	

      //Do the same for west rate files
      rateFile.open(TString::Format("../Asymmetry/BinByBinComparison/%s_Wrun%i_anaChD.dat", prefix.c_str(),rn));
      rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
      rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
    
      bin = 0; 
    
      while ( rateFile >> en >> fg >> bg && bin<nBins)  { 
	westSpectrum[nRun][bin] = (fg - bg);
	westSpectrumErr[nRun][bin] = sqrt( westSpectrum[nRun][bin]/time );
	++bin;
      }
    
      rateFile.close();
  

      nRun++;
    
    }

    std::cout << "\nMoving to Kurie Fits...\n\n";
    TFile *fout = new TFile(TString::Format("%s/FinalEndpoints/%s_octet_%i_Erecon.root",getenv("ENDPOINT_ANALYSIS"),prefix.c_str(),octet),"RECREATE");  
    std::cout << "Made output rootfile...\n\n";

    TH1D *specE = new TH1D("specE","East Spectrum", nBins, minRange, maxRange);
    TH1D *specW = new TH1D("specW","West Spectrum", nBins, minRange, maxRange);  

    // loop over each bin for each run and make weighted average

    for ( int bin=1; bin<=nBins; ++bin ) {

      double numerE = 0.;
      double denomE = 0.;
      double numerW = 0.;
      double denomW = 0.;
    
      for ( UInt_t i=0; i<runs.size(); ++i ) {
	double weight = eastSpectrum[i][bin-1]>0. ? 1./(eastSpectrumErr[i][bin-1]*eastSpectrumErr[i][bin-1]) : 1.;
	numerE += eastSpectrum[i][bin-1]*weight;
	denomE += weight;
      
	weight = westSpectrum[i][bin-1]>0. ? 1./(westSpectrumErr[i][bin-1]*westSpectrumErr[i][bin-1]) : 1.;
	numerW += westSpectrum[i][bin-1]*weight;
	denomW += weight;
      }
    
      specE->SetBinContent(bin, numerE>0. ? numerE/denomE : 0.);
      specE->SetBinError(bin, numerE>0. ? TMath::Sqrt(1./denomE) : 0.01);
      specW->SetBinContent(bin, numerW>0. ? numerW/denomW : 0.);
      specW->SetBinError(bin, numerW>0. ? TMath::Sqrt(1./denomW) : 0.01);
    
    }

    TGraphErrors kurieE, kurieW; //These will hold the kurie plots

    //Now do some Kurie Fitting and write endpoints to file

    // open file for writing out all endpoints for each detector in a txt file
    std::ofstream epfile(TString::Format("%s/FinalEndpoints/%s_finalEndpoints_Octet_%i.txt",
					 getenv("ENDPOINT_ANALYSIS"),prefix.c_str(),octet).Data());

    KurieFitter kf;

    kf.FitSpectrum(specE, 300., 600., 1.);
    epfile << "East_Endpoint: " << kf.returnK0() << " +/- " << kf.returnK0err() << "\n";

    kurieE = kf.returnKuriePlot();
    kurieE.SetName("kurieE");
    kurieE.Write();

    kf.FitSpectrum(specW, 300., 600., 1.);
    epfile << "West_Endpoint: " << kf.returnK0() << " +/- " << kf.returnK0err();
    epfile.close();

    kurieW = kf.returnKuriePlot();
    kurieW.SetName("kurieW");
    kurieW.Write();

  

    fout->Write();
    delete fout;
  }



  return 0;
  
}

  
std::vector <int>  readOctetFile(int octet) {

  std::vector <int> runs;
  
  char fileName[200];
  sprintf(fileName,"%s/All_Octets/octet_list_%i.dat",getenv("OCTET_LIST"),octet);
  std::ifstream infile(fileName);
  std::string runTypeHold;
  int runNumberHold;
  int numRuns = 0;
  // Populate the map with the runType and runNumber from this octet
  while (infile >> runTypeHold >> runNumberHold) {

    if (runTypeHold=="A2" || runTypeHold=="A5" || runTypeHold=="A7" || runTypeHold=="A10" || 
	runTypeHold=="B2" || runTypeHold=="B5" || runTypeHold=="B7" || runTypeHold=="B10" )  {
     
      runs.push_back(runNumberHold);

    }
    numRuns++;
  }

  infile.close();
 
  std::cout << "Read in octet file for octet " << octet << "\n";
  return runs;

};
