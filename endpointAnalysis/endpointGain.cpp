#include "BetaDecayTools.hh"
#include "DataTree.hh"
#include <TString.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TLeaf.h>

#include <iostream>
#include <fstream>


std::vector <Int_t> badOct {7,60,61,62,63,64,65,66};


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

std::vector <int>  readOctetFileForBGruns(int octet) {

  std::vector <int> runs;
  
  char fileName[200];
  sprintf(fileName,"%s/All_Octets/octet_list_%i.dat",getenv("OCTET_LIST"),octet);
  std::ifstream infile(fileName);
  std::string runTypeHold;
  int runNumberHold;
  int numRuns = 0;
  // Populate the map with the runType and runNumber from this octet
  while (infile >> runTypeHold >> runNumberHold) {

    if (runTypeHold=="A1" || runTypeHold=="A4" || runTypeHold=="A9" || runTypeHold=="A12" || 
	runTypeHold=="B1" || runTypeHold=="B4" || runTypeHold=="B9" || runTypeHold=="B12" )  {
     
      runs.push_back(runNumberHold);

    }
    numRuns++;
  }

  infile.close();
 
  std::cout << "Read in octet file for octet " << octet << "\n";
  return runs;

};

std::vector < std::vector <Double_t> > readPMTbackgroundRates(Int_t octet) {

  std::vector < std::vector <Double_t> > PMT;

  TString dir = TString::Format("%s/reference_rates/backgroundRatesByAnaChoice/",getenv("ANALYSIS_CODE"));

  TString filename = ( dir + TString("pmtRefSpectra_Octets_") +
		       ( octet<=59 ? TString("0-59") : TString("60-121") ) + TString(".txt") );

  std::ifstream infile(filename.Data());

  Double_t totalRefTime[8];

  std::string label;
  infile >> label >> totalRefTime[0] >> totalRefTime[1] >> totalRefTime[2] >> totalRefTime[3]
	 >> totalRefTime[4] >> totalRefTime[5] >> totalRefTime[6] >> totalRefTime[7];

  double binMid=0.;
  std::vector <Double_t> pmt(8, 0.);

  while ( infile >> binMid >> pmt[0] >> pmt[1] >> pmt[2] >> pmt[3] 
	  >> pmt[4] >> pmt[5] >> pmt[6] >> pmt[7] ) {

    for ( int i=0; i<8; ++i ) {
      pmt[i] = pmt[i]/totalRefTime[i];
      std::cout << pmt[i] << "  ";
    }
    std::cout << std::endl;

    PMT.push_back(pmt);
  }
  
  return PMT;

}



int main(int argc, char *argv[]) {

  
  if (argc!=2) {
    std::cout << "Usage: ./endpointGain.exe [octet]\n";
    //std::cout << "The code will produce comparisons for every octet in the range given,\non an octet-by-octet basis, and as a whole, using the Super-Sum\n";
    exit(0);
  }
  

  int octet = atoi(argv[1]);

  if ( std::find(badOct.begin(), badOct.end(),octet) != badOct.end() ) {

    std::cout << "Bad Octet... \n";

    std::ofstream gainFile(TString::Format("%s/EndpointGain/endpointGain_octet-%i.dat", getenv("ENDPOINT_ANALYSIS"),octet));

    for ( int i=0; i<8; ++i )  gainFile << 1. << std::endl;

    gainFile.close();

    return 0;
  }
  

  int nBins = 100;
  double minRange = 0.;
  double maxRange = 1000.;
  
  std::vector <int> runs = readOctetFile(octet);
  std::vector <int> bgruns = readOctetFileForBGruns(octet);

  std::vector < std::vector <Double_t> > pmtBackgroundRates = readPMTbackgroundRates(octet);

  

  ////////////////// Begin with data files /////////////////////

  // Vectors for creating the individual runs events and errors
  std::vector < std::vector < std::vector <Double_t> > > pmtSpec;
  std::vector < std::vector < std::vector <Double_t> > > pmtSpecErr;
    
  pmtSpec.resize(runs.size(),std::vector<std::vector<Double_t>>(8,std::vector<Double_t>(nBins,0.)));
  pmtSpecErr.resize(runs.size(),std::vector<std::vector<Double_t>>(8,std::vector<Double_t>(nBins,0.)));

  std::vector < std::vector < std::vector <Double_t> > > bgpmtSpec;
  std::vector < std::vector < std::vector <Double_t> > > bgpmtSpecErr;
    
  bgpmtSpec.resize(runs.size(),std::vector<std::vector<Double_t>>(8,std::vector<Double_t>(nBins,0.)));
  bgpmtSpecErr.resize(runs.size(),std::vector<std::vector<Double_t>>(8,std::vector<Double_t>(nBins,0.)));

  // Now loop over each run to determine the individual spectra, then fill their appropriate bins in the vector

  TH1D *spec[8]; // All 8 PMTs signals
  TH1D *bgspec[8]; // All 8 PMTs bg signals
  TH1D *simspec[8]; // All 8 PMTs signals

  int nRun = 0;
  
  std::vector <Double_t> totalTime(runs.size(),0.); // Holds the runLengths 
  std::vector <Double_t> bgtotalTime(runs.size(),0.); // Holds the runLengths 
  
  //runs.resize(0);
  for ( auto rn : runs ) {

    for ( int i=0; i<8; ++i ) {
      spec[i] = new TH1D(TString::Format("PMT%i",i),TString::Format("PMT %i",i),
		      nBins, minRange, maxRange);
    }
    
    // DataTree structure
    DataTree t;

    // Input ntuple
    char tempIn[500];
    sprintf(tempIn, "%s/replay_pass3_%i.root", getenv("REPLAY_PASS3"),rn);
    
    t.setupInputTree(std::string(tempIn),"pass3");

    unsigned int nevents = t.getEntries();

    t.getEvent(nevents-1);
    totalTime[nRun] = t.Time;

    double r2E = 0.; //position of event squared
    double r2W = 0.;
    
    for (unsigned int n=0 ; n<nevents ; n++ ) {

      t.getEvent(n);

      r2E = t.xE.center*t.xE.center + t.yE.center*t.yE.center;
      r2W = t.xW.center*t.xW.center + t.yW.center*t.yW.center;

      if ( t.PID==1 && t.Side<2 && t.Type==0 && t.Erecon>0. ) {

	if ( t.Side==0 ) {
	  if ( t.xeRC>6 || t.yeRC>6 ) continue; //only look at MWPC signal on East
	  else if ( t.xE.mult<1 || t.yE.mult<1 ) continue;
	}
	else if ( t.Side==1 ) {
	  if ( t.xwRC>6 || t.ywRC>6 ) continue; //only look at MWPC signal on West
	  else if ( t.xW.mult<1 || t.yW.mult<1 ) continue;      
	}

	if ( r2E > 30.*30. || r2W > 30.*30. ) continue;

	if ( t.Side==0 ) {
	  spec[0]->Fill(t.ScintE_bare.e1);
	  spec[1]->Fill(t.ScintE_bare.e2);
	  spec[2]->Fill(t.ScintE_bare.e3);
	  spec[3]->Fill(t.ScintE_bare.e4);
	}

	if ( t.Side==1 ) {
	  spec[4]->Fill(t.ScintW_bare.e1);
	  spec[5]->Fill(t.ScintW_bare.e2);
	  spec[6]->Fill(t.ScintW_bare.e3);
	  spec[7]->Fill(t.ScintW_bare.e4);
	}

      }
    }

    for ( int p=0; p<8; ++p ) {
      for ( int bin=1; bin<=nBins; ++bin ) {
	pmtSpec[nRun][p][bin-1] = (double)spec[p]->GetBinContent(bin)/totalTime[nRun];
	pmtSpecErr[nRun][p][bin-1] = spec[p]->GetBinError(bin)/totalTime[nRun];
      }
    }

    for ( int i=0; i<8; ++i ) { 
      // std::cout << "deleting spec[" << i << "]\n";
      delete spec[i];
    }
    nRun++;
    
  }

  nRun = 0;

  for ( auto rn : bgruns ) {

    for ( int i=0; i<8; ++i ) {
      bgspec[i] = new TH1D(TString::Format("bgPMT%i",i),TString::Format("bg PMT %i",i),
		      nBins, minRange, maxRange);
    }
    
    // DataTree structure
    DataTree t;

    // Input ntuple
    char tempIn[500];
    sprintf(tempIn, "%s/replay_pass3_%i.root", getenv("REPLAY_PASS3"),rn);
    
    t.setupInputTree(std::string(tempIn),"pass3");

    unsigned int nevents = t.getEntries();

    t.getEvent(nevents-1);
    bgtotalTime[nRun] = t.Time;

    double r2E = 0.; //position of event squared
    double r2W = 0.;
    
    for (unsigned int n=0 ; n<nevents ; n++ ) {

      t.getEvent(n);

      r2E = t.xE.center*t.xE.center + t.yE.center*t.yE.center;
      r2W = t.xW.center*t.xW.center + t.yW.center*t.yW.center;

      if ( t.PID==1 && t.Side<2 && t.Type==0 && t.Erecon>0. ) {

	if ( t.Side==0 ) {
	  if ( t.xeRC>6 || t.yeRC>6 ) continue; //only look at MWPC signal on East
	  else if ( t.xE.mult<1 || t.yE.mult<1 ) continue;
	}
	else if ( t.Side==1 ) {
	  if ( t.xwRC>6 || t.ywRC>6 ) continue; //only look at MWPC signal on West
	  else if ( t.xW.mult<1 || t.yW.mult<1 ) continue;      
	}

	if ( r2E > 30.*30. || r2W > 30.*30. ) continue;

	if ( t.Side==0 ) {
	  bgspec[0]->Fill(t.ScintE_bare.e1);
	  bgspec[1]->Fill(t.ScintE_bare.e2);
	  bgspec[2]->Fill(t.ScintE_bare.e3);
	  bgspec[3]->Fill(t.ScintE_bare.e4);
	}

	if ( t.Side==1 ) {
	  bgspec[4]->Fill(t.ScintW_bare.e1);
	  bgspec[5]->Fill(t.ScintW_bare.e2);
	  bgspec[6]->Fill(t.ScintW_bare.e3);
	  bgspec[7]->Fill(t.ScintW_bare.e4);
	}

      }
    }
    std::cout << "Made it through filling hists\n";
    for ( int p=0; p<8; ++p ) {
      for ( int bin=1; bin<=nBins; ++bin ) {
	bgpmtSpec[nRun][p][bin-1] = (double)bgspec[p]->GetBinContent(bin)/bgtotalTime[nRun];
	bgpmtSpecErr[nRun][p][bin-1] = ( bgspec[p]->GetBinContent(bin)>20 ?
					 bgspec[p]->GetBinError(bin)/bgtotalTime[nRun] :
					 TMath::Sqrt(pmtBackgroundRates[bin-1][p]/bgtotalTime[nRun]) );
      }
    }

    for ( int i=0; i<8; ++i ) { 
      // std::cout << "deleting spec[" << i << "]\n";
      delete bgspec[i];
    }
    nRun++;
    
  }

  

  ////////////// Now for simulation //////////////////
  // Vectors for creating the individual runs events and errors
  std::vector < std::vector < std::vector <Double_t> > > simSpec;
  std::vector < std::vector < std::vector <Double_t> > > simSpecErr;
    
  simSpec.resize(runs.size(),std::vector<std::vector<Double_t>>(8,std::vector<Double_t>(nBins,0.)));
  simSpecErr.resize(runs.size(),std::vector<std::vector<Double_t>>(8,std::vector<Double_t>(nBins,0.)));

  

  // Now loop over each run to determine the individual spectra, then fill their appropriate bins in the vector

  nRun = 0;
  
  for ( auto rn : runs ) {

    for ( int i=0; i<8; ++i ) {
      simspec[i] = new TH1D(TString::Format("SIM%i",i),TString::Format("SIM %i",i),
			 nBins, minRange, maxRange);
    }
    
    TFile *f = new TFile(TString::Format("%s/beta/revCalSim_%i_Beta.root",getenv("REVCALSIM"),rn), "READ");
    TTree *Tin = (TTree*)f->Get("revCalSim");

    std::cout << "Reading from " << TString::Format("%s/beta_highStatistics/revCalSim_%i_Beta.root",getenv("REVCALSIM"),rn).Data() << "\n";
    
    Double_t e0,e1,e2,e3,e4,e5,e6,e7;
    MWPC xE, yE, xW, yW;
    int PID, Side, Type;
    Double_t Erecon;
    
    
    Tin->SetBranchAddress("PID", &PID);
    Tin->SetBranchAddress("type", &Type);
    Tin->SetBranchAddress("side", &Side); 
    Tin->SetBranchAddress("Erecon",&Erecon);
    Tin->SetBranchAddress("xE",&xE);
    Tin->SetBranchAddress("yE",&yE);
    Tin->SetBranchAddress("xW",&xW);
    Tin->SetBranchAddress("yW",&yW);
    Tin->GetBranch("PMT")->GetLeaf("Evis0")->SetAddress(&e0);
    Tin->GetBranch("PMT")->GetLeaf("Evis1")->SetAddress(&e1);
    Tin->GetBranch("PMT")->GetLeaf("Evis2")->SetAddress(&e2);
    Tin->GetBranch("PMT")->GetLeaf("Evis3")->SetAddress(&e3);
    Tin->GetBranch("PMT")->GetLeaf("Evis4")->SetAddress(&e4);
    Tin->GetBranch("PMT")->GetLeaf("Evis5")->SetAddress(&e5);
    Tin->GetBranch("PMT")->GetLeaf("Evis6")->SetAddress(&e6);
    Tin->GetBranch("PMT")->GetLeaf("Evis7")->SetAddress(&e7);
    /*
      Tin->GetBranch("xE")->GetLeaf("center")->SetAddress(&EmwpcX);
      Tin->GetBranch("yE")->GetLeaf("center")->SetAddress(&EmwpcY);
      Tin->GetBranch("xW")->GetLeaf("center")->SetAddress(&WmwpcX);
      Tin->GetBranch("yW")->GetLeaf("center")->SetAddress(&WmwpcY);
      Tin->GetBranch("xE")->GetLeaf("nClipped")->SetAddress(&xE_nClipped);
      Tin->GetBranch("yE")->GetLeaf("nClipped")->SetAddress(&yE_nClipped);
      Tin->GetBranch("xW")->GetLeaf("nClipped")->SetAddress(&xW_nClipped);
      Tin->GetBranch("yW")->GetLeaf("nClipped")->SetAddress(&yW_nClipped);
      Tin->GetBranch("xE")->GetLeaf("mult")->SetAddress(&xE_mult);
      Tin->GetBranch("yE")->GetLeaf("mult")->SetAddress(&yE_mult);
      Tin->GetBranch("xW")->GetLeaf("mult")->SetAddress(&xW_mult);
      Tin->GetBranch("yW")->GetLeaf("mult")->SetAddress(&yW_mult);*/
    
    
    double r2E = 0.; //position of event squared
    double r2W = 0.;

    UInt_t nevents = Tin->GetEntriesFast();
    
    for (unsigned int n=0 ; n<nevents ; n++ ) {
      
      Tin->GetEvent(n);
      
      r2E = xE.center*xE.center + yE.center*yE.center;
      r2W = xW.center*xW.center + yW.center*yW.center;

      if ( PID==1 && Side<2 && Type==0 && Erecon>0. ) {

	if ( Side==0 && ( xE.mult<1 || yE.mult<1 ) ) continue;
	else if ( Side==1 && ( xW.mult<1 || yW.mult<1 ) ) continue;

	if ( r2E > 30.*30. || r2W > 30.*30. ) continue;

	
	if ( Side==0 ) {
	  simspec[0]->Fill(e0);
	  simspec[1]->Fill(e1);
	  simspec[2]->Fill(e2);
	  simspec[3]->Fill(e3);
	}

	if ( Side==1 ) {
	  simspec[4]->Fill(e4);
	  simspec[5]->Fill(e5);
	  simspec[6]->Fill(e6);
	  simspec[7]->Fill(e7);
	}

      }
    }

    for ( int p=0; p<8; ++p ) {
      for ( int bin=1; bin<=nBins; ++bin ) {
	simSpec[nRun][p][bin-1] = (double)simspec[p]->GetBinContent(bin)/totalTime[nRun];
	simSpecErr[nRun][p][bin-1] = simspec[p]->GetBinError(bin)/totalTime[nRun];
      }
    }

    for ( int i=0; i<8; ++i ) { 
      // std::cout << "deleting spec[" << i << "]\n";
      delete simspec[i];
    }
    nRun++;
    delete f;
    
  }

  //Now we take the weighted average over the runs in the octet
  
  TFile *fout = new TFile(TString::Format("%s/EndpointGain/endpointGain_octet-%i.root",
					  getenv("ENDPOINT_ANALYSIS"),octet),"RECREATE");
  //TFile *fout = new TFile(TString::Format("endpointGain_octet-%i.root",
  //					  octet),"RECREATE");
  
  std::cout << "Made output rootfile...\n\n";

  // Data //

  for ( int i=0; i<8; ++i ) {

    spec[i] = new TH1D(TString::Format("PMT%i",i),TString::Format("PMT %i",i),
		    nBins, minRange, maxRange);

  }
 
  for ( int p=0; p<8; ++p ) {
    for ( int bin=1; bin<=nBins; ++bin ) {
      
      double numer = 0.;
      double denom = 0.;
      
      for ( UInt_t i=0; i<runs.size(); ++i ) {
	
	//First background subtract each rate (keeping the error on the rate as just the counting
	// error
	if (i==0) std::cout << bin << ": " <<  pmtSpec[i][p][bin-1] << " - " << bgpmtSpec[i][p][bin-1] << " = ";
	pmtSpec[i][p][bin-1] -= bgpmtSpec[i][p][bin-1];//pmtBackgroundRates[bin-1][p];	
	pmtSpecErr[i][p][bin-1] = TMath::Sqrt( 
					      TMath::Power(bgpmtSpecErr[i][p][bin-1],2) + 
					      TMath::Power(pmtSpecErr[i][p][bin-1],2) );
	if (i==0) std::cout << pmtSpec[i][p][bin-1] << std::endl;

	double weight = pmtSpecErr[i][p][bin-1]>0. ? 1./(pmtSpecErr[i][p][bin-1]*pmtSpecErr[i][p][bin-1]) : 0.;
	numer += pmtSpec[i][p][bin-1]*weight;
	denom += weight;
      }

      spec[p]->SetBinContent(bin, denom>0. ? numer/denom : 0.);
      spec[p]->SetBinError(bin, denom>0. ? TMath::Sqrt(1./denom) : 0. );
    }
  }

  // Sim //
  for ( int i=0; i<8; ++i ) {

    simspec[i] = new TH1D(TString::Format("SIM%i",i),TString::Format("SIM %i",i),
		    nBins, minRange, maxRange);

  }
 
  for ( int p=0; p<8; ++p ) {
    for ( int bin=1; bin<=nBins; ++bin ) {

      double numer = 0.;
      double denom = 0.;
      
      for ( UInt_t i=0; i<runs.size(); ++i ) {
	double weight = simSpec[i][p][bin-1]>0. ? 1./(simSpecErr[i][p][bin-1]*simSpecErr[i][p][bin-1]) : 0.;
	numer += simSpec[i][p][bin-1]*weight;
	denom += weight;
      }

      simspec[p]->SetBinContent(bin, denom>0. ? numer/denom : 0.);
      simspec[p]->SetBinError(bin, denom>0. ? TMath::Sqrt(1./denom) : 0. );
    }
  }

  
  
  
  /////////////////////////// Kurie Fitting ////////////////////

  // First I'm going to Kurie fit the simulated endpoint, and then
  // I will iterate until the data until the endpoint matches this endpoint

  std::vector <Double_t> gain(8,0.);

  TGraphErrors kurie[8];
  TGraphErrors simkurie[8];

  KurieFitter kf, simkf;

  for ( int i=0; i<8; ++i ) {
    
    simkf.FitSpectrum(simspec[i],300.,500.,1.); //Fit the simulated spectrum to get relative endpoint

    kf.setActualW0( simkf.returnW0() ); // Setting the comparison endpoint to the extracted ep from sim
    kf.IterativeKurie(spec[i],300.,500.,1.,1.e-7);
   

    kurie[i] = kf.returnKuriePlot();
    kurie[i].SetName(TString::Format("data%i",i));
    kurie[i].Write();
    
    simkurie[i] = simkf.returnKuriePlot();
    simkurie[i].SetName(TString::Format("sim%i",i));
    simkurie[i].Write();

    gain[i] = kf.returnAlpha();

  }

  
  /*///// Getting the gains... 

  std::vector <Double_t> alpha_data(8,0.);
  std::vector <Double_t> alpha_sim(8,0.);
  std::vector <Double_t> gain(8,0.);

  TGraphErrors kurie[8];
  TGraphErrors simkurie[8];

  KurieFitter kf, simkf;

  for ( int i=0; i<8; ++i ) {
    
    kf.IterativeKurie(spec[i],300.,500.,1.1,1.e-6);
    simkf.IterativeKurie(simspec[i],300.,500.,1.1,1.e-6);

    kurie[i] = kf.returnKuriePlot();
    kurie[i].SetName(TString::Format("data%i",i));
    kurie[i].Write();
    
    simkurie[i] = simkf.returnKuriePlot();
    simkurie[i].SetName(TString::Format("sim%i",i));
    simkurie[i].Write();

    alpha_data[i] = kf.returnAlpha();
    alpha_sim[i] = simkf.returnAlpha();

    gain[i] = alpha_sim[i]>0. ? alpha_data[i]/alpha_sim[i] : 1.;

    }*/

  /// Write out pmt endpoint gain corrections
  std::ofstream gainFile(TString::Format("%s/EndpointGain/endpointGain_octet-%i.dat",
  					 getenv("ENDPOINT_ANALYSIS"),octet));
  // std::ofstream gainFile(TString::Format("endpointGain_octet-%i.dat",octet));

  for ( auto g : gain ) gainFile << g << std::endl;

  gainFile.close();
  

  fout->Write();
  delete fout;

  return 0;
  
}

  
