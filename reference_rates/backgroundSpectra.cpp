/* 
    Produce the reference background spectra, as well as files which hold
    reference errors for low count bins

    Also makes background reference spectra for each individual PMT
    for analysis choice D (Type 0 only) for use in endpoint gain matching
*/


#include "MBUtils.hh"
#include "DataTree.hh"
#include "calibrationTools.hh"

#include <vector>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>

#include <TRandom3.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH1D.h>
#include <TChain.h>
#include <TMath.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TString.h>
#include <TStyle.h>

const bool useEvis = false;           // Uses Evis rather than Erecon

const bool useRCclasses = true;      // If this is true, we only use "good" response class 
                                     // events as defined by C. Swank (triangular MWPC responses)

std::map<Int_t,std::string> runType; //List of all runs in octet and their types
std::vector <Int_t> bgRuns_SFon; //list of bg runs with spin flipper on 
std::vector <Int_t> bgRuns_SFoff; //list of bg runs with spin flipper off
Double_t totalTimeON=0.;
Double_t totalTimeOFF=0.;
Double_t totalBLINDTimeON[2]={0.,0.};
Double_t totalBLINDTimeOFF[2]={0.,0.};

std::vector <Double_t> pmtTotalTime(8,0.); // Holds total time for each PMT, since some aren't always functioning

TString anaChoice[10] = {"A","B","C","D","E","F","G","H","J","K"};

std::vector < std::vector < std::vector <Double_t> > >  sfON(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.))); 
std::vector < std::vector < std::vector <Double_t> > > sfON_err(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.)));

std::vector < std::vector < std::vector <Double_t> > >  sfOFF(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.))); 
std::vector < std::vector < std::vector <Double_t> > > sfOFF_err(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.))); 

//Vector for individual PMT background
std::vector < std::vector < Double_t > > pmtBackground(120, std::vector<Double_t>(8,0.) ); 


std::vector <Int_t> badOct = {7,60,61,62,63,64,65,66}; 


std::vector <Int_t> getPMTQuality(Int_t runNumber) {
  //Read in PMT quality file
  std::cout << "Reading in PMT Quality file ...\n";
  std::vector <Int_t>  pmtQuality (8,0);
  Char_t temp[200];
  sprintf(temp,"%s/residuals/PMT_runQuality_master.dat",getenv("ANALYSIS_CODE")); 
  std::ifstream pmt;
  std::cout << temp << std::endl;
  pmt.open(temp);
  Int_t run_hold;
  while (pmt >> run_hold >> pmtQuality[0] >> pmtQuality[1] >> pmtQuality[2]
	 >> pmtQuality[3] >> pmtQuality[4] >> pmtQuality[5]
	 >> pmtQuality[6] >> pmtQuality[7]) {
    if (run_hold==runNumber) break;
    if (pmt.fail()) break;
  }
  pmt.close();
  if (run_hold!=runNumber) {
    std::cout << "Run not found in PMT quality file!" << std::endl;
    exit(0);
  }
  return pmtQuality;
}

void writeRatesToFile(int octMin, int octMax) {

  TString fn_base = TString::Format("backgroundRatesByAnaChoice/ReferenceSpectra_Octets-%i-%i_",octMin,octMax); 
  std::ofstream sf_ON;
  std::ofstream sf_OFF;

  std::ofstream pmt(TString::Format("backgroundRatesByAnaChoice/pmtRefSpectra_Octets_%i-%i.txt",octMin,octMax).Data());

  for (int anaCh=0; anaCh<10; anaCh++) {
    sf_ON.open(TString::Format("%ssfON-AnaCh-%s.txt",fn_base.Data(),anaChoice[anaCh].Data()).Data());
    sf_OFF.open(TString::Format("%ssfOFF-AnaCh-%s.txt",fn_base.Data(),anaChoice[anaCh].Data()).Data());
    
    sf_ON << std::setprecision(9);
    sf_OFF << std::setprecision(9);
    
    sf_ON << "totalTime\t" << totalTimeON << std::endl 
 	  << "totalTimeE\t" << totalBLINDTimeON[0] << std::endl
	  << "totalTimeW\t" << totalBLINDTimeON[1] << std::endl;

    sf_OFF << "totalTime\t" << totalTimeOFF << std::endl
	   << "totalTimeE\t" << totalBLINDTimeOFF[0] << std::endl
	   << "totalTimeW\t" << totalBLINDTimeOFF[1] << std::endl;

    pmt << "totalTimes:\t";
    for ( int p=0; p<8; ++p ) pmt << "\t\t" << pmtTotalTime[p];
    pmt << std::endl;

    for (int bin=0; bin<120; bin++) {

      Double_t binMidPoint = (double)bin*10.+5.;

      
      sf_ON << binMidPoint << "\t\t" << sfON[anaCh][0][bin] << "\t\t" << sfON_err[anaCh][0][bin] << "\t\t" 
	   << sfON[anaCh][1][bin] << "\t\t" << sfON_err[anaCh][1][bin] << "\n";

      sf_OFF << binMidPoint << "\t\t" << sfOFF[anaCh][0][bin] << "\t\t" << sfOFF_err[anaCh][0][bin] << "\t\t" 
	   << sfOFF[anaCh][1][bin] << "\t\t" << sfOFF_err[anaCh][1][bin] << "\n";
      
      pmt << binMidPoint << "\t\t";
      for ( auto p : pmtBackground[bin] ) pmt << p << "\t\t";
      pmt << std::endl;
    
    }
    
    sf_ON.close();
    sf_OFF.close();
    pmt.close();
      
  }
  
  std::cout << "Wrote All Reference rates to file!\n";

};


void readOctetFileForBGruns(int octet) {
  char fileName[200];
  sprintf(fileName,"%s/All_Octets/octet_list_%i.dat",getenv("OCTET_LIST"),octet);
  std::ifstream infile(fileName);
  std::string runTypeHold;
  int runNumberHold;
  
  // Populate the map with the runType and runNumber from this octet
  while (infile >> runTypeHold >> runNumberHold) {

    runType[runNumberHold] = runTypeHold;

    if ( runTypeHold=="A1" || runTypeHold=="A12" || runTypeHold=="B4" || runTypeHold=="B9" ) {
      bgRuns_SFoff.push_back(runNumberHold);
    }
    if ( runTypeHold=="B1" || runTypeHold=="B12" || runTypeHold=="A4" || runTypeHold=="A9" )  {
      bgRuns_SFon.push_back(runNumberHold);
    }

  }

  infile.close();
 
  std::cout << "Read in octet file for octet " << octet << "\n";
};

std::string getRunTypeFromOctetFile(int run) {
  for (auto const& rt : runType) {
    if (rt.first == run) return rt.second;
  }
  
  return "BAD";
  
};



void doBackgroundSpectra (int octetMin, int octetMax) 
{

  gStyle->SetOptStat(0);
  double fiducialCut = 50.; //mm
 
  std::cout << "Calculating BG spectra..." << std::endl;

  //Load all the bg run numbers
  
  for ( int i=octetMin ; i<=octetMax ; i++ ) {

    if ( std::find(badOct.begin(), badOct.end(),i) != badOct.end() ) continue;  //Checking if octet should be ignored for data quality reasons
    readOctetFileForBGruns(i);

  }
   
  char temp[200];

  //Root file to output to
  TFile *outfile = new TFile(TString::Format("BackgroundSpectra_%i-%i%s.root",octetMin,octetMax,(useEvis?"_Evis":"")),"RECREATE");

  //Make all the pertinent histograms

  TH1D* histOFF[3][2]; //We only go to type 2, then we separate them later and fill other histograms

  histOFF[0][0] = new TH1D("Type0E_sfOFF","Type0E",120,0.,1200.);
  histOFF[0][1] = new TH1D("Type0W_sfOFF","Type0W",120,0.,1200.);
  histOFF[1][0] = new TH1D("Type1E_sfOFF","Type1E",120,0.,1200.);
  histOFF[1][1] = new TH1D("Type1W_sfOFF","Type1W",120,0.,1200.);
  histOFF[2][0] = new TH1D("Type23E_sfOFF","Type23E",120,0.,1200.);
  histOFF[2][1] = new TH1D("Type23W_sfOFF","Type23W",120,0.,1200.);
  
  // Separated type 2
  TH1D* histOFF2[2];
  histOFF2[0] = new TH1D("Type2E_sfOFF","Type2E",120,0.,1200.);
  histOFF2[1] = new TH1D("Type2W_sfOFF","Type2W",120,0.,1200.);
  // Separated type 3
  TH1D* histOFF3[2];
  histOFF3[0] = new TH1D("Type3E_sfOFF","Type3E",120,0.,1200.);
  histOFF3[1] = new TH1D("Type3W_sfOFF","Type3W",120,0.,1200.);

  TH1D* histON[3][2]; //We only go to type 2, then we separate them later and fill other histograms

  histON[0][0] = new TH1D("Type0E_sfON","Type0E",120,0.,1200.);
  histON[0][1] = new TH1D("Type0W_sfON","Type0W",120,0.,1200.);
  histON[1][0] = new TH1D("Type1E_sfON","Type1E",120,0.,1200.);
  histON[1][1] = new TH1D("Type1W_sfON","Type1W",120,0.,1200.);
  histON[2][0] = new TH1D("Type23E_sfON","Type23E",120,0.,1200.);
  histON[2][1] = new TH1D("Type23W_sfON","Type23W",120,0.,1200.);
  
  // Separated type 2
  TH1D* histON2[2];
  histON2[0] = new TH1D("Type2E_sfON","Type2E",120,0.,1200.);
  histON2[1] = new TH1D("Type2W_sfON","Type2W",120,0.,1200.);
  // Separated type 3
  TH1D* histON3[2];
  histON3[0] = new TH1D("Type3E_sfON","Type3E",120,0.,1200.);
  histON3[1] = new TH1D("Type3W_sfON","Type3W",120,0.,1200.);

  // Individual PMT histograms for Type0 
  TH1D *pmt[8];
  for ( UInt_t i=0; i<8; ++i ) {
    pmt[i] = new TH1D(TString::Format("pmt%i",i), TString::Format("pmt %i",i),120, 0., 1200.);
  }


  BackscatterSeparator sep; // This is the object which takes care of 2/3 separation
  
  
  //Process SF off runs first
  for ( auto rn : bgRuns_SFoff ) {

    std::cout << "Processing run " << rn << "... ";

    

    sep.LoadCutCurve(rn); // Load the backscatter separator curve

    // DataTree structure
    DataTree t;

    // Input ntuple
    char tempIn[500];
    sprintf(tempIn, "%s/replay_pass3_%i.root", getenv("REPLAY_PASS3"),rn);
    //sprintf(tempIn, "%s/replay_pass3_%i.root", getenv("SPEC_REPLAY_FILES"),rn);
    
    t.setupInputTree(std::string(tempIn),"pass3");


    unsigned int nevents = t.getEntries();

    t.getEvent(nevents-1);
    totalTimeOFF += t.Time;
    totalBLINDTimeOFF[0] += t.TimeE;
    totalBLINDTimeOFF[1] += t.TimeW;

    std::vector <Int_t> pmtQuality = getPMTQuality(rn);

    for ( int p=0; p<8; ++p ) {
      if ( pmtQuality[p]) pmtTotalTime[p]+=t.Time;
    }

    double r2E = 0.; //position of event squared
    double r2W = 0.;
    
    for (unsigned int n=0 ; n<nevents ; n++ ) {
      
      
      t.getEvent(n);

      
      if ( t.PID==1 && t.Side<2 && t.Type<4 && t.Erecon>0. ) {

        if ( t.Type!=0 ) {
          if ( t.xeRC>6 || t.yeRC>6 || t.xwRC>6 || t.ywRC>6 ) continue;
	  else if ( t.xE.mult<1 || t.yE.mult<1 || t.xW.mult<1 || t.yW.mult<1 ) continue;
	}
        else {
          if ( t.Side==0 ) {
            if ( t.xeRC>6 || t.yeRC>6 ) continue; //only look at MWPC signal on East
	    else if ( t.xE.mult<1 || t.yE.mult<1 ) continue;
          }
          else if ( t.Side==1 ) {
            if ( t.xwRC>6 || t.ywRC>6 ) continue; //only look at MWPC signal on West
	    else if ( t.xW.mult<1 || t.yW.mult<1 ) continue;      
          }
	}
	
	
	//Type 2/3 separation... ADD IN AFTER DOING MWPC CAL
	if (t.Erecon>0. && t.Type==2) {
	  
	  if (t.Side==0) {
	    t.Type = sep.separate23(t.EMWPC_E);
	    t.Side = t.Type==2 ? 1 : 0;
	  }
	  else if (t.Side==1) {
	    t.Type = sep.separate23(t.EMWPC_W);
	    t.Side = t.Type==2 ? 0 : 1;
	  }
	  
	}
	
	
	//***********************************************************************************************************************
	// Filling rate histograms with "good" events to calculate the corrections
	
	
	r2E = t.xE.center*t.xE.center + t.yE.center*t.yE.center;
	r2W = t.xW.center*t.xW.center + t.yW.center*t.yW.center;

	if ( r2E<(fiducialCut*fiducialCut) && r2W<(fiducialCut*fiducialCut) )	  {
		
	  //Type 0
	  if (t.Type==0) {
	    histOFF[0][t.Side]->Fill( !useEvis ? t.Erecon : ( t.Side==0 ? t.EvisE : t.EvisW ) );
	    
	    if ( t.Side==0 ) {
	      if ( pmtQuality[0] ) pmt[0]->Fill(t.ScintE.e1);
	      if ( pmtQuality[1] ) pmt[1]->Fill(t.ScintE.e2);
	      if ( pmtQuality[2] ) pmt[2]->Fill(t.ScintE.e3);
	      if ( pmtQuality[3] ) pmt[3]->Fill(t.ScintE.e4);
	    }
	    else if (t.Side==1) {
	      if ( pmtQuality[4] ) pmt[4]->Fill(t.ScintW.e1);
	      if ( pmtQuality[5] ) pmt[5]->Fill(t.ScintW.e2);
	      if ( pmtQuality[6] ) pmt[6]->Fill(t.ScintW.e3);
	      if ( pmtQuality[7] ) pmt[7]->Fill(t.ScintW.e4);
	    }
	  }
	  
	  //Type 1
	  if (t.Type==1) histOFF[1][t.Side]->Fill(t.Erecon);
	
	  //Type 23... This puts them back where they would go if they weren't separated
	  if (t.Type==2 || t.Type==3) {
	    if (t.Side==0) { 
	      //histOFF[2][0]->Fill(t.Erecon);
	      if (t.Type==3) histOFF[2][0]->Fill(t.Erecon);
	      else histOFF[2][1]->Fill(t.Erecon);
	    }
	    else if (t.Side==1) {
	      //histOFF[2][1]->Fill(t.Erecon);
	      if (t.Type==3) histOFF[2][1]->Fill(t.Erecon);
	      else histOFF[2][0]->Fill(t.Erecon);
	    }
	  }
	
	  //Type 2
	  if (t.Type==2) histOFF2[t.Side]->Fill(t.Erecon);
	
	  //Type 3
	  if (t.Type==3) histOFF3[t.Side]->Fill(t.Erecon);

	}
      }
      
    }

    std::cout << "Finished Run " << rn << std::endl;
  }

  //Process SF ON runs now
  for ( auto rn : bgRuns_SFon ) {

    std::cout << "Processing run " << rn << "... ";

    sep.LoadCutCurve(rn); // Load the backscatter separator curve
    
    // DataTree structure
    DataTree t;

    // Input ntuple
    char tempIn[500];
    sprintf(tempIn, "%s/replay_pass3_%i.root", getenv("REPLAY_PASS3"),rn);
    
    t.setupInputTree(std::string(tempIn),"pass3");


    unsigned int nevents = t.getEntries();

    t.getEvent(nevents-1);
    totalTimeON += t.Time;
    totalBLINDTimeON[0] += t.TimeE;
    totalBLINDTimeON[1] += t.TimeW;

    std::vector <Int_t> pmtQuality = getPMTQuality(rn);

    for ( int p=0; p<8; ++p ) {
      if ( pmtQuality[p]) pmtTotalTime[p]+=t.Time;
    }

    double r2E = 0.; //position of event squared
    double r2W = 0.;
    
    for (unsigned int n=0 ; n<nevents ; n++ ) {
      
      
      t.getEvent(n);

      if ( t.PID==1 && t.Side<2 && t.Type<4 && t.Erecon>0. ) {

        if ( t.Type!=0 ) {
          if ( t.xeRC>6 || t.yeRC>6 || t.xwRC>6 || t.ywRC>6 ) continue;
	  else if ( t.xE.mult<1 || t.yE.mult<1 || t.xW.mult<1 || t.yW.mult<1 ) continue;
	}
        else {
          if ( t.Side==0 ) {
            if ( t.xeRC>6 || t.yeRC>6 ) continue; //only look at MWPC signal on East
	    else if ( t.xE.mult<1 || t.yE.mult<1 ) continue;
          }
          else if ( t.Side==1 ) {
            if ( t.xwRC>6 || t.ywRC>6 ) continue; //only look at MWPC signal on West
	    else if ( t.xW.mult<1 || t.yW.mult<1 ) continue;      
          }
	}
	
	
	//Type 2/3 separation... ADD IN AFTER DOING MWPC CAL
	if (t.Erecon>0. && t.Type==2) {
	  
	  if (t.Side==0) {
	    t.Type = sep.separate23(t.EMWPC_E);
	    t.Side = t.Type==2 ? 1 : 0;
	  }
	  else if (t.Side==1) {
	    t.Type = sep.separate23(t.EMWPC_W);
	    t.Side = t.Type==2 ? 0 : 1;
	  }
	  
	}
	
	
	//***********************************************************************************************************************
	// Filling rate histograms with "good" events to calculate the corrections
	
	
	r2E = t.xE.center*t.xE.center + t.yE.center*t.yE.center;
	r2W = t.xW.center*t.xW.center + t.yW.center*t.yW.center;

	if ( r2E<(fiducialCut*fiducialCut) && r2W<(fiducialCut*fiducialCut) )	  {
		
	  //Type 0
	  if (t.Type==0) {
	    histON[0][t.Side]->Fill( !useEvis ? t.Erecon : ( t.Side==0 ? t.EvisE : t.EvisW ) );

	    if ( t.Side==0 ) {
	      if ( pmtQuality[0] ) pmt[0]->Fill(t.ScintE.e1);
	      if ( pmtQuality[1] ) pmt[1]->Fill(t.ScintE.e2);
	      if ( pmtQuality[2] ) pmt[2]->Fill(t.ScintE.e3);
	      if ( pmtQuality[3] ) pmt[3]->Fill(t.ScintE.e4);
	    }
	    else if (t.Side==1) {
	      if ( pmtQuality[4] ) pmt[4]->Fill(t.ScintW.e1);
	      if ( pmtQuality[5] ) pmt[5]->Fill(t.ScintW.e2);
	      if ( pmtQuality[6] ) pmt[6]->Fill(t.ScintW.e3);
	      if ( pmtQuality[7] ) pmt[7]->Fill(t.ScintW.e4);
	    }
	  }
	
	  //Type 1
	  if (t.Type==1) histON[1][t.Side]->Fill(t.Erecon);
	
	  //Type 23... This puts them back where they would go if they weren't separated
	  if (t.Type==2 || t.Type==3) {
	    if (t.Side==0) { 
	      //histON[2][0]->Fill(t.Erecon);
	      if (t.Type==3) histON[2][0]->Fill(t.Erecon);
	      else histON[2][1]->Fill(t.Erecon);
	    }
	    else if (t.Side==1) {
	      //histON[2][1]->Fill(t.Erecon);
	      if (t.Type==3) histON[2][1]->Fill(t.Erecon);
	      else histON[2][0]->Fill(t.Erecon);
	    }
	  }
	
	  //Type 2
	  if (t.Type==2) histON2[t.Side]->Fill(t.Erecon);
	
	  //Type 3
	  if (t.Type==3) histON3[t.Side]->Fill(t.Erecon);
  
	}
      }
      
    }

    
    std::cout << "Finished Run " << rn << std::endl;
  }

  //Fill vectors with bin contents
   
  for (int anaCh = 1; anaCh<11 ; anaCh++) {
    
    int type_low=-1, type_high=-1;
    bool sep23 = false;
    bool sep2 = false;
    bool sep3 = false;
    
    if (anaCh==1) { type_low=0; type_high=2; }                             //All types, 2/3 not separated
    else if (anaCh==3 || anaCh==5) { type_low=0; type_high=1; sep23 = sep2 = sep3 = true;} // All event types, 2/3 separated
    else if (anaCh==2) { type_low=0; type_high=1;}                         // Type 0 and 1
    else if (anaCh==4) { type_low=0; type_high=0;}                         // Type 0
    else if (anaCh==6) { type_low=1; type_high=1;}                         // Type 1
    else if (anaCh==7) { type_low=type_high=2; }                           // Type 2/3 not separated
    else if (anaCh==8) { sep23 = sep2 = sep3 = true; }                 // Type 2/3, separated
    else if (anaCh==9) { sep23 = sep2 = true;}                   // Type 2
    else if (anaCh==10) { sep23 = sep3 = true;}                  // Type 3
    
    for (unsigned int side=0; side<2; side++) {
      for (unsigned int bin=1; bin<=120; bin++) {
	for (int type=type_low; type<=type_high; type++) {	 
	  
	  
	  sfON[anaCh-1][side][bin-1] += type>-1 ? histON[type][side]->GetBinContent(bin) : 0.;
	  sfON_err[anaCh-1][side][bin-1] += type>-1 ? power(histON[type][side]->GetBinError(bin),2) : 0.;
	  
	  sfOFF[anaCh-1][side][bin-1] += type>-1 ? histOFF[type][side]->GetBinContent(bin) : 0.;
	  sfOFF_err[anaCh-1][side][bin-1] += type>-1 ? power(histOFF[type][side]->GetBinError(bin),2) : 0.;

	}
	  
	if ( sep23 ) {
	  
	  sfON[anaCh-1][side][bin-1] += sep2 ? histON2[side]->GetBinContent(bin) : 0.;
	  sfON_err[anaCh-1][side][bin-1] += sep2 ? power(histON2[side]->GetBinError(bin),2) : 0.;
	  sfON[anaCh-1][side][bin-1] += sep3 ? histON3[side]->GetBinContent(bin) : 0.;
	  sfON_err[anaCh-1][side][bin-1] += sep3 ? power(histON3[side]->GetBinError(bin),2) : 0.;
	  
	  sfOFF[anaCh-1][side][bin-1] += sep2 ? histOFF2[side]->GetBinContent(bin) : 0.;
	  sfOFF_err[anaCh-1][side][bin-1] += sep2 ? power(histOFF2[side]->GetBinError(bin),2) : 0.;
	  sfOFF[anaCh-1][side][bin-1] += sep3 ? histOFF3[side]->GetBinContent(bin) : 0.;
	  sfOFF_err[anaCh-1][side][bin-1] += sep3 ? power(histOFF3[side]->GetBinError(bin),2) : 0.;
	  
	}	  
	  
	sfON_err[anaCh-1][side][bin-1] = sqrt(sfON_err[anaCh-1][side][bin-1]);
	sfOFF_err[anaCh-1][side][bin-1] = sqrt(sfOFF_err[anaCh-1][side][bin-1]);
	  
	//std::cout << anaChoice_A2[side][bin] << " " << anaChoice_A2_err[side][bin] << std::endl;
	
	
	// Now divide by the total time associated with that spin state's runs
	//sfON[anaCh-1][side][bin-1] /= totalTimeON;
	//sfON_err[anaCh-1][side][bin-1] /= totalTimeON;
	
	//sfOFF[anaCh-1][side][bin-1] /= totalTimeOFF;
	//sfOFF_err[anaCh-1][side][bin-1] /= totalTimeOFF;
      }
    }
  }

  // Filling individual PMT background vector
  for ( UInt_t b=0; b<120; ++b ) {
    for ( UInt_t p=0; p<8; ++p ) {
      pmtBackground[b][p] = pmt[p]->GetBinContent(b+1); //Note that bin 0 is underflow in TH1
    }
  }

  // Scale by 10^3mHz / (10 keV per bin) / total Time for drawing histograms
  // with nice units
  for (int s = 0; s<2; s++) {
    for (int t=0; t<3; t++) {

      histOFF[t][s]->SetYTitle("mHz/keV");
      histON[t][s]->SetYTitle("mHz/keV");
      
      histOFF[t][s]->Scale(1.e2/totalTimeOFF);
      histON[t][s]->Scale(1.e2/totalTimeON);
    }

    histON2[s]->SetYTitle("mHz/keV");
    histON3[s]->SetYTitle("mHz/keV");
    histOFF2[s]->SetYTitle("mHz/keV");
    histOFF3[s]->SetYTitle("mHz/keV");

    histON2[s]->Scale(1.e2/totalTimeON);
    histON3[s]->Scale(1.e2/totalTimeON);
    histOFF2[s]->Scale(1.e2/totalTimeOFF);
    histOFF3[s]->Scale(1.e2/totalTimeOFF);
  }
  
  for ( UInt_t p=0; p<8; ++p ) {
    
    pmt[p]->SetYTitle("mHz/keV");
    pmt[p]->SetXTitle("E_{vis} (keV)");
    pmt[p]->Scale( 1.e2/(pmtTotalTime[p]) );

  }
  
  
  

  outfile->Write();
  outfile->Close();
  std::cout << std::endl;


  /*for (int side = 0; side<2; side++) {
    for (int type=0; type<3; type++) {
      delete histOFF[type][side];
      delete histON[type][side];
    }
    delete histOFF2[side]; delete histOFF3[side]; delete histON2[side]; delete histON3[side];
    }*/
  

    
};
  

int main(int argc, char *argv[]) {

  int octetMin = atoi(argv[1]);
  int octetMax = atoi(argv[2]);

  doBackgroundSpectra(octetMin, octetMax);
  writeRatesToFile(octetMin, octetMax);

  //tests
  /*UInt_t XePeriod = getXeRunPeriod(atoi(argv[1]));
  vector < vector <Double_t> > triggerFunc = getTriggerFunctionParams(XePeriod,7);
  Double_t triggProbE = triggerProbability(triggerFunc[0],25.);
  Double_t triggProbW = triggerProbability(triggerFunc[1],25.);
  cout << triggProbE << " " << triggProbW << std::endl;*/


}
  
