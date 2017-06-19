/* 
    Produce the reference background spectra, as well as files which hold
    reference errors for low count bins
*/


#include "MBUtils.hh"
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

const bool useEvis = true;           // Uses Evis rather than Erecon

using namespace std;

std::map<Int_t,std::string> runType; //List of all runs in octet and their types
std::vector <Int_t> fgRuns_SFon; //list of fg runs with spin flipper on 
std::vector <Int_t> fgRuns_SFoff; //list of fg runs with spin flipper off
Double_t totalTimeON=0.;
Double_t totalTimeOFF=0.;
Double_t totalBLINDTimeON[2]={0.,0.};
Double_t totalBLINDTimeOFF[2]={0.,0.};

TString anaChoice[10] = {"A","B","C","D","E","F","G","H","J","K"};

std::vector < std::vector < std::vector <Double_t> > >  sfON(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.))); 
std::vector < std::vector < std::vector <Double_t> > > sfON_err(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.)));

std::vector < std::vector < std::vector <Double_t> > >  sfOFF(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.))); 
std::vector < std::vector < std::vector <Double_t> > > sfOFF_err(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.))); 


std::vector <Int_t> badOct = {7,9,59,60,61,62,63,64,65,66,70,92}; 

int separate23(int side, double mwpcEn) {
  int type = 2;
  if (side==0) {
    type = ( mwpcEn>4.14 ) ? 3 : 2;
  }
  
  if (side==1) {
    type = ( mwpcEn>4.14 ) ? 3 : 2;
  }
  return type;
};

void writeRatesToFile(int octMin, int octMax) {
  
  TString fn_base = TString::Format("foregroundRatesByAnaChoice/ReferenceSpectra_Octets-%i-%i_",octMin,octMax);
  ofstream sf_ON;
  ofstream sf_OFF;
  
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

    
    for (int bin=0; bin<120; bin++) {

      Double_t binMidPoint = (double)bin*10.+5.;

      sf_ON << binMidPoint << "\t\t" << sfON[anaCh][0][bin] << "\t\t" << sfON_err[anaCh][0][bin] << "\t\t" 
	   << sfON[anaCh][1][bin] << "\t\t" << sfON_err[anaCh][1][bin] << "\n";

      sf_OFF << binMidPoint << "\t\t" << sfOFF[anaCh][0][bin] << "\t\t" << sfOFF_err[anaCh][0][bin] << "\t\t" 
	   << sfOFF[anaCh][1][bin] << "\t\t" << sfOFF_err[anaCh][1][bin] << "\n";
    
    }
    
    sf_ON.close();
    sf_OFF.close();
  }
  std::cout << "Wrote All Reference rates to file!\n";

};


void readOctetFileForFGruns(int octet) {
  char fileName[200];
  sprintf(fileName,"%s/All_Octets/octet_list_%i.dat",getenv("OCTET_LIST"),octet);
  std::ifstream infile(fileName);
  std::string runTypeHold;
  int runNumberHold;
  
  // Populate the map with the runType and runNumber from this octet
  while (infile >> runTypeHold >> runNumberHold) {

    runType[runNumberHold] = runTypeHold;

    if ( runTypeHold=="A2" || runTypeHold=="A10" || runTypeHold=="B5" || runTypeHold=="B7" ) {
      fgRuns_SFoff.push_back(runNumberHold);
    }
    if ( runTypeHold=="B2" || runTypeHold=="B10" || runTypeHold=="A5" || runTypeHold=="A7" )  {
      fgRuns_SFon.push_back(runNumberHold);
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



void doForegroundSpectra (int octetMin, int octetMax) 
{

  gStyle->SetOptStat(0);
  double fiducialCut = 50.; //mm
 
  cout << "Calculating FG spectra..." << endl;

  //Load all the fg run numbers
  
  for ( int i=octetMin ; i<=octetMax ; i++ ) {

    if ( std::find(badOct.begin(), badOct.end(),i) != badOct.end() ) continue;  //Checking if octet should be ignored for data quality reasons
    readOctetFileForFGruns(i);

  }
   
  char temp[200];

  //Root file to output to
  TFile *outfile = new TFile(TString::Format("SIM_ForegroundSpectra_%i-%i%s.root",octetMin,octetMax,(useEvis?"_Evis":"")),"RECREATE");

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

  
  BackscatterSeparator sep; // This is the object which takes care of 2/3 separation
  

  //Process SF off runs first
  for ( auto rn : fgRuns_SFoff ) {

    sep.LoadCutCurve(rn); // Load the backscatter separator curve
    
    std::string infile;
    TFile *input;
    TTree *Tin;
    
    sprintf(temp,"revCalSim_%i_Beta.root",rn);
    infile = getenv("REVCALSIM")+std::string("/beta/")+std::string(temp);
    input = new TFile(infile.c_str(), "READ");
    Tin = (TTree*)input->Get("revCalSim");

    double MWPCEnergyE=0., MWPCEnergyW=0.;
    double mwpcE[3] = {0.};
    double mwpcW[3] = {0.};
    double EmwpcX=0., EmwpcY=0., WmwpcX=0., WmwpcY=0., TimeE=0., TimeW=0., Time=0., Erecon=0.; //Branch Variables being read in
    int PID, type, side; // basic analysis tags
    int nClipped_EX, nClipped_EY, nClipped_WX, nClipped_WY;
    int xE_mult=0, yE_mult=0, xW_mult=0, yW_mult=0;

    double EdepQE=0., EdepQW=0.;
    
    Tin->SetBranchAddress("PID", &PID);
    Tin->SetBranchAddress("type", &type);
    Tin->SetBranchAddress("side", &side); 
    Tin->SetBranchAddress("Erecon",&Erecon);
    
    Tin->GetBranch("EdepQ")->GetLeaf("EdepQE")->SetAddress(&EdepQE);
    Tin->GetBranch("EdepQ")->GetLeaf("EdepQW")->SetAddress(&EdepQW);

    Tin->GetBranch("xE")->GetLeaf("center")->SetAddress(&EmwpcX);
    Tin->GetBranch("yE")->GetLeaf("center")->SetAddress(&EmwpcY);
    Tin->GetBranch("xW")->GetLeaf("center")->SetAddress(&WmwpcX);
    Tin->GetBranch("yW")->GetLeaf("center")->SetAddress(&WmwpcY);
    Tin->GetBranch("xE")->GetLeaf("nClipped")->SetAddress(&nClipped_EX);
    Tin->GetBranch("yE")->GetLeaf("nClipped")->SetAddress(&nClipped_EY);
    Tin->GetBranch("xW")->GetLeaf("nClipped")->SetAddress(&nClipped_WX);
    Tin->GetBranch("yW")->GetLeaf("nClipped")->SetAddress(&nClipped_WY);
    Tin->GetBranch("xE")->GetLeaf("mult")->SetAddress(&xE_mult);
    Tin->GetBranch("yE")->GetLeaf("mult")->SetAddress(&yE_mult);
    Tin->GetBranch("xW")->GetLeaf("mult")->SetAddress(&xW_mult);
    Tin->GetBranch("yW")->GetLeaf("mult")->SetAddress(&yW_mult);
    Tin->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyE")->SetAddress(&MWPCEnergyE);
    Tin->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyW")->SetAddress(&MWPCEnergyW);
    

    unsigned int nevents = Tin->GetEntriesFast();

    double r2E = 0.; //position of event squared
    double r2W = 0.;

    for (unsigned int n=0 ; n<nevents ; n++ ) {
      
      
      Tin->GetEvent(n);

      if ( PID==1 && side<2 && type<4 && Erecon>0.) {

	if ( type!=0 ) {
	  if ( xE_mult<1 || yE_mult<1 || xW_mult<1 || yW_mult<1 ) continue;
	}
	else {
	  if ( side==0 && ( xE_mult<1 || yE_mult<1 ) ) continue;
	  else if ( side==1 && ( xW_mult<1 || yW_mult<1 ) ) continue;
	}

	//Type 2/3 separation... ADD IN AFTER DOING MWPC CAL
	if (Erecon>0. && type==2) {
	  
	  if (side==0) {
	    type = sep.separate23(MWPCEnergyE);
	    side = type==2 ? 1 : 0;
	  }
	  else if (side==1) {
	    type = sep.separate23(MWPCEnergyW);
	    side = type==2 ? 0 : 1;
	  }
	  
	}
	
	
	//***********************************************************************************************************************
	// Filling rate histograms with "good" events to calculate the corrections
	
	r2E=EmwpcX*EmwpcX+EmwpcY*EmwpcY;
	r2W=WmwpcX*WmwpcX+WmwpcY*WmwpcY;

	if ( r2E<(fiducialCut*fiducialCut) && r2W<(fiducialCut*fiducialCut) )	  {
		
	  //Type 0
	  if (type==0) histOFF[0][side]->Fill( !useEvis ? Erecon : (side==0 ? EdepQE : EdepQW) );
	
	  //Type 1
	  if (type==1) histOFF[1][side]->Fill(Erecon);
	
	  //Type 23
	  if (type==2 || type==3) {
	    if (side==0) { 
	      //histOFF[2][0]->Fill(Erecon);
	      if (type==3) histOFF[2][0]->Fill(Erecon);
	      else histOFF[2][1]->Fill(Erecon);
	    }
	    else if (side==1) {
	      //histOFF[2][1]->Fill(Erecon);
	      if (type==3) histOFF[2][1]->Fill(Erecon);
	      else histOFF[2][0]->Fill(Erecon);
	    }
	  }
	
	  //Type 2
	  if (type==2) histOFF2[side]->Fill(Erecon);
	
	  //Type 3
	  if (type==3) histOFF3[side]->Fill(Erecon);

	}
      }
      
    }

    input->Close();
    if (input) delete input;
    cout << "Finished Run " << rn << endl;
  }

  //Process SF ON runs now
  for ( auto rn : fgRuns_SFon ) {

    sep.LoadCutCurve(rn); // Load the backscatter separator curve

    sprintf(temp,"revCalSim_%i_Beta.root",rn);
    std::string infile = getenv("REVCALSIM")+std::string("/beta/")+std::string(temp);
    TFile *input = new TFile(infile.c_str(), "READ");
    TTree *Tin = (TTree*)input->Get("revCalSim");


    double MWPCEnergyE=0., MWPCEnergyW=0.;
    double mwpcE[3] = {0.};
    double mwpcW[3] = {0.};
    double EmwpcX=0., EmwpcY=0., WmwpcX=0., WmwpcY=0., TimeE=0., TimeW=0., Time=0., Erecon=0.; //Branch Variables being read in
    int PID, type, side; // basic analysis tags
    int nClipped_EX, nClipped_EY, nClipped_WX, nClipped_WY;
    int xE_mult=0, yE_mult=0, xW_mult=0, yW_mult=0;

    double EdepQE=0., EdepQW=0.;
    
    Tin->SetBranchAddress("PID", &PID);
    Tin->SetBranchAddress("type", &type);
    Tin->SetBranchAddress("side", &side); 
    Tin->SetBranchAddress("Erecon",&Erecon);
    
    Tin->GetBranch("EdepQ")->GetLeaf("EdepQE")->SetAddress(&EdepQE);
    Tin->GetBranch("EdepQ")->GetLeaf("EdepQW")->SetAddress(&EdepQW);
    
    Tin->GetBranch("xE")->GetLeaf("center")->SetAddress(&EmwpcX);
    Tin->GetBranch("yE")->GetLeaf("center")->SetAddress(&EmwpcY);
    Tin->GetBranch("xW")->GetLeaf("center")->SetAddress(&WmwpcX);
    Tin->GetBranch("yW")->GetLeaf("center")->SetAddress(&WmwpcY);
    Tin->GetBranch("xE")->GetLeaf("nClipped")->SetAddress(&nClipped_EX);
    Tin->GetBranch("yE")->GetLeaf("nClipped")->SetAddress(&nClipped_EY);
    Tin->GetBranch("xW")->GetLeaf("nClipped")->SetAddress(&nClipped_WX);
    Tin->GetBranch("yW")->GetLeaf("nClipped")->SetAddress(&nClipped_WY);
    Tin->GetBranch("xE")->GetLeaf("mult")->SetAddress(&xE_mult);
    Tin->GetBranch("yE")->GetLeaf("mult")->SetAddress(&yE_mult);
    Tin->GetBranch("xW")->GetLeaf("mult")->SetAddress(&xW_mult);
    Tin->GetBranch("yW")->GetLeaf("mult")->SetAddress(&yW_mult);
    Tin->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyE")->SetAddress(&MWPCEnergyE);
    Tin->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyW")->SetAddress(&MWPCEnergyW);


    unsigned int nevents = Tin->GetEntriesFast();


    double r2E = 0.; //position of event squared
    double r2W = 0.;

    for (unsigned int n=0 ; n<nevents ; n++ ) {
      
      
      Tin->GetEvent(n);

      if ( PID==1 && side<2 && type<4 && Erecon>0.) {

	if ( type!=0 ) {
	  if ( xE_mult<1 || yE_mult<1 || xW_mult<1 || yW_mult<1 ) continue;
	}
	else {
	  if ( side==0 && ( xE_mult<1 || yE_mult<1 ) ) continue;
	  else if ( side==1 && ( xW_mult<1 || yW_mult<1 ) ) continue;
	}

	//Type 2/3 separation... ADD IN AFTER DOING MWPC CAL
	if (Erecon>0. && type==2) {
	  
	  if (side==0) {
	    type = sep.separate23(MWPCEnergyE);
	    side = type==2 ? 1 : 0;
	  }
	  else if (side==1) {
	    type = sep.separate23(MWPCEnergyW);
	    side = type==2 ? 0 : 1;
	  }
	  
	}
	
	
	//***********************************************************************************************************************
	// Filling rate histograms with "good" events to calculate the corrections
	
	r2E=EmwpcX*EmwpcX+EmwpcY*EmwpcY;
	r2W=WmwpcX*WmwpcX+WmwpcY*WmwpcY;

	if ( r2E<(fiducialCut*fiducialCut) && r2W<(fiducialCut*fiducialCut) )	  {
		
	  //Type 0
	  if (type==0) histON[0][side]->Fill( !useEvis ? Erecon : (side==0 ? EdepQE : EdepQW) );
	
	  //Type 1
	  if (type==1) histON[1][side]->Fill(Erecon);
	
	  //Type 23
	  if (type==2 || type==3) {
	    if (side==0) { 
	      //histON[2][0]->Fill(Erecon);
	      if (type==3) histON[2][0]->Fill(Erecon);
	      else histON[2][1]->Fill(Erecon);
	    }
	    else if (side==1) {
	      //histON[2][1]->Fill(Erecon);
	      if (type==3) histON[2][1]->Fill(Erecon);
	      else histON[2][0]->Fill(Erecon);
	    }
	  }
	
	  //Type 2
	  if (type==2) histON2[side]->Fill(Erecon);
	
	  //Type 3
	  if (type==3) histON3[side]->Fill(Erecon);

	}
      }
      
    }

    input->Close();
    if (input) delete input;
    cout << "Finished Run " << rn << endl;
  }

  //Fill vectors with bin contents
   
  for (int anaCh = 1; anaCh<11 ; anaCh++) {
    
    int type_low=-1, type_high=-1;
    bool sep23 = false;
    bool sep2 = false;
    bool sep3 = false;
    
    if (anaCh==1) { type_low=0; type_high=2; }                             //All types, 2/3 not separated
    else if (anaCh==3 || anaCh==5) { type_low=0; type_high=1; sep23=true;} // All event types, 2/3 separated
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

  // Scale by 10^3mHz / (10 keV per bin) / total Time for drawing histograms
  // with nice units
  /*for (int s = 0; s<2; s++) {
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
    }*/
  
  
  

  outfile->Write();
  outfile->Close();
  cout << endl;


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

  doForegroundSpectra(octetMin, octetMax);
  //writeRatesToFile(octetMin, octetMax);

  //tests
  /*UInt_t XePeriod = getXeRunPeriod(atoi(argv[1]));
  vector < vector <Double_t> > triggerFunc = getTriggerFunctionParams(XePeriod,7);
  Double_t triggProbE = triggerProbability(triggerFunc[0],25.);
  Double_t triggProbW = triggerProbability(triggerFunc[1],25.);
  cout << triggProbE << " " << triggProbW << endl;*/


}
  
