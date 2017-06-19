 /* Code to take a run number, retrieve it's runperiod, and construct the 
weighted spectra which would be seen as a reconstructed energy on one side
of the detector. Also applies the trigger functions */

#include "DeltaExpProcessor.h"
#include "positionMapHandler.hh"
#include "sourcePeaks.h"
#include "runInfo.h"
#include "calibrationTools.hh"
#include "TriggerMap.hh"
#include "BetaSpectrum.hh"
#include "MBUtils.hh"
#include "MWPCPositionResponse.hh"


#include <vector>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>

#include <TRandom3.h>
#include <TTree.h>
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


using namespace std;

std::map<Int_t,std::string> runType; //List of all runs in octet and their types
std::vector <Int_t> betaRuns; //list of beta decay runs


// MonteCarlo Super Ratios
std::vector < std::vector <Double_t> > oct_A_SR(10,std::vector<Double_t>(120,0.)); 
std::vector < std::vector <Double_t> > quartA_A_SR(10,std::vector<Double_t>(120,0.));
std::vector < std::vector <Double_t> > quartB_A_SR(10,std::vector<Double_t>(120,0.));
std::vector < std::vector <Double_t> > pairA1_A_SR(10,std::vector<Double_t>(120,0.));
std::vector < std::vector <Double_t> > pairA2_A_SR(10,std::vector<Double_t>(120,0.));
std::vector < std::vector <Double_t> > pairB1_A_SR(10,std::vector<Double_t>(120,0.));
std::vector < std::vector <Double_t> > pairB2_A_SR(10,std::vector<Double_t>(120,0.));
std::vector < std::vector <Double_t> > oct_A_SR_err(10,std::vector<Double_t>(120,0.)); 
std::vector < std::vector <Double_t> > quartA_A_SR_err(10,std::vector<Double_t>(120,0.));
std::vector < std::vector <Double_t> > quartB_A_SR_err(10,std::vector<Double_t>(120,0.));
std::vector < std::vector <Double_t> > pairA1_A_SR_err(10,std::vector<Double_t>(120,0.));
std::vector < std::vector <Double_t> > pairA2_A_SR_err(10,std::vector<Double_t>(120,0.));
std::vector < std::vector <Double_t> > pairB1_A_SR_err(10,std::vector<Double_t>(120,0.));
std::vector < std::vector <Double_t> > pairB2_A_SR_err(10,std::vector<Double_t>(120,0.));


/*void calculate_A_SR() {

  for (unsigned int anaCh = 1; anaCh<11; anaCh++) {
    for (unsigned int bin=0; bin<120; bin++) {
      
      if (oct_sfOFF[anaCh-1][0][bin]>0. && oct_sfOFF[anaCh-1][1][bin]>0. && oct_sfON[anaCh-1][0][bin]>0. && oct_sfON[anaCh-1][1][bin]>0.) { 
	Double_t R = oct_sfOFF[anaCh-1][1][bin]*oct_sfON[anaCh-1][0][bin]/(oct_sfON[anaCh-1][1][bin]*oct_sfOFF[anaCh-1][0][bin]);
	Double_t deltaR = sqrt(R*R*(power(oct_sfOFF_err[0][bin]/oct_sfOFF[anaCh-1][0][bin],2)+power(oct_sfON_err[1][bin]/oct_sfON[anaCh-1][1][bin],2)+
				    power(oct_sfOFF_err[1][bin]/oct_sfOFF[anaCh-1][1][bin],2)+power(oct_sfON_err[0][bin]/oct_sfON[anaCh-1][0][bin],2)));
	oct_A_SR[anaCh-1][bin] = (1.-sqrt(R))/(1+sqrt(R));
	oct_A_SR_err[anaCh-1][bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
      }
      else {
	oct_A_SR[anaCh-1][bin] = 0.;
	oct_A_SR_err[anaCh-1][bin] = 0.;
      }
      
      if (quartA_sfOFF[anaCh-1][0][bin]>0. && quartA_sfOFF[anaCh-1][1][bin]>0. && quartA_sfON[anaCh-1][0][bin]>0. && quartA_sfON[anaCh-1][1][bin]>0.) { 
	Double_t R = quartA_sfOFF[anaCh-1][1][bin]*quartA_sfON[anaCh-1][0][bin]/(quartA_sfON[anaCh-1][1][bin]*quartA_sfOFF[anaCh-1][0][bin]);
	Double_t deltaR = sqrt(R*R*(power(quartA_sfOFF_err[0][bin]/quartA_sfOFF[anaCh-1][0][bin],2)+power(quartA_sfON_err[1][bin]/quartA_sfON[anaCh-1][1][bin],2)+
				    power(quartA_sfOFF_err[1][bin]/quartA_sfOFF[anaCh-1][1][bin],2)+power(quartA_sfON_err[0][bin]/quartA_sfON[anaCh-1][0][bin],2)));
	quartA_A_SR[anaCh-1][bin] = (1.-sqrt(R))/(1+sqrt(R));
	quartA_A_SR_err[anaCh-1][bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
      }
      else {
	quartA_A_SR[anaCh-1][bin] = 0.;
	quartA_A_SR_err[anaCh-1][bin] = 0.;
      }

      if (quartB_sfOFF[anaCh-1][0][bin]>0. && quartB_sfOFF[anaCh-1][1][bin]>0. && quartB_sfON[anaCh-1][0][bin]>0. && quartB_sfON[anaCh-1][1][bin]>0.) { 
	Double_t R = quartB_sfOFF[anaCh-1][1][bin]*quartB_sfON[anaCh-1][0][bin]/(quartB_sfON[anaCh-1][1][bin]*quartB_sfOFF[anaCh-1][0][bin]);
	Double_t deltaR = sqrt(R*R*(power(quartB_sfOFF_err[0][bin]/quartB_sfOFF[anaCh-1][0][bin],2)+power(quartB_sfON_err[1][bin]/quartB_sfON[anaCh-1][1][bin],2)+
				    power(quartB_sfOFF_err[1][bin]/quartB_sfOFF[anaCh-1][1][bin],2)+power(quartB_sfON_err[0][bin]/quartB_sfON[anaCh-1][0][bin],2)));
	quartB_A_SR[anaCh-1][bin] = (1.-sqrt(R))/(1+sqrt(R));
	quartB_A_SR_err[anaCh-1][bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
      }
      else {
	quartB_A_SR[anaCh-1][bin] = 0.;
	quartB_A_SR_err[anaCh-1][bin] = 0.;
      }

      if (pairA1_sfOFF[anaCh-1][0][bin]>0. && pairA1_sfOFF[anaCh-1][1][bin]>0. && pairA1_sfON[anaCh-1][0][bin]>0. && pairA1_sfON[anaCh-1][1][bin]>0.) { 
	Double_t R = pairA1_sfOFF[anaCh-1][1][bin]*pairA1_sfON[anaCh-1][0][bin]/(pairA1_sfON[anaCh-1][1][bin]*pairA1_sfOFF[anaCh-1][0][bin]);
	Double_t deltaR = sqrt(R*R*(power(pairA1_sfOFF_err[0][bin]/pairA1_sfOFF[anaCh-1][0][bin],2)+power(pairA1_sfON_err[1][bin]/pairA1_sfON[anaCh-1][1][bin],2)+
				    power(pairA1_sfOFF_err[1][bin]/pairA1_sfOFF[anaCh-1][1][bin],2)+power(pairA1_sfON_err[0][bin]/pairA1_sfON[anaCh-1][0][bin],2)));
	pairA1_A_SR[anaCh-1][bin] = (1.-sqrt(R))/(1+sqrt(R));
	pairA1_A_SR_err[anaCh-1][bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
      }
      else {
	pairA1_A_SR[anaCh-1][bin] = 0.;
	pairA1_A_SR_err[anaCh-1][bin] = 0.;
      }

      if (pairA2_sfOFF[anaCh-1][0][bin]>0. && pairA2_sfOFF[anaCh-1][1][bin]>0. && pairA2_sfON[anaCh-1][0][bin]>0. && pairA2_sfON[anaCh-1][1][bin]>0.) { 
	Double_t R = pairA2_sfOFF[anaCh-1][1][bin]*pairA2_sfON[anaCh-1][0][bin]/(pairA2_sfON[anaCh-1][1][bin]*pairA2_sfOFF[anaCh-1][0][bin]);
	Double_t deltaR = sqrt(R*R*(power(pairA2_sfOFF_err[0][bin]/pairA2_sfOFF[anaCh-1][0][bin],2)+power(pairA2_sfON_err[1][bin]/pairA2_sfON[anaCh-1][1][bin],2)+
				    power(pairA2_sfOFF_err[1][bin]/pairA2_sfOFF[anaCh-1][1][bin],2)+power(pairA2_sfON_err[0][bin]/pairA2_sfON[anaCh-1][0][bin],2)));
	pairA2_A_SR[anaCh-1][bin] = (1.-sqrt(R))/(1+sqrt(R));
	pairA2_A_SR_err[anaCh-1][bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
      }
      else {
	pairA2_A_SR[anaCh-1][bin] = 0.;
	pairA2_A_SR_err[anaCh-1][bin] = 0.;
      }

      if (pairB1_sfOFF[anaCh-1][0][bin]>0. && pairB1_sfOFF[anaCh-1][1][bin]>0. && pairB1_sfON[anaCh-1][0][bin]>0. && pairB1_sfON[anaCh-1][1][bin]>0.) { 
	Double_t R = pairB1_sfOFF[anaCh-1][1][bin]*pairB1_sfON[anaCh-1][0][bin]/(pairB1_sfON[anaCh-1][1][bin]*pairB1_sfOFF[anaCh-1][0][bin]);
	Double_t deltaR = sqrt(R*R*(power(pairB1_sfOFF_err[0][bin]/pairB1_sfOFF[anaCh-1][0][bin],2)+power(pairB1_sfON_err[1][bin]/pairB1_sfON[anaCh-1][1][bin],2)+
				    power(pairB1_sfOFF_err[1][bin]/pairB1_sfOFF[anaCh-1][1][bin],2)+power(pairB1_sfON_err[0][bin]/pairB1_sfON[anaCh-1][0][bin],2)));
	pairB1_A_SR[anaCh-1][bin] = (1.-sqrt(R))/(1+sqrt(R));
	pairB1_A_SR_err[anaCh-1][bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
      }
      else {
	pairB1_A_SR[anaCh-1][bin] = 0.;
	pairB1_A_SR_err[anaCh-1][bin] = 0.;
      }

      if (pairB2_sfOFF[anaCh-1][0][bin]>0. && pairB2_sfOFF[anaCh-1][1][bin]>0. && pairB2_sfON[anaCh-1][0][bin]>0. && pairB2_sfON[anaCh-1][1][bin]>0.) { 
	Double_t R = pairB2_sfOFF[anaCh-1][1][bin]*pairB2_sfON[anaCh-1][0][bin]/(pairB2_sfON[anaCh-1][1][bin]*pairB2_sfOFF[anaCh-1][0][bin]);
	Double_t deltaR = sqrt(R*R*(power(pairB2_sfOFF_err[0][bin]/pairB2_sfOFF[anaCh-1][0][bin],2)+power(pairB2_sfON_err[1][bin]/pairB2_sfON[anaCh-1][1][bin],2)+
				    power(pairB2_sfOFF_err[1][bin]/pairB2_sfOFF[anaCh-1][1][bin],2)+power(pairB2_sfON_err[0][bin]/pairB2_sfON[anaCh-1][0][bin],2)));
	pairB2_A_SR[anaCh-1][bin] = (1.-sqrt(R))/(1+sqrt(R));
	pairB2_A_SR_err[anaCh-1][bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
      }
      else {
	pairB2_A_SR[anaCh-1][bin] = 0.;
	pairB2_A_SR_err[anaCh-1][bin] = 0.;
      }
    }
  }
  };*/

void writeCorrectionFactorsToFile(Int_t octet) {
  TString fn_base = TString::Format("DeltaExp_OctetByOctetCorrections/ThOverProc_Octet-%i_Analysis-",octet); 
  std::ofstream oct; 
  std::ofstream quartA;
  std::ofstream quartB; 
  std::ofstream pairA1; 
  std::ofstream pairA2; 
  std::ofstream pairB1; 
  std::ofstream pairB2; 

  for (int i=1; i<11; i++) {
    oct.open(TString::Format("%s%i.txt",fn_base.Data(),i));
    //quartA.open(TString::Format("%s%i_quartA.txt",fn_base.Data(),i));
    //quartB.open(TString::Format("%s%i_quartB.txt",fn_base.Data(),i));
    //pairA1.open(TString::Format("%s%i_pairA1.txt",fn_base.Data(),i));
    //pairA2.open(TString::Format("%s%i_pairA2.txt",fn_base.Data(),i));
    //pairB1.open(TString::Format("%s%i_pairB1.txt",fn_base.Data(),i));
    //pairB2.open(TString::Format("%s%i_pairB2.txt",fn_base.Data(),i));

    for (int bin=0; bin<120; bin++) {
      Double_t binMidPoint = (double)bin*10.+5.;
      Double_t A_theory = A0_PDG*asymmetryCorrectionFactor(binMidPoint)*beta(binMidPoint)/2.;
      oct << binMidPoint << "\t" << ( fabs(oct_A_SR[i-1][bin] )>0.000001 ? A_theory/oct_A_SR[i-1][bin] : 1.) << "\t" 
	  << ( fabs(oct_A_SR[i-1][bin] )>0.000001 ? fabs(A_theory/power(oct_A_SR[i-1][bin],2)*oct_A_SR_err[i-1][bin]) : 1.) << "\n";
      /*quartA << binMidPoint << "\t" << A_theory/quartA_A_SR[i-1][bin] << "\t" << A_theory/power(quartA_A_SR[i-1][bin],2)*quartA_A_SR_err[i-1][bin] << "\n";
      quartB << binMidPoint << "\t" << A_theory/quartB_A_SR[i-1][bin] << "\t" << A_theory/power(quartB_A_SR[i-1][bin],2)*quartB_A_SR_err[i-1][bin] << "\n";
      pairA1 << binMidPoint << "\t" << A_theory/pairA1_A_SR[i-1][bin] << "\t" << A_theory/power(pairA1_A_SR[i-1][bin],2)*pairA1_A_SR_err[i-1][bin] << "\n";
      pairA2 << binMidPoint << "\t" << A_theory/pairA2_A_SR[i-1][bin] << "\t" << A_theory/power(pairA2_A_SR[i-1][bin],2)*pairA2_A_SR_err[i-1][bin] << "\n";
      pairB1 << binMidPoint << "\t" << A_theory/pairB1_A_SR[i-1][bin] << "\t" << A_theory/power(pairB1_A_SR[i-1][bin],2)*pairB1_A_SR_err[i-1][bin] << "\n";
      pairB2 << binMidPoint << "\t" << A_theory/pairB2_A_SR[i-1][bin] << "\t" << A_theory/power(pairB2_A_SR[i-1][bin],2)*pairB2_A_SR_err[i-1][bin] << "\n";
      */
    }
    oct.close();
    quartA.close();
    quartB.close();
    pairA1.close();
    pairA2.close();
    pairB1.close();
    pairB2.close();
    
  }
  std::cout << "Wrote All Corrections to file!\n";

};


void readOctetFile(int octet) {
  char fileName[200];
  sprintf(fileName,"%s/All_Octets/octet_list_%i.dat",getenv("OCTET_LIST"),octet);
  std::ifstream infile(fileName);
  std::string runTypeHold;
  int runNumberHold;
  int numRuns = 0;
  // Populate the map with the runType and runNumber from this octet
  while (infile >> runTypeHold >> runNumberHold) {

    runType[runNumberHold] = runTypeHold;

    if (runTypeHold=="A2" || runTypeHold=="A5" || runTypeHold=="A7" || runTypeHold=="A10" || 
	runTypeHold=="B2" || runTypeHold=="B5" || runTypeHold=="B7" || runTypeHold=="B10" )  {
     
      betaRuns.push_back(runNumberHold);
      std::cout << betaRuns[numRuns] << std::endl;

    }
    numRuns++;
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


int getPolarization(int run) {

  std::ifstream infile(TString::Format("%s/masterBetaRunList.txt",getenv("OCTET_LIST")).Data());

  int rn = 0;
  std::string type = "";
  bool runFound = false;

  while ( infile >> rn >> type ) {
    if ( rn == run ) { runFound = true; infile.close(); continue; }
  }

  if ( !runFound ) throw "You didn't pick a run in an octet...";
  //Flipper OFF is p = -1
  if ( type=="A1" || type=="A2" || type=="A10" || type=="A12" || type=="B4" || type=="B5" || type=="B7" || type=="B9" ) return -1;
  else if ( type=="B1" || type=="B2" || type=="B10" || type=="B12" || type=="A4" || type=="A5" || type=="A7" || type=="A9" ) return 1;
  
  else {
    std::cout << "You chose a Depol Run.. no polarization\n";
    return 0;
  }
};


vector <vector <double> > returnSourcePosition (Int_t runNumber, string src) {
  Char_t temp[500];
  sprintf(temp,"%s/source_list_%i.dat",getenv("SOURCE_LIST"),runNumber);
  ifstream file(temp);
  cout << src << endl;
  int num = 0;
  file >> num;
  cout << num << endl;
  int srcNum = 0;
  string src_check;
  for (int i=0; i<num;srcNum++,i++) {
    file >> src_check;
    cout << src_check << endl;
    if (src_check==src) break;   
  }
  cout << "The source Number is: " << srcNum << endl;
  if (srcNum==num) {
    cout << "Didn't find source in that run\n"; exit(0);
  }
  file.close();
  
  sprintf(temp,"%s/source_positions_%i.dat",getenv("SOURCE_POSITIONS"),runNumber);
  file.open(temp);
  
  vector < vector < double > > srcPos;
  srcPos.resize(2,vector <double> (3,0.));
  
  for (int i=0; i<srcNum+1; i++) {
    for (int j=0; j<2; j++) {
      for (int jj=0; jj<3; jj++) {
	file >> srcPos[j][jj];
      }
    }
  }
  return srcPos;
}

vector <Int_t> getEreconPMTQuality(Int_t runNumber) {
  //Read in PMT quality file
  cout << "Reading in PMT Quality file ...\n";
  vector <Int_t>  pmtQuality (8,0);
  Char_t temp[200];
  sprintf(temp,"%s/residuals/PMT_runQuality_master.dat",getenv("ANALYSIS_CODE")); 
  ifstream pmt;
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
    cout << "Run not found in PMT quality file!" << endl;
    exit(0);
  }
  return pmtQuality;
}

vector < Double_t > GetAlphaValues(Int_t runPeriod)
{
  Char_t temp[500];
  vector < Double_t > alphas (8,0.);
  sprintf(temp,"%s/simulation_comparison/nPE_per_keV/nPE_per_keV_%i.dat",getenv("ANALYSIS_CODE"),runPeriod);
  ifstream infile;
  infile.open(temp);
  Int_t i = 0;

  while (infile >> alphas[i]) { std::cout << alphas[i] << std::endl; i++; }
  return alphas;
}


std::vector < std::vector <Double_t> > loadPMTpedestals(Int_t runNumber) {

  Char_t temp[500];
  std::vector < std::vector < Double_t > > peds (8,std::vector <Double_t> (2,0.));
  if (runNumber>20000) sprintf(temp,"%s/PMT_pedestals_%i.dat",getenv("PEDESTALS"),runNumber);
  else sprintf(temp,"%s/pedestal_widths_%i.dat",getenv("PEDESTALS"),runNumber);
  ifstream infile;
  infile.open(temp);

  Int_t i = 0;
  Int_t run;

  while (infile >> run >> peds[i][0] >> peds[i][1]) { std::cout << "Pedestal " << i << ": " << peds[i][0] << " " << peds[i][1] << std::endl; i++; }
  return peds;

};

std::vector <Double_t> loadGainFactors(Int_t runNumber) {

  // Read gain corrections file                                                                                                                
  char tempFileGain[500];
  sprintf(tempFileGain, "%s/gain_bismuth_%i.dat",getenv("GAIN_BISMUTH"), runNumber);
  std::cout << "... Reading: " << tempFileGain << std::endl;

  std::vector <Double_t> gainCorrection(8,1.);
  std::vector <Double_t> fitMean(8,1.);
  ifstream fileGain(tempFileGain);
  for (int i=0; i<8; i++) {
    fileGain >> fitMean[i] >> gainCorrection[i];
  }
  std::cout << "...   PMT E1: " << gainCorrection[0] << std::endl;
  std::cout << "...   PMT E2: " << gainCorrection[1] << std::endl;
  std::cout << "...   PMT E3: " << gainCorrection[2] << std::endl;
  std::cout << "...   PMT E4: " << gainCorrection[3] << std::endl;
  std::cout << "...   PMT W1: " << gainCorrection[4] << std::endl;
  std::cout << "...   PMT W2: " << gainCorrection[5] << std::endl;
  std::cout << "...   PMT W3: " << gainCorrection[6] << std::endl;
  std::cout << "...   PMT W4: " << gainCorrection[7] << std::endl;

  return gainCorrection;
};


void SetUpTree(TTree *tree) {
  tree->Branch("PID", &PID, "PID/I");
  tree->Branch("side", &side, "side/I");
  tree->Branch("type", &type, "type/I");
  tree->Branch("Erecon", &Erecon,"Erecon/D");
  
  tree->Branch("Evis",&evis,"EvisE/D:EvisW");
  tree->Branch("Edep",&edep,"EdepE/D:EdepW");
  tree->Branch("EdepQ",&edepQ,"EdepQE/D:EdepQW");
  tree->Branch("primKE",&primKE,"primKE/D");
  tree->Branch("primTheta",&primTheta,"primTheta/D");
  tree->Branch("AsymWeight",&AsymWeight,"AsymWeight/D");
  
  tree->Branch("time",&Time,"timeE/D:timeW");
  tree->Branch("MWPCEnergy",&mwpcE,"MWPCEnergyE/D:MWPCEnergyW");
  tree->Branch("MWPCPos",&mwpc_pos,"MWPCPosE[3]/D:MWPCPosW[3]");
  tree->Branch("nClipped_EX",&nClipped_EX,"nClipped_EX/I");
  tree->Branch("nClipped_EY",&nClipped_EY,"nClipped_EY/I");
  tree->Branch("nClipped_WX",&nClipped_WX,"nClipped_WX/I");
  tree->Branch("nClipped_WY",&nClipped_WY,"nClipped_WY/I");
  tree->Branch("ScintPos",&scint_pos,"ScintPosE[3]/D:ScintPosW[3]");
  tree->Branch("ScintPosAdjusted",&scint_pos_adj,"ScintPosAdjE[3]/D:ScintPosAdjW[3]");
  tree->Branch("PMT",&pmt,"Evis0/D:Evis1:Evis2:Evis3:Evis4:Evis5:Evis6:Evis7:etaEvis0/D:etaEvis1:etaEvis2:etaEvis3:etaEvis4:etaEvis5:etaEvis6:etaEvis7:nPE0/D:nPE1:nPE2:nPE3:nPE4:nPE5:nPE6:nPE7");
  
  tree->Branch("xE",&xE,"center/D:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/D:height");
  tree->Branch("yE",&yE,"center/D:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/D:height");
  tree->Branch("xW",&xW,"center/D:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/D:height");
  tree->Branch("yW",&yW,"center/D:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/D:height");
  

}
  

void calcDeltaExp (int octet) 
{

 
  cout << "Calculating Delta Exp for octet " << octet << endl;

  readOctetFile(octet);

  std::vector < std::vector < std::vector <Double_t> > >  A2(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.))); 
  std::vector < std::vector < std::vector <Double_t> > > A2_err(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.))); 
  std::vector < std::vector < std::vector <Double_t> > > A5(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.)));  
  std::vector < std::vector < std::vector <Double_t> > > A5_err(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.))); 
  std::vector < std::vector < std::vector <Double_t> > > A7(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.)));  
  std::vector < std::vector < std::vector <Double_t> > > A7_err(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.))); 
  std::vector < std::vector < std::vector <Double_t> > > A10(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.)));  
  std::vector < std::vector < std::vector <Double_t> > > A10_err(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.))); 
  std::vector < std::vector < std::vector <Double_t> > > B2(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.)));  
  std::vector < std::vector < std::vector <Double_t> > > B2_err(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.))); 
  std::vector < std::vector < std::vector <Double_t> > > B5(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.)));  
  std::vector < std::vector < std::vector <Double_t> > > B5_err(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.))); 
  std::vector < std::vector < std::vector <Double_t> > > B7(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.)));  
  std::vector < std::vector < std::vector <Double_t> > > B7_err(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.))); 
  std::vector < std::vector < std::vector <Double_t> > > B10(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.)));  
  std::vector < std::vector < std::vector <Double_t> > > B10_err(10,std::vector < std::vector <Double_t> > (2, std::vector<Double_t>(120,0.))); 


  TRandom3 *event = new TRandom3(0);
  UInt_t evtON = (int)(event->Rndm()*97000000.);
  UInt_t evtOFF = (int)(event->Rndm()*97000000.);
  UInt_t evt = 0;
  delete event;

  for (auto& runNumber : betaRuns) {
  
    UInt_t BetaEvents = 24000000; //Number of Beta events for each beta run... This is pure events, not triggered events
    

    std::cout << "Processing " << BetaEvents << " events for run " << runNumber << "...\n";



    ///////////////////////// SETTING GAIN OF FIRST and second DYNYODE
    Double_t g_d1 = 4.;
    Double_t g_d2 = 4.;
    Double_t g_d = 16.;
    Double_t g_rest = 12500.;

    /////// Loading other run dependent quantities
    vector <Int_t> pmtQuality = getEreconPMTQuality(runNumber); // Get the quality of the PMTs for that run
    UInt_t calibrationPeriod = getSrcRunPeriod(runNumber); // retrieve the calibration period for this run
    UInt_t XePeriod = getXeRunPeriod(runNumber); // Get the proper Xe run period for the Trigger functions
    //GetPositionMap(XePeriod);
    PositionMap posmap(5.0, 50.); //Load position map with 5 mm bins
    posmap.readPositionMap(XePeriod,"endpoint");
    vector <Double_t> alpha = GetAlphaValues(calibrationPeriod); // fill vector with the alpha (nPE/keV) values for this run period


    /////////////// Load trigger functions //////////////////
    TriggerFunctions trigger(runNumber);

    std::vector < std::vector <Double_t> > pedestals = loadPMTpedestals(runNumber);
    std::vector < Double_t > PMTgain = loadGainFactors(runNumber);

    //Load the simulated relationship between EQ and Etrue
    EreconParameterization eRecon(runNumber);

    Int_t pol = getPolarization(runNumber);

    LinearityCurve linCurve(calibrationPeriod,false);

    //Decide which simulation to use...
    std::string simLocation;
    TChain *chain = new TChain("anaTree");
  
    if (runNumber<20000) simLocation = string(getenv("SIM_2011_2012"));
    else if (runNumber>=21087 && runNumber<21679) simLocation = string(getenv("SIM_2012_2013_ISOBUTANE"));
    else simLocation = string(getenv("SIM_2012_2013"));


    std::cout << "Using simulation from " << simLocation << "...\n";

    
    //Read in simulated data and put in a TChain
    
    int numFiles = 2000 ; 
    Char_t temp[200];
    
    for (int i=0; i<numFiles; i++) {
      
      if (pol==-1) sprintf(temp,"%s/Beta_polE/analyzed_%i.root",simLocation.c_str(),i);
      if (pol==1) sprintf(temp,"%s/Beta_polW/analyzed_%i.root",simLocation.c_str(),i);
      
      chain->AddFile(temp);
      
    }
    
    
    // Set the addresses of the information read in from the simulation file
    chain->SetBranchAddress("MWPCEnergy",&mwpcE);
    chain->SetBranchAddress("time",&Time);
    chain->SetBranchAddress("Edep",&edep);
    chain->SetBranchAddress("EdepQ",&edepQ);
    chain->SetBranchAddress("MWPCPos",&mwpc_pos);
    chain->SetBranchAddress("ScintPos",&scint_pos);
    chain->SetBranchAddress("primKE",&primKE);
    chain->SetBranchAddress("primTheta",&primTheta);
    chain->SetBranchAddress("Cath_EX",Cath_EX);
    chain->SetBranchAddress("Cath_EY",Cath_EY);
    chain->SetBranchAddress("Cath_WX",Cath_WX);
    chain->SetBranchAddress("Cath_WY",Cath_WY);
    
    
    //These are for feeding in Xuan's simulations... this needs to be updated so that I can pass a flag and change these on the fly
    //chain->SetBranchAddress("PrimaryParticleSpecies",&primaryID);
    //chain->SetBranchAddress("PrimaryParticleSpecies",&primaryID);
    //chain->SetBranchAddress("mwpcEnergy",&mwpcE);
    //chain->SetBranchAddress("scintTimeToHit",&Time);
    //chain->SetBranchAddress("scintillatorEdep",&edep);
    //chain->SetBranchAddress("scintillatorEdepQuenched",&edepQ);
    //chain->SetBranchAddress("MWPCPos",&mwpc_pos);
    //chain->SetBranchAddress("ScintPos",&scint_pos);
    //chain->SetBranchAddress("primaryKE",&Eprim);
    
    
   
    //Set random number generator
    TRandom3 *seed = new TRandom3(4*runNumber); // seed generator, always same for given run for reproducibility
    TRandom3 *rand0 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );
    TRandom3 *rand1 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );
    TRandom3 *rand2 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );
    TRandom3 *rand3 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );
    
    //Initialize random numbers for 8 pmt trigger probabilities
    TRandom3 *randPMTE1 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );
    TRandom3 *randPMTE2 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );
    TRandom3 *randPMTE3 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );
    TRandom3 *randPMTE4 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );
    TRandom3 *randPMTW1 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );
    TRandom3 *randPMTW2 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );
    TRandom3 *randPMTW3 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );
    TRandom3 *randPMTW4 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );

    std::vector <Double_t> triggRandVec(4,0.);

    // Wirechamber information
    bool EastScintTrigger, WestScintTrigger, EMWPCTrigger, WMWPCTrigger; //Trigger booleans
    Double_t MWPCAnodeThreshold=0.2; // keV dep in the wirechamber.. 

    //Get total number of events in TChain
    UInt_t nevents = chain->GetEntries();
    cout << "events = " << nevents << endl;

    std::string runType = getRunTypeFromOctetFile(runNumber);
    
    if ( runType=="BAD" ) { std::cout << "BAD RUNTYPE FOUND FOR RUN\n"; exit(0); }
    
    /*Int_t evtStart = 0;
    
    if ( runType=="A2" || runType=="A5" ) evtStart = 0;
    if ( runType=="A7" || runType=="A10" ) evtStart = 25000000;
    if ( runType=="B2" || runType=="B5" ) evtStart = 50000000;
    if ( runType=="B7" || runType=="B10" ) evtStart = 75000000;
    */
  
    
    UInt_t evtTally = 0; //To keep track of the number of events 

    // To start where the last sim of it's run type left off
    if ( runType=="A2" || runType=="A10" || runType=="B5" || runType=="B7" ) evt = evtOFF;
    else evt = evtON;

    vector < vector <Int_t> > gridPoint;

    BackscatterSeparator sep;
    sep.LoadCutCurve(runNumber);

    //make evt type histograms

    TFile *specFile = new TFile(TString::Format("spectra/spectra_%i.root",runNumber),"RECREATE");
    
    TH1D* hist[3][2]; //We only go to type 2, then we separate them later and fill other histograms

    hist[0][0] = new TH1D("Type0E","Type0E",120,0.,1200.);
    hist[0][1] = new TH1D("Type0W","Type0W",120,0.,1200.);
    hist[1][0] = new TH1D("Type1E","Type1E",120,0.,1200.);
    hist[1][1] = new TH1D("Type1W","Type1W",120,0.,1200.);
    hist[2][0] = new TH1D("Type23E","Type23E",120,0.,1200.);
    hist[2][1] = new TH1D("Type23W","Type23W",120,0.,1200.);
   
    // Separated type 2
    TH1D* hist2[2];
    hist2[0] = new TH1D("Type2E","Type2E",120,0.,1200.);
    hist2[1] = new TH1D("Type2W","Type2W",120,0.,1200.);
    // Separated type 3
    TH1D* hist3[2];
    hist3[0] = new TH1D("Type3E","Type3E",120,0.,1200.);
    hist3[1] = new TH1D("Type3W","Type3W",120,0.,1200.);

    //Create simulation output file
    //Char_t outputfile[500];
    //sprintf(outputfile,"revCalSim_%i_%s.root",runNumber,"Beta");
    //sprintf(outputfile,"revCalSim_%i_%s.root",runNumber,source.c_str());
    //TFile *outfile = new TFile(outputfile, "RECREATE");
    
    //Setup the output tree
    //    TTree *tree = new TTree("revCalSim", "revCalSim");
    //SetUpTree(tree); //Setup the output tree and branches
    

    //Read in events and determine evt type based on triggers
    while (evtTally<=BetaEvents) {

      if (evt>=nevents) evt=0; //Wrapping the events back to the beginning
      EastScintTrigger = WestScintTrigger = EMWPCTrigger = WMWPCTrigger = false; //Resetting triggers each event

      chain->GetEvent(evt);

      if (evtTally%100000==0) {std::cout << " Event Number " << evtTally << " in run " << runNumber << "  (" << evt << ")" << std::endl;}//cout << "filled event " << evt << endl;


      // Go ahead and increment the event counters
      evt++;
      evtTally++;

      //Checking that the event occurs within the fiducial volume in the simulation to minimize
      // contamination from edge effects and interactions with detector walls
      Double_t fidCut = 50.;
      
    
      scint_pos_adj.ScintPosAdjE[0] = -scint_pos.ScintPosE[0]*sqrt(0.6)*10.;
      scint_pos_adj.ScintPosAdjE[1] = scint_pos.ScintPosE[1]*sqrt(0.6)*10.;
      scint_pos_adj.ScintPosAdjW[0] = scint_pos.ScintPosW[0]*sqrt(0.6)*10.;
      scint_pos_adj.ScintPosAdjW[1] = scint_pos.ScintPosW[1]*sqrt(0.6)*10.;
      scint_pos_adj.ScintPosAdjE[2] = scint_pos.ScintPosE[2]*10.;
      scint_pos_adj.ScintPosAdjW[2] = scint_pos.ScintPosW[2]*10.;

      std::vector <double> posex(3,0.);
      std::vector <double> poswx(3,0.);
      std::vector <double> posey(3,0.);
      std::vector <double> poswy(3,0.);

      double dCath_EX[16]{0.};
      double dCath_EY[16]{0.};
      double dCath_WX[16]{0.};
      double dCath_WY[16]{0.};

      // Different def of wires in data
      for ( int i=0; i<16; ++i ) {
	dCath_EX[15-i] = (double)Cath_EX[i];
	dCath_EY[15-i] = (double)Cath_EY[i];
	dCath_WX[i] = (double)Cath_WX[i];
	dCath_WY[15-i] = (double)Cath_WY[i];
      }
      
      MWPCCathodeHandler cathResp(dCath_EX,dCath_EY,dCath_WX,dCath_WY);
      cathResp.loadCathodeModelParams(runNumber);
      
      cathResp.findAllPositions(true,false);

      posex = cathResp.getPosEX();
      posey = cathResp.getPosEY();
      poswx = cathResp.getPosWX();
      poswy = cathResp.getPosWY();
      

      //std::cout << mwpcAdjE[0] << "\t" << mwpcAdjE[1] << mwpcAdjW[0] << "\t" << mwpcAdjW[1] <<"\n"; 
      
      nClipped_EX = nClipped_EY = nClipped_WX = nClipped_WY = 0;
     
      nClipped_EX = cathResp.getnClippedEX();
      nClipped_EY = cathResp.getnClippedEY();
      nClipped_WX = cathResp.getnClippedWX();
      nClipped_WY = cathResp.getnClippedWY();

      xE.center = posex[0] * sqrt(0.6);
      yE.center = posey[0] * sqrt(0.6);
      xW.center = poswx[0] * sqrt(0.6);
      yW.center = poswy[0] * sqrt(0.6);

      xE.width = posex[1] * sqrt(0.6);
      yE.width = posey[1] * sqrt(0.6);
      xW.width = poswx[1] * sqrt(0.6);
      yW.width = poswy[1] * sqrt(0.6);
      
      xE.height = posex[2];
      yE.height = posey[2];
      xW.height = poswx[2];
      yW.height = poswy[2];

      xE.mult = cathResp.getMultEX();
      yE.mult = cathResp.getMultEY();
      xW.mult = cathResp.getMultWX();
      yW.mult = cathResp.getMultWY();

      xE.nClipped = cathResp.getnClippedEX();
      yE.nClipped = cathResp.getnClippedEY();
      xW.nClipped = cathResp.getnClippedWX();
      yW.nClipped = cathResp.getnClippedWY();

      xE.maxWire = cathResp.getMaxWireEX();
      yE.maxWire = cathResp.getMaxWireEY();
      xW.maxWire = cathResp.getMaxWireWX();
      yW.maxWire = cathResp.getMaxWireWY();

      xE.maxValue = dCath_EX[xE.maxWire];
      yE.maxValue = dCath_EY[yE.maxWire];
      xW.maxValue = dCath_WX[xW.maxWire];
      yW.maxValue = dCath_WY[yW.maxWire];

      xE.rawCenter = cathResp.getWirePosEX(xE.maxWire);
      yE.rawCenter = cathResp.getWirePosEY(yE.maxWire);
      xW.rawCenter = cathResp.getWirePosWX(xW.maxWire);
      yW.rawCenter = cathResp.getWirePosWY(yW.maxWire);


      std::vector <Double_t> eta = posmap.getInterpolatedEta( xE.center,yE.center,
				      xW.center,yW.center );

      //for (UInt_t iii=0; iii<eta.size(); iii++) std::cout << eta[iii] << std::endl;
    
      
      //MWPC triggers
      if (mwpcE.MWPCEnergyE>MWPCAnodeThreshold) EMWPCTrigger=true;
      if (mwpcE.MWPCEnergyW>MWPCAnodeThreshold) WMWPCTrigger=true;


      //If there is no wirechamber trigger, skip the event
      if ( !(EMWPCTrigger || WMWPCTrigger) ) continue;
      
      Double_t pmtEnergyLowerLimit = 1.; //To put a hard cut on the weight
    
      std::vector<Double_t> ADCvecE(4,0.);
      std::vector<Double_t> ADCvecW(4,0.);

      //East Side smeared PMT energies
    
    
      for (UInt_t p=0; p<4; p++) {
	if ( edep.EdepE>0. ) { //Check to make sure that there is light to see in the scintillator
	
	  if (eta[p]>0.) {

	    pmt.etaEvis[p] = (1./(alpha[p]*g_d*g_rest)) * (rand3->Poisson(g_rest*rand2->Poisson(g_d*rand1->Poisson(alpha[p]*eta[p]*edepQ.EdepQE))));
	    Double_t ADC = linCurve.applyInverseLinCurve(p, pmt.etaEvis[p]);// + rand0->Gaus(0.,pedestals[p][1]); //Take into account non-zero width of the pedestal
	    ADCvecE[p] = ADC;
	    pmt.etaEvis[p] = linCurve.applyLinCurve(p, ADC);
	    pmt.Evis[p] = pmt.etaEvis[p]/eta[p];

	  }

	  else { //To avoid dividing by zero.. these events won't be used in analysis since they are outside the fiducial cut
	  
	    pmt.etaEvis[p] = (1./(alpha[p]*g_d*g_rest)) * (rand3->Poisson(g_rest*rand2->Poisson(g_d*rand1->Poisson(alpha[p]*edepQ.EdepQE))));
	    Double_t ADC = linCurve.applyInverseLinCurve(p, pmt.etaEvis[p]);// + rand0->Gaus(0.,pedestals[p][1]); //Take into account non-zero width of the pedestal
	    ADCvecE[p] = ADC;
	    pmt.etaEvis[p] = linCurve.applyLinCurve(p, ADC);
	    pmt.Evis[p] = pmt.etaEvis[p];

	  }
	
	  pmt.nPE[p] = alpha[p]*pmt.etaEvis[p];
	}
	// If eQ is 0...
	else {
	  pmt.etaEvis[p] = 0.;
	  pmt.Evis[p] = 0.;
	  pmt.nPE[p] = 0.;
	}
      }
  
      
      //Calculate the weighted energy on a side
      Double_t numer=0., denom=0.;
      for (UInt_t p=0;p<4;p++) {
	numer += (pmtQuality[p] ) ? pmt.nPE[p] : 0.;
	denom += (pmtQuality[p] ) ? eta[p]*alpha[p] : 0.;
      }

      //Now we apply the trigger probability
      Double_t totalEnE = denom>0. ? numer/denom : 0.;
      evis.EvisE = totalEnE;

      triggRandVec[0] = randPMTE1->Rndm();
      triggRandVec[1] = randPMTE2->Rndm();
      triggRandVec[2] = randPMTE3->Rndm();
      triggRandVec[3] = randPMTE4->Rndm();
    
      if (  trigger.decideEastTrigger(ADCvecE,triggRandVec) ) EastScintTrigger = true; 
      
      //West Side
      for ( UInt_t p=4; p<8; p++ ) {
	if ( !(p==5 && runNumber>16983 && runNumber<17249)  &&  edep.EdepW>0. ) { //Check to make sure that there is light to see in the scintillator and that run isn't one where PMTW2 was dead
	
	  if (eta[p]>0.) {

	    pmt.etaEvis[p] = (1./(alpha[p]*g_d*g_rest)) * (rand3->Poisson(g_rest*rand2->Poisson(g_d*rand1->Poisson(alpha[p]*eta[p]*edepQ.EdepQW))));
	    Double_t ADC = linCurve.applyInverseLinCurve(p, pmt.etaEvis[p]);// + rand0->Gaus(0.,pedestals[p][1]); //Take into account non-zero width of the pedestal
	    ADCvecW[p-4] = ADC;
	    //cout << ADCvecW[p-4] << endl;
	    pmt.etaEvis[p] = linCurve.applyLinCurve(p, ADC);	  
	    pmt.Evis[p] = pmt.etaEvis[p]/eta[p];

	  }

	  else { //To avoid dividing by zero.. these events won't be used in analysis since they are outside the fiducial cut

	    pmt.etaEvis[p] = (1./(alpha[p]*g_d*g_rest)) * (rand3->Poisson(g_rest*rand2->Poisson(g_d*rand1->Poisson(alpha[p]*edepQ.EdepQW))));
	    Double_t ADC = linCurve.applyInverseLinCurve(p, pmt.etaEvis[p]);// + rand0->Gaus(0.,pedestals[p][1]); //Take into account non-zero width of the pedestal
	    ADCvecW[p-4] = ADC;
	    pmt.etaEvis[p] = linCurve.applyLinCurve(p, ADC);
	    pmt.Evis[p] = pmt.etaEvis[p];
	  
	  }
	
	  pmt.nPE[p] = alpha[p]*pmt.etaEvis[p];
	}
	// If PMT is dead and EQ=0...
	else {
	  pmt.etaEvis[p] = 0.;
	  pmt.Evis[p] = 0.;
	  pmt.nPE[p] = 0.;
	}
      }
      
      
      //Calculate the total weighted energy
      numer=denom=0.;
      for (UInt_t p=4;p<8;p++) {
	numer += ( pmtQuality[p] ) ? pmt.nPE[p] : 0.;
	denom += ( pmtQuality[p] ) ? eta[p]*alpha[p] : 0.;
      }
      //Now we apply the trigger probability
      Double_t totalEnW = denom>0. ? numer/denom : 0.;
      evis.EvisW = totalEnW;

      triggRandVec[0] = randPMTW1->Rndm();
      triggRandVec[1] = randPMTW2->Rndm();
      triggRandVec[2] = randPMTW3->Rndm();
      triggRandVec[3] = randPMTW4->Rndm();
    
      if ( trigger.decideWestTrigger(ADCvecW,triggRandVec) ) WestScintTrigger = true;

      //if (totalEnW<25.) std::cout << totalEnW << " " << WestScintTrigger << "\n";

            
      //Fill proper total event histogram based on event type
      PID=6; //This is an unidentified particle
      type=4; //this is a lost event
      side=2; //This means there are no scintillator triggers

      //Type 0 East
      if (EastScintTrigger && EMWPCTrigger && !WestScintTrigger && !WMWPCTrigger) {
	PID=1;
	type=0;
	side=0;
	//finalEn[0]->Fill(evis.EvisE);
	//cout << "Type 0 East E = " << totalEnE << endl;
      }
      //Type 0 West
      else if (WestScintTrigger && WMWPCTrigger && !EastScintTrigger && !EMWPCTrigger) {
	PID=1;
	type=0;
	side=1;
	//finalEn[1]->Fill(totalEnW);
      }
      //Type 1 
      else if (EastScintTrigger && EMWPCTrigger && WestScintTrigger && WMWPCTrigger) {
	PID=1;
	type=1;
	//East
	if (Time.timeE<Time.timeW) {
	  //finalEn[2]->Fill(totalEnE);
	  side=0;
	}
	//West
	else if (Time.timeE>Time.timeW) {
	  //finalEn[3]->Fill(totalEnW);
	  side=1;
	}
      }
      //Type 2/3 East
      else if (EastScintTrigger && EMWPCTrigger && !WestScintTrigger && WMWPCTrigger) {
	PID=1;
	type=2;
	side=0;
	//finalEn[4]->Fill(totalEnE);
	//cout << "Type 2/3 East E = " << totalEnE << endl;
      }
      //Type 2/3 West
      else if (!EastScintTrigger && EMWPCTrigger && WestScintTrigger && WMWPCTrigger) {
	PID=1;
	type=2;
	side=1;
	//finalEn[5]->Fill(totalEnW);
	//cout << "Type 2/3 East W = " << totalEnW << endl;
      }   
      //Gamma events and missed events (Type 4)
      else {
	if (!WMWPCTrigger && !EMWPCTrigger) {
	  if (EastScintTrigger && !WestScintTrigger) {
	    PID=0;
	    type=0;
	    side=0;
	  }
	  else if (!EastScintTrigger && WestScintTrigger) {
	    PID=0;
	    type=0;
	    side=1;
	  }
	  else if (EastScintTrigger && WestScintTrigger) {
	    PID=0;
	    type=1;
	    if (Time.timeE<Time.timeW) {
	      side=0;
	    }
	    else {
	      side=1;
	    }
	  }
	  else {
	    PID=6;
	    type=4;
	    side=2;
	  }
	}
	else {
	  PID=1;
	  type=4;
	  side = (WMWPCTrigger && EMWPCTrigger) ? 2 : (WMWPCTrigger ? 1 : 0); //Side==2 means the event triggered both wirechambers, but neither scint
	}
      }
    
      //Calculate Erecon
      Erecon = -1.;
      Int_t typeIndex = (type==0 || type==4) ? 0:(type==1 ? 1:2); //for retrieving the parameters from EQ2Etrue
      if (side==0) {
	Double_t totalEvis = type==1 ? (evis.EvisE+evis.EvisW):evis.EvisE;
	if (evis.EvisE>0. && totalEvis>0.) {
	  Erecon = eRecon.getErecon(0,typeIndex,totalEvis);
	}
	else Erecon=-1.;
      }
      if (side==1) {
	Double_t totalEvis = type==1 ? (evis.EvisE+evis.EvisW):evis.EvisW;
	if (evis.EvisW>0. && totalEvis>0.) {
	  Erecon = eRecon.getErecon(1,typeIndex,totalEvis);
	}
	else Erecon=-1.;
      }    

      //Write to tree so Xuan can also process it
      //if (PID>=0) tree->Fill();
	    

      //***********************************************************************************************************************
      // Filling rate histograms with "good" events to calculate the corrections

      if ( PID==1 && side<2 && type<4 && Erecon>0. ) {


	if ( type!=0 ) {
	  if (xE.mult<1 || yE.mult<1 || xW.mult<1 || yW.mult<1)  continue;
	}
	else {
	  if ( side==0 ) {
	    if ( xE.mult<1 || yE.mult<1 )  continue;
	  }
	  else if ( side==1 ) {
	    if ( xW.mult<1 || yW.mult<1 )  continue;
	  }
	}		
	
	
	
	// These are the reconstructed position using the cathode segments
	double r2E = xE.center*xE.center + yE.center*yE.center;
	double r2W = xW.center*xW.center + yW.center*yW.center;

	if ( r2E<(fidCut*fidCut) && r2W<(fidCut*fidCut ) ) {
	
	
	  //Type 0
	  if (type==0) hist[0][side]->Fill(Erecon);
	
	  //Type 1
	  if (type==1) hist[1][side]->Fill(Erecon);
	
	  //Type 23
	  if (type==2) { 
	    
	    hist[2][side]->Fill(Erecon);
	      
	    //Type 2/3 separation
	    if (side==0) {
	      type = sep.separate23(mwpcE.MWPCEnergyE);
	      side = type==2 ? 1 : 0;
	    }
	    else if (side==1) {
	      type = sep.separate23(mwpcE.MWPCEnergyW);
	      side = type==2 ? 0 : 1;
	    }
	    
	    //Separated Type 2
	    if (type==2) hist2[side]->Fill(Erecon); 
	    //Separated Type 3
	    if (type==3) hist3[side]->Fill(Erecon);
	  }

	}
      }

    }

    
    cout << endl;

    //setting the events
    if ( runType=="A2" || runType=="A10" || runType=="B5" || runType=="B7" ) evtOFF = evt;
    else evtON = evt;


    delete chain; delete seed; delete rand0; delete rand1; delete rand2; delete rand3;
    delete randPMTE1; delete randPMTE2; delete randPMTE3; delete randPMTE4; delete randPMTW1; delete randPMTW2; delete randPMTW3; delete randPMTW4;

    //    outfile->Write();
    //outfile->Close();
    

    //Fill all the vectors according to their event types and analysis choices

    

    for (int anaCh = 1; anaCh<11 ; anaCh++) {

      int type_low=-1, type_high=-1;
      bool sep23 = false;
      bool sep2 = false;
      bool sep3 = false;

      if (anaCh==1) { type_low=0; type_high=2; }                             // "A" All types, 2/3 not separated
      else if (anaCh==3 || anaCh==5) { type_low=0; type_high=1; sep23 = sep2 = sep3 = true;} // "C & E" All event types, 2/3 separated
      else if (anaCh==2) { type_low=0; type_high=1;}                         // "B" Type 0 and 1
      else if (anaCh==4) { type_low=0; type_high=0;}                         // "D" Type 0
      else if (anaCh==6) { type_low=1; type_high=1;}                         // "F" Type 1
      else if (anaCh==7) { type_low=type_high=2; }                           // "G" Type 2/3 not separated
      else if (anaCh==8) { sep23 = sep2 = sep3 = true; }                 // "H" Type 2/3, separated
      else if (anaCh==9) { sep23 = sep2 = true;}                   // "J" Type 2
      else if (anaCh==10) { sep23 = sep3 = true;}                  // "K" Type 3
      
      for (unsigned int side=0; side<2; side++) {
	for (unsigned int bin=1; bin<=120; bin++) {
	  for (int type=type_low; type<=type_high; type++) {	 

	    if ( runType == "A2" ) {
	      A2[anaCh-1][side][bin-1] += type>-1 ? hist[type][side]->GetBinContent(bin) : 0.;
	      A2_err[anaCh-1][side][bin-1] += type>-1 ? power(hist[type][side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "A5" ) {
	      A5[anaCh-1][side][bin-1] += type>-1 ? hist[type][side]->GetBinContent(bin) : 0.;
	      A5_err[anaCh-1][side][bin-1] += type>-1 ? power(hist[type][side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "A7" ) {
	      A7[anaCh-1][side][bin-1] += type>-1 ? hist[type][side]->GetBinContent(bin) : 0.;
	      A7_err[anaCh-1][side][bin-1] += type>-1 ? power(hist[type][side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "A10" ) {
	      A10[anaCh-1][side][bin-1] += type>-1 ? hist[type][side]->GetBinContent(bin) : 0.;
	      A10_err[anaCh-1][side][bin-1] += type>-1 ? power(hist[type][side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "B2" ) {
	      B2[anaCh-1][side][bin-1] += type>-1 ? hist[type][side]->GetBinContent(bin) : 0.;
	      B2_err[anaCh-1][side][bin-1] += type>-1 ? power(hist[type][side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "B5" ) {
	      B5[anaCh-1][side][bin-1] += type>-1 ? hist[type][side]->GetBinContent(bin) : 0.;
	      B5_err[anaCh-1][side][bin-1] += type>-1 ? power(hist[type][side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "B7" ) {
	      B7[anaCh-1][side][bin-1] += type>-1 ? hist[type][side]->GetBinContent(bin) : 0.;
	      B7_err[anaCh-1][side][bin-1] += type>-1 ? power(hist[type][side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "B10" ) {
	      B10[anaCh-1][side][bin-1] += type>-1 ? hist[type][side]->GetBinContent(bin) : 0.;
	      B10_err[anaCh-1][side][bin-1] += type>-1 ? power(hist[type][side]->GetBinError(bin),2) : 0.;
	    }

	  }
	  
	  if ( sep23 ) {
	    if ( runType == "A2" ) {
	      A2[anaCh-1][side][bin-1] += sep2 ? hist2[side]->GetBinContent(bin) : 0.;
	      A2_err[anaCh-1][side][bin-1] += sep2 ? power(hist2[side]->GetBinError(bin),2) : 0.;
	      A2[anaCh-1][side][bin-1] += sep3 ? hist3[side]->GetBinContent(bin) : 0.;
	      A2_err[anaCh-1][side][bin-1] += sep3 ? power(hist3[side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "A5" ) {
	      A5[anaCh-1][side][bin-1] += sep2 ? hist2[side]->GetBinContent(bin) : 0.;
	      A5_err[anaCh-1][side][bin-1] += sep2 ? power(hist2[side]->GetBinError(bin),2) : 0.;
	      A5[anaCh-1][side][bin-1] += sep3 ? hist3[side]->GetBinContent(bin) : 0.;
	      A5_err[anaCh-1][side][bin-1] += sep3 ? power(hist3[side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "A7" ) {
	      A7[anaCh-1][side][bin-1] += sep2 ? hist2[side]->GetBinContent(bin) : 0.;
	      A7_err[anaCh-1][side][bin-1] += sep2 ? power(hist2[side]->GetBinError(bin),2) : 0.;
	      A7[anaCh-1][side][bin-1] += sep3 ? hist3[side]->GetBinContent(bin) : 0.;
	      A7_err[anaCh-1][side][bin-1] += sep3 ? power(hist3[side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "A10" ) {
	      A10[anaCh-1][side][bin-1] += sep2 ? hist2[side]->GetBinContent(bin) : 0.;
	      A10_err[anaCh-1][side][bin-1] += sep2 ? power(hist2[side]->GetBinError(bin),2) : 0.;
	      A10[anaCh-1][side][bin-1] += sep3 ? hist3[side]->GetBinContent(bin) : 0.;
	      A10_err[anaCh-1][side][bin-1] += sep3 ? power(hist3[side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "B2" ) {
	      B2[anaCh-1][side][bin-1] += sep2 ? hist2[side]->GetBinContent(bin) : 0.;
	      B2_err[anaCh-1][side][bin-1] += sep2 ? power(hist2[side]->GetBinError(bin),2) : 0.;
	      B2[anaCh-1][side][bin-1] += sep3 ? hist3[side]->GetBinContent(bin) : 0.;
	      B2_err[anaCh-1][side][bin-1] += sep3 ? power(hist3[side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "B5" ) {
	      B5[anaCh-1][side][bin-1] += sep2 ? hist2[side]->GetBinContent(bin) : 0.;
	      B5_err[anaCh-1][side][bin-1] += sep2 ? power(hist2[side]->GetBinError(bin),2) : 0.;
	      B5[anaCh-1][side][bin-1] += sep3 ? hist3[side]->GetBinContent(bin) : 0.;
	      B5_err[anaCh-1][side][bin-1] += sep3 ? power(hist3[side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "B7" ) {
	      B7[anaCh-1][side][bin-1] += sep2 ? hist2[side]->GetBinContent(bin) : 0.;
	      B7_err[anaCh-1][side][bin-1] += sep2 ? power(hist2[side]->GetBinError(bin),2) : 0.;
	      B7[anaCh-1][side][bin-1] += sep3 ? hist3[side]->GetBinContent(bin) : 0.;
	      B7_err[anaCh-1][side][bin-1] += sep3 ? power(hist3[side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "B10" ) {
	      B10[anaCh-1][side][bin-1] += sep2 ? hist2[side]->GetBinContent(bin) : 0.;
	      B10_err[anaCh-1][side][bin-1] += sep2 ? power(hist2[side]->GetBinError(bin),2) : 0.;
	      B10[anaCh-1][side][bin-1] += sep3 ? hist3[side]->GetBinContent(bin) : 0.;
	      B10_err[anaCh-1][side][bin-1] += sep3 ? power(hist3[side]->GetBinError(bin),2) : 0.;
	    }
	  }
	    
	  
	  if ( runType == "A2"  ) A2_err[anaCh-1][side][bin-1] = sqrt(A2_err[anaCh-1][side][bin-1]);
	  if ( runType == "A5"  ) A5_err[anaCh-1][side][bin-1] = sqrt(A5_err[anaCh-1][side][bin-1]);
	  if ( runType == "A7"  ) A7_err[anaCh-1][side][bin-1] = sqrt(A7_err[anaCh-1][side][bin-1]);
	  if ( runType == "A10" ) A10_err[anaCh-1][side][bin-1] = sqrt(A10_err[anaCh-1][side][bin-1]);
	  if ( runType == "B2"  ) B2_err[anaCh-1][side][bin-1] = sqrt(B2_err[anaCh-1][side][bin-1]);
	  if ( runType == "B5"  ) B5_err[anaCh-1][side][bin-1] = sqrt(B5_err[anaCh-1][side][bin-1]);
	  if ( runType == "B7"  ) B7_err[anaCh-1][side][bin-1] = sqrt(B7_err[anaCh-1][side][bin-1]);
	  if ( runType == "B10" ) B10_err[anaCh-1][side][bin-1]=sqrt(B10_err[anaCh-1][side][bin-1]);
	  //std::cout << anaChoice_A2[side][bin] << " " << anaChoice_A2_err[side][bin] << std::endl;
	}
      }
    }
   

    specFile->Write();
    delete specFile;

    /*  for (int side = 0; side<2; side++) {
	  for (int type=0; type<3; type++) {
	    if (hist[type][side]) delete hist[type][side];
	  }
	  if (hist2[side]) delete hist2[side]; 
5A	  if (hist3[side]) delete hist3[side];
}*/
   
  }
    
    
  // Now we calculate all of the asymmetries bin-by-bin for all anaChoices
  for (int anaCh = 1; anaCh<11 ; anaCh++) {
    
    //Octet Asymmetry
    for (unsigned int bin=0; bin<120; bin++) {
      
      double R = 0.;
      double deltaR = 0.;
      double sfON[2]={0.};
      double sfOFF[2]={0.};
      double sfON_err[2]={0.};
      double sfOFF_err[2]={0.};
      
      for (unsigned int side=0; side<2; side++) {
	double weightsum=0.;
	
	if (A2[anaCh-1][side][bin]>0. && A10[anaCh-1][side][bin]>0. && B5[anaCh-1][side][bin]>0. && B7[anaCh-1][side][bin]>0. ) { 

	  weightsum = 1./A2[anaCh-1][side][bin] + 1./A10[anaCh-1][side][bin] + 1./B5[anaCh-1][side][bin] + 1./B7[anaCh-1][side][bin];
	
	}
	else weightsum = 0.;

	sfOFF[side] = weightsum>0. ? 4./weightsum : 0.;
	sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
	
	if (A5[anaCh-1][side][bin]>0. && A7[anaCh-1][side][bin]>0. && B2[anaCh-1][side][bin]>0. && B10[anaCh-1][side][bin]>0. ) { 
	  //std::cout << "made it in here\n";
	  weightsum = 1./A5[anaCh-1][side][bin] + 1./A7[anaCh-1][side][bin] + 1./B2[anaCh-1][side][bin] + 1./B10[anaCh-1][side][bin];
	
	}
	else weightsum = 0.;
	//if (anaCh==1 && bin<80) std::cout << weightsum << std::endl;
	sfON[side] = weightsum>0. ? 4./weightsum : 0.;
	sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
		
      }

      //      if (anaCh==1 && bin<80) std::cout <<  sfOFF[0] << " " <<  sfOFF[1] << " " << sfON[0] << " " <<  sfON[1] << std::endl;


      if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
	R = sfOFF[1]*sfON[0]/(sfON[1]*sfOFF[0]);
	deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
	oct_A_SR[anaCh-1][bin] = (1.-sqrt(R))/(1+sqrt(R));
	oct_A_SR_err[anaCh-1][bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
      }
      else {
	oct_A_SR[anaCh-1][bin] = 0.;
	oct_A_SR_err[anaCh-1][bin] = 0.;
      }
    }

    /*//QuartetA Asymmetry
    for (unsigned int bin=0; bin<120; bin++) {
      
      double R = 0.;
      double deltaR = 0.;
      double sfON[2]={0.};
      double sfOFF[2]={0.};
      double sfON_err[2]={0.};
      double sfOFF_err[2]={0.};
      
      for (unsigned int side=0; side<2; side++) {
	double weightsum=0.;
	
	sfOFF[side] = (A2_err[anaCh-1][side][bin]>0.?power(1./A2_err[anaCh-1][side][bin],2)*A2[anaCh-1][side][bin]:0.) + (A10_err[anaCh-1][side][bin]>0.?power(1./A10_err[anaCh-1][side][bin],2)*A10[anaCh-1][side][bin]:0.) + (B5_err[anaCh-1][side][bin]>0.?power(1./B5_err[anaCh-1][side][bin],2)*B5[anaCh-1][side][bin]:0.) + (B7_err[anaCh-1][side][bin]>0.?power(1./B7_err[anaCh-1][side][bin],2)*B7[anaCh-1][side][bin]:0.);
	weightsum = (A2_err[anaCh-1][side][bin]>0.?power(1./A2_err[anaCh-1][side][bin],2):0.) + (A10_err[anaCh-1][side][bin]>0.?power(1./A10_err[anaCh-1][side][bin],2):0.) + (B5_err[anaCh-1][side][bin]>0.?power(1./B5_err[anaCh-1][side][bin],2):0.) + (B7_err[anaCh-1][side][bin]>0.?power(1./B7_err[anaCh-1][side][bin],2):0.);
	
	sfOFF[side] = weightsum>0. ? sfOFF[side]/weightsum : 0.;
	sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
	weightsum=0.;
	sfON[side] = (A5_err[anaCh-1][side][bin]>0.?power(1./A5_err[anaCh-1][side][bin],2)*A5[anaCh-1][side][bin]:0.) + (A7_err[anaCh-1][side][bin]>0.?power(1./A7_err[anaCh-1][side][bin],2)*A7[anaCh-1][side][bin]:0.) + (B2_err[anaCh-1][side][bin]>0.?power(1./B2_err[anaCh-1][side][bin],2)*B2[anaCh-1][side][bin]:0.) + (B10_err[anaCh-1][side][bin]>0.?power(1./B10_err[anaCh-1][side][bin],2)*B10[anaCh-1][side][bin]:0.);
	weightsum = (A5_err[anaCh-1][side][bin]>0.?power(1./A5_err[anaCh-1][side][bin],2):0.) + (A7_err[anaCh-1][side][bin]>0.?power(1./A7_err[anaCh-1][side][bin],2):0.) + (B2_err[anaCh-1][side][bin]>0.?power(1./B2_err[anaCh-1][side][bin],2):0.) + (B10_err[anaCh-1][side][bin]>0.?power(1./B10_err[anaCh-1][side][bin],2):0.);
	
	
	sfON[side] = weightsum>0. ? sfON[side]/weightsum : 0.;
	sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
	
      }
      if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
	R = sfOFF[1]*sfON[0]/(sfON[1]*sfOFF[0]);
	deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
	oct_A_SR[anaCh-1][bin] = (1.-sqrt(R))/(1+sqrt(R));
	oct_A_SR_err[anaCh-1][bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
      }
      else {
	oct_A_SR[anaCh-1][bin] = 0.;
	oct_A_SR_err[anaCh-1][bin] = 0.;
      }
      }*/


    
  }
  
};

int main(int argc, char *argv[]) {

  int octet = atoi(argv[1]);

  calcDeltaExp(octet);
  writeCorrectionFactorsToFile(octet);

  //tests
  /*UInt_t XePeriod = getXeRunPeriod(atoi(argv[1]));
  vector < vector <Double_t> > triggerFunc = getTriggerFunctionParams(XePeriod,7);
  Double_t triggProbE = triggerProbability(triggerFunc[0],25.);
  Double_t triggProbW = triggerProbability(triggerFunc[1],25.);
  cout << triggProbE << " " << triggProbW << endl;*/


}
  
