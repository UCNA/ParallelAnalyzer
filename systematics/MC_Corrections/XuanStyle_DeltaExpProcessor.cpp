 /* Code to take a run number, retrieve it's runperiod, and construct the 
weighted spectra which would be seen as a reconstructed energy on one side
of the detector. Also applies the trigger functions */

#include "DeltaExpProcessor.h"
#include "posMapReader.h"
#include "positionMapHandler.hh"
#include "sourcePeaks.h"
#include "runInfo.h"
#include "calibrationTools.hh"
#include "TriggerMap.hh"
#include "../Asymmetry/SQLinterface.hh"
#include "BetaSpectrum.hh"
#include "MBUtils.hh"

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

void writeCorrectionFactorsToFile(Int_t octet) {
  TString fn_base = TString::Format("lowStats_corrections/ThOverProc_Octet-%i_Analysis-",octet); 
  //  TString fn_base = TString::Format("DeltaExp_OctetByOctetCorrections/ThOverProc_Octet-%i_Analysis-",octet); 
  ofstream oct; 
  ofstream quartA;
  ofstream quartB; 
  ofstream pairA1; 
  ofstream pairA2; 
  ofstream pairB1; 
  ofstream pairB2; 

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
 
  std::string dbAddress = std::string(getenv("UCNADBADDRESS"));
  std::string dbname = std::string(getenv("UCNADB"));
  std::string dbUser = std::string(getenv("UCNADBUSER"));
  std::string dbPass = std::string(getenv("UCNADBPASS"));
  
  char cmd[200];
  sprintf(cmd,"SELECT flipper FROM run WHERE run_number=%i;",run);

  SQLdatabase *db = new SQLdatabase(dbname, dbAddress, dbUser, dbPass);
  db->fetchQuery(cmd);
  std::string flipperStatus = db->returnQueryEntry();
  delete db;

  std::cout << flipperStatus << std::endl;
  if (flipperStatus=="On") return 1;
  else if (flipperStatus=="Off") return -1;
  else {
    std::cout <<  "Polarization isn't applicaple or you chose a Depol Run";
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

vector <Int_t> getPMTQuality(Int_t runNumber) {
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

  
//Get the conversion from EQ2Etrue
std::vector < std::vector < std::vector <double> > > getEQ2EtrueParams(int runNumber) {
  ifstream infile;
  std::string basePath = getenv("ANALYSIS_CODE"); 
  if (runNumber<16000) basePath+=std::string("/simulation_comparison/EQ2EtrueConversion/2011-2012_EQ2EtrueFitParams.dat");
  else if (runNumber<20000) basePath+=std::string("/simulation_comparison/EQ2EtrueConversion/2011-2012_EQ2EtrueFitParams.dat");
  else if (runNumber<21628 && runNumber>21087) basePath+=std::string("/simulation_comparison/EQ2EtrueConversion/2012-2013_isobutane_EQ2EtrueFitParams.dat");
  else if (runNumber<24000) basePath+=std::string("/simulation_comparison/EQ2EtrueConversion/2012-2013_EQ2EtrueFitParams.dat");
  else {
    std::cout << "Bad runNumber passed to getEQ2EtrueParams\n";
    exit(0);
  }
  infile.open(basePath.c_str());
  std::vector < std::vector < std::vector < double > > > params;
  params.resize(2,std::vector < std::vector < double > > (3, std::vector < double > (6,0.)));

  char holdType[10];
  int side=0, type=0;
  while (infile >> holdType >> params[side][type][0] >> params[side][type][1] >> params[side][type][2] >> params[side][type][3] >> params[side][type][4] >> params[side][type][5]) { 
    std::cout << holdType << " " << params[side][type][0] << " " << params[side][type][1] << " " << params[side][type][2] << " " << params[side][type][3] << " " << params[side][type][4] << " " << params[side][type][5] << std::endl;
    type+=1;
    if (type==3) {type=0; side=1;}
  }
  return params;
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
  tree->Branch("ScintPos",&scint_pos,"ScintPosE[3]/D:ScintPosW[3]");
  tree->Branch("ScintPosAdjusted",&scint_pos_adj,"ScintPosAdjE[3]/D:ScintPosAdjW[3]");
  tree->Branch("PMT",&pmt,"Evis0/D:Evis1:Evis2:Evis3:Evis4:Evis5:Evis6:Evis7:etaEvis0/D:etaEvis1:etaEvis2:etaEvis3:etaEvis4:etaEvis5:etaEvis6:etaEvis7:nPE0/D:nPE1:nPE2:nPE3:nPE4:nPE5:nPE6:nPE7");
  
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


  
  

  for (auto& runNumber : betaRuns) {
  
   
    

    std::cout << "Processing delta_exp for run " << runNumber << "...\n";

    Double_t fidCut = 50.;

    
    

    std::string runType = getRunTypeFromOctetFile(runNumber);
    
    if ( runType=="BAD" ) { std::cout << "BAD RUNTYPE FOUND FOR RUN\n"; exit(0); }
    

    //make evt type histograms
    
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


    TFile *f = new TFile(TString::Format("%s/beta/revCalSim_%i_Beta.root",getenv("REVCALSIM"),runNumber),"READ");
    TTree *tree = (TTree*)f->Get("revCalSim");

    tree->SetBranchAddress("PID", &PID);
    tree->SetBranchAddress("side", &side);
    tree->SetBranchAddress("type", &type);
    tree->SetBranchAddress("Erecon", &Erecon);
    
    tree->SetBranchAddress("Evis",&evis);
    tree->SetBranchAddress("Edep",&edep);
    tree->SetBranchAddress("EdepQ",&edepQ);
    tree->SetBranchAddress("primKE",&primKE);
    tree->SetBranchAddress("primTheta",&primTheta);
    tree->SetBranchAddress("AsymWeight",&AsymWeight);
    
    tree->SetBranchAddress("time",&Time);
    tree->SetBranchAddress("MWPCEnergy",&mwpcE);
    tree->SetBranchAddress("MWPCPos",&mwpc_pos);
    tree->SetBranchAddress("ScintPos",&scint_pos);
    tree->SetBranchAddress("ScintPosAdjusted",&scint_pos_adj);
    tree->SetBranchAddress("PMT",&pmt);   

    int nClipped_EX, nClipped_EY, nClipped_WX, nClipped_WY;
    tree->SetBranchAddress("nClipped_EX",&nClipped_EX);
    tree->SetBranchAddress("nClipped_EY",&nClipped_EY);
    tree->SetBranchAddress("nClipped_WX",&nClipped_WX);
    tree->SetBranchAddress("nClipped_WY",&nClipped_WY);
    
    Int_t nEvents = tree->GetEntriesFast();

    for (int i=0; i<nEvents; i++) {
      
      tree->GetEvent(i);


      if (PID!=1) continue;  

      if (type>3) continue;

      if (side>1) continue;

      if (Erecon<=0.) continue;

      //Cut clipped events
      if ( type!=0 ) {
	if (nClipped_EX>0 || nClipped_EY>0 || nClipped_WX>0 || nClipped_WY>0) continue;
      }
      else {
	if ( side==0 && ( nClipped_EX>0 || nClipped_EY>0 ) ) continue;
	else if ( side==1 && ( nClipped_WX>0 || nClipped_WY>0 ) ) continue;
      }

      if ( sqrt(mwpc_pos.MWPCPosE[0]*mwpc_pos.MWPCPosE[0] + mwpc_pos.MWPCPosE[1]*mwpc_pos.MWPCPosE[1])*sqrt(0.6)*10.>fidCut
      	   || sqrt(mwpc_pos.MWPCPosW[0]*mwpc_pos.MWPCPosW[0] + mwpc_pos.MWPCPosW[1]*mwpc_pos.MWPCPosW[1])*sqrt(0.6)*10.>fidCut ) continue;
      //if ( sqrt(scint_pos.ScintPosE[0]*scint_pos.ScintPosE[0]+scint_pos.ScintPosE[1]+scint_pos.ScintPosE[1])*sqrt(0.6)*10.>fidCut
      //	   || sqrt(scint_pos.ScintPosW[0]*scint_pos.ScintPosW[0]+scint_pos.ScintPosW[1]*scint_pos.ScintPosW[1])*sqrt(0.6)*10.>fidCut ) continue;

      
      
      //      std::cout << i << " " << PID << " " << side << " " << type << " " << Erecon << std::endl;


      // This separates 2/3 if true. Only make it false if they have already been separated in processed file
      if (true) {
	if (type==2) {

	  if (side==0) {
	    type = separate23(side,mwpcE.MWPCEnergyE);
	    side = type==2 ? 1 : 0;
	  }
	  else if (side==1) {
	    type = separate23(side,mwpcE.MWPCEnergyW);
	    side = type==2 ? 0 : 1;
	  }
	  	  
	}
      }

      //Type 0                                                                                                                             
      if (type==0) hist[0][side]->Fill(Erecon);

      //Type 1                                                                                                                             
      if (type==1) hist[1][side]->Fill(Erecon);
      
      

      //Type 23                                                                                                                            
      if (type==2 || type==3) {
	if (side==0) {
	  //hist[2][0]->Fill(Erecon);
	  if (type==3) hist[2][0]->Fill(Erecon);
	  else hist[2][1]->Fill(Erecon);
	}
	else if (side==1) {
	  if (type==3) hist[2][1]->Fill(Erecon);
	  else hist[2][0]->Fill(Erecon);
	}
      }

      //Type 2                                                                                                                             
      if (type==2) hist2[side]->Fill(Erecon);

      //Type 3                                                                                                                             
      if (type==3) hist3[side]->Fill(Erecon);
      
    }

    f->Close();
    if (f) delete f;
  

   

    //Fill all the vectors according to their event types and analysis choices

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

	  
	  A2_err[anaCh-1][side][bin-1] = sqrt(A2_err[anaCh-1][side][bin-1]);
	  A5_err[anaCh-1][side][bin-1] = sqrt(A5_err[anaCh-1][side][bin-1]);
	  A7_err[anaCh-1][side][bin-1] = sqrt(A7_err[anaCh-1][side][bin-1]);
	  A10_err[anaCh-1][side][bin-1] = sqrt(A10_err[anaCh-1][side][bin-1]);
	  B2_err[anaCh-1][side][bin-1] = sqrt(B2_err[anaCh-1][side][bin-1]);
	  B5_err[anaCh-1][side][bin-1] = sqrt(B5_err[anaCh-1][side][bin-1]);
	  B7_err[anaCh-1][side][bin-1] = sqrt(B7_err[anaCh-1][side][bin-1]);
	  B10_err[anaCh-1][side][bin-1]=sqrt(B10_err[anaCh-1][side][bin-1]);
	  //std::cout << anaChoice_A2[side][bin] << " " << anaChoice_A2_err[side][bin] << std::endl;
	}
      }

    }
   
    for (int side = 0; side<2; side++) {
	  for (int type=0; type<3; type++) {
	    delete hist[type][side];
	  }
	  delete hist2[side]; delete hist3[side];
    }
   
  }
    
    
  // Now we calculate all of the asymmetries bin-by-bin for all anaChoices
  for (int anaCh = 1; anaCh<11 ; anaCh++) {
    
    //Octet Asymmetry
    for (unsigned int bin=0; bin<120; bin++) {

      //      if (anaCh==1 && bin<80) std::cout << A2[anaCh-1][0][bin] << std::endl;
	  
      
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
  
