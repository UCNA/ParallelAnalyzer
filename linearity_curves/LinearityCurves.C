#include <vector>
#include <algorithm>
#include <map>
#include "../include/sourcePeaks.h"

unsigned int find_vec_location_int(std::vector<Int_t> vec, Int_t val)
{
  UInt_t size = vec.size();
  //for (int i=0;i<size;i++){cout << vec[i] << endl;}

  std::vector<Int_t>::iterator it = std::find(vec.begin(),vec.end(),val);
  if (it!=vec.end()) {
    return it - vec.begin();
  }
  else cout << "Can't locate run " << val << " in PMT Quality list" << endl; exit(0);
  
}

unsigned int find_vec_location_double(std::vector<Double_t> vec, Double_t val)
{
  UInt_t size = vec.size();
  //for (int i=0;i<size;i++){cout << vec[i] << endl;}

  std::vector<Double_t>::iterator it = std::find(vec.begin(),vec.end(),val);
  if (it!=vec.end()) {
    return it - vec.begin();
  }
  else cout << "Can't locate " << val << " in your vector!" << endl; exit(0);
  
}

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

bool isSourceInFidCut(Int_t run, string src, Double_t fidCut, Int_t side) {

  std::vector < std::vector <Double_t> > pos = returnSourcePosition(run,src);

  return ( fidCut*fidCut > ( pos[side][0]*pos[side][0] + pos[side][1]*pos[side][1] ) ) ? true : false;

}

void LinearityCurves(Int_t runPeriod, bool useTanh=false)
{

  bool quadratic = true;
  
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  // Style options
  gROOT->SetStyle("Plain");
  //gStyle->SetOptStat(11);
  gStyle->SetOptStat(0);
  gStyle->SetStatFontSize(0.020);
  gStyle->SetOptFit(0); // 1111
  gStyle->SetOptTitle(1);
  gStyle->SetTitleFontSize(0.05);
  //gStyle->SetTitleSize(0.04,"XY");
  //gStyle->SetTitleX(0.17);
  //gStyle->SetTitleAlign(13);
  gStyle->SetTitleOffset(1.20, "x");
  gStyle->SetTitleOffset(1.10, "y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetNdivisions(510,"X");
  gStyle->SetNdivisions(510,"Y");
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.12);

  //Run Range for this envelope
  //Int_t runLow = 17359;
  //Int_t runHigh = 19959;
  Int_t calibrationPeriod = runPeriod;

  // Read data file
  char temp[500];
  sprintf(temp, "%s/residuals/source_runs_RunPeriod_%i.dat",getenv("ANALYSIS_CODE"),calibrationPeriod);
  ifstream filein(temp);
  sprintf(temp, "%s/residuals/source_runs_peakErrors_RunPeriod_%i.dat",getenv("ANALYSIS_CODE"),calibrationPeriod);
  ifstream fileinErr(temp);
  sprintf(temp, "%s/residuals/SIM_source_runs_RunPeriod_%i.dat",getenv("ANALYSIS_CODE"),calibrationPeriod);
  ifstream simFileIn(temp);
  sprintf(temp, "%s/residuals/SIM_source_runs_peakErrors_RunPeriod_%i.dat",getenv("ANALYSIS_CODE"),calibrationPeriod);
  ifstream simFileInErr(temp);

  Int_t i = 0;
  Int_t run[500], run_hold=0;
  string sourceName[500], sourceName_hold="";
  Double_t adcE1[500], adcE2[500], adcE3[500], adcE4[500];
  Double_t adcW1[500], adcW2[500], adcW3[500], adcW4[500];
  Double_t adcE1_err[500], adcE2_err[500], adcE3_err[500], adcE4_err[500];
  Double_t adcW1_err[500], adcW2_err[500], adcW3_err[500], adcW4_err[500];
  Double_t EqE1[500], EqE2[500], EqE3[500], EqE4[500];
  Double_t EqW1[500], EqW2[500], EqW3[500], EqW4[500];
  Double_t EqE1_err[500], EqE2_err[500], EqE3_err[500], EqE4_err[500];
  Double_t EqW1_err[500], EqW2_err[500], EqW3_err[500], EqW4_err[500];
  Double_t resE1[500], resE2[500], resE3[500], resE4[500];
  Double_t resW1[500], resW2[500], resW3[500], resW4[500];
  Double_t err[500];

  TString Eq_hold[8], ADC_hold[8]; //this is a string to take care of a -nan issue with bad PMT values
  
  //Load the peak information
  while (filein >> run[i] >> sourceName[i]
	 >> ADC_hold[0] >> ADC_hold[1] >> ADC_hold[2] >> ADC_hold[3]
	 >> ADC_hold[4] >> ADC_hold[5] >> ADC_hold[6] >> ADC_hold[7]) {

    adcE1[i] = TMath::IsNaN(atof(ADC_hold[0].Data())) || atof(ADC_hold[0].Data())<5.0 ? -10000. : atof(ADC_hold[0].Data()); // if these are bad, we make them -10000
    adcE2[i] = TMath::IsNaN(atof(ADC_hold[1].Data())) || atof(ADC_hold[1].Data())<5.0 ? -10000. : atof(ADC_hold[1].Data()); // so that they don't affect the fits right away and have
    adcE3[i] = TMath::IsNaN(atof(ADC_hold[2].Data())) || atof(ADC_hold[2].Data())<5.0 ? -10000. : atof(ADC_hold[2].Data()); // a chance to recover on the next go around
    adcE4[i] = TMath::IsNaN(atof(ADC_hold[3].Data())) || atof(ADC_hold[3].Data())<5.0 ? -10000. : atof(ADC_hold[3].Data());
    adcW1[i] = TMath::IsNaN(atof(ADC_hold[4].Data())) || atof(ADC_hold[4].Data())<5.0 ? -10000. : atof(ADC_hold[4].Data());
    adcW2[i] = TMath::IsNaN(atof(ADC_hold[5].Data())) || atof(ADC_hold[5].Data())<5.0 ? -10000. : atof(ADC_hold[5].Data());
    adcW3[i] = TMath::IsNaN(atof(ADC_hold[6].Data())) || atof(ADC_hold[6].Data())<5.0 ? -10000. : atof(ADC_hold[6].Data());
    adcW4[i] = TMath::IsNaN(atof(ADC_hold[7].Data())) || atof(ADC_hold[7].Data())<5.0 ? -10000. : atof(ADC_hold[7].Data());

    while (!(run_hold==run[i] && sourceName_hold==sourceName[i])) {
      simFileIn >> run_hold >> sourceName_hold >> Eq_hold[0] >> Eq_hold[1] >> Eq_hold[2]
		>> Eq_hold[3] >> Eq_hold[4] >> Eq_hold[5] >> Eq_hold[6] >> Eq_hold[7];
    }
    //cout << run_hold << " " << sourceName_hold << " ";
    if (run_hold==run[i] && sourceName_hold==sourceName[i]) {
      
      EqE1[i] = TMath::IsNaN(atof(Eq_hold[0].Data())) || atof(Eq_hold[0].Data())<5.0 ? 0. : atof(Eq_hold[0].Data()); // If these are bad, we make them zero so we can see the issue
      EqE2[i] = TMath::IsNaN(atof(Eq_hold[1].Data())) || atof(Eq_hold[1].Data())<5.0 ? 0. : atof(Eq_hold[1].Data());
      EqE3[i] = TMath::IsNaN(atof(Eq_hold[2].Data())) || atof(Eq_hold[2].Data())<5.0 ? 0. : atof(Eq_hold[2].Data());
      EqE4[i] = TMath::IsNaN(atof(Eq_hold[3].Data())) || atof(Eq_hold[3].Data())<5.0 ? 0. : atof(Eq_hold[3].Data());
      EqW1[i] = TMath::IsNaN(atof(Eq_hold[4].Data())) || atof(Eq_hold[4].Data())<5.0 ? 0. : atof(Eq_hold[4].Data());
      EqW2[i] = TMath::IsNaN(atof(Eq_hold[5].Data())) || atof(Eq_hold[5].Data())<5.0 ? 0. : atof(Eq_hold[5].Data());
      EqW3[i] = TMath::IsNaN(atof(Eq_hold[6].Data())) || atof(Eq_hold[6].Data())<5.0 ? 0. : atof(Eq_hold[6].Data());
      EqW4[i] = TMath::IsNaN(atof(Eq_hold[7].Data())) || atof(Eq_hold[7].Data())<5.0 ? 0. : atof(Eq_hold[7].Data());


      //simFileIn >> EqE1[i] >> EqE2[i] >> EqE3[i] >> EqE4[i]
      //	>> EqW1[i] >> EqW2[i] >> EqW3[i] >> EqW4[i];
      //cout << EqE1[i] << " " << EqE2[i] << " " << EqE3[i] << " " << EqE4[i] << " "
      //   << EqW1[i] << " " << EqW2[i] << " " << EqW3[i] << " " << EqW4[i] << endl;
    }
    else {
      cout << "The simulation and data source peak files don't match. You need to re-make the source calibration files for both sim and data\n";
      exit(0);
    }

    if (filein.fail()) break;
    i++;
  }


  //Load the peak errors
  Int_t j = 0;
  run_hold = 0;
  sourceName_hold = "";
  string w2;

  while (fileinErr >> run[j] >> sourceName[j]
	 >> ADC_hold[0] >> ADC_hold[1] >> ADC_hold[2] >> ADC_hold[3]
	 >> ADC_hold[4] >> ADC_hold[5] >> ADC_hold[6] >> ADC_hold[7]) {

    adcE1_err[j] = TMath::IsNaN(atof(ADC_hold[0].Data())) ? 100. : atof(ADC_hold[0].Data());
    adcE2_err[j] = TMath::IsNaN(atof(ADC_hold[1].Data())) ? 100. : atof(ADC_hold[1].Data());
    adcE3_err[j] = TMath::IsNaN(atof(ADC_hold[2].Data())) ? 100. : atof(ADC_hold[2].Data());
    adcE4_err[j] = TMath::IsNaN(atof(ADC_hold[3].Data())) ? 100. : atof(ADC_hold[3].Data());
    adcW1_err[j] = TMath::IsNaN(atof(ADC_hold[4].Data())) ? 100. : atof(ADC_hold[4].Data());
    adcW2_err[j] = TMath::IsNaN(atof(ADC_hold[5].Data())) ? 100. : atof(ADC_hold[5].Data());
    adcW3_err[j] = TMath::IsNaN(atof(ADC_hold[6].Data())) ? 100. : atof(ADC_hold[6].Data());
    adcW4_err[j] = TMath::IsNaN(atof(ADC_hold[7].Data())) ? 100. : atof(ADC_hold[7].Data());
      
    while (!(run_hold==run[j] && sourceName_hold==sourceName[j])) {
      simFileInErr >> run_hold >> sourceName_hold >> Eq_hold[0] >> Eq_hold[1] >> Eq_hold[2]
		>> Eq_hold[3] >> Eq_hold[4] >> Eq_hold[5] >> Eq_hold[6] >> Eq_hold[7];
    }
    cout << run_hold << " " << sourceName_hold << " ";
    if (run_hold==run[j] && sourceName_hold==sourceName[j]) { 
      EqE1_err[j] = TMath::IsNaN(atof(Eq_hold[0].Data())) ? 100. : atof(Eq_hold[0].Data());
      EqE2_err[j] = TMath::IsNaN(atof(Eq_hold[1].Data())) ? 100. : atof(Eq_hold[1].Data());
      EqE3_err[j] = TMath::IsNaN(atof(Eq_hold[2].Data())) ? 100. : atof(Eq_hold[2].Data());
      EqE4_err[j] = TMath::IsNaN(atof(Eq_hold[3].Data())) ? 100. : atof(Eq_hold[3].Data());
      EqW1_err[j] = TMath::IsNaN(atof(Eq_hold[4].Data())) ? 100. : atof(Eq_hold[4].Data());
      EqW2_err[j] = TMath::IsNaN(atof(Eq_hold[5].Data())) ? 100. : atof(Eq_hold[5].Data());
      EqW3_err[j] = TMath::IsNaN(atof(Eq_hold[6].Data())) ? 100. : atof(Eq_hold[6].Data());
      EqW4_err[j] = TMath::IsNaN(atof(Eq_hold[7].Data())) ? 100. : atof(Eq_hold[7].Data());
      //simFileInErr >> EqE1_err[j] >> EqE2_err[j] >> EqE3_err[j] >> EqE4_err[j]
      //	>> EqW1_err[j] >> EqW2_err[j] >> EqW3_err[j] >> EqW4_err[j];
      //cout << EqE1_err[j] << " " << EqE2_err[j] << " " << EqE3_err[j] << " " << EqE4_err[j] << " "
      //   << EqW1_err[j] << " " << EqW2_err[j] << " " << EqW3_err[j] << " " << EqW4_err[j] << endl;
    }
    else {
      cout << "The simulation and data source peak files don't match. You need to re-make the source calibration files for both sim and data\n";
      exit(0);
    }

    if (fileinErr.fail()) break;
    j++;
  }
  
  std::cout << i << " " << j << std::endl;

  if (i!=j) {
    std::cout << "Source peaks and source errors don't match up!!\n";
    exit(0);
  }
  Int_t num = i;
  cout << "Number of data points: " << num << endl;

  

  //Read in PMT quality file to save whether or not to use a PMT for each run
  vector<vector<int> > pmtQuality;
  vector <int> pmthold(8,0);
  vector<Int_t> pmtRun;
  sprintf(temp,"%s/residuals/PMT_runQuality_SrcPeriod_%i.dat",getenv("ANALYSIS_CODE"),calibrationPeriod);
  ifstream pmt;
  pmt.open(temp);
  Int_t EPMT1,EPMT2,EPMT3,EPMT4,WPMT1,WPMT2,WPMT3,WPMT4;
  Int_t numRuns=0;
  while (pmt >> run_hold >> pmthold[0] >> pmthold[1] >> pmthold[2]
	 >> pmthold[3] >> pmthold[4] >> pmthold[5]
	 >> pmthold[6] >> pmthold[7]) {
    pmtRun.push_back(run_hold);
    pmtQuality.push_back(pmthold);

//>> EPMT1 >> EPMT2 >> EPMT3 >> EPMT4 >> WPMT1 >> WPMT2 >> WPMT3 >> WPMT4;
    //cout << run_hold << endl;
    
    numRuns++;
    if (pmt.fail()) break;
  }
  pmt.close();
  for (int i=0;i<pmtRun.size();i++) {
    cout << pmtRun[i] << " " << pmtQuality[i][0] << " " << pmtQuality[i][1] << " " << pmtQuality[i][2] << 
      " " << pmtQuality[i][3] << " " << pmtQuality[i][4] << " " << pmtQuality[i][5] << " " << pmtQuality[i][6] << 
      " " << pmtQuality[i][7] << endl;
  }

  
  // Filling new vectors for each PMT with ADC values only when the PMT is 
  // usable
  vector<int> runE1,runE2,runE3,runE4,runW1,runW2,runW3,runW4;
  vector<Double_t> EQE1,EQE2,EQE3,EQE4,EQW1,EQW2,EQW3,EQW4;
  vector<Double_t> EQE1_err,EQE2_err,EQE3_err,EQE4_err,EQW1_err,EQW2_err,EQW3_err,EQW4_err;
  vector<string> nameE1,nameE2,nameE3,nameE4,nameW1,nameW2,nameW3,nameW4;
  vector<Double_t> ADCE1, ADCE2, ADCE3, ADCE4;
  vector<Double_t> ADCW1, ADCW2, ADCW3, ADCW4;
  vector<Double_t> ADCE1_err, ADCE2_err, ADCE3_err, ADCE4_err;
  vector<Double_t> ADCW1_err, ADCW2_err, ADCW3_err, ADCW4_err;
  vector<Double_t> ResE1, ResE2, ResE3, ResE4;
  vector<Double_t> ResW1, ResW2, ResW3, ResW4;

  Double_t fiducialCut = 35.; //This is the cut for the center of the source position,
                              // but only for inclusion in the linearity curve

  // Fill run, adc, and eQ vectors
  for (Int_t i=0; i<num;i++) {
    UInt_t runPos = find_vec_location_int(pmtRun,run[i]);
    cout << "Found run " <<  run[i] << " in PMT Quality\n";
    /*int src_hold;
    if (sourceName[i]=="Ce") src_hold = 0;
    else if (sourceName[i]=="Sn") src_hold=1;
    else if (sourceName[i]=="Bi2") src_hold=2;
    else if (sourceName[i]=="Bi1") src_hold=3;*/

    //Sources to ignore in fit
    if (sourceName[i]=="Cd" || sourceName[i]=="In") continue; //|| sourceName[i]=="Bi2" 

    
    if ( isSourceInFidCut(run[i],sourceName[i].substr(0,2), fiducialCut, 0) ) {

      //if (pmtQuality[runPos][0]) {
      runE1.push_back(run[i]);
      EQE1.push_back(EqE1[i]);
      EQE1_err.push_back(EqE1_err[i]);
      nameE1.push_back(sourceName[i]);
      ADCE1.push_back(adcE1[i]);
      ADCE1_err.push_back(adcE1_err[i]);
      //}

      //if (pmtQuality[runPos][1]) {
      runE2.push_back(run[i]);
      EQE2.push_back(EqE2[i]);
      EQE2_err.push_back(EqE2_err[i]);
      nameE2.push_back(sourceName[i]);
      ADCE2.push_back(adcE2[i]);
      ADCE2_err.push_back(adcE2_err[i]);
      //}
      //if (pmtQuality[runPos][2]) {
      runE3.push_back(run[i]);
      EQE3.push_back(EqE3[i]);
      EQE3_err.push_back(EqE3_err[i]);
      nameE3.push_back(sourceName[i]);
      ADCE3.push_back(adcE3[i]);
      ADCE3_err.push_back(adcE3_err[i]);
      //}
      //if (pmtQuality[runPos][3]) {
      runE4.push_back(run[i]);
      EQE4.push_back(EqE4[i]);
      EQE4_err.push_back(EqE4_err[i]);
      nameE4.push_back(sourceName[i]);
      ADCE4.push_back(adcE4[i]);
      ADCE4_err.push_back(adcE4_err[i]);
      //}
    }

    if ( isSourceInFidCut(run[i],sourceName[i].substr(0,2), fiducialCut, 1) ) {

      //if (pmtQuality[runPos][4]) 
      if (run[i]!=17361 && run[i]!=17360) {
      runW1.push_back(run[i]);
      EQW1.push_back(EqW1[i]);
      EQW1_err.push_back(EqW1_err[i]);
      nameW1.push_back(sourceName[i]);
      ADCW1.push_back(adcW1[i]);
      ADCW1_err.push_back(adcW1_err[i]);
      }
      if (run[i]<16983 || run[i]>17249) { //pmtQuality[runPos][5]) {
	runW2.push_back(run[i]);
	EQW2.push_back(EqW2[i]);
	EQW2_err.push_back(EqW2_err[i]);
	nameW2.push_back(sourceName[i]);
	ADCW2.push_back(adcW2[i]);
	ADCW2_err.push_back(adcW2_err[i]);
      }
      //if (pmtQuality[runPos][6]) {
      runW3.push_back(run[i]);
      EQW3.push_back(EqW3[i]);
      EQW3_err.push_back(EqW3_err[i]);
      nameW3.push_back(sourceName[i]);
      ADCW3.push_back(adcW3[i]);
      ADCW3_err.push_back(adcW3_err[i]);
      //}
      //if (pmtQuality[runPos][7]) {
      runW4.push_back(run[i]);
      EQW4.push_back(EqW4[i]);
      EQW4_err.push_back(EqW4_err[i]);
      nameW4.push_back(sourceName[i]);
      ADCW4.push_back(adcW4[i]);
      ADCW4_err.push_back(adcW4_err[i]);
      //}
    }

  }


  //Resize res vectors
  ResE1.resize(runE1.size());
  ResE2.resize(runE2.size());
  ResE3.resize(runE3.size());
  ResE4.resize(runE4.size());
  ResW1.resize(runW1.size());
  ResW2.resize(runW2.size());
  ResW3.resize(runW3.size());
  ResW4.resize(runW4.size());
  
  //Checking that there are enough sources to make a linearity curve
  /*if (std::find(EQE1.begin(),EQE1.end(),peakCe)==EQE1.end() || std::find(EQE1.begin(),EQE1.end(),peakSn)==EQE1.end() || std::find(EQE1.begin(),EQE1.end(),peakBiHigh)==EQE1.end()) { 
    cout << "Not enough sources to construct quadratic linearity curve\n";
    num=0;
    }*/
  if (std::find(nameE1.begin(),nameE1.end(),"Ce")==nameE1.end() || std::find(nameE1.begin(),nameE1.end(),"Sn")==nameE1.end() || (std::find(nameE1.begin(),nameE1.end(),"Bi1")==nameE1.end() && std::find(nameE1.begin(),nameE1.end(),"Bi2")==nameE1.end())) { 
    cout << "Not enough sources to construct quadratic linearity curve\n";
    num=0;
  }
    
  //Finding the lowest ADC/EQ peak in each PMT for fitting purposes
  std::vector <Double_t> lowestEQ(8,10000.);
  std::vector <Double_t> lowestADC(8,10000.);
  for (Int_t i=0; i<ADCE1.size(); i++) { if (ADCE1[i] < lowestADC[0] && ADCE1[i]>0. ) { lowestADC[0] = ADCE1[i]; lowestEQ[0] = EQE1[i]; } }
  for (Int_t i=0; i<ADCE2.size(); i++) { if (ADCE2[i] < lowestADC[1] && ADCE2[i]>0. ) { lowestADC[1] = ADCE2[i]; lowestEQ[1] = EQE2[i]; } }
  for (Int_t i=0; i<ADCE3.size(); i++) { if (ADCE3[i] < lowestADC[2] && ADCE3[i]>0. ) { lowestADC[2] = ADCE3[i]; lowestEQ[2] = EQE3[i]; } }
  for (Int_t i=0; i<ADCE4.size(); i++) { if (ADCE4[i] < lowestADC[3] && ADCE4[i]>0. ) { lowestADC[3] = ADCE4[i]; lowestEQ[3] = EQE4[i]; } }
  for (Int_t i=0; i<ADCW1.size(); i++) { if (ADCW1[i] < lowestADC[4] && ADCW1[i]>0. ) { lowestADC[4] = ADCW1[i]; lowestEQ[4] = EQW1[i]; } }
  for (Int_t i=0; i<ADCW2.size(); i++) { if (ADCW2[i] < lowestADC[5] && ADCW2[i]>0. ) { lowestADC[5] = ADCW2[i]; lowestEQ[5] = EQW2[i]; } }
  for (Int_t i=0; i<ADCW3.size(); i++) { if (ADCW3[i] < lowestADC[6] && ADCW3[i]>0. ) { lowestADC[6] = ADCW3[i]; lowestEQ[6] = EQW3[i]; } }
  for (Int_t i=0; i<ADCW4.size(); i++) { if (ADCW4[i] < lowestADC[7] && ADCW4[i]>0. ) { lowestADC[7] = ADCW4[i]; lowestEQ[7] = EQW4[i]; } }
  

  //File to hold the linearity curves for each calibration run period
  sprintf(temp,"%s/lin_curves_srcCal_Period_%i.dat",getenv("LINEARITY_CURVES"),calibrationPeriod);
  ofstream linCurves(temp);

  // Fit function
  TF1 *fitADC = new TF1("fitADC", "([0] + [1]*x + [2]*x*x)", 0., 2500.0);
  fitADC->SetParameter(0, 0.0);
  fitADC->SetParLimits(0, -20., 20.);
  //fitADC->FixParameter(0, 0.0);
  fitADC->SetParameter(1, 1.0);
  fitADC->SetParameter(2, 0.0);
  if ( calibrationPeriod<13) fitADC->SetParLimits(2, -0.00005, 0.00005);
  else fitADC->SetParLimits(2, -0.000075, 0.000075);
  if (!quadratic) fitADC->FixParameter(2, 0.0);
  //fitADC->FixParameter(0, 0.0);

  /*TF1 *fitADC = new TF1("fitADC", "([0] + [1]*x + [2]*x*x)*(0.5+0.5*TMath::TanH((x-[4])/[5]))+([3]*x)*(0.5-0.5*TMath::TanH((x-[4])/[5]))", 0., 2500.0);
  fitADC->SetParameter(0, 0.0);
  fitADC->SetParameter(1, 1.0);
  fitADC->SetParameter(2, 0.0);
  fitADC->FixParameter(0, 0.0);

  fitADC->SetParameter(3, 1.0);
  fitADC->SetParameter(4, 50.);
  fitADC->SetParameter(5, 1.);
  //fitADC->FixParameter(5, 75.0);
  //fitADC->FixParameter(2, 0.);
  //fitADC->SetParLimits(5, 1., 20.);
  
  
  //fitADC->SetParLimits(3, 0.1, 2.);
  */

  fitADC->SetNpx(100000);
  fitADC->SetLineColor(2);
  fitADC->SetLineWidth(1);
  //fitADC->SetLineStyle(10);

  Double_t offset=0., slope=0., quad=0., cubic=0., spreadOfTanH=0., meanOfTanH=0., lowEnSlope=0.; //Fit Parameters
  Double_t shiftOfTanH=0.;

  Double_t x1_text, y1_text;

  //Creating canvases and TPads for all of the calibrations curves
  TCanvas *c1 = new TCanvas("c1", "c1", 1400, 1000);

  TPad *E1 = new TPad("E1","East 1", 0.0, 0.65, 0.25, 1.0);
  TPad *E1res = new TPad("E1res","East 1 residuals", 0.0, 0.5, 0.25, 0.65);
  TPad *E2 = new TPad("E2","East 2", 0.25, 0.65, 0.5, 1.0);
  TPad *E2res = new TPad("E2res","East 2 residuals", 0.25, 0.5, 0.5, 0.65);
  TPad *E3 = new TPad("E3","East 3", 0.5, 0.65, 0.75, 1.0);
  TPad *E3res = new TPad("E3res","East 3 residuals", 0.5, 0.5, 0.75, 0.65);
  TPad *E4 = new TPad("E1","East 1", 0.75, 0.65, 1.0, 1.0);
  TPad *E4res = new TPad("E4res","East 4 residuals", 0.75, 0.5, 1.0, 0.65);
  E1->Draw();
  E1res->Draw();
  E2->Draw();
  E2res->Draw();
  E3->Draw();
  E3res->Draw();
  E4->Draw();
  E4res->Draw();

  TPad *W1 = new TPad("W1","West 1", 0.0, 0.15, 0.25, 0.5);
  TPad *W1res = new TPad("W1res","West 1 residuals", 0.0, 0.0, 0.25, 0.15);
  TPad *W2 = new TPad("W2","West 2", 0.25, 0.15, 0.5, 0.5);
  TPad *W2res = new TPad("W2res","West 2 residuals", 0.25, 0.0, 0.5, 0.15);
  TPad *W3 = new TPad("W3","West 3", 0.5, 0.15, 0.75, 0.5);
  TPad *W3res = new TPad("W3res","West 3 residuals", 0.5, 0.0, 0.75, 0.15);
  TPad *W4 = new TPad("W1","West 1", 0.75, 0.15, 1.0, 0.5);
  TPad *W4res = new TPad("W4res","West 4 residuals", 0.75, 0.0, 1.0, 0.15);
  W1->Draw();
  W1res->Draw();
  W2->Draw();
  W2res->Draw();
  W3->Draw();
  W3res->Draw();
  W4->Draw();
  W4res->Draw();

  // East 1

  sprintf(temp,"%s/residuals/residuals_East_runPeriod_%i_PMTE1.dat",getenv("ANALYSIS_CODE"),calibrationPeriod);
  ofstream oFileE1(temp);
  vector <Double_t> fitEQ_E1(runE1.size(),0);

  

  if (runE1.size()>0 && std::find(nameE1.begin(),nameE1.end(),"Ce")!=nameE1.end() && std::find(nameE1.begin(),nameE1.end(),"Sn")!=nameE1.end() && std::find(nameE1.begin(),nameE1.end(),"Bi1")!=nameE1.end()  )  {
 
    //Calculating the average of the lower Ce ADC values and assigning a value to lowFitThreshold
    Double_t maxADC = 0., maxEQ = 0.;
    for (UInt_t ii=0; ii<runE1.size(); ii++) {
      if (ADCE1[ii]>maxADC) maxADC = ADCE1[ii];
      if (EQE1[ii]>maxEQ) maxEQ = EQE1[ii];
    }

    E1->cd();

    TGraphErrors *grE1 = new TGraphErrors(runE1.size(),&ADCE1[0],&EQE1[0],&ADCE1_err[0],&EQE1_err[0]);
    grE1->SetTitle("East PMT 1");
    grE1->SetMarkerColor(1);
    grE1->SetLineColor(1);
    grE1->SetMarkerStyle(20);
    grE1->SetMarkerSize(0.75);
    grE1->GetXaxis()->SetTitle("g(t)#upointADC [channels]");
    grE1->GetXaxis()->SetTitleOffset(1.2);
    grE1->GetXaxis()->CenterTitle();
    grE1->GetYaxis()->SetTitle("#eta(x,y)#upointE_{Q} [keV]");
    grE1->GetYaxis()->SetTitleOffset(1.6);
    grE1->GetYaxis()->CenterTitle();
    grE1->GetXaxis()->SetLimits(0.0,maxADC+40.);
    grE1->SetMinimum(0.0);
    grE1->SetMaximum(maxEQ+40.);
    grE1->Draw("AP");


    grE1->Fit("fitADC", "","", 0., maxADC+40.);

    offset = fitADC->GetParameter(0);
    slope = fitADC->GetParameter(1);
    quad = fitADC->GetParameter(2);

    shiftOfTanH = offset>=0. ? lowestADC[0]/2. : (lowestADC[0]-offset/slope)/2.;
    lowEnSlope = (quad*shiftOfTanH*shiftOfTanH + slope*shiftOfTanH + offset)/shiftOfTanH;
    spreadOfTanH = shiftOfTanH/3.;

    if (useTanh) linCurves << offset << " " << slope << " " << quad  << " " << lowEnSlope << " " << shiftOfTanH << " " << spreadOfTanH << endl;
    else linCurves << offset << " " << slope << " " << quad << endl;
  
    x1_text = 1200;
    y1_text = 100;

    /*TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("East PMT 1");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();*/

    // Calculate residuals in [keV]
    
    for (int j=0; j<runE1.size(); j++) {
    
      fitEQ_E1[j]    = fitADC->Eval(ADCE1[j]);//offset + slope*ADCE1[j] + quad*ADCE1[j]*ADCE1[j] + cubic*ADCE1[j]*ADCE1[j]*ADCE1[j];
      if (nameE1[j]=="Ce") {
	ResE1[j] = 100.*(fitEQ_E1[j] - EQE1[j])/fitEQ_E1[j];
	oFileE1 << "Ce_East" << " " << runE1[j] << " " << ResE1[j] << endl;
	if (sqrt(ResE1[j]*ResE1[j])>2.5) cout << "Ce_East" << " PMT1 " << runE1[j] << " " << ResE1[j] << endl;
	//ResE1[j] = (fitEQ - peakCe)/peakCe * 100.;
      }
      else if (nameE1[j]=="Sn") {
	ResE1[j] = 100.*(fitEQ_E1[j] - EQE1[j])/fitEQ_E1[j];
	//ResE1[j] = (fitEQ - peakSn)/peakSn * 100.;
	oFileE1 << "Sn_East" << " " << runE1[j] << " " << ResE1[j] << endl;
	if (sqrt(ResE1[j]*ResE1[j])>2.5) cout << "Sn_East" << " " << "PMT1 " << runE1[j] << " " << ResE1[j] << endl;
      }
      else if (nameE1[j]=="Bi1") {
	ResE1[j] = 100.*(fitEQ_E1[j] - EQE1[j])/fitEQ_E1[j];
	//ResE1[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE1 << "Bi1_East" << " " << runE1[j] << " " << ResE1[j] << endl;
	if (sqrt(ResE1[j]*ResE1[j])>2.5) cout << "Bi1_East" << " " << "PMT1 " << runE1[j] << " " << ResE1[j] << endl;
      }
      else if (nameE1[j]=="Bi2") {
	ResE1[j] = 100.*(fitEQ_E1[j] - EQE1[j])/fitEQ_E1[j];
	//ResE1[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE1 << "Bi2_East" << " " << runE1[j] << " " << ResE1[j] << endl;
	if (sqrt(ResE1[j]*ResE1[j])>2.5) cout << "Bi2_East" << " " << "PMT1 " << runE1[j] << " " << ResE1[j] << endl;	
      }
      else if (nameE1[j]=="In") {
	ResE1[j] = 100.*(fitEQ_E1[j] - EQE1[j])/fitEQ_E1[j];
	//ResE1[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE1 << "In_East" << " " << runE1[j] << " " << ResE1[j] << endl;
	if (sqrt(ResE1[j]*ResE1[j])>2.5) cout << "In_East" << " " << "PMT1 " << runE1[j] << " " << ResE1[j] << endl;	
      }
    }
  

    E1res->cd();

    TGraphErrors *grE1r = new TGraphErrors(runE1.size(),&ADCE1[0],&ResE1[0],err,err);
    grE1r->SetTitle("");
    grE1r->SetMarkerColor(1);
    grE1r->SetLineColor(1);
    grE1r->SetMarkerStyle(20);
    grE1r->SetMarkerSize(0.75);
    //grE1r->GetXaxis()->SetTitle("g(t)#upointADC [channels]");
    grE1r->GetXaxis()->SetTitleOffset(1.2);
    grE1r->GetXaxis()->CenterTitle();
    grE1r->GetYaxis()->SetTitle("% Residuals");
    grE1r->GetYaxis()->SetTitleSize(0.1);
    grE1r->GetYaxis()->SetLabelSize(0.07);
    grE1r->GetXaxis()->SetLabelSize(0);
    grE1r->GetYaxis()->SetTitleOffset(0.6);
    grE1r->GetYaxis()->CenterTitle();
    grE1r->GetXaxis()->SetLimits(0.0,maxADC+40.);
    grE1r->SetMinimum(-10.0);
    grE1r->SetMaximum( 10.0);
    grE1r->Draw("AP");

    const Int_t n = 2;
    Double_t x[n] = {0, maxADC+40.};
    Double_t y[n] = {0.0, 0.0};

    TGraph *gr1 = new TGraph(n,x,y);
    gr1->Draw("Same");
    gr1->SetLineWidth(2);
    gr1->SetLineColor(1);
    gr1->SetLineStyle(2);

    x1_text = 1200;
    y1_text = -40;

    /*TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("East PMT 1");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();*/
  }

  else { 
    if (useTanh) linCurves << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << endl;
    else linCurves << 0. << " " << 0. << " " << 0. << endl;
    oFileE1 << "PMT NOT USABLE";}
  oFileE1.close();

  // East 2

  sprintf(temp,"../residuals/residuals_East_runPeriod_%i_PMTE2.dat",calibrationPeriod);
  ofstream oFileE2(temp);
  vector <Double_t> fitEQ_E2(runE2.size(),0);


  if (runE2.size()>0 && std::find(nameE2.begin(),nameE2.end(),"Ce")!=nameE2.end() && std::find(nameE2.begin(),nameE2.end(),"Sn")!=nameE2.end() && std::find(nameE2.begin(),nameE2.end(),"Bi1")!=nameE2.end()  )  {

    //Calculating the average of the lower Ce ADC values and assigning a value to lowFitThreshold
    Double_t maxADC = 0., maxEQ=0.;
    for (UInt_t ii=0; ii<runE2.size(); ii++) {
      if (ADCE2[ii]>maxADC) maxADC = ADCE2[ii];
      if (EQE2[ii]>maxEQ) maxEQ = EQE2[ii];
    }
    

    E2->cd();

    TGraphErrors *grE2 = new TGraphErrors(runE2.size(),&ADCE2[0],&EQE2[0],&ADCE2_err[0],&EQE2_err[0]);
    grE2->SetTitle("East PMT 2");
    grE2->SetMarkerColor(1);
    grE2->SetLineColor(1);
    grE2->SetMarkerStyle(20);
    grE2->SetMarkerSize(0.75);
    grE2->GetXaxis()->SetTitle("g(t)#upointADC [channels]");
    grE2->GetXaxis()->SetTitleOffset(1.2);
    grE2->GetXaxis()->CenterTitle();
    grE2->GetYaxis()->SetTitle("#eta(x,y)#upointE_{Q} [keV]");
    grE2->GetYaxis()->SetTitleOffset(1.6);
    grE2->GetYaxis()->CenterTitle();
    grE2->GetXaxis()->SetLimits(0.0,maxADC+40.);
    grE2->SetMinimum(0.0);
    grE2->SetMaximum(maxEQ+40.);
    grE2->Draw("AP");

    grE2->Fit("fitADC", "","", 0., maxADC+40.);

    offset = fitADC->GetParameter(0);
    slope = fitADC->GetParameter(1);
    quad = fitADC->GetParameter(2);

    shiftOfTanH = offset>=0. ? lowestADC[1]/2. : (lowestADC[1]-offset/slope)/2.;
    lowEnSlope = (quad*shiftOfTanH*shiftOfTanH + slope*shiftOfTanH + offset)/shiftOfTanH;
    spreadOfTanH = shiftOfTanH/3.;

  
    if (useTanh) linCurves << offset << " " << slope << " " << quad  << " " << lowEnSlope << " " << shiftOfTanH << " " << spreadOfTanH << endl;
    else linCurves << offset << " " << slope << " " << quad << endl;

    x1_text = 1200;
    y1_text = 100;

    /*TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("East PMT 2");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();*/

    // Calculate residuals in [keV]
    
    for (int j=0; j<runE2.size(); j++) {
      fitEQ_E2[j]    = fitADC->Eval(ADCE2[j]);//offset + slope*ADCE2[j] + quad*ADCE2[j]*ADCE2[j] + cubic*ADCE2[j]*ADCE2[j]*ADCE2[j];
      if (nameE2[j]=="Ce") {
	ResE2[j] = 100.*(fitEQ_E2[j] - EQE2[j])/fitEQ_E2[j];
	//ResE2[j] = (fitEQ - peakCe)/peakCe * 100.;
	oFileE2 << "Ce_East" << " " << runE2[j] << " " << ResE2[j] << endl;
	if (sqrt(ResE2[j]*ResE2[j])>2.5) cout<< "Ce_East" << " " << "PMT2 " << runE2[j] << " " << ResE2[j] << endl;  
      }
      else if (nameE2[j]=="Sn") {
	ResE2[j] = 100.*(fitEQ_E2[j] - EQE2[j])/fitEQ_E2[j];
	//ResE2[j] = (fitEQ - peakSn)/peakSn * 100.;
	oFileE2 << "Sn_East" << " " << runE2[j] << " " << ResE2[j] << endl;
	if (sqrt(ResE2[j]*ResE2[j])>2.5) cout<< "Sn_East" << " " << "PMT2 " << runE2[j] << " " << ResE2[j] << endl;  
      }
      else if (nameE2[j]=="Bi1") {
	ResE2[j] = 100.*(fitEQ_E2[j] - EQE2[j])/fitEQ_E2[j];
	//ResE2[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE2 << "Bi1_East" << " " << runE2[j] << " " << ResE2[j] << endl;
	if (sqrt(ResE2[j]*ResE2[j])>2.5) cout<< "Bi1_East" << " " << "PMT2 " << runE2[j] << " " << ResE2[j] << endl;  
      }
      else if (nameE2[j]=="Bi2") {
	ResE2[j] = 100.*(fitEQ_E2[j] - EQE2[j])/fitEQ_E2[j];
	//ResE2[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE2 << "Bi2_East" << " " << runE2[j] << " " << ResE2[j] << endl;
	if (sqrt(ResE2[j]*ResE2[j])>2.5) cout<< "Bi2_East" << " " << "PMT2 " << runE2[j] << " " << ResE2[j] << endl;  
      }
      else if (nameE2[j]=="In") {
	ResE2[j] = 100.*(fitEQ_E2[j] - EQE2[j])/fitEQ_E2[j];
	//ResE2[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE2 << "In_East" << " " << runE2[j] << " " << ResE2[j] << endl;
	if (sqrt(ResE2[j]*ResE2[j])>2.5) cout<< "In_East" << " " << "PMT2 " << runE2[j] << " " << ResE2[j] << endl;  
      }
    }

    // East 2 residuals
    E2res->cd();

    TGraphErrors *grE2r = new TGraphErrors(runE2.size(),&ADCE2[0],&ResE2[0],err,err);
    grE2r->SetTitle("");
    grE2r->SetMarkerColor(1);
    grE2r->SetLineColor(1);
    grE2r->SetMarkerStyle(20);
    grE2r->SetMarkerSize(0.75);
    //grE2r->GetXaxis()->SetTitle("g(t)#upointADC [channels]");
    grE2r->GetXaxis()->SetTitleOffset(1.2);
    grE2r->GetXaxis()->CenterTitle();
    grE2r->GetYaxis()->SetTitle("% Residuals");
    grE2r->GetYaxis()->SetTitleSize(0.1);
    grE2r->GetYaxis()->SetLabelSize(0.07);
    grE2r->GetXaxis()->SetLabelSize(0);
    grE2r->GetYaxis()->SetTitleOffset(0.6);
    grE2r->GetYaxis()->CenterTitle();
    grE2r->GetXaxis()->SetLimits(0.0,maxADC+40.);
    grE2r->SetMinimum(-10.0);
    grE2r->SetMaximum( 10.0);
    grE2r->Draw("AP");

    const Int_t n = 2;
    Double_t x[n] = {0, maxADC+40.};
    Double_t y[n] = {0.0, 0.0};

    TGraph *gr1 = new TGraph(n,x,y);
    gr1->Draw("Same");
    gr1->SetLineWidth(2);
    gr1->SetLineColor(1);
    gr1->SetLineStyle(2);

    x1_text = 1200;
    y1_text = -40;

    /*TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("East PMT 2");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();*/
  }

  else {
    if (useTanh) linCurves << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << endl;
    else linCurves << 0. << " " << 0. << " " << 0. << endl;
    oFileE2 << "PMT NOT USABLE";}
  oFileE2.close();

  // East 3

  sprintf(temp,"../residuals/residuals_East_runPeriod_%i_PMTE3.dat",calibrationPeriod);
  ofstream oFileE3(temp);
  vector <Double_t> fitEQ_E3(runE3.size(),0);

  if (runE3.size()>0 && std::find(nameE3.begin(),nameE3.end(),"Ce")!=nameE3.end() && std::find(nameE3.begin(),nameE3.end(),"Sn")!=nameE3.end() && std::find(nameE3.begin(),nameE3.end(),"Bi1")!=nameE3.end()  )  {

    //Calculating the average of the lower Ce ADC values and assigning a value to lowFitThreshold
    Double_t maxADC = 0., maxEQ = 0.;
    for (UInt_t ii=0; ii<runE3.size(); ii++) {
      if (ADCE3[ii]>maxADC) maxADC = ADCE3[ii];
      if (EQE3[ii]>maxEQ) maxEQ = EQE3[ii];
    }

    E3->cd();

    TGraphErrors *grE3 = new TGraphErrors(runE3.size(),&ADCE3[0],&EQE3[0],&ADCE3_err[0],&EQE3_err[0]);
    grE3->SetTitle("East PMT 3");
    grE3->SetMarkerColor(1);
    grE3->SetLineColor(1);
    grE3->SetMarkerStyle(20);
    grE3->SetMarkerSize(0.75);
    grE3->GetXaxis()->SetTitle("g(t)#upointADC [channels]");
    grE3->GetXaxis()->SetTitleOffset(1.2);
    grE3->GetXaxis()->CenterTitle();
    grE3->GetYaxis()->SetTitle("#eta(x,y)#upointE_{Q} [keV]");
    grE3->GetYaxis()->SetTitleOffset(1.6);
    grE3->GetYaxis()->CenterTitle();
    grE3->GetXaxis()->SetLimits(0.0,maxADC+40.);
    grE3->SetMinimum(0.0);
    grE3->SetMaximum(maxEQ+40.);
    grE3->Draw("AP");

    //fitADC->SetParLimits(0, 0., 20.);

    grE3->Fit("fitADC", "","", 0., maxADC+40.);

    //fitADC->SetParLimits(0, -20., 20.);

    offset = fitADC->GetParameter(0);
    slope = fitADC->GetParameter(1);
    quad = fitADC->GetParameter(2);

    shiftOfTanH = offset>=0. ? lowestADC[2]/2. : (lowestADC[2]-offset/slope)/2.;
    lowEnSlope = (quad*shiftOfTanH*shiftOfTanH + slope*shiftOfTanH + offset)/shiftOfTanH;
    spreadOfTanH = shiftOfTanH/3.;

    if (useTanh) linCurves << offset << " " << slope << " " << quad  << " " << lowEnSlope << " " << shiftOfTanH << " " << spreadOfTanH << endl;
    else linCurves << offset << " " << slope << " " << quad << endl;

    /*grE3->Fit("fitADC", "","",0.,maxADC+40.);

    offset = fitADC->GetParameter(0);
    slope = fitADC->GetParameter(1);
    quad = fitADC->GetParameter(2);
    
    linCurves << offset << " " << slope << " " << quad << << " " << lowEnSlope << " " << shiftOfTanH << " " << spreadOfTanH << endl;
    */

    x1_text = 1200;
    y1_text = 100;

    /*TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("East PMT 3");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();*/

    // Calculate residuals in [keV]
    
    for (int j=0; j<runE3.size(); j++) {
      fitEQ_E3[j]    = fitADC->Eval(ADCE3[j]);//offset + slope*ADCE3[j] + quad*ADCE3[j]*ADCE3[j] + cubic*ADCE3[j]*ADCE3[j]*ADCE3[j];
      if (nameE3[j]=="Ce") {
	ResE3[j] = 100.*(fitEQ_E3[j] - EQE3[j])/fitEQ_E3[j];
	//ResE3[j] = (fitEQ - peakCe)/peakCe * 100.;
	oFileE3 << "Ce_East" << " " << runE3[j] << " " << ResE3[j] << endl;
	if (sqrt(ResE3[j]*ResE3[j])>2.5) cout<< "Ce_East" << " " << "PMT3 " << runE3[j] << " " << ResE3[j] << endl;  
      }
      else if (nameE3[j]=="Sn") {
	ResE3[j] =100.*(fitEQ_E3[j] - EQE3[j])/fitEQ_E3[j];
	//ResE3[j] = (fitEQ - peakSn)/peakSn * 100.;
	oFileE3 << "Sn_East" << " " << runE3[j] << " " << ResE3[j] << endl;
	if (sqrt(ResE3[j]*ResE3[j])>2.5) cout << "Sn_East" << " " << "PMT3 " << runE3[j] << " " << ResE3[j] << endl;  
      }
      else if (nameE3[j]=="Bi1") {
	ResE3[j] = 100.*(fitEQ_E3[j] - EQE3[j])/fitEQ_E3[j];
	//ResE3[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE3 << "Bi1_East" << " " << runE3[j] << " " << ResE3[j] << endl;
	if (sqrt(ResE3[j]*ResE3[j])>2.5) cout << "Bi1_East" << " " << "PMT3 " << runE3[j] << " " << ResE3[j] << endl;  
      }
      else if (nameE3[j]=="Bi2") {
	ResE3[j] = 100.*(fitEQ_E3[j] - EQE3[j])/fitEQ_E3[j];
	//ResE3[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE3 << "Bi2_East" << " " << runE3[j] << " " << ResE3[j] << endl;
	if (sqrt(ResE3[j]*ResE3[j])>2.5) cout << "Bi2_East" << " " << "PMT3 " << runE3[j] << " " << ResE3[j] << endl;  
      }
      else if (nameE3[j]=="In") {
	ResE3[j] = 100.*(fitEQ_E3[j] - EQE3[j])/fitEQ_E3[j];
	//ResE3[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE3 << "In_East" << " " << runE3[j] << " " << ResE3[j] << endl;
	if (sqrt(ResE3[j]*ResE3[j])>2.5) cout << "In_East" << " " << "PMT3 " << runE3[j] << " " << ResE3[j] << endl;  
      }
    }

    // East 3 residuals
    E3res->cd();

    TGraphErrors *grE3r = new TGraphErrors(runE3.size(),&ADCE3[0],&ResE3[0],err,err);
    grE3r->SetTitle("");
    grE3r->SetMarkerColor(1);
    grE3r->SetLineColor(1);
    grE3r->SetMarkerStyle(20);
    grE3r->SetMarkerSize(0.75);
    //grE3r->GetXaxis()->SetTitle("g(t)#upointADC [channels]");
    grE3r->GetXaxis()->SetTitleOffset(1.2);
    grE3r->GetXaxis()->CenterTitle();
    grE3r->GetYaxis()->SetTitle("% Residuals");
    grE3r->GetYaxis()->SetTitleSize(0.1);
    grE3r->GetYaxis()->SetLabelSize(0.07);
    grE3r->GetXaxis()->SetLabelSize(0);
    grE3r->GetYaxis()->SetTitleOffset(0.6);
    grE3r->GetYaxis()->CenterTitle();
    grE3r->GetXaxis()->SetLimits(0.0,maxADC+40.);
    grE3r->SetMinimum(-10.0);
    grE3r->SetMaximum( 10.0);
    grE3r->Draw("AP");

    static const Int_t n = 2;
    Double_t x[n] = {0, maxADC+40.};
    Double_t y[n] = {0.0, 0.0};

    TGraph *gr1 = new TGraph(n,x,y);
    gr1->Draw("Same");
    gr1->SetLineWidth(2);
    gr1->SetLineColor(1);
    gr1->SetLineStyle(2);

    x1_text = 1200;
    y1_text = -40;

    /*TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("East PMT 3");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();*/
  }

  else { 
    if (useTanh) linCurves << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << endl;
    else linCurves << 0. << " " << 0. << " " << 0. << endl;
    oFileE3 << "PMT NOT USABLE";}
  oFileE3.close();

  // East 4

  sprintf(temp,"%s/residuals/residuals_East_runPeriod_%i_PMTE4.dat",getenv("ANALYSIS_CODE"), calibrationPeriod);
  ofstream oFileE4(temp);
  vector <Double_t> fitEQ_E4(runE4.size(),0);

  if (runE4.size()>0 && std::find(nameE4.begin(),nameE4.end(),"Ce")!=nameE4.end() && std::find(nameE4.begin(),nameE4.end(),"Sn")!=nameE4.end() && std::find(nameE4.begin(),nameE4.end(),"Bi1")!=nameE4.end()  )  {
    //Calculating the average of the lower Ce ADC values and assigning a value to lowFitThreshold
    Double_t maxADC = 0., maxEQ = 0.;
    for (UInt_t ii=0; ii<runE4.size(); ii++) {
      if (ADCE4[ii]>maxADC) maxADC = ADCE4[ii];
      if (EQE4[ii]>maxEQ) maxEQ = EQE4[ii];
    }

    E4->cd();

    TGraphErrors *grE4 = new TGraphErrors(runE4.size(),&ADCE4[0],&EQE4[0],&ADCE4_err[0],&EQE4_err[0]);
    grE4->SetTitle("East PMT 4");
    grE4->SetMarkerColor(1);
    grE4->SetLineColor(1);
    grE4->SetMarkerStyle(20);
    grE4->SetMarkerSize(0.75);
    grE4->GetXaxis()->SetTitle("g(t)#upointADC [channels]");
    grE4->GetXaxis()->SetTitleOffset(1.2);
    grE4->GetXaxis()->CenterTitle();
    grE4->GetYaxis()->SetTitle("#eta(x,y)#upointE_{Q} [keV]");
    grE4->GetYaxis()->SetTitleOffset(1.6);
    grE4->GetYaxis()->CenterTitle();
    grE4->GetXaxis()->SetLimits(0.0,maxADC+40.);
    grE4->SetMinimum(0.0);
    grE4->SetMaximum(maxEQ+40.);
    grE4->Draw("AP");

    grE4->Fit("fitADC", "","", 0., maxADC+40.);

    offset = fitADC->GetParameter(0);
    slope = fitADC->GetParameter(1);
    quad = fitADC->GetParameter(2);
   
    shiftOfTanH = offset>=0. ? lowestADC[3]/2. : (lowestADC[3]-offset/slope)/2.;
    lowEnSlope = (quad*shiftOfTanH*shiftOfTanH + slope*shiftOfTanH + offset)/shiftOfTanH;
    spreadOfTanH = shiftOfTanH/3.;

    if (useTanh) linCurves << offset << " " << slope << " " << quad  << " " << lowEnSlope << " " << shiftOfTanH << " " << spreadOfTanH << endl;
    else linCurves << offset << " " << slope << " " << quad << endl;

    /*grE4->Fit("fitADC", "","",0.,maxADC+40.);

    offset = fitADC->GetParameter(0);
    slope = fitADC->GetParameter(1);
    quad = fitADC->GetParameter(2);
    
    linCurves << offset << " " << slope << " " << quad << " " << lowEnSlope << " " << shiftOfTanH << " " << spreadOfTanH << endl;
    */

    x1_text = 1200;
    y1_text = 100;

    /*TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("East PMT 4");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();*/

    // Calculate residuals in [keV]
    
    for (int j=0; j<runE4.size(); j++) {
      fitEQ_E4[j]    = fitADC->Eval(ADCE4[j]);//offset + slope*ADCE4[j] + quad*ADCE4[j]*ADCE4[j] + cubic*ADCE4[j]*ADCE4[j]*ADCE4[j];
      if (nameE4[j]=="Ce") {
	ResE4[j] = 100.*(fitEQ_E4[j] - EQE4[j]) / fitEQ_E4[j];
	//ResE4[j] = (fitEQ - peakCe)/peakCe * 100.;
	oFileE4 << "Ce_East" << " " << runE4[j] << " " << ResE4[j] << endl;
	if (sqrt(ResE4[j]*ResE4[j])>2.5) cout << "Ce_East" << " " << "PMT4 " << runE4[j] << " " << ResE4[j] << endl;  
      }
      else if (nameE4[j]=="Sn") {
	ResE4[j] = 100.*(fitEQ_E4[j] - EQE4[j]) / fitEQ_E4[j];
	//ResE4[j] = (fitEQ - peakSn)/peakSn * 100.;
	oFileE4 << "Sn_East" << " " << runE4[j] << " " << ResE4[j] << endl;
	if (sqrt(ResE4[j]*ResE4[j])>2.5) cout<< "Sn_East" << " " << "PMT4 " << runE4[j] << " " << ResE4[j] << endl;  
      }
      else if (nameE4[j]=="Bi1") {
	ResE4[j] = 100.*(fitEQ_E4[j] - EQE4[j]) / fitEQ_E4[j];
	//ResE4[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE4 << "Bi1_East" << " " << runE4[j] << " " << ResE4[j] << endl;
	if (sqrt(ResE4[j]*ResE4[j])>2.5) cout<< "Bi1_East" << " " << "PMT4 " << runE4[j] << " " << ResE4[j] << endl;  
      }
      else if (nameE4[j]=="Bi2") {
	ResE4[j] = 100.*(fitEQ_E4[j] - EQE4[j]) / fitEQ_E4[j];
	//ResE4[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE4 << "Bi2_East" << " " << runE4[j] << " " << ResE4[j] << endl;
	if (sqrt(ResE4[j]*ResE4[j])>2.5) cout<< "Bi2_East" << " " << "PMT4 " << runE4[j] << " " << ResE4[j] << endl;  
      }
      else if (nameE4[j]=="In") {
	ResE4[j] = 100.*(fitEQ_E4[j] - EQE4[j]) / fitEQ_E4[j];
	//ResE4[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE4 << "In_East" << " " << runE4[j] << " " << ResE4[j] << endl;
	if (sqrt(ResE4[j]*ResE4[j])>2.5) cout<< "In_East" << " " << "PMT4 " << runE4[j] << " " << ResE4[j] << endl;  
      }
    }

    // East 4 residuals
    E4res->cd();

    TGraphErrors *grE4r = new TGraphErrors(runE4.size(),&ADCE4[0],&ResE4[0],err,err);
    grE4r->SetTitle("");
    grE4r->SetMarkerColor(1);
    grE4r->SetLineColor(1);
    grE4r->SetMarkerStyle(20);
    grE4r->SetMarkerSize(0.75);
    //grE4r->GetXaxis()->SetTitle("g(t)#upointADC [channels]");
    grE4r->GetXaxis()->SetTitleOffset(1.2);
    grE4r->GetXaxis()->CenterTitle();
    grE4r->GetYaxis()->SetTitle("% Residuals");
    grE4r->GetYaxis()->SetTitleSize(0.1);
    grE4r->GetYaxis()->SetLabelSize(0.07);
    grE4r->GetXaxis()->SetLabelSize(0);
    grE4r->GetYaxis()->SetTitleOffset(0.6);
    grE4r->GetYaxis()->CenterTitle();
    grE4r->GetXaxis()->SetLimits(0.0,maxADC+40.);
    grE4r->SetMinimum(-10.0);
    grE4r->SetMaximum( 10.0);
    grE4r->Draw("AP");

    static const Int_t n = 2;
    Double_t x[n] = {0, maxADC+40.};
    Double_t y[n] = {0.0, 0.0};

    TGraph *gr1 = new TGraph(n,x,y);
    gr1->Draw("Same");
    gr1->SetLineWidth(2);
    gr1->SetLineColor(1);
    gr1->SetLineStyle(2);

    x1_text = 1200;
    y1_text = -40;

    /*TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("East PMT 4");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();*/
  }

  else { 
    if (useTanh) linCurves << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << endl;
    else linCurves << 0. << " " << 0. << " " << 0. << endl;
    oFileE4 << "PMT NOT USABLE";}
  oFileE4.close();

 
  ////////////////////////////////////////////////////////////////////

  // West 1

  sprintf(temp,"%s/residuals/residuals_West_runPeriod_%i_PMTW1.dat",getenv("ANALYSIS_CODE"),calibrationPeriod);
  ofstream oFileW1(temp);
  vector <Double_t> fitEQ_W1(runW1.size(),0);

  if (runW1.size()>0 && std::find(nameW1.begin(),nameW1.end(),"Ce")!=nameW1.end() && std::find(nameW1.begin(),nameW1.end(),"Sn")!=nameW1.end() && std::find(nameW1.begin(),nameW1.end(),"Bi1")!=nameW1.end()  )  {
    //Calculating the average of the lower Ce ADC values and assigning a value to lowFitThreshold
    Double_t maxADC = 0., maxEQ = 0.;
    for (UInt_t ii=0; ii<runW1.size(); ii++) {
      if (ADCW1[ii]>maxADC) maxADC = ADCW1[ii];
      if (EQW1[ii]>maxEQ) maxEQ = EQW1[ii];
    }

    W1->cd();

    TGraphErrors *grW1 = new TGraphErrors(runW1.size(),&ADCW1[0],&EQW1[0],&ADCW1_err[0],&EQW1_err[0]);
    grW1->SetTitle("West PMT 1");
    grW1->SetMarkerColor(1);
    grW1->SetLineColor(1);
    grW1->SetMarkerStyle(20);
    grW1->SetMarkerSize(0.75);
    grW1->GetXaxis()->SetTitle("g(t)#upointADC [channels]");
    grW1->GetXaxis()->SetTitleOffset(1.2);
    grW1->GetXaxis()->CenterTitle();
    grW1->GetYaxis()->SetTitle("#eta(x,y)#upointE_{Q} [keV]");
    grW1->GetYaxis()->SetTitleOffset(1.6);
    grW1->GetYaxis()->CenterTitle();
    grW1->GetXaxis()->SetLimits(0.0,maxADC+40.);
    grW1->SetMinimum(0.0);
    grW1->SetMaximum(maxEQ+40.);
    grW1->Draw("AP");

    grW1->Fit("fitADC", "","", 0., maxADC+40.);

    offset = fitADC->GetParameter(0);
    slope = fitADC->GetParameter(1);
    quad = fitADC->GetParameter(2);
    
    shiftOfTanH = offset>=0. ? lowestADC[4]/2. : (lowestADC[4]-offset/slope)/2.;
    lowEnSlope = (quad*shiftOfTanH*shiftOfTanH + slope*shiftOfTanH + offset)/shiftOfTanH;
    spreadOfTanH = shiftOfTanH/3.;

    if (useTanh) linCurves << offset << " " << slope << " " << quad  << " " << lowEnSlope << " " << shiftOfTanH << " " << spreadOfTanH << endl;
    else linCurves << offset << " " << slope << " " << quad << endl;

    /*grW1->Fit("fitADC", "","",0.,maxADC+40.);

    offset = fitADC->GetParameter(0);
    slope = fitADC->GetParameter(1);
    quad = fitADC->GetParameter(2);
    
    linCurves << offset << " " << slope << " " << quad <<  " " << lowEnSlope << " " << shiftOfTanH << " " << spreadOfTanH << endl;
    */ 

    x1_text = 1200;
    y1_text = 100;

    /*TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("West PMT 1");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();*/
  
    // Calculate residuals in [keV]
    //Double_t fitEQ_W1[num];
    for (int j=0; j<runW1.size(); j++) {
      fitEQ_W1[j]    = fitADC->Eval(ADCW1[j]);//offset + slope*ADCW1[j] + quad*ADCW1[j]*ADCW1[j] + cubic*ADCW1[j]*ADCW1[j]*ADCW1[j];
      if (nameW1[j]=="Ce") {
	ResW1[j] = 100.*(fitEQ_W1[j] - EQW1[j])/fitEQ_W1[j];
	oFileW1 << "Ce_West" << " " << runW1[j] << " " << ResW1[j] << endl;
	if (sqrt(ResW1[j]*ResW1[j])>2.5) cout<< "Ce_West" << " " << "PMT1 " << runW1[j] << " " << ResW1[j] << endl;  
      }
      else if (nameW1[j]=="Sn") {
	ResW1[j] = 100.*(fitEQ_W1[j] - EQW1[j])/fitEQ_W1[j];
	oFileW1 << "Sn_West" << " " <<  runW1[j] << " " << ResW1[j] << endl;
	if (sqrt(ResW1[j]*ResW1[j])>2.5) cout<< "Sn_West" << " " << "PMT1 " << runW1[j] << " " << ResW1[j] << endl;  
      }
      else if (nameW1[j]=="Bi1") {
	ResW1[j] =100.*(fitEQ_W1[j] - EQW1[j])/fitEQ_W1[j];
	oFileW1 << "Bi1_West" << " " << runW1[j] << " " << ResW1[j] << endl;
	if (sqrt(ResW1[j]*ResW1[j])>2.5) cout<< "Bi1_West" << " " << "PMT1 " << runW1[j] << " " << ResW1[j] << endl;  
      }
      else if (nameW1[j]=="Bi2") {
	ResW1[j] = 100.*(fitEQ_W1[j] - EQW1[j])/fitEQ_W1[j];
	oFileW1 << "Bi2_West" << " " << runW1[j] << " " << ResW1[j] << endl;
	if (sqrt(ResW1[j]*ResW1[j])>2.5) cout<< "Bi2_West" << " " << "PMT1 " << runW1[j] << " " << ResW1[j] << endl;  
      }
      else if (nameW1[j]=="In") {
	ResW1[j] = 100.*(fitEQ_W1[j] - EQW1[j])/fitEQ_W1[j];
	oFileW1 << "In_West" << " " << runW1[j] << " " << ResW1[j] << endl;
	if (sqrt(ResW1[j]*ResW1[j])>2.5) cout<< "In_West" << " " << "PMT1 " << runW1[j] << " " << ResW1[j] << endl;  
      }
    }

    // West 1 residuals
    W1res->cd();

    TGraphErrors *grW1r = new TGraphErrors(runW1.size(),&ADCW1[0],&ResW1[0],err,err);
    grW1r->SetTitle("");
    grW1r->SetMarkerColor(1);
    grW1r->SetLineColor(1);
    grW1r->SetMarkerStyle(20);
    grW1r->SetMarkerSize(0.75);
    //grW1r->GetXaxis()->SetTitle("g(t)#upointADC [channels]");
    grW1r->GetXaxis()->SetTitleOffset(1.2);
    grW1r->GetXaxis()->CenterTitle();
    grW1r->GetYaxis()->SetTitle("% Residuals");
    grW1r->GetYaxis()->SetTitleSize(0.1);
    grW1r->GetYaxis()->SetLabelSize(0.07);
    grW1r->GetXaxis()->SetLabelSize(0);
    grW1r->GetYaxis()->SetTitleOffset(0.6);
    grW1r->GetYaxis()->CenterTitle();
    grW1r->GetXaxis()->SetLimits(0.0,maxADC+40.);
    grW1r->SetMinimum(-10.0);
    grW1r->SetMaximum( 10.0);
    grW1r->Draw("AP");

    static const Int_t n = 2;
    Double_t x[n] = {0, maxADC+40.};
    Double_t y[n] = {0.0, 0.0};

    TGraph *gr1 = new TGraph(n,x,y);
    gr1->Draw("Same");
    gr1->SetLineWidth(2);
    gr1->SetLineColor(1);
    gr1->SetLineStyle(2);

    x1_text = 1200;
    y1_text = -40;

    /*TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("West PMT 1");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();*/
  }

  else {
    if (useTanh) linCurves << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << endl;
    else linCurves << 0. << " " << 0. << " " << 0. << endl;
    oFileW1 << "PMT NOT USABLE";}
  oFileW1.close();

  // West 2

  sprintf(temp,"%s/residuals/residuals_West_runPeriod_%i_PMTW2.dat",getenv("ANALYSIS_CODE"),calibrationPeriod);
  ofstream oFileW2(temp);
  vector <Double_t> fitEQ_W2(runW2.size(),0);


  if (runW2.size()>0 && std::find(nameW2.begin(),nameW2.end(),"Ce")!=nameW2.end() && std::find(nameW2.begin(),nameW2.end(),"Sn")!=nameW2.end() && std::find(nameW2.begin(),nameW2.end(),"Bi1")!=nameW2.end()  )  {
  
    //Calculating the average of the lower Ce ADC values and assigning a value to lowFitThreshold
    Double_t maxADC = 0., maxEQ = 0.;
    for (UInt_t ii=0; ii<runW2.size(); ii++) {
      if (ADCW2[ii]>maxADC) maxADC = ADCW2[ii];
      if (EQW2[ii]>maxEQ) maxEQ = EQW2[ii];
    }

    W2->cd();

    TGraphErrors *grW2 = new TGraphErrors(runW2.size(),&ADCW2[0],&EQW2[0],&ADCW2_err[0],&EQW2_err[0]);
    grW2->SetTitle("West PMT 2");
    grW2->SetMarkerColor(1);
    grW2->SetLineColor(1);
    grW2->SetMarkerStyle(20);
    grW2->SetMarkerSize(0.75);
    grW2->GetXaxis()->SetTitle("g(t)#upointADC [channels]");
    grW2->GetXaxis()->SetTitleOffset(1.2);
    grW2->GetXaxis()->CenterTitle();
    grW2->GetYaxis()->SetTitle("#eta(x,y)#upointE_{Q} [keV]");
    grW2->GetYaxis()->SetTitleOffset(1.6);
    grW2->GetYaxis()->CenterTitle();
    grW2->GetXaxis()->SetLimits(0.0,maxADC+40.);
    grW2->SetMinimum(0.0);
    grW2->SetMaximum(maxEQ+40.);
    grW2->Draw("AP");

    grW2->Fit("fitADC", "","", 0., maxADC+40.);

    offset = fitADC->GetParameter(0);
    slope = fitADC->GetParameter(1);
    quad = fitADC->GetParameter(2);
    
    shiftOfTanH = offset>=0. ? lowestADC[5]/2. : (lowestADC[5]-offset/slope)/2.;
    lowEnSlope = (quad*shiftOfTanH*shiftOfTanH + slope*shiftOfTanH + offset)/shiftOfTanH;
    spreadOfTanH = shiftOfTanH/3.;

    if (useTanh) linCurves << offset << " " << slope << " " << quad  << " " << lowEnSlope << " " << shiftOfTanH << " " << spreadOfTanH << endl;
    else linCurves << offset << " " << slope << " " << quad << endl;


    /*grW2->Fit("fitADC", "","",0.,maxADC+40.);

    offset = fitADC->GetParameter(0);
    slope = fitADC->GetParameter(1);
    quad = fitADC->GetParameter(2);
    
    linCurves << offset << " " << slope << " " << quad << " " << lowEnSlope << " " << shiftOfTanH << " " << spreadOfTanH << endl;
    */ 

    x1_text = 1200;
    y1_text = 100;

    /*TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("West PMT 2");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();*/

    // Calculate residuals in [keV]
    //Double_t fitEQ_W2[num];
    for (int j=0; j<runW2.size(); j++) {
      fitEQ_W2[j]    = fitADC->Eval(ADCW2[j]);//offset + slope*ADCW2[j] + quad*ADCW2[j]*ADCW2[j] + cubic*ADCW2[j]*ADCW2[j]*ADCW2[j];
      if (nameW2[j]=="Ce") {
	ResW2[j] = 100.*(fitEQ_W2[j] - EQW2[j])/fitEQ_W2[j];
	oFileW2 << "Ce_West" << " " << runW2[j] << " " << ResW2[j] << endl;
	if (sqrt(ResW2[j]*ResW2[j])>2.5) cout<< "Ce_West" << " " << "PMT2 " << runW2[j] << " " << ResW2[j] << endl;  
      }
      else if (nameW2[j]=="Sn") {
	ResW2[j] = 100.*(fitEQ_W2[j] - EQW2[j])/fitEQ_W2[j];
	oFileW2 << "Sn_West" << " " << runW2[j] << " " << ResW2[j] << endl;
	if (sqrt(ResW2[j]*ResW2[j])>2.5) cout<< "Sn_West" << " " << "PMT2 " << runW2[j] << " " << ResW2[j] << endl;  
      }
      else if (nameW2[j]=="Bi1") {
	ResW2[j] = 100.*(fitEQ_W2[j] - EQW2[j])/fitEQ_W2[j];
	oFileW2 << "Bi1_West" << " " << runW2[j] << " " << ResW2[j] << endl;
	if (sqrt(ResW2[j]*ResW2[j])>2.5) cout<< "Bi1_West" << " " << "PMT2 " << runW2[j] << " " << ResW2[j] << endl;  
      }
      else if (nameW2[j]=="Bi2") {
	ResW2[j] = 100.*(fitEQ_W2[j] - EQW2[j])/fitEQ_W2[j];
	oFileW2 << "Bi2_West" << " " << runW2[j] << " " << ResW2[j] << endl;
	if (sqrt(ResW2[j]*ResW2[j])>2.5) cout<< "Bi2_West" << " " << "PMT2 " << runW2[j] << " " << ResW2[j] << endl;  
      }
      else if (nameW2[j]=="In") {
	ResW2[j] = 100.*(fitEQ_W2[j] - EQW2[j])/fitEQ_W2[j];
	oFileW2 << "In_West" << " " << runW2[j] << " " << ResW2[j] << endl;
	if (sqrt(ResW2[j]*ResW2[j])>2.5) cout<< "In_West" << " " << "PMT2 " << runW2[j] << " " << ResW2[j] << endl;  
      }
    }

    // West 2 residuals
    W2res->cd();

    TGraphErrors *grW2r = new TGraphErrors(runW2.size(),&ADCW2[0],&ResW2[0],err,err);
    grW2r->SetTitle("");
    grW2r->SetMarkerColor(1);
    grW2r->SetLineColor(1);
    grW2r->SetMarkerStyle(20);
    grW2r->SetMarkerSize(0.75);
    //grW2r->GetXaxis()->SetTitle("g(t)#upointADC [channels]");
    grW2r->GetXaxis()->SetTitleOffset(1.2);
    grW2r->GetXaxis()->CenterTitle();
    grW2r->GetYaxis()->SetTitle("% Residuals");
    grW2r->GetYaxis()->SetTitleSize(0.1);
    grW2r->GetYaxis()->SetLabelSize(0.07);
    grW2r->GetXaxis()->SetLabelSize(0);
    grW2r->GetYaxis()->SetTitleOffset(0.6);
    grW2r->GetYaxis()->CenterTitle();
    grW2r->GetXaxis()->SetLimits(0.0,maxADC+40.);
    grW2r->SetMinimum(-10.0);
    grW2r->SetMaximum( 10.0);
    grW2r->Draw("AP");

    static const Int_t n = 2;
    Double_t x[n] = {0, maxADC+40.};
    Double_t y[n] = {0.0, 0.0};

    TGraph *gr1 = new TGraph(n,x,y);
    gr1->Draw("Same");
    gr1->SetLineWidth(2);
    gr1->SetLineColor(1);
    gr1->SetLineStyle(2);

    x1_text = 1200;
    y1_text = -40;

    /*TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("West PMT 2");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();*/
  }

  else { 
    if (useTanh) linCurves << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << endl;
    else linCurves << 0. << " " << 0. << " " << 0. << endl;
    oFileW2 << "PMT NOT USABLE";}
  oFileW2.close();

  // West 3

  sprintf(temp,"%s/residuals/residuals_West_runPeriod_%i_PMTW3.dat",getenv("ANALYSIS_CODE"),calibrationPeriod);
  ofstream oFileW3(temp);
  vector <Double_t> fitEQ_W3(runW3.size(),0);

  if (runW3.size()>0 && std::find(nameW3.begin(),nameW3.end(),"Ce")!=nameW3.end() && std::find(nameW3.begin(),nameW3.end(),"Sn")!=nameW3.end() && std::find(nameW3.begin(),nameW3.end(),"Bi1")!=nameW3.end()  )  {

    //Calculating the average of the lower Ce ADC values and assigning a value to lowFitThreshold
    Double_t maxADC = 0., maxEQ = 0.;
    for (UInt_t ii=0; ii<runW3.size(); ii++) {
      if (ADCW3[ii]>maxADC) maxADC = ADCW3[ii];
      if (EQW3[ii]>maxEQ) maxEQ = EQW3[ii];
    }
    
    W3->cd();

    TGraphErrors *grW3 = new TGraphErrors(runW3.size(),&ADCW3[0],&EQW3[0],&ADCW3_err[0],&EQW3_err[0]);
    grW3->SetTitle("West PMT 3");
    grW3->SetMarkerColor(1);
    grW3->SetLineColor(1);
    grW3->SetMarkerStyle(20);
    grW3->SetMarkerSize(0.75);
    grW3->GetXaxis()->SetTitle("g(t)#upointADC [channels]");
    grW3->GetXaxis()->SetTitleOffset(1.2);
    grW3->GetXaxis()->CenterTitle();
    grW3->GetYaxis()->SetTitle("#eta(x,y)#upointE_{Q} [keV]");
    grW3->GetYaxis()->SetTitleOffset(1.6);
    grW3->GetYaxis()->CenterTitle();
    grW3->GetXaxis()->SetLimits(0.0,maxADC+40.);
    grW3->SetMinimum(0.0);
    grW3->SetMaximum(maxEQ+40.);
    grW3->Draw("AP");

    grW3->Fit("fitADC", "","", 0., maxADC+40.);

    offset = fitADC->GetParameter(0);
    slope = fitADC->GetParameter(1);
    quad = fitADC->GetParameter(2);
    
    shiftOfTanH = offset>=0. ? lowestADC[6]/2. : (lowestADC[6]-offset/slope)/2.;
    lowEnSlope = (quad*shiftOfTanH*shiftOfTanH + slope*shiftOfTanH + offset)/shiftOfTanH;
    spreadOfTanH = shiftOfTanH/3.;

    if (useTanh) linCurves << offset << " " << slope << " " << quad  << " " << lowEnSlope << " " << shiftOfTanH << " " << spreadOfTanH << endl;
    else linCurves << offset << " " << slope << " " << quad << endl;


    /*grW3->Fit("fitADC", "","",0.,maxADC+40.);
  
    offset = fitADC->GetParameter(0);
    slope = fitADC->GetParameter(1);
    quad = fitADC->GetParameter(2);
    
    linCurves << offset << " " << slope << " " << quad << " " << lowEnSlope << " " << shiftOfTanH << " " << spreadOfTanH << endl;
    */

    x1_text = 1200;
    y1_text = 100;

    /*TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("West PMT 3");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();*/

    // Calculate residuals in [keV]
    //Double_t fitEQ_W3[num];
    for (int j=0; j<runW3.size(); j++) {
      fitEQ_W3[j]    = fitADC->Eval(ADCW3[j]);//offset + slope*ADCW3[j] + quad*ADCW3[j]*ADCW3[j] + cubic*ADCW3[j]*ADCW3[j]*ADCW3[j];
      if (nameW3[j]=="Ce") {
	ResW3[j] = 100.*(fitEQ_W3[j] - EQW3[j])/fitEQ_W3[j];
	oFileW3 << "Ce_West" << " " << runW3[j] << " " << ResW3[j] << endl;
	if (sqrt(ResW3[j]*ResW3[j])>2.5) cout<< "Ce_West" << " " << "PMT3 " << runW3[j] << " " << ResW3[j] << endl;  
      }
      else if (nameW3[j]=="Sn") {
	ResW3[j] = 100.*(fitEQ_W3[j] - EQW3[j])/fitEQ_W3[j];
	oFileW3 << "Sn_West" << " " << runW3[j] << " " << ResW3[j] << endl;
	if (sqrt(ResW3[j]*ResW3[j])>2.5) cout<< "Sn_West" << " " << "PMT3 " << runW3[j] << " " << ResW3[j] << endl;  
      }
      else if (nameW3[j]=="Bi1") {
	ResW3[j] = 100.*(fitEQ_W3[j] - EQW3[j])/fitEQ_W3[j];
	oFileW3 << "Bi1_West" << " " << runW3[j] << " " << ResW3[j] << endl;
	if (sqrt(ResW3[j]*ResW3[j])>2.5) cout<< "Bi1_West" << " " << "PMT3 " << runW3[j] << " " << ResW3[j] << endl;  
      }
      else if (nameW3[j]=="Bi2") {
	ResW3[j] = 100.*(fitEQ_W3[j] - EQW3[j])/fitEQ_W3[j];
	oFileW3 << "Bi2_West" << " " << runW3[j] << " " << ResW3[j] << endl;
	if (sqrt(ResW3[j]*ResW3[j])>2.5) cout<< "Bi2_West" << " " << "PMT3 " << runW3[j] << " " << ResW3[j] << endl;  
      }
      else if (nameW3[j]=="In") {
	ResW3[j] = 100.*(fitEQ_W3[j] - EQW3[j])/fitEQ_W3[j];
	oFileW3 << "In_West" << " " << runW3[j] << " " << ResW3[j] << endl;
	if (sqrt(ResW3[j]*ResW3[j])>2.5) cout<< "In_West" << " " << "PMT3 " << runW3[j] << " " << ResW3[j] << endl;  
      }
    }

    // West 3 residuals
    W3res->cd();

    TGraphErrors *grW3r = new TGraphErrors(runW3.size(),&ADCW3[0],&ResW3[0],err,err);
    grW3r->SetTitle("");
    grW3r->SetMarkerColor(1);
    grW3r->SetLineColor(1);
    grW3r->SetMarkerStyle(20);
    grW3r->SetMarkerSize(0.75);
    //grW3r->GetXaxis()->SetTitle("g(t)#upointADC [channels]");
    grW3r->GetXaxis()->SetTitleOffset(1.2);
    grW3r->GetXaxis()->CenterTitle();
    grW3r->GetYaxis()->SetTitle("% Residuals");
    grW3r->GetYaxis()->SetTitleSize(0.1);
    grW3r->GetYaxis()->SetLabelSize(0.07);
    grW3r->GetXaxis()->SetLabelSize(0);
    grW3r->GetYaxis()->SetTitleOffset(0.6);
    grW3r->GetYaxis()->CenterTitle();
    grW3r->GetXaxis()->SetLimits(0.0,maxADC+40.);
    grW3r->SetMinimum(-10.0);
    grW3r->SetMaximum( 10.0);
    grW3r->Draw("AP");

    static const Int_t n = 2;
    Double_t x[n] = {0, maxADC+40.};
    Double_t y[n] = {0.0, 0.0};

    TGraph *gr1 = new TGraph(n,x,y);
    gr1->Draw("Same");
    gr1->SetLineWidth(2);
    gr1->SetLineColor(1);
    gr1->SetLineStyle(2);

    x1_text = 1200;
    y1_text = -40;

    /*TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("West PMT 3");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();*/
  }

  else { 
    if (useTanh) linCurves << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << endl;
    else linCurves << 0. << " " << 0. << " " << 0. << endl;
    oFileW3 << "PMT NOT USABLE";}
  oFileW3.close();

  // West 4

  sprintf(temp,"%s/residuals/residuals_West_runPeriod_%i_PMTW4.dat",getenv("ANALYSIS_CODE"),calibrationPeriod);
  ofstream oFileW4(temp);
  vector <Double_t> fitEQ_W4(runW4.size(),0);

  if (runW4.size()>0 && std::find(nameW4.begin(),nameW4.end(),"Ce")!=nameW4.end() && std::find(nameW4.begin(),nameW4.end(),"Sn")!=nameW4.end() && std::find(nameW4.begin(),nameW4.end(),"Bi1")!=nameW4.end()  )  {
    //Calculating the average of the lower Ce ADC values and assigning a value to lowFitThreshold
    Double_t maxADC = 0., maxEQ = 0.;
    for (UInt_t ii=0; ii<runW4.size(); ii++) {
      if (ADCW4[ii]>maxADC) maxADC = ADCW4[ii];
      if (EQW4[ii]>maxEQ) maxEQ = EQW4[ii];
    }


    W4->cd();

    TGraphErrors *grW4 = new TGraphErrors(runW4.size(),&ADCW4[0],&EQW4[0],&ADCW4_err[0],&EQW4_err[0]);
    grW4->SetTitle("West PMT 4");
    grW4->SetMarkerColor(1);
    grW4->SetLineColor(1);
    grW4->SetMarkerStyle(20);
    grW4->SetMarkerSize(0.75);
    grW4->GetXaxis()->SetTitle("g(t)#upointADC [channels]");
    grW4->GetXaxis()->SetTitleOffset(1.2);
    grW4->GetXaxis()->CenterTitle();
    grW4->GetYaxis()->SetTitle("#eta(x,y)#upointE_{Q} [keV]");
    grW4->GetYaxis()->SetTitleOffset(1.6);
    grW4->GetYaxis()->CenterTitle();
    grW4->GetXaxis()->SetLimits(0.0,maxADC+40.);
    grW4->SetMinimum(0.0);
    grW4->SetMaximum(maxEQ+40.);
    grW4->Draw("AP");

    grW4->Fit("fitADC", "","", 0., maxADC+40.);

    offset = fitADC->GetParameter(0);
    slope = fitADC->GetParameter(1);
    quad = fitADC->GetParameter(2);

    shiftOfTanH = offset>=0. ? lowestADC[7]/2. : (lowestADC[7]-offset/slope)/2.;
    lowEnSlope = (quad*shiftOfTanH*shiftOfTanH + slope*shiftOfTanH + offset)/shiftOfTanH;
    spreadOfTanH = shiftOfTanH/3.;

    if (useTanh) linCurves << offset << " " << slope << " " << quad  << " " << lowEnSlope << " " << shiftOfTanH << " " << spreadOfTanH << endl;
    else linCurves << offset << " " << slope << " " << quad << endl;

    //fitADC->SetRange(0., 3500.);
    /*grW4->Fit("fitADC", "","",0.,maxADC+40.);

    offset = fitADC->GetParameter(0);
    slope = fitADC->GetParameter(1);
    quad = fitADC->GetParameter(2);
    
    linCurves << offset << " " << slope << " " << quad << " " << lowEnSlope << " " << shiftOfTanH << " " << spreadOfTanH << endl;
    */

    linCurves.close();
  
    x1_text = 1200;
    y1_text = 100;

    /*TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("West PMT 4");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();*/

    // Calculate residuals in [keV]
    //Double_t fitEQ_W4[num];
    for (int j=0; j<runW4.size(); j++) {
      fitEQ_W4[j]    = fitADC->Eval(ADCW4[j]);//offset + slope*ADCW4[j] + quad*ADCW4[j]*ADCW4[j] + cubic*ADCW4[j]*ADCW4[j]*ADCW4[j];
      if (nameW4[j]=="Ce") {
	ResW4[j] = 100.*(fitEQ_W4[j] - EQW4[j])/fitEQ_W4[j];
	oFileW4 << "Ce_West" << " " << runW4[j] << " " << ResW4[j] << endl;
	if (sqrt(ResW4[j]*ResW4[j])>2.5) cout<< "Ce_West" << " " << "PMT4 " << runW4[j] << " " << ResW4[j] << endl;  
      }
      else if (nameW4[j]=="Sn") {
	ResW4[j] = 100.*(fitEQ_W4[j] - EQW4[j])/fitEQ_W4[j];
	oFileW4 << "Sn_West" << " " << runW4[j] << " " << ResW4[j] << endl;
	if (sqrt(ResW4[j]*ResW4[j])>2.5) cout<< "Sn_West" << " " << "PMT4 " << runW4[j] << " " << ResW4[j] << endl;  
      }
      else if (nameW4[j]=="Bi1"){
	ResW4[j] = 100.*(fitEQ_W4[j] - EQW4[j])/fitEQ_W4[j];
	oFileW4 << "Bi1_West" << " " << runW4[j] << " " << ResW4[j] << endl;
	if (sqrt(ResW4[j]*ResW4[j])>2.5) cout<< "Bi1_West" << " " << "PMT4 " << runW4[j] << " " << ResW4[j] << endl;  
      }
      else if (nameW4[j]=="Bi2"){
	ResW4[j] = 100.*(fitEQ_W4[j] - EQW4[j])/fitEQ_W4[j];
	oFileW4 << "Bi2_West" << " " << runW4[j] << " " << ResW4[j] << endl;
	if (sqrt(ResW4[j]*ResW4[j])>2.5) cout<< "Bi2_West" << " " << "PMT4 " << runW4[j] << " " << ResW4[j] << endl;  
      }
      else if (nameW4[j]=="In"){
	ResW4[j] = 100.*(fitEQ_W4[j] - EQW4[j])/fitEQ_W4[j];
	oFileW4 << "In_West" << " " << runW4[j] << " " << ResW4[j] << endl;
	if (sqrt(ResW4[j]*ResW4[j])>2.5) cout<< "In_West" << " " << "PMT4 " << runW4[j] << " " << ResW4[j] << endl;  
      }
    }

    // West 4 residuals
    W4res->cd();

    TGraphErrors *grW4r = new TGraphErrors(runW4.size(),&ADCW4[0],&ResW4[0],err,err);
    grW4r->SetTitle("");
    grW4r->SetMarkerColor(1);
    grW4r->SetLineColor(1);
    grW4r->SetMarkerStyle(20);
    grW4r->SetMarkerSize(0.75);
    //grW4r->GetXaxis()->SetTitle("g(t)#upointADC [channels]");
    grW4r->GetXaxis()->SetTitleOffset(1.2);
    grW4r->GetXaxis()->CenterTitle();
    grW4r->GetYaxis()->SetTitle("% Residuals");
    grW4r->GetYaxis()->SetTitleSize(0.1);
    grW4r->GetYaxis()->SetLabelSize(0.07);
    grW4r->GetXaxis()->SetLabelSize(0);
    grW4r->GetYaxis()->SetTitleOffset(0.6);
    grW4r->GetYaxis()->CenterTitle();
    grW4r->GetXaxis()->SetLimits(0.0,maxADC+40.);
    grW4r->SetMinimum(-10.0);
    grW4r->SetMaximum( 10.0);
    grW4r->Draw("AP");

    static const Int_t n = 2;
    Double_t x[n] = {0, maxADC+40.};
    Double_t y[n] = {0.0, 0.0};

    TGraph *gr1 = new TGraph(n,x,y);
    gr1->Draw("Same");
    gr1->SetLineWidth(2);
    gr1->SetLineColor(1);
    gr1->SetLineStyle(2);

    x1_text = 1200;
    y1_text = -40;

    /*TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("West PMT 4");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();*/
  }
  
  else { 
    if (useTanh) linCurves << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << endl;
    else linCurves << 0. << " " << 0. << " " << 0. << endl;
    oFileW4 << "PMT NOT USABLE";
  }
  oFileW4.close();
  linCurves.close();

  TString pdfFile = TString::Format("linCurves_SrcPeriod_%i.pdf",runPeriod);
  c1->Print(pdfFile);


}
