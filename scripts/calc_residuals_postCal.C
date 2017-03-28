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


void calc_residuals_postCal(Int_t runPeriod)
{
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
  sprintf(temp, "%s/residuals/source_runs_Erecon_RunPeriod_%i.dat",getenv("ANALYSIS_CODE"),calibrationPeriod);
  ifstream filein_data(temp);
  sprintf(temp, "%s/residuals/SIM_source_runs_Erecon_RunPeriod_%i.dat",getenv("ANALYSIS_CODE"),calibrationPeriod);
  ifstream filein_sim(temp);

  Int_t i = 0;
  Int_t run[500], simRun=0; 
  string sourceName[500], simSourceName="";
  Double_t dataEreconE[500], dataSigmaE[500], simEreconE[500], simSigmaE[500];
  Double_t dataEreconW[500], dataSigmaW[500], simEreconW[500], simSigmaW[500];
  //Double_t res[500];
  //Double_t err[500];
  while (filein_data >> run[i] >> sourceName[i]
	 >> dataEreconE[i] >> dataSigmaE[i] >> dataEreconW[i] >> dataSigmaW[i]) {
    filein_sim >> simRun >> simSourceName;
    if (simRun==run[i] && simSourceName==sourceName[i]) {
      filein_sim >> simEreconE[i] >> simSigmaE[i] >> simEreconW[i] >> simSigmaW[i]; 
    }
    else {
      cout << "Data and sim files don't match. Rerun MakeSourceCalibrationFiles!\n"; 
      exit(0);
    }
    if (filein_data.fail()) break;
    i++;
  }
  filein_data.close();
  filein_sim.close();
  Int_t num = i;
  cout << "Number of data peaks: " << num << endl;


  ///////////////////////////////////////////////////////////////////


  sprintf(temp,"%s/residuals/residuals_Erecon_runPeriod_%i.dat",getenv("ANALYSIS_CODE"),calibrationPeriod);
  ofstream oFile(temp);

  vector <double> resE(num,0);
  vector <double> resW(num,0);
  vector <double> x(num,0);
  
  for (int j=0; j<num; j++) {
  
    if (sourceName[j]=="Ce") {
      resE[j] = simEreconE[j] - dataEreconE[j];
      resW[j] = simEreconW[j] - dataEreconW[j];
      oFile << "Ce" << " " << (int) run[j] << " " << resE[j] << " " << resW[j] << endl;
    }
    if (sourceName[j]=="In") {
      resE[j] = simEreconE[j] - dataEreconE[j];
      resW[j] = simEreconW[j] - dataEreconW[j];
      oFile << "In" << " " << (int) run[j] << " " << resE[j] << " " << resW[j] << endl;
    }
    else if (sourceName[j]=="Sn") {
      resE[j] = simEreconE[j] - dataEreconE[j];
      resW[j] = simEreconW[j] - dataEreconW[j];
      oFile << "Sn" << " " << (int) run[j] << " " << resE[j] << " " << resW[j] << endl;
    }
    else if (sourceName[j]=="Bi2") {
      resE[j] = simEreconE[j] - dataEreconE[j];
      resW[j] = simEreconW[j] - dataEreconW[j];
      oFile << "Bi2" << " " << (int) run[j] << " " << resE[j] << " " << resW[j] << endl;
    }
    else if (sourceName[j]=="Bi1") {
      resE[j] = simEreconE[j] - dataEreconE[j];
      resW[j] = simEreconW[j] - dataEreconW[j];
      oFile << "Bi1" << " " << (int) run[j] << " " << resE[j] << " " << resW[j] << endl;
   }

  }
  oFile.close();
  
  
}
