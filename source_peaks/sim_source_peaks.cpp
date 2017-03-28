#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include <sstream>

// ROOT libraries
#include "TRandom3.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TH1D.h>

#include "sourcePeaks.h"

#include "peaks.hh"

using namespace std;

struct PMT {
  double Evis[8];
  double etaEvis[8];
  double nPE[8]; 
};

int main(int argc, char *argv[])
{
  //cout.setf(ios::fixed, ios::floatfield);
  //cout.precision(12);

  // Run number integer
  istringstream ss(argv[1]);
  int runNumber;
  ss >> runNumber;
  cout << "runNumber = " << runNumber << endl;

  // Read source list
  int nSources = 0;
  string sourceName[3];
  string sourceNameLong[3];
  char tempList[500];
  sprintf(tempList, "%s/source_list_%s.dat",getenv("SOURCE_LIST"), argv[1]);
  //cout << tempList << endl;
  ifstream fileList(tempList);
  if (fileList.is_open()) fileList >> nSources;
  cout << " ... Number of sources: " << nSources << endl;
  for (int n=0; n<nSources; n++) {
    fileList >> sourceName[n];
    cout << "  " << sourceName[n] << endl;
  }

  //Checking if one of the sources is Bi so we can add in that we will fit the second Bi peak as well
  bool useLowBiPeak=false;
  int BiPeakIndex = 0;
  for (int n=0; n<nSources; n++) {
    if (sourceName[n]=="Bi") {
      useLowBiPeak=true;
      BiPeakIndex=n;
      continue;
    }
  }
      

  //Adding in here that we are only looking for runs with Ce, Sn, In, or Bi in them
  bool correctSource = false;
  bool useSource[3] = {false,false,false};
  for (int n=0; n<nSources; n++) {
    if (sourceName[n]=="Ce") {
      sourceNameLong[n] = "Ce139";
      correctSource = true;
      useSource[n]=true;
      continue;
    }
    if (sourceName[n]=="Sn") {
      sourceNameLong[n] = "Sn113";
      correctSource = true;
      useSource[n]=true;
      continue;
    }
    if (sourceName[n]=="Bi") {
      sourceNameLong[n] = "Bi207";
      correctSource = true;
      useSource[n]=true;
      continue;
    }
    if (sourceName[n]=="In") {
      string side = getIndiumSide(runNumber);
      if (side=="West") sourceNameLong[n] = "In114W";
      else if (side=="East") sourceNameLong[n] = "In114E";
      else {cout << "No Indium side determined\n"; continue;}
      correctSource = true;
      useSource[n]=true;
      continue; 
    }
  }
  
  if (!correctSource) {
    cout << "Source run with no Ce, Sn, In, or Bi\n";
    exit(0);
  }


  // Open output ntuple
  char tempOut[500];
  sprintf(tempOut, "%s/source_peaks/source_peaks_%s.root",getenv("REVCALSIM"), argv[1]);
  TFile *fileOut = new TFile(tempOut,"RECREATE");

  // Output histograms
  int nBin = 406;

  TH1D *his[3][8];
  his[0][0] = new TH1D("his1_E0", "", nBin,-30.0,2000.0);
  his[0][1] = new TH1D("his1_E1", "", nBin,-30.0,2000.0);
  his[0][2] = new TH1D("his1_E2", "", nBin,-30.0,2000.0);
  his[0][3] = new TH1D("his1_E3", "", nBin,-30.0,2000.0);
  his[0][4] = new TH1D("his1_W0", "", nBin,-30.0,2000.0);
  his[0][5] = new TH1D("his1_W1", "", nBin,-30.0,2000.0);
  his[0][6] = new TH1D("his1_W2", "", nBin,-30.0,2000.0);
  his[0][7] = new TH1D("his1_W3", "", nBin,-30.0,2000.0);

  his[1][0] = new TH1D("his2_E0", "", nBin,-30.0,2000.0);
  his[1][1] = new TH1D("his2_E1", "", nBin,-30.0,2000.0);
  his[1][2] = new TH1D("his2_E2", "", nBin,-30.0,2000.0);
  his[1][3] = new TH1D("his2_E3", "", nBin,-30.0,2000.0);
  his[1][4] = new TH1D("his2_W0", "", nBin,-30.0,2000.0);
  his[1][5] = new TH1D("his2_W1", "", nBin,-30.0,2000.0);
  his[1][6] = new TH1D("his2_W2", "", nBin,-30.0,2000.0);
  his[1][7] = new TH1D("his2_W3", "", nBin,-30.0,2000.0);

  his[2][0] = new TH1D("his3_E0", "", nBin,-30.0,2000.0);
  his[2][1] = new TH1D("his3_E1", "", nBin,-30.0,2000.0);
  his[2][2] = new TH1D("his3_E2", "", nBin,-30.0,2000.0);
  his[2][3] = new TH1D("his3_E3", "", nBin,-30.0,2000.0);
  his[2][4] = new TH1D("his3_W0", "", nBin,-30.0,2000.0);
  his[2][5] = new TH1D("his3_W1", "", nBin,-30.0,2000.0);
  his[2][6] = new TH1D("his3_W2", "", nBin,-30.0,2000.0);
  his[2][7] = new TH1D("his3_W3", "", nBin,-30.0,2000.0);


  for (Int_t src=0; src<nSources; src++) {
    if (!useSource[src]) continue;
    
    // Open input ntuple
    char tempIn[500];
    sprintf(tempIn, "%s/sources/revCalSim_%i_%s.root",getenv("REVCALSIM"), runNumber,sourceNameLong[src].c_str() );
    cout << tempIn << endl;
    TFile *fileIn = new TFile(tempIn, "READ");
    TTree *Tin = (TTree*)(fileIn->Get("revCalSim"));

    PMT pmt;
    Int_t PID, type, side;
    // Variables
    Tin->SetBranchAddress("PMT",&pmt);

    //Tin->SetBranchAddress("EvisE", &EvisE);
    //Tin->SetBranchAddress("EvisW", &EvisW);

    Tin->SetBranchAddress("PID",  &PID);
    Tin->SetBranchAddress("type", &type);
    Tin->SetBranchAddress("side", &side);

    int nEvents = Tin->GetEntries();
    cout << "Processing " << argv[1] << " ... " << endl;
    cout << "... nEvents = " << nEvents << endl;

  // Loop over events
    for (int i=0; i<nEvents; i++) {
      Tin->GetEvent(i);

      // Use Type 0 events
      if (type != 0 || PID!=1) continue;
      
      if (side == 0) {
	for (int p=0; p<4; p++) {
	  his[src][p]->Fill(pmt.etaEvis[p]);
	}
      }
      if (side == 1) {
	for (int p=4; p<8; p++) {
	  his[src][p]->Fill(pmt.etaEvis[p]);
	}	
      }
    } 
    fileIn->Close();
  }
  

  // Find maximum bin
  double maxBin[3][8]={0.};
  double maxCounts[3][8]={0.};
  for (int n=0; n<nSources; n++) {
    for (int j=0; j<8; j++) {
      
      maxBin[n][j] = his[n][j]->GetMaximumBin();
      maxCounts[n][j] = his[n][j]->GetBinContent(maxBin[n][j]);

    }
  }

  // Define histogram fit ranges
  double xLow[3][8]={0.}, xHigh[3][8]={0.};
  for (int n=0; n<nSources; n++) {
    
    for (int j=0; j<8; j++) {
      for (int i=maxBin[n][j]; i<nBin; i++) {
        if (his[n][j]->GetBinContent(i+1) < 0.33*maxCounts[n][j]) {
          xHigh[n][j] = his[n][j]->GetBinCenter(i+1);
	  //Check to make sure the value isn't too close to the maximum bin...
          if ((i-maxBin[n][j])<5) xHigh[n][j] = his[n][j]->GetXaxis()->GetBinCenter(maxBin[n][j])+250.;
	  break;
        }
	if (i==(nBin-1)) xHigh[n][j] = his[n][j]->GetBinCenter(i);
      }
      if (sourceName[n]!="Bi") {
	for (int i=maxBin[n][j]; i>0; i--) {
	  if (his[n][j]->GetBinContent(i-1) < 0.33*maxCounts[n][j]) {
	    xLow[n][j] = his[n][j]->GetBinCenter(i-1);
	    //if ((maxBin[n][j]-i)<5) xLow[n][j] = binCenterMax[n][j]-200.;
	    break;
	  }
	}
      }
      else {
	for (int i=maxBin[n][j]*0.5; i>0; i--) {
	  if (his[n][j]->GetBinContent(i-1) < 0.33*0.5*maxCounts[n][j]) {
	    xLow[n][j] = his[n][j]->GetBinCenter(i-1);
	    //if ((maxBin[n][j]-i)<5) xLow[n][j] = binCenterMax[n][j]-200.;
	    break;
	  }
	}
      }
    }
    
  }

  double fitMean[3][8]={0.};
  double lowBiFitMean[8]={0.};
  double fitSigma[3][8]={0.};
  double lowBiFitSigma[8]={0.};

  for (int n=0; n<nSources; n++) {

    if (useSource[n]) {

      for (int j=0; j<8; j++) {

	if (sourceName[n]!="Bi") {

	  //Double_t rangeLow = 5.;
	  //Double_t rangeHigh = 4096.;
	  SinglePeakHist sing(his[n][j], xLow[n][j], xHigh[n][j]);

	  if (sing.isGoodFit()) {
	    fitMean[n][j] = sing.ReturnMean();
	    fitSigma[n][j] = sing.ReturnSigma();
	  }

	  else  {
	    cout << "Run " << runNumber << " can't converge on " << sourceName[n] << " peak in PMT " << j << ". Trying one more time......" << endl;
	    sing.FitHist(maxBin[n][j], 40., his[n][j]->GetBinContent(maxBin[n][j]));

	    if (sing.isGoodFit()) { 
	      fitMean[n][j] = sing.ReturnMean();
	      fitSigma[n][j] = sing.ReturnSigma();
	    }

	    else cout << "RUN " << runNumber << " CAN'T CONVERGE ON " << sourceName[n] << " PEAK IN PMT " << j  << endl;
	  }
	}

	else {

	  //Double_t rangeLow = 5.;
	  //Double_t rangeHigh = 4096.;
	  DoublePeakHist doub(his[n][j], xLow[n][j], xHigh[n][j]);

	  if (doub.isGoodFit()) {	    
	    fitMean[n][j] = doub.ReturnMean1();
	    lowBiFitMean[j] = doub.ReturnMean2();
	    fitSigma[n][j] = doub.ReturnSigma1();
	    lowBiFitSigma[j] = doub.ReturnSigma2();
	  }

	  else  {
	    cout << "Run " << runNumber << " can't converge on " << sourceName[n] << " peak in PMT " << j << ". Trying one more time......" << endl;
	    doub.FitHist(maxBin[n][j], 40., his[n][j]->GetBinContent(maxBin[n][j]), 0.5*maxBin[n][j], 40., 0.5*his[n][j]->GetBinContent(maxBin[n][j]));

	    if (doub.isGoodFit()) {
	      fitMean[n][j] = doub.ReturnMean1();
	      lowBiFitMean[j] = doub.ReturnMean2();
	      fitSigma[n][j] = doub.ReturnSigma1();
	      lowBiFitSigma[j] = doub.ReturnSigma2();
	    }

	    else cout << "RUN " << runNumber << " CAN'T CONVERGE ON " << sourceName[n] << " PEAK IN PMT " << j  << endl;
	    
	  }
	}

      }
    }
  }


  // Write results to file
  char tempResults[500];
  sprintf(tempResults, "%s/source_peaks/source_peaks_%s.dat",getenv("REVCALSIM"), argv[1]);
  ofstream outResultsMean(tempResults);
  sprintf(tempResults, "%s/source_peaks/source_widths_%s.dat",getenv("REVCALSIM"), argv[1]);
  ofstream outResultsSigma(tempResults);

  for (int n=0; n<nSources; n++) {
    if (useSource[n]) {
      if (sourceName[n]=="Bi") sourceName[n]=sourceName[n]+"1";
      outResultsMean << runNumber << " "
		     << sourceName[n] << " "
		     << fitMean[n][0] << " "
		     << fitMean[n][1] << " "
		     << fitMean[n][2] << " "
		     << fitMean[n][3] << " "
		     << fitMean[n][4] << " "
		     << fitMean[n][5] << " "
		     << fitMean[n][6] << " "
		     << fitMean[n][7] << endl;
      outResultsSigma << runNumber << " "
		     << sourceName[n] << " "
		     << fitSigma[n][0] << " "
		     << fitSigma[n][1] << " "
		     << fitSigma[n][2] << " "
		     << fitSigma[n][3] << " "
		     << fitSigma[n][4] << " "
		     << fitSigma[n][5] << " "
		     << fitSigma[n][6] << " "
		     << fitSigma[n][7] << endl;
    }
  }

  if (useLowBiPeak) {
    outResultsMean << runNumber << " "
		   << "Bi2" << " "
		   << lowBiFitMean[0] << " "
		   << lowBiFitMean[1] << " "
		   << lowBiFitMean[2] << " "
		   << lowBiFitMean[3] << " "
		   << lowBiFitMean[4] << " "
		   << lowBiFitMean[5] << " "
		   << lowBiFitMean[6] << " "
		   << lowBiFitMean[7] << endl;
    outResultsSigma << runNumber << " "
		   << "Bi2" << " "
		   << lowBiFitSigma[0] << " "
		   << lowBiFitSigma[1] << " "
		   << lowBiFitSigma[2] << " "
		   << lowBiFitSigma[3] << " "
		   << lowBiFitSigma[4] << " "
		   << lowBiFitSigma[5] << " "
		   << lowBiFitSigma[6] << " "
		   << lowBiFitSigma[7] << endl;
  }
  outResultsMean.close();
  outResultsSigma.close();

  // Write output ntuple
  fileOut->Write();
  fileOut->Close();

  return 0;
}
