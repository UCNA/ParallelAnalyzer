// This will fit the simulation peaks in Erecon for the weighted
// energy over all good PMTs


#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

// ROOT libraries
#include "TRandom3.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TH1D.h>

#include "runInfo.h"
#include "sourcePeaks.h"
#include "positionMapHandler.hh"
#include "peaks.hh"
//#include "calibrationTools.hh"

using namespace std;

struct ScintPosAdjusted {
  double ScintPosAdjE[3];
  double ScintPosAdjW[3];
};

struct PMT {
  double Evis[8];
  double etaEvis[8];
  double nPE[8]; 
};

//Returns the average eta value from the position map
vector < vector < Double_t > > getMeanEtaForSources(Int_t run) 
{
  Char_t temp[500];
  vector < vector < Double_t > > eta (3,vector <Double_t> (8,0.));
  sprintf(temp,"%s/meanEta_%i.dat",getenv("SOURCE_POSITIONS"),run);
  ifstream infile;
  infile.open(temp);
  Int_t i = 0;
  std::string name, rn;

  while (infile >> rn >> name >> eta[i][0] >> eta[i][1] >> eta[i][2] >> eta[i][3] >> eta[i][4] >> eta[i][5] >> eta[i][6] >> eta[i][7] ) {
    //if (name==src) break;
    //std::cout << eta[i][0] << eta[i][1] << eta[i][2] << eta[i][3] << eta[i][4] << eta[i][5] << eta[i][6] << eta[i][7] << std::endl;
    i++;
  }
  //if (name!=src) {std::cout << "Didn't find " << src << " in this run!\n"; exit(0);}
  return eta;
}

int main(int argc, char *argv[])
{
  //cout.setf(ios::fixed, ios::floatfield);
  //cout.precision(12);

  // Run number integer
  int runNumber = atoi(argv[1]);
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

  vector < vector <Double_t> > meanEta = getMeanEtaForSources(runNumber); // fill vector with the calculated mean position correction for the source

  
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
  sprintf(tempOut, "%s/source_peaks/source_peaks_%s_Erecon.root",getenv("REVCALSIM"), argv[1]);
  TFile *fileOut = new TFile(tempOut,"RECREATE");


  //Set up all the output histograms

  Int_t nBin = 600;
  TH1D *hisEvisTot[3][2];
  hisEvisTot[0][0] = new TH1D("hisEvis1E", "source1 East", nBin,-50.0,2400.0);
  hisEvisTot[0][1] = new TH1D("hisEvis1W", "source1 West", nBin,-50.0,2400.0);
  hisEvisTot[1][0] = new TH1D("hisEvis2E", "source2 East", nBin,-50.0,2400.0);
  hisEvisTot[1][1] = new TH1D("hisEvis2W", "source2 West", nBin,-50.0,2400.0);
  hisEvisTot[2][0] = new TH1D("hisEvis3E", "source3 East", nBin,-50.0,2400.0);
  hisEvisTot[2][1] = new TH1D("hisEvis3W", "source3 West", nBin,-50.0,2400.0);

  TH1D *hisEvis[3][8];
  hisEvis[0][0] = new TH1D("hisEvis1_E0", "", nBin,-50.0,2400.0);
  hisEvis[0][1] = new TH1D("hisEvis1_E1", "", nBin,-50.0,2400.0);
  hisEvis[0][2] = new TH1D("hisEvis1_E2", "", nBin,-50.0,2400.0);
  hisEvis[0][3] = new TH1D("hisEvis1_E3", "", nBin,-50.0,2400.0);
  hisEvis[0][4] = new TH1D("hisEvis1_W0", "", nBin,-50.0,2400.0);
  hisEvis[0][5] = new TH1D("hisEvis1_W1", "", nBin,-50.0,2400.0);
  hisEvis[0][6] = new TH1D("hisEvis1_W2", "", nBin,-50.0,2400.0);
  hisEvis[0][7] = new TH1D("hisEvis1_W3", "", nBin,-50.0,2400.0);

  hisEvis[1][0] = new TH1D("hisEvis2_E0", "", nBin,-50.0,2400.0);
  hisEvis[1][1] = new TH1D("hisEvis2_E1", "", nBin,-50.0,2400.0);
  hisEvis[1][2] = new TH1D("hisEvis2_E2", "", nBin,-50.0,2400.0);
  hisEvis[1][3] = new TH1D("hisEvis2_E3", "", nBin,-50.0,2400.0);
  hisEvis[1][4] = new TH1D("hisEvis2_W0", "", nBin,-50.0,2400.0);
  hisEvis[1][5] = new TH1D("hisEvis2_W1", "", nBin,-50.0,2400.0);
  hisEvis[1][6] = new TH1D("hisEvis2_W2", "", nBin,-50.0,2400.0);
  hisEvis[1][7] = new TH1D("hisEvis2_W3", "", nBin,-50.0,2400.0);

  hisEvis[2][0] = new TH1D("hisEvis3_E0", "", nBin,-50.0,2400.0);
  hisEvis[2][1] = new TH1D("hisEvis3_E1", "", nBin,-50.0,2400.0);
  hisEvis[2][2] = new TH1D("hisEvis3_E2", "", nBin,-50.0,2400.0);
  hisEvis[2][3] = new TH1D("hisEvis3_E3", "", nBin,-50.0,2400.0);
  hisEvis[2][4] = new TH1D("hisEvis3_W0", "", nBin,-50.0,2400.0);
  hisEvis[2][5] = new TH1D("hisEvis3_W1", "", nBin,-50.0,2400.0);
  hisEvis[2][6] = new TH1D("hisEvis3_W2", "", nBin,-50.0,2400.0);
  hisEvis[2][7] = new TH1D("hisEvis3_W3", "", nBin,-50.0,2400.0);
  

  TH1D *hisEreconTot[3][2];
  hisEreconTot[0][0] = new TH1D("hisErecon1E", "source1 East", nBin,0.0,2400.0);
  hisEreconTot[0][1] = new TH1D("hisErecon1W", "source1 West", nBin,0.0,2400.0);
  hisEreconTot[1][0] = new TH1D("hisErecon2E", "source2 East", nBin,0.0,2400.0);
  hisEreconTot[1][1] = new TH1D("hisErecon2W", "source2 West", nBin,0.0,2400.0);
  hisEreconTot[2][0] = new TH1D("hisErecon3E", "source3 East", nBin,0.0,2400.0);
  hisEreconTot[2][1] = new TH1D("hisErecon3W", "source3 West", nBin,0.0,2400.0);

  

  // Load calibration info
  PositionMap posmap(5.0);
  posmap.readPositionMap(getXeRunPeriod(runNumber));
  
  vector < Double_t > eta;
  std::vector < std::vector <Int_t> > numDataPoints(3, std::vector <Int_t>(2,0)); //Holds the number of data points for each side for each source
  std::vector < std::vector <Double_t> > aveEta(3, std::vector <Double_t>(8,0.)); //Holds the average value of eta for the the data being read in for each source and each PMT

  for (Int_t src=0; src<nSources; src++) {
    if (!useSource[src]) continue;
    
    // Open input ntuple
    char tempIn[500];
    sprintf(tempIn, "%s/sources/revCalSim_%i_%s.root",getenv("REVCALSIM"), runNumber,sourceNameLong[src].c_str() );
    cout << tempIn << endl;
    TFile *fileIn = new TFile(tempIn, "READ");
    TTree *Tin = (TTree*)(fileIn->Get("revCalSim"));

    PMT pmt;
    Double_t Erecon;
    Int_t PID, type, side;
    Double_t Evis[2];
    ScintPosAdjusted pos;
    Int_t nClipped_EX=0, nClipped_EY=0, nClipped_WX=0, nClipped_WY=0;
    // Variables
    Tin->SetBranchAddress("PMT",&pmt);
    Tin->SetBranchAddress("Erecon", &Erecon); 
    Tin->SetBranchAddress("Evis", Evis);
    Tin->SetBranchAddress("ScintPosAdjusted", &pos);
    Tin->SetBranchAddress("PID",  &PID);
    Tin->SetBranchAddress("type", &type);
    Tin->SetBranchAddress("side", &side);
    Tin->SetBranchAddress("nClipped_EX", &nClipped_EX);
    Tin->SetBranchAddress("nClipped_EY", &nClipped_EY);
    Tin->SetBranchAddress("nClipped_WX", &nClipped_WX);
    Tin->SetBranchAddress("nClipped_WY", &nClipped_WY);

    int nEvents = Tin->GetEntries();
    cout << "Processing " << argv[1] << " ... " << endl;
    cout << "... nEvents = " << nEvents << endl;

  // Loop over events
    for (int i=0; i<nEvents; i++) {
      Tin->GetEvent(i);
      
      // Use Type 0 events
      if (type != 0 || PID!=1) continue;

      // Cut the clipped events
      if ( nClipped_EX>0 || nClipped_EY>0 || nClipped_WX>0 || nClipped_WY>0 ) continue; 

      eta = posmap.getInterpolatedEta(pos.ScintPosAdjE[0], pos.ScintPosAdjE[1], pos.ScintPosAdjW[0], pos.ScintPosAdjW[1]);
      
      if (side == 0) {

	numDataPoints[src][0]+=1;
	aveEta[src][0] += eta[0];
	aveEta[src][1] += eta[1];
	aveEta[src][2] += eta[2];
	aveEta[src][3] += eta[3];
	
	if (Erecon>=0.) hisEreconTot[src][0]->Fill(Erecon);
	hisEvisTot[src][0]->Fill(Evis[0]);
	
	for (int p=0; p<4; p++) {
	  hisEvis[src][p]->Fill(pmt.Evis[p]);
	}
      }
      
      if (side == 1) {

	numDataPoints[src][1]+=1;
	aveEta[src][4] += eta[4];
	aveEta[src][5] += eta[5];
	aveEta[src][6] += eta[6];
	aveEta[src][7] += eta[7];
	
	if (Erecon>=0.) hisEreconTot[src][1]->Fill(Erecon);
	hisEvisTot[src][1]->Fill(Evis[1]);
	
	for (int p=4; p<8; p++) {
	  hisEvis[src][p]->Fill(pmt.Evis[p]);
	}	
      }
    } 
    fileIn->Close();
  }
  

  //Finish average eta calculation
  for (int src=0; src<nSources; src++) {
    for (int pmt=0; pmt<8; pmt++) {

      aveEta[src][pmt] = pmt<4 ? aveEta[src][pmt]/(double)numDataPoints[src][0] : aveEta[src][pmt]/(double)numDataPoints[src][1];
      
    }
  }


  
  // Find maximum bin
  double maxBin[3][8]={0.};
  double maxCounts[3][8]={0.};
  for (int n=0; n<nSources; n++) {
    for (int j=0; j<8; j++) {
      
      maxBin[n][j] = hisEvis[n][j]->GetMaximumBin();
      maxCounts[n][j] = hisEvis[n][j]->GetBinContent(maxBin[n][j]);

    }
  }

  // Define histogram fit ranges
  double xLow[3][8]={0.}, xHigh[3][8]={0.};
  for (int n=0; n<nSources; n++) {
    
    for (int j=0; j<8; j++) {
      for (int i=maxBin[n][j]; i<nBin; i++) {
        if (hisEvis[n][j]->GetBinContent(i+1) < 0.25*maxCounts[n][j]) {
          xHigh[n][j] = hisEvis[n][j]->GetBinCenter(i+1);
	  //Check to make sure the value isn't too close to the maximum bin...
          if ((i-maxBin[n][j])<5) xHigh[n][j] = hisEvis[n][j]->GetXaxis()->GetBinCenter(maxBin[n][j])+250.;
	  break;
        }
	if (i==(nBin-1)) xHigh[n][j] = hisEvis[n][j]->GetBinCenter(i);
      }
      if (sourceName[n]!="Bi") {
	for (int i=maxBin[n][j]; i>0; i--) {
	  if (hisEvis[n][j]->GetBinContent(i-1) < 0.25*maxCounts[n][j]) {
	    xLow[n][j] = hisEvis[n][j]->GetBinCenter(i-1);
	    //if ((maxBin[n][j]-i)<5) xLow[n][j] = binCenterMax[n][j]-200.;
	    break;
	  }
	}
      }
      else {
	for (int i=maxBin[n][j]*0.5; i>0; i--) {
	  if (hisEvis[n][j]->GetBinContent(i-1) < 0.2*0.5*maxCounts[n][j]) {
	    xLow[n][j] = hisEvis[n][j]->GetBinCenter(i-1);
	    //if ((maxBin[n][j]-i)<5) xLow[n][j] = binCenterMax[n][j]-200.;
	    break;
	  }
	}
      }
    }
    
  }

  double fitMean[3][8]={0.};
  double fitMeanError[3][8]={0.};
  double lowBiFitMean[8]={0.};
  double lowBiFitMeanError[8]={0.};
  double fitSigma[3][8]={0.};
  double lowBiFitSigma[8]={0.};

  for (int n=0; n<nSources; n++) {

    if (useSource[n]) {

      for (int j=0; j<8; j++) {

	if (sourceName[n]!="Bi") {

	  //Double_t rangeLow = 5.;
	  //Double_t rangeHigh = 4096.;
	  SinglePeakHist sing(hisEvis[n][j], xLow[n][j], xHigh[n][j]);

	  if (sing.isGoodFit()) {
	    fitMean[n][j] = sing.ReturnMean();
	    fitMeanError[n][j] = sing.ReturnMeanError();
	    fitSigma[n][j] = sing.ReturnSigma();
	  }

	  else  {
	    cout << "Run " << runNumber << " can't converge on " << sourceName[n] << " peak in PMT " << j << ". Trying one more time......" << endl;
	    sing.SetRangeMin(xLow[n][j]);
	    sing.SetRangeMax(xHigh[n][j]);
	    sing.FitHist(maxBin[n][j], 40., hisEvis[n][j]->GetBinContent(maxBin[n][j]));

	    if (sing.isGoodFit()) { 
	      fitMean[n][j] = sing.ReturnMean();
	      fitMeanError[n][j] = sing.ReturnMeanError();
	      fitSigma[n][j] = sing.ReturnSigma();
	    }

	    else cout << "RUN " << runNumber << " CAN'T CONVERGE ON " << sourceName[n] << " PEAK IN PMT " << j  << endl;
	  }
	}

	else {

	  //Double_t rangeLow = 5.;
	  //Double_t rangeHigh = 4096.;
	  DoublePeakHist doub(hisEvis[n][j], xLow[n][j], xHigh[n][j]);

	  if (doub.isGoodFit()) {	    
	    fitMean[n][j] = doub.ReturnMean1();
	    fitMeanError[n][j] = doub.ReturnMean1Error();
	    lowBiFitMean[j] = doub.ReturnMean2();
	    lowBiFitMeanError[j] = doub.ReturnMean2Error();
	    fitSigma[n][j] = doub.ReturnSigma1();
	    lowBiFitSigma[j] = doub.ReturnSigma2();
	  }

	  else  {
	    cout << "Run " << runNumber << " can't converge on " << sourceName[n] << " peak in PMT " << j << ". Trying one more time......" << endl;
	    doub.SetRangeMin(xLow[n][j]);
	    doub.SetRangeMax(xHigh[n][j]);
	    doub.FitHist(maxBin[n][j], 40., hisEvis[n][j]->GetBinContent(maxBin[n][j]), 0.5*maxBin[n][j], 40., 0.5*hisEvis[n][j]->GetBinContent(maxBin[n][j]));

	    if (doub.isGoodFit()) {
	      fitMean[n][j] = doub.ReturnMean1();
	      fitMeanError[n][j] = doub.ReturnMean1Error();
	      lowBiFitMean[j] = doub.ReturnMean2();
	      lowBiFitMeanError[j] = doub.ReturnMean2Error();
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
  sprintf(tempResults, "%s/source_peaks/source_peaks_%s_Evis.dat",getenv("REVCALSIM"), argv[1]);
  ofstream outResultsMean(tempResults);
  sprintf(tempResults, "%s/source_peaks/source_peaks_errors_%s_Evis.dat",getenv("REVCALSIM"), argv[1]);
  ofstream outResultsMeanError(tempResults);
  sprintf(tempResults, "%s/source_peaks/source_widths_%s_Evis.dat",getenv("REVCALSIM"), argv[1]);
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
      outResultsMeanError << runNumber << " "
		     << sourceName[n] << " "
		     << fitMeanError[n][0] << " "
		     << fitMeanError[n][1] << " "
		     << fitMeanError[n][2] << " "
		     << fitMeanError[n][3] << " "
		     << fitMeanError[n][4] << " "
		     << fitMeanError[n][5] << " "
		     << fitMeanError[n][6] << " "
		     << fitMeanError[n][7] << endl;
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
    outResultsMeanError << runNumber << " "
		   << "Bi2" << " "
		   << lowBiFitMeanError[0] << " "
		   << lowBiFitMeanError[1] << " "
		   << lowBiFitMeanError[2] << " "
		   << lowBiFitMeanError[3] << " "
		   << lowBiFitMeanError[4] << " "
		   << lowBiFitMeanError[5] << " "
		   << lowBiFitMeanError[6] << " "
		   << lowBiFitMeanError[7] << endl;
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
  outResultsMeanError.close();
  outResultsSigma.close();


  //Now for the total Evis peaks and widths

  // Find maximum bin
  double maxBinEvisTot[3][2]={0.};
  double maxCountsEvisTot[3][2]={0.};
  for (int n=0; n<nSources; n++) {
    for (int j=0; j<2; j++) {
      maxBinEvisTot[n][j] = hisEvisTot[n][j]->GetMaximumBin();
      maxCountsEvisTot[n][j] = hisEvisTot[n][j]->GetBinContent(maxBinEvisTot[n][j]);
    }
  }


  // Define histogram fit ranges
  double xLowEvisTot[3][2]={0.}, xHighEvisTot[3][2]={0.};
  for (int n=0; n<nSources; n++) { 
    for (int j=0; j<2; j++) {
      for (int i=maxBinEvisTot[n][j]; i<nBin; i++) {
	if (hisEvisTot[n][j]->GetBinContent(i+1) < 0.33*maxCountsEvisTot[n][j]) {
	  xHighEvisTot[n][j] = hisEvisTot[n][j]->GetBinCenter(i+1);
	  //Check to make sure the value isn't too close to the maximum bin...
	  if ((i-maxBinEvisTot[n][j])<5) xHighEvisTot[n][j] = hisEvisTot[n][j]->GetXaxis()->GetBinCenter(maxBinEvisTot[n][j])+250.;
	  break;
	}
	if (i==(nBin-1)) xHighEvisTot[n][j] = hisEvisTot[n][j]->GetBinCenter(i);
      }
      if (sourceName[n].substr(0,2)!="Bi") {
	for (int i=maxBinEvisTot[n][j]; i>0; i--) {
	  if (hisEvisTot[n][j]->GetBinContent(i-1) < 0.33*maxCountsEvisTot[n][j]) {
	    xLowEvisTot[n][j] = hisEvisTot[n][j]->GetBinCenter(i-1);
	    //if ((maxBinEvisTot[n][j]-i)<5) xLowEvisTot[n][j] = binCenterMax[n][j]-200.;
	    break;
	  }
	}
      }
      else {
	for (int i=maxBinEvisTot[n][j]*0.5; i>0; i--) {
	  if (hisEvisTot[n][j]->GetBinContent(i-1) < 0.33*0.5*maxCountsEvisTot[n][j]) {
	    xLowEvisTot[n][j] = hisEvisTot[n][j]->GetBinCenter(i-1);
	    //if ((maxBinEvisTot[n][j]-i)<5) xLowEvisTot[n][j] = binCenterMax[n][j]-200.;
	    break;
	  }
	}
      }
    }
  }

  double fitMeanEvisTot[3][2]={0.};
  double lowBiFitMeanEvisTot[2]={0.};
  double fitSigmaEvisTot[3][2]={0.};
  double lowBiFitSigmaEvisTot[2]={0.};

  for (int n=0; n<nSources; n++) {

    if (useSource[n]) {
      
      for (int j=0; j<2; j++) {

	if (sourceName[n].substr(0,2)!="Bi") {
	
	  
	  //Double_t rangeLowEvisTot = 5.;
	  //Double_t rangeHighEvisTot = 4096.;
	  SinglePeakHist sing(hisEvisTot[n][j], xLowEvisTot[n][j], xHighEvisTot[n][j]);
	  
	  if (sing.isGoodFit()) {
	    fitMeanEvisTot[n][j] = sing.ReturnMean();
	    fitSigmaEvisTot[n][j] = sing.ReturnSigma();
	  }
	  
	  else  {
	    cout << "Run " << runNumber << " can't converge on " << sourceName[n] << " EvisTot peak.... Trying one more time......" << endl;
	    sing.SetRangeMin(xLowEvisTot[n][j]);
	    sing.SetRangeMax(xHighEvisTot[n][j]);
	    sing.FitHist(maxBinEvisTot[n][j], 40., hisEvisTot[n][j]->GetBinContent(maxBinEvisTot[n][j]));
	    
	    if (sing.isGoodFit()) { 
	      fitMeanEvisTot[n][j] = sing.ReturnMean();
	      fitSigmaEvisTot[n][j] = sing.ReturnSigma();
	    }
	    
	    else cout << "RUN " << runNumber << " CAN'T CONVERGE ON " << sourceName[n] << " EvisTot PEAK "<<  endl;
	  }
	}
	
	else {
	  
	  
	  //DoublePeakHist doub(hisEvisTot[n][j], xLowEvisTot[n][j], xHighEvisTot[n][j]);
	  SinglePeakHist singBi1(hisEvisTot[n][j], 800., xHighEvisTot[n][j]);
	  SinglePeakHist singBi2(hisEvisTot[n][j], 340., 530.);
	  
	  if (singBi1.isGoodFit()) {	    
	    fitMeanEvisTot[n][j] = singBi1.ReturnMean();	 
	    fitSigmaEvisTot[n][j] = singBi1.ReturnSigma();
	  }
	  else  {
	    cout << "Run " << runNumber << " can't converge on " << sourceName[n] << " EvisTot peak.... Trying one more time......" << endl;
	    singBi1.SetRangeMin(775.);
	    singBi1.SetRangeMax(xHighEvisTot[n][j]);
	    singBi1.FitHist(910., 40., hisEvisTot[n][j]->GetBinContent(maxBinEvisTot[n][j]));
	    
	    if (singBi1.isGoodFit()) {
	      fitMeanEvisTot[n][j] = singBi1.ReturnMean();
	      fitSigmaEvisTot[n][j] = singBi1.ReturnSigma();
	    }
	    
	    else cout << "RUN " << runNumber << " CAN'T CONVERGE ON " << sourceName[n] << " EvisTot PEAK " << endl;
	    
	  }
	  
	  if (singBi2.isGoodFit()) {	    
	    lowBiFitMeanEvisTot[j] = singBi2.ReturnMean();	 
	    lowBiFitSigmaEvisTot[j] = singBi2.ReturnSigma();
	  }
	  else  {
	    cout << "Run " << runNumber << " can't converge on " << sourceName[n] << " EvisTot peak.... Trying one more time......" << endl;
	    singBi2.SetRangeMin(330.);
	    singBi2.SetRangeMax(530.);
	    singBi2.FitHist(440., 40., 0.4*hisEvisTot[n][j]->GetBinContent(maxBinEvisTot[n][j]));
	    
	    if (singBi2.isGoodFit()) {
	      lowBiFitMeanEvisTot[j] = singBi2.ReturnMean();
	      lowBiFitSigmaEvisTot[j] = singBi2.ReturnSigma();
	    }
	    
	    else cout << "RUN " << runNumber << " CAN'T CONVERGE ON " << sourceName[n] << " EvisTot PEAK " << endl;
	  }
	}
      }
      
    }
  }
  
  
  sprintf(tempResults, "%s/source_peaks/source_peaks_%s_EvisTot.dat",getenv("REVCALSIM"), argv[1]);
  ofstream outResultsEvisTot(tempResults);
  
  for (int n=0; n<nSources; n++) {
    if (useSource[n]) {
      if (sourceName[n]=="Bi") sourceName[n]=sourceName[n]+"1";
     
      outResultsEvisTot << runNumber << " "
		       << sourceName[n] << " "
		       << fitMeanEvisTot[n][0] << " "
		       << fitSigmaEvisTot[n][0] << " "
		       << fitMeanEvisTot[n][1] << " " 
		       << fitSigmaEvisTot[n][1] << std::endl;
      
    }
  }
  
  if (useLowBiPeak) {
    outResultsEvisTot << runNumber << " "
		     << "Bi2" << " "
		     << lowBiFitMeanEvisTot[0] << " "
		     << lowBiFitSigmaEvisTot[0] << " "
		     << lowBiFitMeanEvisTot[1] << " " 
		     << lowBiFitSigmaEvisTot[1] << std::endl;
  }
  outResultsEvisTot.close();




  //Now for the Erecon peaks and widths

  // Find maximum bin
  double maxBinEreconTot[3][2]={0.};
  double maxCountsEreconTot[3][2]={0.};
  for (int n=0; n<nSources; n++) {
    for (int j=0; j<2; j++) {
      maxBinEreconTot[n][j] = hisEreconTot[n][j]->GetMaximumBin();
      maxCountsEreconTot[n][j] = hisEreconTot[n][j]->GetBinContent(maxBinEreconTot[n][j]);
    }
  }


  // Define histogram fit ranges
  double xLowEreconTot[3][2]={0.}, xHighEreconTot[3][2]={0.};
  for (int n=0; n<nSources; n++) { 
    for (int j=0; j<2; j++) {
      for (int i=maxBinEreconTot[n][j]; i<nBin; i++) {
	if (hisEreconTot[n][j]->GetBinContent(i+1) < 0.33*maxCountsEreconTot[n][j]) {
	  xHighEreconTot[n][j] = hisEreconTot[n][j]->GetBinCenter(i+1);
	  //Check to make sure the value isn't too close to the maximum bin...
	  if ((i-maxBinEreconTot[n][j])<5) xHighEreconTot[n][j] = hisEreconTot[n][j]->GetXaxis()->GetBinCenter(maxBinEreconTot[n][j])+250.;
	  break;
	}
	if (i==(nBin-1)) xHighEreconTot[n][j] = hisEreconTot[n][j]->GetBinCenter(i);
      }
      if (sourceName[n].substr(0,2)!="Bi") {
	for (int i=maxBinEreconTot[n][j]; i>0; i--) {
	  if (hisEreconTot[n][j]->GetBinContent(i-1) < 0.33*maxCountsEreconTot[n][j]) {
	    xLowEreconTot[n][j] = hisEreconTot[n][j]->GetBinCenter(i-1);
	    //if ((maxBinEreconTot[n][j]-i)<5) xLowEreconTot[n][j] = binCenterMax[n][j]-200.;
	    break;
	  }
	}
      }
      else {
	for (int i=maxBinEreconTot[n][j]*0.5; i>0; i--) {
	  if (hisEreconTot[n][j]->GetBinContent(i-1) < 0.33*0.5*maxCountsEreconTot[n][j]) {
	    xLowEreconTot[n][j] = hisEreconTot[n][j]->GetBinCenter(i-1);
	    //if ((maxBinEreconTot[n][j]-i)<5) xLowEreconTot[n][j] = binCenterMax[n][j]-200.;
	    break;
	  }
	}
      }
    }
  }

  double fitMeanEreconTot[3][2]={0.};
  double lowBiFitMeanEreconTot[2]={0.};
  double fitSigmaEreconTot[3][2]={0.};
  double lowBiFitSigmaEreconTot[2]={0.};

  for (int n=0; n<nSources; n++) {

    if (useSource[n]) {
      
      for (int j=0; j<2; j++) {

	if (sourceName[n].substr(0,2)!="Bi") {
	
	  
	  //Double_t rangeLowEreconTot = 5.;
	  //Double_t rangeHighEreconTot = 4096.;
	  SinglePeakHist sing(hisEreconTot[n][j], xLowEreconTot[n][j], xHighEreconTot[n][j]);
	  
	  if (sing.isGoodFit()) {
	    fitMeanEreconTot[n][j] = sing.ReturnMean();
	    fitSigmaEreconTot[n][j] = sing.ReturnSigma();
	  }
	  
	  else  {
	    cout << "Run " << runNumber << " can't converge on " << sourceName[n] << " EreconTot peak.... Trying one more time......" << endl;
	    sing.SetRangeMin(xLowEreconTot[n][j]);
	    sing.SetRangeMax(xHighEreconTot[n][j]);
	    sing.FitHist(maxBinEreconTot[n][j], 40., hisEreconTot[n][j]->GetBinContent(maxBinEreconTot[n][j]));
	    
	    if (sing.isGoodFit()) { 
	      fitMeanEreconTot[n][j] = sing.ReturnMean();
	      fitSigmaEreconTot[n][j] = sing.ReturnSigma();
	    }
	    
	    else cout << "RUN " << runNumber << " CAN'T CONVERGE ON " << sourceName[n] << " EreconTot PEAK "<<  endl;
	  }
	}
	
	else {
	  
	  
	  //DoublePeakHist doub(hisEreconTot[n][j], xLowEreconTot[n][j], xHighEreconTot[n][j]);
	  SinglePeakHist singBi1(hisEreconTot[n][j], 850., xHighEreconTot[n][j]);
	  SinglePeakHist singBi2(hisEreconTot[n][j], 390., 580.);
	  
	  if (singBi1.isGoodFit()) {	    
	    fitMeanEreconTot[n][j] = singBi1.ReturnMean();	 
	    fitSigmaEreconTot[n][j] = singBi1.ReturnSigma();
	  }
	  else  {
	    cout << "Run " << runNumber << " can't converge on " << sourceName[n] << " EreconTot peak.... Trying one more time......" << endl;
	    singBi1.SetRangeMin(850.);
	    singBi1.SetRangeMax(990.);
	    singBi1.FitHist(960., 40., hisEreconTot[n][j]->GetBinContent(maxBinEreconTot[n][j]));
	    
	    if (singBi1.isGoodFit()) {
	      fitMeanEreconTot[n][j] = singBi1.ReturnMean();
	      fitSigmaEreconTot[n][j] = singBi1.ReturnSigma();
	    }
	    
	    else cout << "RUN " << runNumber << " CAN'T CONVERGE ON " << sourceName[n] << " EreconTot PEAK " << endl;
	    
	  }
	  
	  if (singBi2.isGoodFit()) {	    
	    lowBiFitMeanEreconTot[j] = singBi2.ReturnMean();	 
	    lowBiFitSigmaEreconTot[j] = singBi2.ReturnSigma();
	  }
	  else  {
	    cout << "Run " << runNumber << " can't converge on " << sourceName[n] << " EreconTot peak.... Trying one more time......" << endl;
	    singBi2.SetRangeMin(390.);
	    singBi2.SetRangeMax(580.);
	    singBi2.FitHist(480., 40., 0.4*hisEreconTot[n][j]->GetBinContent(maxBinEreconTot[n][j]));
	    
	    if (singBi2.isGoodFit()) {
	      lowBiFitMeanEreconTot[j] = singBi2.ReturnMean();
	      lowBiFitSigmaEreconTot[j] = singBi2.ReturnSigma();
	    }
	    
	    else cout << "RUN " << runNumber << " CAN'T CONVERGE ON " << sourceName[n] << " EreconTot PEAK " << endl;
	  }
	}
      }
      
    }
  }
  
  
  sprintf(tempResults, "%s/source_peaks/source_peaks_%s_EreconTot.dat",getenv("REVCALSIM"), argv[1]);
  ofstream outResultsEreconTot(tempResults);
  
  for (int n=0; n<nSources; n++) {
    if (useSource[n]) {
      if (sourceName[n]=="Bi") sourceName[n]=sourceName[n]+"1";
     
      outResultsEreconTot << runNumber << " "
		       << sourceName[n] << " "
		       << fitMeanEreconTot[n][0] << " "
		       << fitSigmaEreconTot[n][0] << " "
		       << fitMeanEreconTot[n][1] << " " 
		       << fitSigmaEreconTot[n][1] << std::endl;
      
    }
  }
  
  if (useLowBiPeak) {
    outResultsEreconTot << runNumber << " "
		     << "Bi2" << " "
		     << lowBiFitMeanEreconTot[0] << " "
		     << lowBiFitSigmaEreconTot[0] << " "
		     << lowBiFitMeanEreconTot[1] << " " 
		     << lowBiFitSigmaEreconTot[1] << std::endl;
  }
  outResultsEreconTot.close();
  

  // NOW I CALCULATE THE ADJUSTED Evis value for calibration purposes
  //std::vector < std::vector < Double_t > > ADCbar(nSources, std::vector <Double_t> (8,0.));
  for (int n=0; n<nSources; n++) {
    for (int p=0; p<8; p++) {
      //ADCbar[n][p] = linearityCurve.applyInverseLinCurve(p,fitMean[n][p]*aveEta[n][p<4?0:1]);
      std::cout << "Source " << n << " PMT " << p << " aveEta :" << meanEta[n][p] << std::endl;
    }
  }
 
  sprintf(tempResults, "%s/source_peaks/source_peaks_%s_etaEvis.dat",getenv("REVCALSIM"), argv[1]);
  ofstream outResultsEtaEvis(tempResults);
  sprintf(tempResults, "%s/source_peaks/source_peaks_errors_%s_etaEvis.dat",getenv("REVCALSIM"), argv[1]);
  ofstream outResultsEtaEvisError(tempResults);

  for (int n=0; n<nSources; n++) {
    if (useSource[n]) {
      if (sourceName[n]=="Bi") sourceName[n]=sourceName[n]+"1";
      outResultsEtaEvis << runNumber << " "
			<< sourceName[n] << " "
			<< fitMean[n][0]*meanEta[n][0] << " "
			<< fitMean[n][1]*meanEta[n][1] << " "
			<< fitMean[n][2]*meanEta[n][2] << " "
			<< fitMean[n][3]*meanEta[n][3] << " "
			<< fitMean[n][4]*meanEta[n][4] << " "
			<< fitMean[n][5]*meanEta[n][5] << " "
			<< fitMean[n][6]*meanEta[n][6] << " "
			<< fitMean[n][7]*meanEta[n][7] << " " << endl;
      outResultsEtaEvisError << runNumber << " "
			<< sourceName[n] << " "
			<< fitMeanError[n][0]*meanEta[n][0] << " "
			<< fitMeanError[n][1]*meanEta[n][1] << " "
			<< fitMeanError[n][2]*meanEta[n][2] << " "
			<< fitMeanError[n][3]*meanEta[n][3] << " "
			<< fitMeanError[n][4]*meanEta[n][4] << " "
			<< fitMeanError[n][5]*meanEta[n][5] << " "
			<< fitMeanError[n][6]*meanEta[n][6] << " "
			<< fitMeanError[n][7]*meanEta[n][7] << " " << endl;
    }
  }

  if (useLowBiPeak) {
    outResultsEtaEvis << runNumber << " "
		      << "Bi2" << " "
		      << lowBiFitMean[0]*meanEta[BiPeakIndex][0] << " "
		      << lowBiFitMean[1]*meanEta[BiPeakIndex][1] << " "
		      << lowBiFitMean[2]*meanEta[BiPeakIndex][2] << " "
		      << lowBiFitMean[3]*meanEta[BiPeakIndex][3] << " "
		      << lowBiFitMean[4]*meanEta[BiPeakIndex][4] << " "
		      << lowBiFitMean[5]*meanEta[BiPeakIndex][5] << " "
		      << lowBiFitMean[6]*meanEta[BiPeakIndex][6] << " "
		      << lowBiFitMean[7]*meanEta[BiPeakIndex][7] << " " << endl;
    outResultsEtaEvisError << runNumber << " "
		      << "Bi2" << " "
		      << lowBiFitMeanError[0]*meanEta[BiPeakIndex][0] << " "
		      << lowBiFitMeanError[1]*meanEta[BiPeakIndex][1] << " "
		      << lowBiFitMeanError[2]*meanEta[BiPeakIndex][2] << " "
		      << lowBiFitMeanError[3]*meanEta[BiPeakIndex][3] << " "
		      << lowBiFitMeanError[4]*meanEta[BiPeakIndex][4] << " "
		      << lowBiFitMeanError[5]*meanEta[BiPeakIndex][5] << " "
		      << lowBiFitMeanError[6]*meanEta[BiPeakIndex][6] << " "
		      << lowBiFitMeanError[7]*meanEta[BiPeakIndex][7] << " " << endl;
		   
  }
  outResultsEtaEvis.close();

  // Write output ntuple
  fileOut->Write();
  fileOut->Close();

  return 0;

}
  
