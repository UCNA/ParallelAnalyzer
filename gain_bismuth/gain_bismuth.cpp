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

#include "fullTreeVariables.h"
#include "MWPCGeometry.h"
#include "pedestals.h"
#include "cuts.h"
#include "basic_reconstruction.h"
#include "runInfo.h"
#include "DataTree.hh"
#include "peaks.hh"

using namespace std;

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
};


int main(int argc, char *argv[])
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  // Run number integer
  istringstream ss(argv[1]);
  int runNumber;
  ss >> runNumber;
  cout << "runNumber = " << runNumber << endl;


  //First check that if the run isn't the reference run, the reference has been run at least
  unsigned int referenceRun = getGainReferenceRun(runNumber);
  std::cout << "Reference Run = " << referenceRun << std::endl;
  char tempIn[500];

  if ((int)referenceRun!=runNumber) {
    sprintf(tempIn,"%s/gain_bismuth_%i.dat",getenv("GAIN_BISMUTH"),referenceRun);
    ifstream refGainFile(tempIn);
    if (refGainFile.good()) {
      refGainFile.close();
    }
    else {
      cout << "The reference run hasn't been run yet! Run reference runs then try again!\n";
      exit(0);
    }
  }

  //Struct to hold the scintillator information (struct is in DataTree.hh)
  Scint ScintE, ScintW;

  // Open output ntuple
  char tempOut[500];
  sprintf(tempOut, "%s/gain_bismuth_%s.root",getenv("GAIN_BISMUTH"),argv[1]);
  TFile *fileOut = new TFile(tempOut,"RECREATE");


  // Output histograms
  int nBin = 200;
  TH1D *his[8];
  his[0] = new TH1D("hisE0", "", nBin,0.0,4000.0);
  his[1] = new TH1D("hisE1", "", nBin,0.0,4000.0);
  his[2] = new TH1D("hisE2", "", nBin,0.0,4000.0);
  his[3] = new TH1D("hisE3", "", nBin,0.0,4000.0);
  his[4] = new TH1D("hisW0", "", nBin,0.0,4000.0);
  his[5] = new TH1D("hisW1", "", nBin,0.0,4000.0);
  his[6] = new TH1D("hisW2", "", nBin,0.0,4000.0);
  his[7] = new TH1D("hisW3", "", nBin,0.0,4000.0);

  // Open input ntuple
  
  sprintf(tempIn, "%s/replay_pass1_%s.root",getenv("REPLAY_PASS1"), argv[1]);
  TFile *fileIn = new TFile(tempIn, "READ");
  TTree *Tin = (TTree*)(fileIn->Get("pass1"));
  int nEvents = Tin->GetEntries();
  cout << "Processing " << argv[1] << " ... " << endl;
  cout << "... nEvents = " << nEvents << endl;

  // Input variables
  Tin->SetBranchAddress("ScintE", &ScintE);
  Tin->SetBranchAddress("ScintW", &ScintW);
  
  Tin->SetBranchAddress("PID", &PID);
  Tin->SetBranchAddress("Type", &type);
  Tin->SetBranchAddress("Side", &side);

  //count Bi pulser events for each PMT to scale the "fuzz" level below
  std::vector <Int_t> nBiEvents(8,0);

  // Loop over events
  for (int i=0; i<nEvents; i++) {
    Tin->GetEvent(i);

    // Select Bi pulser events
    if (PID != 4) continue;

    // Fill PMT histograms
    his[0]->Fill(ScintE.q1);
    his[1]->Fill(ScintE.q2);
    his[2]->Fill(ScintE.q3);
    his[3]->Fill(ScintE.q4);
    his[4]->Fill(ScintW.q1);
    his[5]->Fill(ScintW.q2);
    his[6]->Fill(ScintW.q3);
    his[7]->Fill(ScintW.q4);

    if (ScintE.q1>500.) nBiEvents[0]++;
    if (ScintE.q2>500.) nBiEvents[1]++;
    if (ScintE.q3>500.) nBiEvents[2]++;
    if (ScintE.q4>500.) nBiEvents[3]++;
    if (ScintW.q1>500.) nBiEvents[4]++;
    if (ScintW.q2>500.) nBiEvents[5]++;
    if (ScintW.q3>500.) nBiEvents[6]++;
    if (ScintW.q4>500.) nBiEvents[7]++;
  }

  // Find maximum bin
  double maxBin[8];
  double maxCounts[8];
  for (int i=0; i<8; i++) {
    maxCounts[i] = -1.0;
  }

  double binCenter[8];
  double binCenterMax[8];
  double binCounts[8];

  // counting from the end bin to avoid weird high peaks at lower energies
  
  for (int j=0; j<8; j++) {
    //int counter = 0;
    int maxCounter = 0;
    int minCounter = 0;
    cout << nBiEvents[j] << endl;
    Int_t fuzzCounts = 0.004*nBiEvents[j];

    for (int i=nBin-1; i>0; i--) {
      binCenter[j] = his[j]->GetBinCenter(i);
      binCounts[j] = his[j]->GetBinContent(i);      

      if ( binCounts[j] > fuzzCounts && his[j]->GetBinContent(i+1) > 0 ) { // use bins above the fuzz and ignore the overflow bin (which isn't the last bin after pedestal subtraction) 
	if ( binCounts[j] >= maxCounts[j] ) {
	  maxCounts[j] = binCounts[j];
	  maxBin[j] = i;
	  binCenterMax[j] = binCenter[j];
	  
	  //if ( minCounter < 3 ) maxCounter++; //increment this when we have a new max. This is to make sure that we actually climb a peak
	  maxCounter++; //increment this when we have a new max. This is to make sure that we actually climb a peak
	  minCounter = 0;

	}
        
	else minCounter++;
	
	if ( minCounter>10 && maxCounter>4 ) break; //Making sure we only get the high energy peak
      }
    }
  }
  

  /*if (runNumber<21274) {
    for (int i=1; i<nBin-1; i++) {
      for (int j=0; j<8; j++) {
	binCenter[j] = his[j]->GetBinCenter(i);
	binCounts[j] = his[j]->GetBinContent(i);
	if (binCenter[j] > 500. && binCenter[j] < 3500. && binCounts[j] >= maxCounts[j]) {
	  maxCounts[j] = binCounts[j];
	  maxBin[j] = i;
	  binCenterMax[j] = binCenter[j];
	}
      }
    }
  }

  else { // taking care of the odd Bi pulser shape in EPMT4
    for (int i=nBin-1; i>0; i--) {
      for (int j=0; j<8; j++) {
	binCenter[j] = his[j]->GetBinCenter(i);
	binCounts[j] = his[j]->GetBinContent(i);
	if (binCenter[j] > 2250. && binCenter[j] < 3500. && binCounts[j] >= maxCounts[j]) {
	  maxCounts[j] = binCounts[j];
	  maxBin[j] = i;
	  binCenterMax[j] = binCenter[j];
	}
      }
    }
    }*/

  

  // Define histogram fit ranges
  double xLow[8], xHigh[8];
  //if (runNumber<21274) {
    for (int n=0; n<8; n++) {
      for (int i=maxBin[n]; i<nBin; i++) {
	if ( his[n]->GetBinContent(i+1) < 0.4*maxCounts[n] ) {
	  xHigh[n] =  (his[n]->GetXaxis()->GetBinCenter(i+1) - binCenterMax[n]) > 225. ? his[n]->GetBinCenter(i+1) : (binCenterMax[n] + 300.);
	  break;
	}
      }
      for (int i=maxBin[n]; i>0; i--) {
	if (his[n]->GetBinContent(i-1) < 0.55*maxCounts[n]) {
	  xLow[n] = ( binCenterMax[n] - his[n]->GetXaxis()->GetBinCenter(i-1) ) > 200. ? his[n]->GetBinCenter(i-1) : (binCenterMax[n] - 250.);
	  break;
	}
      }
    }
    //}
    /*else { // taking care of the odd Bi pulser shape in EPMT4
    for (int n=0; n<8; n++) {
      for (int i=maxBin[n]; i<nBin; i++) {
	if (his[n]->GetBinContent(i+1) < 0.4*maxCounts[n]) {
	  xHigh[n] = his[n]->GetBinCenter(i+1);
	  break;
	}
      }
      for (int i=maxBin[n]; i>0; i--) {
	if (his[n]->GetBinContent(i-1) < 0.4*maxCounts[n]) {
	  xLow[n] = his[n]->GetBinCenter(i-1);
	  break;
	}
      }
    }
    xHigh[3] = binCenterMax[3]+400.;
    xLow[3] = binCenterMax[3]-300.;
    }*/

  // Fit parameters
  double fitMean[8];
  for (int m=0; m<8; m++) {
    fitMean[m] = -1.0;
  }


  SinglePeakHist *fitter;

  // Fit for East PMT #0
  fitter = new SinglePeakHist(his[0],xLow[0],xHigh[0]);
  fitMean[0] = fitter->ReturnMean();
  delete fitter;

  cout << "fitMean[0] = " << fitMean[0] << endl;

  // Fit for East PMT #1
  fitter = new SinglePeakHist(his[1],xLow[1],xHigh[1]);
  fitMean[1] = fitter->ReturnMean();
  delete fitter;

  cout << "fitMean[1] = " << fitMean[1] << endl;

  // Fit for East PMT #2
  fitter = new SinglePeakHist(his[2],xLow[2],xHigh[2]);
  fitMean[2] = fitter->ReturnMean();
  delete fitter;

  cout << "fitMean[2] = " << fitMean[2] << endl;

  // Fit for East PMT #3
  fitter = new SinglePeakHist(his[3],xLow[3],xHigh[3]);
  fitMean[3] = fitter->ReturnMean();
  delete fitter;

  cout << "fitMean[3] = " << fitMean[3] << endl;


  // Fit for West PMT #0
  fitter = new SinglePeakHist(his[4],xLow[4],xHigh[4]);
  fitMean[4] = fitter->ReturnMean();
  delete fitter;

  cout << "fitMean[4] = " << fitMean[4] << endl;

  // Fit for West PMT #1
  fitter = new SinglePeakHist(his[5],xLow[5],xHigh[5]);
  fitMean[5] = fitter->ReturnMean();
  delete fitter;

  cout << "fitMean[5] = " << fitMean[5] << endl;

  // Fit for West PMT #2
  fitter = new SinglePeakHist(his[6],xLow[6],xHigh[6]);
  fitMean[6] = fitter->ReturnMean();
  delete fitter;

  cout << "fitMean[6] = " << fitMean[6] << endl;

  // Fit for West PMT #3
  fitter = new SinglePeakHist(his[7],xLow[7],xHigh[7]);
  fitMean[7] = fitter->ReturnMean();
  delete fitter;

  cout << "fitMean[7] = " << fitMean[7] << endl;

  /* TF1 *gaussian_fit_E0 = new TF1("gaussian_fit_E0",
                                 "[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",
                                 xLow[0], xHigh[0]);
  gaussian_fit_E0->SetParameter(0, maxCounts[0]);
  gaussian_fit_E0->SetParameter(1, binCenterMax[0]);
  gaussian_fit_E0->SetParameter(2, 200.);
  gaussian_fit_E0->SetLineColor(2);
  his[0]->Fit("gaussian_fit_E0", "LR");
  fitMean[0] = gaussian_fit_E0->GetParameter(1);
  cout << "fitMean[0] = " << fitMean[0] << endl;

  // Fit for East PMT #1
  TF1 *gaussian_fit_E1 = new TF1("gaussian_fit_E1",
                                 "[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",
                                 xLow[1], xHigh[1]);
  gaussian_fit_E1->SetParameter(0, maxCounts[1]);
  gaussian_fit_E1->SetParameter(1, binCenterMax[1]);
  gaussian_fit_E1->SetParameter(2, 200.);
  gaussian_fit_E1->SetLineColor(2);
  his[1]->Fit("gaussian_fit_E1", "LR");
  fitMean[1] = gaussian_fit_E1->GetParameter(1);
  cout << "fitMean[1] = " << fitMean[1] << endl;

  // Fit for East PMT #2
  TF1 *gaussian_fit_E2 = new TF1("gaussian_fit_E2",
                                 "[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",
                                 xLow[2], xHigh[2]);
  gaussian_fit_E2->SetParameter(0, maxCounts[2]);
  gaussian_fit_E2->SetParameter(1, binCenterMax[2]);
  gaussian_fit_E2->SetParameter(2, 200.);
  gaussian_fit_E2->SetLineColor(2);
  his[2]->Fit("gaussian_fit_E2", "LR");
  fitMean[2] = gaussian_fit_E2->GetParameter(1);
  cout << "fitMean[2] = " << fitMean[2] << endl;

  // Fit for East PMT #3
  TF1 *gaussian_fit_E3 = new TF1("gaussian_fit_E3",
                                 "[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",
                                 xLow[3], xHigh[3]);
  gaussian_fit_E3->SetParameter(0, maxCounts[3]);
  gaussian_fit_E3->SetParameter(1, binCenterMax[3]);
  gaussian_fit_E3->SetParameter(2, 200.);
  gaussian_fit_E3->SetLineColor(2);
  his[3]->Fit("gaussian_fit_E3", "LR");
  fitMean[3] = gaussian_fit_E3->GetParameter(1);
  cout << "fitMean[3] = " << fitMean[3] << endl;

  // Fit for West PMT #0
  TF1 *gaussian_fit_W0 = new TF1("gaussian_fit_W0",
                                 "[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",
                                 xLow[4], xHigh[4]);
  gaussian_fit_W0->SetParameter(0, maxCounts[4]);
  gaussian_fit_W0->SetParameter(1, binCenterMax[4]);
  gaussian_fit_W0->SetParameter(2, 200.);
  gaussian_fit_W0->SetLineColor(2);
  his[4]->Fit("gaussian_fit_W0", "LR");
  fitMean[4] = gaussian_fit_W0->GetParameter(1);
  cout << "fitMean[4] = " << fitMean[4] << endl;

  // Fit for West PMT #1
  TF1 *gaussian_fit_W1 = new TF1("gaussian_fit_W1",
                                 "[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",
                                 xLow[5], xHigh[5]);
  gaussian_fit_W1->SetParameter(0, maxCounts[5]);
  gaussian_fit_W1->SetParameter(1, binCenterMax[5]);
  gaussian_fit_W1->SetParameter(2, 200.);
  gaussian_fit_W1->SetLineColor(2);
  his[5]->Fit("gaussian_fit_W1", "LR");
  fitMean[5] = gaussian_fit_W1->GetParameter(1);
  cout << "fitMean[5] = " << fitMean[5] << endl;

  // Fit for West PMT #2
  TF1 *gaussian_fit_W2 = new TF1("gaussian_fit_W2",
                                 "[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",
                                 xLow[6], xHigh[6]);
  gaussian_fit_W2->SetParameter(0, maxCounts[6]);
  gaussian_fit_W2->SetParameter(1, binCenterMax[6]);
  gaussian_fit_W2->SetParameter(2, 200.);
  gaussian_fit_W2->SetLineColor(2);
  his[6]->Fit("gaussian_fit_W2", "LR");
  fitMean[6] = gaussian_fit_W2->GetParameter(1);
  cout << "fitMean[6] = " << fitMean[6] << endl;

  // Fit for West PMT #3
  TF1 *gaussian_fit_W3 = new TF1("gaussian_fit_W3",
                                 "[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",
                                 xLow[7], xHigh[7]);
  gaussian_fit_W3->SetParameter(0, maxCounts[7]);
  gaussian_fit_W3->SetParameter(1, binCenterMax[7]);
  gaussian_fit_W3->SetParameter(2, 200.);
  gaussian_fit_W3->SetLineColor(2);
  his[7]->Fit("gaussian_fit_W3", "LR");
  fitMean[7] = gaussian_fit_W3->GetParameter(1);
  cout << "fitMean[7] = " << fitMean[7] << endl;*/

  //Get PMT quality status
  std::vector<Int_t> pmtquality = getPMTQuality(runNumber);

  // Calculate gain correction factors
  double referenceMean[8];
  double gainCorrection[8];

  if ((int)referenceRun!=runNumber) {
    sprintf(tempIn,"%s/gain_bismuth_%i.dat",getenv("GAIN_BISMUTH"),referenceRun);
    ifstream refGainFile(tempIn);
    double numHold;
    int ii=0;
    while (refGainFile >> referenceMean[ii] >> numHold) ii++;
    refGainFile.close();
      
    for (int n=0; n<8; n++) {
      gainCorrection[n] = fitMean[n]>0. ? (referenceMean[n] / fitMean[n]) : 1.;
    }    
  }
  else {
    for (int n=0; n<8; n++) {
      gainCorrection[n] = 1.;
    }
  }

  // Close input ntuple
  fileIn->Close();

  // Write output ntuple
  fileOut->Write();
  fileOut->Close();

  // Write fit results to file
  char tempFit[500];
  sprintf(tempFit, "%s/gain_bismuth_%s.dat", getenv("GAIN_BISMUTH"),argv[1]);
  ofstream outFit(tempFit);

  for (int n=0; n<8; n++) {
    //outFit << gainCorrection[n] << endl;
    outFit << fitMean[n] << " " << gainCorrection[n] << endl;
  }

  return 0;
}
