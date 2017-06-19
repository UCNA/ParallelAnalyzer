//////////////////////////////////////////////////////////
// Applies Bi pulser gain monitoring corrections to 
// Scintillator PMT signals
//////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>

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
#include "DataTree.hh"

#include "replay_pass2.h"


using namespace std;

bool OnlyReplayBadFiles = false;

//Read in beam drops
std::vector < std::vector < Double_t > > readBeamDrops(Int_t runNumber) {
  
  std::vector < std::vector < Double_t > > bd; // beam drop times (low,high)
  
  ifstream ifile(TString::Format("%s/beamCuts_%i.dat",getenv("CUTS"),runNumber));
  
  std::string holdTxt;
  std::string rateBefore, rateAfter;
  Int_t numCuts;
  
  ifile >> holdTxt >> rateBefore >> holdTxt
	>> holdTxt >> rateAfter  >> holdTxt
	>> holdTxt >> numCuts;
  
  Double_t tlow, tup;
  for ( Int_t i=0; i<numCuts; ++i ) {
    std::vector < Double_t > pair;
    ifile >> tlow >> tup;

    pair.push_back(tlow);
    pair.push_back(tup);
    
    bd.push_back(pair);
  }
  
  ifile.close();
  
  return bd;
};


vector <Int_t> getPMTQuality(Int_t runNumber) {
  //Read in PMT quality file
  cout << "Reading in PMT Quality file ...\n";
  vector <Int_t>  pmtQuality (8,0);
  Char_t temp[200];
  
  // NOTE: for now i'm using EreconQuality until I can see if the bad EPMT4 crap 
  // can be salvaged...
  sprintf(temp,"%s/residuals/PMT_EreconQuality_master.dat",getenv("ANALYSIS_CODE")); 
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


  
  // the second argument is a boolean. If it is given and it is false, then don't use the beam cuts
  bool useBeamCuts = true;
  if ( argc == 3 ) {
    if ( std::string(argv[2])==std::string("false") || std::string(argv[2])==std::string("0") ) useBeamCuts = false;
  }

  cout << "Run " << argv[1] << " ..." << endl;
  cout << "... Applying Bi pulser gain corrections and cutting beam drops and bursts if applicable ..." << endl;

  // Check if file is corrupt
  char tempOut[500];
  sprintf(tempOut, "%s/replay_pass2_%s.root",getenv("REPLAY_PASS2"), argv[1]);

  if ( OnlyReplayBadFiles ) {
     
    if ( checkIfReplayFileIsGood(std::string(tempOut)) == 1 ) return 1;
  
    else {
      std::ofstream badRuns("badRuns.txt", std::fstream::app);
      badRuns << argv[1] << "\n";
      badRuns.close();
    }
  }

  // Read gain corrections file
  char tempFileGain[500];
  sprintf(tempFileGain, "%s/gain_bismuth_%s.dat",getenv("GAIN_BISMUTH"), argv[1]);
  cout << "... Reading: " << tempFileGain << endl;

  //Read in PMT Quality
  std::vector<Int_t> pmtquality = getPMTQuality(atoi(argv[1])); 

  double fitMean[8], gainCorrection[8];
  ifstream fileGain(tempFileGain);
  for (int i=0; i<8; i++) {
    fileGain >> fitMean[i] >> gainCorrection[i];
    gainCorrection[i] = pmtquality[i] ? gainCorrection[i] : 1.; //This sets the gain to 1 if the bi pulser was bad or something was wrong with the pmt
  }
  cout << "...   PMT E1: " << gainCorrection[0] << endl;
  cout << "...   PMT E2: " << gainCorrection[1] << endl;
  cout << "...   PMT E3: " << gainCorrection[2] << endl;
  cout << "...   PMT E4: " << gainCorrection[3] << endl;
  cout << "...   PMT W1: " << gainCorrection[4] << endl;
  cout << "...   PMT W2: " << gainCorrection[5] << endl;
  cout << "...   PMT W3: " << gainCorrection[6] << endl;
  cout << "...   PMT W4: " << gainCorrection[7] << endl;

  ///////////////// Read in the beam drops /////////////////////////
  
  std::vector < std::vector < Double_t > > beamDropTimes = readBeamDrops(atoi(argv[1]));
  UInt_t numBeamDrops = beamDropTimes.size();
  std::vector < Double_t > deltaTime;
  for ( auto t : beamDropTimes ) deltaTime.push_back( t[1] - t[0] );
  
  
  Double_t beamBurstCut = 0.05; //s

  //sprintf(tempOut, "replay_pass2_%s.root", argv[1]);
  DataTree *t = new DataTree();
  t->makeOutputTree(std::string(tempOut),"pass2");

  char tempIn[500];
  sprintf(tempIn, "%s/replay_pass1_%s.root", getenv("REPLAY_PASS1"),argv[1]);
  //sprintf(tempIn, "../replay_pass1/replay_pass1_%s.root", argv[1]);
  t->setupInputTree(std::string(tempIn),"pass1");
  
  int nEvents = t->getEntries();
  cout << "... Processing nEvents = " << nEvents << endl;
  
  // First I need to calculate the scale factor between the blinded and real time
  // so I can subtract off the proper time from events after beam cuts
  
  t->getEvent(nEvents-1);
  Double_t scaleE = t->TimeE/t->Time;
  Double_t scaleW = t->TimeW/t->Time;
  
  Double_t lastGoodTimeE, lastGoodTimeW, lastGoodTime;
  lastGoodTimeE = lastGoodTimeW = lastGoodTime = 0.;

  Double_t runLengthBlindE = 0.;
  Double_t runLengthBlindW = 0.; 
  Double_t runLengthTrue = 0.;
  Double_t old_runLengthBlindE = 0.;
  Double_t old_runLengthBlindW = 0.; 
  Double_t old_runLengthTrue = 0.;
  
  // Loop over events
  for (Int_t i=0; i<nEvents; i++) {
    t->getEvent(i);
    
    if ( i%10000==0 ) std::cout << i << std::endl;

    // Apply gain correction factors
    t->ScintE.q1 = t->ScintE.q1*gainCorrection[0];
    t->ScintE.q2 = t->ScintE.q2*gainCorrection[1];
    t->ScintE.q3 = t->ScintE.q3*gainCorrection[2];
    t->ScintE.q4 = t->ScintE.q4*gainCorrection[3];
    
    t->ScintW.q1 = t->ScintW.q1*gainCorrection[4];
    t->ScintW.q2 = t->ScintW.q2*gainCorrection[5];
    t->ScintW.q3 = t->ScintW.q3*gainCorrection[6];
    t->ScintW.q4 = t->ScintW.q4*gainCorrection[7];

    t->ScintE_bare.q1 = t->ScintE_bare.q1*gainCorrection[0];
    t->ScintE_bare.q2 = t->ScintE_bare.q2*gainCorrection[1];
    t->ScintE_bare.q3 = t->ScintE_bare.q3*gainCorrection[2];
    t->ScintE_bare.q4 = t->ScintE_bare.q4*gainCorrection[3];
    
    t->ScintW_bare.q1 = t->ScintW_bare.q1*gainCorrection[4];
    t->ScintW_bare.q2 = t->ScintW_bare.q2*gainCorrection[5];
    t->ScintW_bare.q3 = t->ScintW_bare.q3*gainCorrection[6];
    t->ScintW_bare.q4 = t->ScintW_bare.q4*gainCorrection[7];
    
    if ( useBeamCuts ) {
 
      if ( t->Tof < beamBurstCut ) t->badTimeFlag = 1;
      
      for ( UInt_t j=0; j<numBeamDrops; ++j ) {
	// If event was in a beam drop, set the bad time flag and set the event time to 
	// the last good event time. This will create a spike at every beam drop...
	if ( t->Time > beamDropTimes[j][0] && t->Time < beamDropTimes[j][1] ) {
	  t->badTimeFlag = 1;
	  t->Time = lastGoodTime;
	  t->TimeE = lastGoodTimeE;
	  t->TimeW = lastGoodTimeW;
	}
	
	// If it isn't in a drop, subtract all previos drops from event times...
	else if ( t->Time > beamDropTimes[j][0] ) {	
	  t->Time = t->Time - deltaTime[j];
	  t->TimeE = t->TimeE - deltaTime[j]*scaleE;
	  t->TimeW = t->TimeW - deltaTime[j]*scaleW;
	}      
      }
      
      if ( t->badTimeFlag==0 ) {
	lastGoodTimeE = t->TimeE;
	lastGoodTimeW = t->TimeW;
	lastGoodTime = t->Time;
      }
      
      if ( i == (nEvents-1) ) {
	runLengthBlindE = t->TimeE;
	runLengthBlindW = t->TimeW;
	runLengthTrue = t->Time;
	old_runLengthBlindE = t->oldTimeE;
	old_runLengthBlindW = t->oldTimeW;
	old_runLengthTrue = t->oldTime;
      }
    }

    t->fillOutputTree();
  }


  if ( useBeamCuts ) {
    //Now I want to create and store a few pertinent values in a file for later...
    char tempFile[200];
    sprintf(tempFile,"%s/runInfo_%s.dat",getenv("RUN_INFO_FILES"),argv[1]);
    ofstream runInfo(tempFile);
    std::cout << "Writing Info to " << tempFile << std::endl;
    
    runInfo << "RunLengthEast\t" << std::setprecision(9) << runLengthBlindE << std::endl;
    runInfo << "RunLengthWest\t" << std::setprecision(9) << runLengthBlindW << std::endl;
    runInfo << "RunLengthTrue\t" << std::setprecision(9) << runLengthTrue << std::endl;
    runInfo << "UCNMon4Integral\t" << std::setprecision(9) << t->UCN_Mon_4_Rate->Integral("width") << std::endl;
    runInfo << "old_RunLengthEast\t" << std::setprecision(9) << old_runLengthBlindE << std::endl;
    runInfo << "old_RunLengthWest\t" << std::setprecision(9) << old_runLengthBlindW << std::endl;
    runInfo << "old_RunLengthTrue\t" << std::setprecision(9) << old_runLengthTrue << std::endl;
    
    runInfo.close();
  }

  // Write output ntuple
  t->writeOutputFile();
  delete t;

  if ( checkIfReplayFileIsGood(std::string(tempOut)) != 1 ) {

    std::ofstream badRuns("badRuns.txt", std::fstream::app);
    badRuns << argv[1] << "\n";
    badRuns.close();

  }

  return 0;
}

