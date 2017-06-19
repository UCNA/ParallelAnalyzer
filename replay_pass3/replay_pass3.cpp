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
#include "posMapReader.h"
#include "positionMapHandler.hh"
#include "calibrationTools.hh"
#include "MWPCPositionResponse.hh"



#include "replay_pass2.h"
#include "replay_pass3.h"

bool OnlyReplayBadFiles = false;

std::vector <Double_t> loadEndpointGain(Int_t runNumber) {

  std::vector <Double_t> gain(8,1.);
  std::ifstream infile(TString::Format("%s/EndpointGain/run-%i_epGain.dat",getenv("ENDPOINT_ANALYSIS"),runNumber));

  if ( infile.is_open() ) {
    for ( auto &g : gain ) infile >> g;
  }

  for ( auto g : gain ) std::cout << g << std::endl;
  
  return gain;

}; // TODO: copy all endpoint gain factors to individual runs in BetaManager.py

std::vector < std::vector <Double_t> > loadPMTpedestals(Int_t runNumber) {

  Char_t temp[500];
  std::vector < std::vector < Double_t > > peds (8,std::vector <Double_t> (2,0.));
  sprintf(temp,"%s/PMT_pedestals_%i.dat",getenv("PEDESTALS"),runNumber);
  ifstream infile;
  infile.open(temp);

  Int_t i = 0;
  Int_t run;

  while (infile >> run >> peds[i][0] >> peds[i][1]) { std::cout << "Pedestal " << i << ": " << peds[i][0] << " " << peds[i][1] << std::endl; i++; }
  return peds;

};


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
};

vector <Int_t> getEreconPMTQuality(Int_t runNumber) {
  //Read in PMT quality file
  cout << "Reading in PMT Quality file ...\n";
  vector <Int_t>  pmtQuality (8,0);
  Char_t temp[200];
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


using namespace std;

int main(int argc, char *argv[])
{

  if ( argc<2 || argc>3 ) {
    std::cout << "USAGE: ./replay_pass3.exe [run] [applyEndpointGain = false]\n\n";
    exit(0);
  }
  
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);
  
  int runNumber = atoi(argv[1]);
  bool applyEndpointGain = false;
  if ( argc==3 && ( TString(argv[2])==TString("true") || atoi(argv[2])==1 ) ) applyEndpointGain = true; 

  
  int nPMT = 8;
  int nParams = 3; //takes a quadratic 

  // Run number integer
  cout << "Run " << runNumber << " ..." << endl;

  char tempOut[500];
  sprintf(tempOut, "%s/replay_pass3_%s.root",getenv("REPLAY_PASS3"), argv[1]);
  //sprintf(tempOut, "replay_pass3_%s.root", argv[1]);

  //Check if the file is good already and quit if it is so that we can 
  // only replay the files that are bad...

  if ( OnlyReplayBadFiles ) {
     
    if ( checkIfReplayFileIsGood(std::string(tempOut)) == 1 ) return 1;
  
    else {
      std::ofstream badRuns("badRuns.txt", std::fstream::app);
      badRuns << argv[1] << "\n";
      badRuns.close();
    }
  }

  // Reading in pedestals file to get cathode pedestals
  // Read pedestals file
  char tempFilePed[500];
  int iRun;
  sprintf(tempFilePed, "%s/pedestals_%s.dat", getenv("PEDESTALS"), argv[1]);
  std::cout << "... Reading: " << tempFilePed << std::endl;

  std::ifstream filePed(tempFilePed);
  for (int i=0; i<8; i++) {
    filePed >> iRun >> pedQadc[i];   
  }
  for (int i=0; i<32; i++) {
    filePed >> iRun >> pedPdc2[i];
  }
  for (int i=0; i<32; i++) {
    filePed >> iRun >> pedPadc[i];
  }

  
  filePed >> iRun >> pedPdc30;
  filePed >> iRun >> pedPdc34;
  


  cout << "... Applying Calibration ..." << endl;

  unsigned int calibrationPeriod = getSrcRunPeriod(runNumber);
  
  LinearityCurve linearityCurve(calibrationPeriod,false); //Get the linearity curve  
  EreconParameterization eRecon(runNumber); //Load the simulated relationship between EQ and Etrue
  WirechamberCal mwpcCal(runNumber);        //Load the Wirechamber Calibration

  std::vector <Double_t> epGain(8,1.); // Loading the endpoint gain factors if they are to be used
  if ( applyEndpointGain ) epGain = loadEndpointGain(runNumber);
  
  std::vector <Int_t> pmtQuality = getEreconPMTQuality(runNumber); //Read in PMT quality file
  std::vector <Double_t> alpha = GetAlphaValues(calibrationPeriod); //Get values for nPE/keV...
    
  PositionMap posmap(5.0,50.); //Reading Scintillator position maps
  posmap.readPositionMap( getXeRunPeriod(runNumber), "endpoint" );

  MWPCPositionMap anodeMap(5., 50.);    // Reading Anode position maps
  anodeMap.readMWPCPositionMap( getXeRunPeriodForMWPCmap(runNumber) ,250.,300.); // Using 250-300 keV because that's the most probable range

  
  DataTree *t = new DataTree(); // DataTree structure for input pass2 and output pass3
  t->makeOutputTree(std::string(tempOut),"pass3");  // Open output ntuple

  // Input ntuple
  char tempIn[500];
  sprintf(tempIn, "%s/replay_pass2_%s.root", getenv("REPLAY_PASS2"),argv[1]);
  t->setupInputTree(std::string(tempIn),"pass2");
 
  int nEvents = t->getEntries();
  cout << "... Processing nEvents = " << nEvents << endl;

  vector < vector <Int_t> > gridPoint;
  vector < Double_t > eta;
  vector < Double_t > old_eta;
  vector < Double_t > gaus_eta;


  // Loop over events
  for (int i=0; i<nEvents; i++) {
    t->getEvent(i);

    if ( i%10000==0 ) std::cout << i << std::endl;

    // Only process event if it's an electron
    
    if ( t->PID==1 ) {

      std::vector <double> posex(3,0.);
      std::vector <double> poswx(3,0.);
      std::vector <double> posey(3,0.);
      std::vector <double> poswy(3,0.);

      MWPCCathodeHandler cathResp(t->Cathodes_Ex,t->Cathodes_Ey,t->Cathodes_Wx,t->Cathodes_Wy,&pedPdc2[16],&pedPdc2[0],&pedPadc[16],&pedPadc[0]);
 
  
      //First do the normal way... weighted average of good events, gaus fit of clipped
      cathResp.findAllPositions(true,false);
      
      posex = cathResp.getPosEX();
      posey = cathResp.getPosEY();
      poswx = cathResp.getPosWX();
      poswy = cathResp.getPosWY();

      t->xE.center = posex[0] * positionProjection;
      t->yE.center = posey[0] * positionProjection;
      t->xW.center = poswx[0] * positionProjection;
      t->yW.center = poswy[0] * positionProjection;

      t->xE.width = posex[1] * positionProjection;
      t->yE.width = posey[1] * positionProjection;
      t->xW.width = poswx[1] * positionProjection;
      t->yW.width = poswy[1] * positionProjection;
      
      t->xE.height = posex[2];
      t->yE.height = posey[2];
      t->xW.height = poswx[2];
      t->yW.height = poswy[2];

      t->xE.mult = cathResp.getMultEX();
      t->yE.mult = cathResp.getMultEY();
      t->xW.mult = cathResp.getMultWX();
      t->yW.mult = cathResp.getMultWY();

      t->xE.nClipped = cathResp.getnClippedEX();
      t->yE.nClipped = cathResp.getnClippedEY();
      t->xW.nClipped = cathResp.getnClippedWX();
      t->yW.nClipped = cathResp.getnClippedWY();

      t->xE.maxWire = cathResp.getMaxWireEX();
      t->yE.maxWire = cathResp.getMaxWireEY();
      t->xW.maxWire = cathResp.getMaxWireWX();
      t->yW.maxWire = cathResp.getMaxWireWY();

      t->xE.maxValue = cathResp.getMaxSignalEX();
      t->yE.maxValue = cathResp.getMaxSignalEY();
      t->xW.maxValue = cathResp.getMaxSignalWX();
      t->yW.maxValue = cathResp.getMaxSignalWY();

      t->xE.cathSum = cathResp.getCathSumEX();
      t->yE.cathSum = cathResp.getCathSumEY();
      t->xW.cathSum = cathResp.getCathSumWX();
      t->yW.cathSum = cathResp.getCathSumWY();

      t->CathSumE = t->xE.cathSum + t->yE.cathSum;
      t->CathSumW = t->xW.cathSum + t->yW.cathSum;

      t->CathMaxE = t->xE.maxValue + t->yE.maxValue;
      t->CathMaxW = t->xW.maxValue + t->yW.maxValue;
	    
      t->xE.rawCenter = cathResp.getWirePosEX(t->xE.maxWire);
      t->yE.rawCenter = cathResp.getWirePosEY(t->yE.maxWire);
      t->xW.rawCenter = cathResp.getWirePosWX(t->xW.maxWire);
      t->yW.rawCenter = cathResp.getWirePosWY(t->yW.maxWire);

      
      //Now do all gaussian fits... 
      //cathResp.loadGainFactors(runNumber);

      cathResp.findAllPositions(true,true);

      posex = cathResp.getPosEX();
      posey = cathResp.getPosEY();
      poswx = cathResp.getPosWX();
      poswy = cathResp.getPosWY();

      t->gaus_xE.center = posex[0] * positionProjection;
      t->gaus_yE.center = posey[0] * positionProjection;
      t->gaus_xW.center = poswx[0] * positionProjection;
      t->gaus_yW.center = poswy[0] * positionProjection;

      t->gaus_xE.width = posex[1] * positionProjection;
      t->gaus_yE.width = posey[1] * positionProjection;
      t->gaus_xW.width = poswx[1] * positionProjection;
      t->gaus_yW.width = poswy[1] * positionProjection;
      
      t->gaus_xE.height = posex[2];
      t->gaus_yE.height = posey[2];
      t->gaus_xW.height = poswx[2];
      t->gaus_yW.height = poswy[2];

      t->gaus_xE.mult = cathResp.getMultEX();
      t->gaus_yE.mult = cathResp.getMultEY();
      t->gaus_xW.mult = cathResp.getMultWX();
      t->gaus_yW.mult = cathResp.getMultWY();

      t->gaus_xE.nClipped = cathResp.getnClippedEX();
      t->gaus_yE.nClipped = cathResp.getnClippedEY();
      t->gaus_xW.nClipped = cathResp.getnClippedWX();
      t->gaus_yW.nClipped = cathResp.getnClippedWY();

      t->gaus_xE.maxWire = cathResp.getMaxWireEX();
      t->gaus_yE.maxWire = cathResp.getMaxWireEY();
      t->gaus_xW.maxWire = cathResp.getMaxWireWX();
      t->gaus_yW.maxWire = cathResp.getMaxWireWY();
      
      t->gaus_xE.maxValue = cathResp.getMaxSignalEX();
      t->gaus_yE.maxValue = cathResp.getMaxSignalEY();
      t->gaus_xW.maxValue = cathResp.getMaxSignalWX();
      t->gaus_yW.maxValue = cathResp.getMaxSignalWY();

      t->gaus_xE.cathSum = cathResp.getCathSumEX();
      t->gaus_yE.cathSum = cathResp.getCathSumEY();
      t->gaus_xW.cathSum = cathResp.getCathSumWX();
      t->gaus_yW.cathSum = cathResp.getCathSumWY();

      t->gaus_xE.rawCenter = cathResp.getWirePosEX(t->xE.maxWire);
      t->gaus_yE.rawCenter = cathResp.getWirePosEY(t->yE.maxWire);
      t->gaus_xW.rawCenter = cathResp.getWirePosWX(t->xW.maxWire);
      t->gaus_yW.rawCenter = cathResp.getWirePosWY(t->yW.maxWire);


      
      // Now for all weighted averages...
      cathResp.purgeGainFactors();
      cathResp.findAllPositions(false,false);

      posex = cathResp.getPosEX();
      posey = cathResp.getPosEY();
      poswx = cathResp.getPosWX();
      poswy = cathResp.getPosWY();

      t->old_xE.center = posex[0] * positionProjection;
      t->old_yE.center = posey[0] * positionProjection;
      t->old_xW.center = poswx[0] * positionProjection;
      t->old_yW.center = poswy[0] * positionProjection;

      t->old_xE.width = posex[1] * positionProjection;
      t->old_yE.width = posey[1] * positionProjection;
      t->old_xW.width = poswx[1] * positionProjection;
      t->old_yW.width = poswy[1] * positionProjection;
      
      t->old_xE.height = posex[2];
      t->old_yE.height = posey[2];
      t->old_xW.height = poswx[2];
      t->old_yW.height = poswy[2];

      t->old_xE.mult = cathResp.getMultEX();
      t->old_yE.mult = cathResp.getMultEY();
      t->old_xW.mult = cathResp.getMultWX();
      t->old_yW.mult = cathResp.getMultWY();

      t->old_xE.nClipped = cathResp.getnClippedEX();
      t->old_yE.nClipped = cathResp.getnClippedEY();
      t->old_xW.nClipped = cathResp.getnClippedWX();
      t->old_yW.nClipped = cathResp.getnClippedWY();

      t->old_xE.maxWire = cathResp.getMaxWireEX();
      t->old_yE.maxWire = cathResp.getMaxWireEY();
      t->old_xW.maxWire = cathResp.getMaxWireWX();
      t->old_yW.maxWire = cathResp.getMaxWireWY();

      t->old_xE.maxValue = cathResp.getMaxSignalEX();
      t->old_yE.maxValue = cathResp.getMaxSignalEY();
      t->old_xW.maxValue = cathResp.getMaxSignalWX();
      t->old_yW.maxValue = cathResp.getMaxSignalWY();

      t->old_xE.cathSum = cathResp.getCathSumEX();
      t->old_yE.cathSum = cathResp.getCathSumEY();
      t->old_xW.cathSum = cathResp.getCathSumWX();
      t->old_yW.cathSum = cathResp.getCathSumWY();

      t->old_xE.rawCenter = cathResp.getWirePosEX(t->xE.maxWire);
      t->old_yE.rawCenter = cathResp.getWirePosEY(t->yE.maxWire);
      t->old_xW.rawCenter = cathResp.getWirePosWX(t->xW.maxWire);
      t->old_yW.rawCenter = cathResp.getWirePosWY(t->yW.maxWire);

      
      /////////////////////////////////////////////////////////////
      /////// Now do the energy reconstruction
      /*std::cout << "Event " << i << std::endl;    
      std::cout << "Optimal: " << t->xE.center << "\t" << t->yE.center << "\t" << t->xW.center << "\t" << t->yW.center << "\n"
		<< "Old: " << t->old_xE.center << "\t" << t->old_yE.center << "\t" << t->old_xW.center << "\t" << t->old_yW.center << "\n"
		<< "Gaus: " << t->gaus_xE.center << "\t" << t->gaus_yE.center << "\t" << t->gaus_xW.center << "\t" << t->gaus_yW.center << "\n\n" ;*/

      eta = posmap.getInterpolatedEta(t->xE.center, t->yE.center, t->xW.center, t->yW.center);
      old_eta = posmap.getInterpolatedEta(t->old_xE.center, t->old_yE.center, t->old_xW.center, t->old_yW.center);
      gaus_eta = posmap.getInterpolatedEta(t->gaus_xE.center, t->gaus_yE.center, t->gaus_xW.center, t->gaus_yW.center);
      
      //First calculate old position reconstruction old_Erecon
      t->ScintE.e1 = linearityCurve.applyLinCurve(0,t->ScintE.q1) * epGain[0];
      t->ScintE.e2 = linearityCurve.applyLinCurve(1,t->ScintE.q2) * epGain[1];
      t->ScintE.e3 = linearityCurve.applyLinCurve(2,t->ScintE.q3) * epGain[2];
      t->ScintE.e4 = linearityCurve.applyLinCurve(3,t->ScintE.q4) * epGain[3];
      
      t->ScintE.e1 = ( old_eta[0]>0. && t->ScintE.e1>0. ) ? t->ScintE.e1 / old_eta[0] : 0.;
      t->ScintE.e2 = ( old_eta[1]>0. && t->ScintE.e2>0. ) ? t->ScintE.e2 / old_eta[1] : 0.;
      t->ScintE.e3 = ( old_eta[2]>0. && t->ScintE.e3>0. ) ? t->ScintE.e3 / old_eta[2] : 0.;
      t->ScintE.e4 = ( old_eta[3]>0. && t->ScintE.e4>0. ) ? t->ScintE.e4 / old_eta[3] : 0.;
      
      t->ScintE.nPE1 = old_eta[0] > 0. ? t->ScintE.e1 * old_eta[0] * alpha[0] : 0.;
      t->ScintE.nPE2 = old_eta[1] > 0. ? t->ScintE.e2 * old_eta[1] * alpha[1] : 0.;
      t->ScintE.nPE3 = old_eta[2] > 0. ? t->ScintE.e3 * old_eta[2] * alpha[2] : 0.;
      t->ScintE.nPE4 = old_eta[3] > 0. ? t->ScintE.e4 * old_eta[3] * alpha[3] : 0.;
      
      t->ScintE.de1 = t->ScintE.nPE1 > 0. ? t->ScintE.e1/sqrt(t->ScintE.nPE1) : 0.;
      t->ScintE.de2 = t->ScintE.nPE2 > 0. ? t->ScintE.e2/sqrt(t->ScintE.nPE2) : 0.;
      t->ScintE.de3 = t->ScintE.nPE3 > 0. ? t->ScintE.e3/sqrt(t->ScintE.nPE3) : 0.;
      t->ScintE.de4 = t->ScintE.nPE4 > 0. ? t->ScintE.e4/sqrt(t->ScintE.nPE4) : 0.;
      
      
      t->ScintW.e1 = linearityCurve.applyLinCurve(4,t->ScintW.q1) * epGain[4];
      t->ScintW.e2 = linearityCurve.applyLinCurve(5,t->ScintW.q2) * epGain[5];
      t->ScintW.e3 = linearityCurve.applyLinCurve(6,t->ScintW.q3) * epGain[6];
      t->ScintW.e4 = linearityCurve.applyLinCurve(7,t->ScintW.q4) * epGain[7];
      
      t->ScintW.e1 = ( old_eta[4]>0. && t->ScintW.e1>0. ) ? t->ScintW.e1 / old_eta[4] : 0.;
      t->ScintW.e2 = ( old_eta[5]>0. && t->ScintW.e2>0. ) ? t->ScintW.e2 / old_eta[5] : 0.;
      t->ScintW.e3 = ( old_eta[6]>0. && t->ScintW.e3>0. ) ? t->ScintW.e3 / old_eta[6] : 0.;
      t->ScintW.e4 = ( old_eta[7]>0. && t->ScintW.e4>0. ) ? t->ScintW.e4 / old_eta[7] : 0.;
      
      t->ScintW.nPE1 = old_eta[4] > 0. ? t->ScintW.e1 * old_eta[4] * alpha[4] : 0.;
      t->ScintW.nPE2 = old_eta[5] > 0. ? t->ScintW.e2 * old_eta[5] * alpha[5] : 0.;
      t->ScintW.nPE3 = old_eta[6] > 0. ? t->ScintW.e3 * old_eta[6] * alpha[6] : 0.;
      t->ScintW.nPE4 = old_eta[7] > 0. ? t->ScintW.e4 * old_eta[7] * alpha[7] : 0.;
      
      t->ScintW.de1 = t->ScintW.nPE1 > 0. ? t->ScintW.e1/sqrt(t->ScintW.nPE1) : 0.;
      t->ScintW.de2 = t->ScintW.nPE2 > 0. ? t->ScintW.e2/sqrt(t->ScintW.nPE2) : 0.;
      t->ScintW.de3 = t->ScintW.nPE3 > 0. ? t->ScintW.e3/sqrt(t->ScintW.nPE3) : 0.;
      t->ScintW.de4 = t->ScintW.nPE4 > 0. ? t->ScintW.e4/sqrt(t->ScintW.nPE4) : 0.;
      
      //std::cout << "Made it here" << std::endl;
      
      
      //Calculate the weighted energy on a side
      
      //EAST
      Double_t numer = ( (pmtQuality[0] && t->ScintE.nPE1>0. ? t->ScintE.nPE1 : 0.) +
			 (pmtQuality[1] && t->ScintE.nPE1>0. ? t->ScintE.nPE2 : 0.) + 
			 (pmtQuality[2] && t->ScintE.nPE1>0. ? t->ScintE.nPE3 : 0.) + 
			 (pmtQuality[3] && t->ScintE.nPE1>0. ? t->ScintE.nPE4 : 0.) );
      
      Double_t denom  = ( (pmtQuality[0] && t->ScintE.nPE1>0. ? alpha[0] * old_eta[0] : 0.) +
			  (pmtQuality[1] && t->ScintE.nPE2>0. ? alpha[1] * old_eta[1] : 0.) +
			  (pmtQuality[2] && t->ScintE.nPE3>0. ? alpha[2] * old_eta[2] : 0.) + 
			  (pmtQuality[3] && t->ScintE.nPE4>0. ? alpha[3] * old_eta[3] : 0.) ); 
      
      t->ScintE.energy = t->EvisE = (denom!=0. ? numer/denom : 0.);
      t->ScintE.denergy = (denom!=0. ? sqrt(t->ScintE.energy/denom) : 0.);
      
      //WEST
      numer = denom = 0.;
      
      numer = ( (pmtQuality[4] && t->ScintW.nPE1>0. ? t->ScintW.nPE1 : 0.) +
		(pmtQuality[5] && t->ScintW.nPE1>0. ? t->ScintW.nPE2 : 0.) + 
		(pmtQuality[6] && t->ScintW.nPE1>0. ? t->ScintW.nPE3 : 0.) + 
		(pmtQuality[7] && t->ScintW.nPE1>0. ? t->ScintW.nPE4 : 0.) );
      
      denom  = ( (pmtQuality[4] && t->ScintW.nPE1>0. ? alpha[4] * old_eta[4] : 0.) +
		 (pmtQuality[5] && t->ScintW.nPE2>0. ? alpha[5] * old_eta[5] : 0.) +
		 (pmtQuality[6] && t->ScintW.nPE3>0. ? alpha[6] * old_eta[6] : 0.) + 
		 (pmtQuality[7] && t->ScintW.nPE4>0. ? alpha[7] * old_eta[7] : 0.) ); 
      
      
      t->ScintW.energy = t->EvisW = (denom!=0. ? numer/denom : 0.);
      t->ScintW.denergy = (denom!=0. ? sqrt(t->ScintW.energy/denom) : 0.);
      
      
      // Determine the reconstructed energy
      
      int typeIndex = t->Type==0 ? 0:(t->Type==1 ? 1:2); //for retrieving the parameters from EQ2Etrue
      
      double totalEvis=0.;
      
      if (t->Side==0) {
	totalEvis = t->Type==1 ? (t->EvisE+t->EvisW):t->EvisE;
	if (t->EvisE>0. && totalEvis>0.) {
	  t->old_Erecon = eRecon.getErecon(0,typeIndex,totalEvis);
	}
	else t->old_Erecon=-1.;
      }
      if (t->Side==1) {
	totalEvis = t->Type==1 ? (t->EvisE+t->EvisW):t->EvisW;
	if (t->EvisW>0. && totalEvis>0.) {
	  t->old_Erecon = eRecon.getErecon(1,typeIndex,totalEvis);
	}
	else t->old_Erecon=-1.;
      }
    

      ////////////////////////////////////////////////////////////////////
      
      //First calculate old position reconstruction gaus_Erecon
      t->ScintE.e1 = linearityCurve.applyLinCurve(0,t->ScintE.q1) * epGain[0];
      t->ScintE.e2 = linearityCurve.applyLinCurve(1,t->ScintE.q2) * epGain[1];
      t->ScintE.e3 = linearityCurve.applyLinCurve(2,t->ScintE.q3) * epGain[2];
      t->ScintE.e4 = linearityCurve.applyLinCurve(3,t->ScintE.q4) * epGain[3];
      
      t->ScintE.e1 = ( gaus_eta[0]>0. && t->ScintE.e1>0. ) ? t->ScintE.e1 / gaus_eta[0] : 0.;
      t->ScintE.e2 = ( gaus_eta[1]>0. && t->ScintE.e2>0. ) ? t->ScintE.e2 / gaus_eta[1] : 0.;
      t->ScintE.e3 = ( gaus_eta[2]>0. && t->ScintE.e3>0. ) ? t->ScintE.e3 / gaus_eta[2] : 0.;
      t->ScintE.e4 = ( gaus_eta[3]>0. && t->ScintE.e4>0. ) ? t->ScintE.e4 / gaus_eta[3] : 0.;
      
      t->ScintE.nPE1 = gaus_eta[0] > 0. ? t->ScintE.e1 * gaus_eta[0] * alpha[0] : 0.;
      t->ScintE.nPE2 = gaus_eta[1] > 0. ? t->ScintE.e2 * gaus_eta[1] * alpha[1] : 0.;
      t->ScintE.nPE3 = gaus_eta[2] > 0. ? t->ScintE.e3 * gaus_eta[2] * alpha[2] : 0.;
      t->ScintE.nPE4 = gaus_eta[3] > 0. ? t->ScintE.e4 * gaus_eta[3] * alpha[3] : 0.;
      
      t->ScintE.de1 = t->ScintE.nPE1 > 0. ? t->ScintE.e1/sqrt(t->ScintE.nPE1) : 0.;
      t->ScintE.de2 = t->ScintE.nPE2 > 0. ? t->ScintE.e2/sqrt(t->ScintE.nPE2) : 0.;
      t->ScintE.de3 = t->ScintE.nPE3 > 0. ? t->ScintE.e3/sqrt(t->ScintE.nPE3) : 0.;
      t->ScintE.de4 = t->ScintE.nPE4 > 0. ? t->ScintE.e4/sqrt(t->ScintE.nPE4) : 0.;
      
      
      t->ScintW.e1 = linearityCurve.applyLinCurve(4,t->ScintW.q1) * epGain[4];
      t->ScintW.e2 = linearityCurve.applyLinCurve(5,t->ScintW.q2) * epGain[5];
      t->ScintW.e3 = linearityCurve.applyLinCurve(6,t->ScintW.q3) * epGain[6];
      t->ScintW.e4 = linearityCurve.applyLinCurve(7,t->ScintW.q4) * epGain[7];
      
      t->ScintW.e1 = ( gaus_eta[4]>0. && t->ScintW.e1>0. ) ? t->ScintW.e1 / gaus_eta[4] : 0.;
      t->ScintW.e2 = ( gaus_eta[5]>0. && t->ScintW.e2>0. ) ? t->ScintW.e2 / gaus_eta[5] : 0.;
      t->ScintW.e3 = ( gaus_eta[6]>0. && t->ScintW.e3>0. ) ? t->ScintW.e3 / gaus_eta[6] : 0.;
      t->ScintW.e4 = ( gaus_eta[7]>0. && t->ScintW.e4>0. ) ? t->ScintW.e4 / gaus_eta[7] : 0.;
      
      t->ScintW.nPE1 = gaus_eta[4] > 0. ? t->ScintW.e1 * gaus_eta[4] * alpha[4] : 0.;
      t->ScintW.nPE2 = gaus_eta[5] > 0. ? t->ScintW.e2 * gaus_eta[5] * alpha[5] : 0.;
      t->ScintW.nPE3 = gaus_eta[6] > 0. ? t->ScintW.e3 * gaus_eta[6] * alpha[6] : 0.;
      t->ScintW.nPE4 = gaus_eta[7] > 0. ? t->ScintW.e4 * gaus_eta[7] * alpha[7] : 0.;
      
      t->ScintW.de1 = t->ScintW.nPE1 > 0. ? t->ScintW.e1/sqrt(t->ScintW.nPE1) : 0.;
      t->ScintW.de2 = t->ScintW.nPE2 > 0. ? t->ScintW.e2/sqrt(t->ScintW.nPE2) : 0.;
      t->ScintW.de3 = t->ScintW.nPE3 > 0. ? t->ScintW.e3/sqrt(t->ScintW.nPE3) : 0.;
      t->ScintW.de4 = t->ScintW.nPE4 > 0. ? t->ScintW.e4/sqrt(t->ScintW.nPE4) : 0.;
      
      //std::cout << "Made it here" << std::endl;
      
      
      //Calculate the weighted energy on a side
      
      //EAST
      numer = 0.;
      numer = ( (pmtQuality[0] && t->ScintE.nPE1>0. ? t->ScintE.nPE1 : 0.) +
		(pmtQuality[1] && t->ScintE.nPE2>0. ? t->ScintE.nPE2 : 0.) + 
		(pmtQuality[2] && t->ScintE.nPE3>0. ? t->ScintE.nPE3 : 0.) + 
		(pmtQuality[3] && t->ScintE.nPE4>0. ? t->ScintE.nPE4 : 0.) );
      
      denom = 0.;
      denom  = ( (pmtQuality[0] && t->ScintE.nPE1>0. ? alpha[0] * gaus_eta[0] : 0.) +
		 (pmtQuality[1] && t->ScintE.nPE2>0. ? alpha[1] * gaus_eta[1] : 0.) +
		 (pmtQuality[2] && t->ScintE.nPE3>0. ? alpha[2] * gaus_eta[2] : 0.) + 
		 (pmtQuality[3] && t->ScintE.nPE4>0. ? alpha[3] * gaus_eta[3] : 0.) ); 
      
      t->ScintE.energy = t->EvisE = (denom!=0. ? numer/denom : 0.);
      t->ScintE.denergy = (denom!=0. ? sqrt(t->ScintE.energy/denom) : 0.);
      
      //WEST
      numer = denom = 0.;
      
      numer = ( (pmtQuality[4] && t->ScintW.nPE1>0. ? t->ScintW.nPE1 : 0.) +
		(pmtQuality[5] && t->ScintW.nPE2>0. ? t->ScintW.nPE2 : 0.) + 
		(pmtQuality[6] && t->ScintW.nPE3>0. ? t->ScintW.nPE3 : 0.) + 
		(pmtQuality[7] && t->ScintW.nPE4>0. ? t->ScintW.nPE4 : 0.) );
      
      denom  = ( (pmtQuality[4] && t->ScintW.nPE1>0. ? alpha[4] * gaus_eta[4] : 0.) +
		 (pmtQuality[5] && t->ScintW.nPE2>0. ? alpha[5] * gaus_eta[5] : 0.) +
		 (pmtQuality[6] && t->ScintW.nPE3>0. ? alpha[6] * gaus_eta[6] : 0.) + 
		 (pmtQuality[7] && t->ScintW.nPE4>0. ? alpha[7] * gaus_eta[7] : 0.) ); 
      
      
      t->ScintW.energy = t->EvisW = (denom!=0. ? numer/denom : 0.);
      t->ScintW.denergy = (denom!=0. ? sqrt(t->ScintW.energy/denom) : 0.);
      
      
      // Determine the reconstructed energy
    
      typeIndex = t->Type==0 ? 0:(t->Type==1 ? 1:2); //for retrieving the parameters from EQ2Etrue
    
      totalEvis=0.;
    
      if (t->Side==0) {
	totalEvis = t->Type==1 ? (t->EvisE+t->EvisW):t->EvisE;
	if (t->EvisE>0. && totalEvis>0.) {
	  t->gaus_Erecon = eRecon.getErecon(0,typeIndex,totalEvis);
	}
	else t->gaus_Erecon=-1.;
      }
      if (t->Side==1) {
	totalEvis = t->Type==1 ? (t->EvisE+t->EvisW):t->EvisW;
	if (t->EvisW>0. && totalEvis>0.) {
	  t->gaus_Erecon = eRecon.getErecon(1,typeIndex,totalEvis);
	}
	else t->gaus_Erecon=-1.;
      }
      
      
      
      /////////////////////////////////////////////////////////////
      // Now for the real Erecon and all of the variables that will be saved to file
      
      t->ScintE.e1 = linearityCurve.applyLinCurve(0,t->ScintE.q1) * epGain[0];
      t->ScintE.e2 = linearityCurve.applyLinCurve(1,t->ScintE.q2) * epGain[1];
      t->ScintE.e3 = linearityCurve.applyLinCurve(2,t->ScintE.q3) * epGain[2];
      t->ScintE.e4 = linearityCurve.applyLinCurve(3,t->ScintE.q4) * epGain[3];
      
      t->ScintE.e1 = ( eta[0]>0. ) ? t->ScintE.e1 / eta[0] : 0.;
      t->ScintE.e2 = ( eta[1]>0. ) ? t->ScintE.e2 / eta[1] : 0.;
      t->ScintE.e3 = ( eta[2]>0. ) ? t->ScintE.e3 / eta[2] : 0.;
      t->ScintE.e4 = ( eta[3]>0. ) ? t->ScintE.e4 / eta[3] : 0.;
      
      t->ScintE.nPE1 = eta[0] > 0. ? t->ScintE.e1 * eta[0] * alpha[0] : 0.;
      t->ScintE.nPE2 = eta[1] > 0. ? t->ScintE.e2 * eta[1] * alpha[1] : 0.;
      t->ScintE.nPE3 = eta[2] > 0. ? t->ScintE.e3 * eta[2] * alpha[2] : 0.;
      t->ScintE.nPE4 = eta[3] > 0. ? t->ScintE.e4 * eta[3] * alpha[3] : 0.;
      
      t->ScintE.de1 = t->ScintE.nPE1 > 0. ? t->ScintE.e1/sqrt(t->ScintE.nPE1) : 0.;
      t->ScintE.de2 = t->ScintE.nPE2 > 0. ? t->ScintE.e2/sqrt(t->ScintE.nPE2) : 0.;
      t->ScintE.de3 = t->ScintE.nPE3 > 0. ? t->ScintE.e3/sqrt(t->ScintE.nPE3) : 0.;
      t->ScintE.de4 = t->ScintE.nPE4 > 0. ? t->ScintE.e4/sqrt(t->ScintE.nPE4) : 0.;
      
      
      t->ScintW.e1 = linearityCurve.applyLinCurve(4,t->ScintW.q1) * epGain[4];
      t->ScintW.e2 = linearityCurve.applyLinCurve(5,t->ScintW.q2) * epGain[5];
      t->ScintW.e3 = linearityCurve.applyLinCurve(6,t->ScintW.q3) * epGain[6];
      t->ScintW.e4 = linearityCurve.applyLinCurve(7,t->ScintW.q4) * epGain[7];
      
      t->ScintW.e1 = ( eta[4]>0. ) ? t->ScintW.e1 / eta[4] : 0.;
      t->ScintW.e2 = ( eta[5]>0. ) ? t->ScintW.e2 / eta[5] : 0.;
      t->ScintW.e3 = ( eta[6]>0. ) ? t->ScintW.e3 / eta[6] : 0.;
      t->ScintW.e4 = ( eta[7]>0. ) ? t->ScintW.e4 / eta[7] : 0.;
      
      t->ScintW.nPE1 = eta[4] > 0. ? t->ScintW.e1 * eta[4] * alpha[4] : 0.;
      t->ScintW.nPE2 = eta[5] > 0. ? t->ScintW.e2 * eta[5] * alpha[5] : 0.;
      t->ScintW.nPE3 = eta[6] > 0. ? t->ScintW.e3 * eta[6] * alpha[6] : 0.;
      t->ScintW.nPE4 = eta[7] > 0. ? t->ScintW.e4 * eta[7] * alpha[7] : 0.;
      
      t->ScintW.de1 = t->ScintW.nPE1 > 0. ? t->ScintW.e1/sqrt(t->ScintW.nPE1) : 0.;
      t->ScintW.de2 = t->ScintW.nPE2 > 0. ? t->ScintW.e2/sqrt(t->ScintW.nPE2) : 0.;
      t->ScintW.de3 = t->ScintW.nPE3 > 0. ? t->ScintW.e3/sqrt(t->ScintW.nPE3) : 0.;
      t->ScintW.de4 = t->ScintW.nPE4 > 0. ? t->ScintW.e4/sqrt(t->ScintW.nPE4) : 0.;
      

      // Fill bare scintillator branch with no endpoint gain
      t->ScintE_bare.q1 = t->ScintE.q1;
      t->ScintE_bare.q2 = t->ScintE.q2;
      t->ScintE_bare.q3 = t->ScintE.q3;
      t->ScintE_bare.q4 = t->ScintE.q4;

      t->ScintW_bare.q1 = t->ScintW.q1;
      t->ScintW_bare.q2 = t->ScintW.q2;
      t->ScintW_bare.q3 = t->ScintW.q3;
      t->ScintW_bare.q4 = t->ScintW.q4;

      t->ScintE_bare.e1 = linearityCurve.applyLinCurve(0,t->ScintE_bare.q1);
      t->ScintE_bare.e2 = linearityCurve.applyLinCurve(1,t->ScintE_bare.q2);
      t->ScintE_bare.e3 = linearityCurve.applyLinCurve(2,t->ScintE_bare.q3);
      t->ScintE_bare.e4 = linearityCurve.applyLinCurve(3,t->ScintE_bare.q4);
      
      t->ScintE_bare.e1 = ( eta[0]>0. ) ? t->ScintE_bare.e1 / eta[0] : 0.;
      t->ScintE_bare.e2 = ( eta[1]>0. ) ? t->ScintE_bare.e2 / eta[1] : 0.;
      t->ScintE_bare.e3 = ( eta[2]>0. ) ? t->ScintE_bare.e3 / eta[2] : 0.;
      t->ScintE_bare.e4 = ( eta[3]>0. ) ? t->ScintE_bare.e4 / eta[3] : 0.;
      
      t->ScintE_bare.nPE1 = eta[0] > 0. ? t->ScintE_bare.e1 * eta[0] * alpha[0] : 0.;
      t->ScintE_bare.nPE2 = eta[1] > 0. ? t->ScintE_bare.e2 * eta[1] * alpha[1] : 0.;
      t->ScintE_bare.nPE3 = eta[2] > 0. ? t->ScintE_bare.e3 * eta[2] * alpha[2] : 0.;
      t->ScintE_bare.nPE4 = eta[3] > 0. ? t->ScintE_bare.e4 * eta[3] * alpha[3] : 0.;
      
      t->ScintE_bare.de1 = t->ScintE_bare.nPE1 > 0. ? t->ScintE_bare.e1/sqrt(t->ScintE_bare.nPE1) : 0.;
      t->ScintE_bare.de2 = t->ScintE_bare.nPE2 > 0. ? t->ScintE_bare.e2/sqrt(t->ScintE_bare.nPE2) : 0.;
      t->ScintE_bare.de3 = t->ScintE_bare.nPE3 > 0. ? t->ScintE_bare.e3/sqrt(t->ScintE_bare.nPE3) : 0.;
      t->ScintE_bare.de4 = t->ScintE_bare.nPE4 > 0. ? t->ScintE_bare.e4/sqrt(t->ScintE_bare.nPE4) : 0.;
      
      
      t->ScintW_bare.e1 = linearityCurve.applyLinCurve(4,t->ScintW_bare.q1);
      t->ScintW_bare.e2 = linearityCurve.applyLinCurve(5,t->ScintW_bare.q2);
      t->ScintW_bare.e3 = linearityCurve.applyLinCurve(6,t->ScintW_bare.q3);
      t->ScintW_bare.e4 = linearityCurve.applyLinCurve(7,t->ScintW_bare.q4);
      
      t->ScintW_bare.e1 = ( eta[4]>0. ) ? t->ScintW_bare.e1 / eta[4] : 0.;
      t->ScintW_bare.e2 = ( eta[5]>0. ) ? t->ScintW_bare.e2 / eta[5] : 0.;
      t->ScintW_bare.e3 = ( eta[6]>0. ) ? t->ScintW_bare.e3 / eta[6] : 0.;
      t->ScintW_bare.e4 = ( eta[7]>0. ) ? t->ScintW_bare.e4 / eta[7] : 0.;
      
      t->ScintW_bare.nPE1 = eta[4] > 0. ? t->ScintW_bare.e1 * eta[4] * alpha[4] : 0.;
      t->ScintW_bare.nPE2 = eta[5] > 0. ? t->ScintW_bare.e2 * eta[5] * alpha[5] : 0.;
      t->ScintW_bare.nPE3 = eta[6] > 0. ? t->ScintW_bare.e3 * eta[6] * alpha[6] : 0.;
      t->ScintW_bare.nPE4 = eta[7] > 0. ? t->ScintW_bare.e4 * eta[7] * alpha[7] : 0.;
      
      t->ScintW_bare.de1 = t->ScintW_bare.nPE1 > 0. ? t->ScintW_bare.e1/sqrt(t->ScintW_bare.nPE1) : 0.;
      t->ScintW_bare.de2 = t->ScintW_bare.nPE2 > 0. ? t->ScintW_bare.e2/sqrt(t->ScintW_bare.nPE2) : 0.;
      t->ScintW_bare.de3 = t->ScintW_bare.nPE3 > 0. ? t->ScintW_bare.e3/sqrt(t->ScintW_bare.nPE3) : 0.;
      t->ScintW_bare.de4 = t->ScintW_bare.nPE4 > 0. ? t->ScintW_bare.e4/sqrt(t->ScintW_bare.nPE4) : 0.;

      //std::cout << "Made it here" << std::endl;
      
      
      //Calculate the weighted energy on a side
      
      //EAST
      numer = 0.;
      numer = ( (pmtQuality[0] ? t->ScintE.nPE1 : 0.) +
		(pmtQuality[1] ? t->ScintE.nPE2 : 0.) + 
		(pmtQuality[2] ? t->ScintE.nPE3 : 0.) + 
		(pmtQuality[3] ? t->ScintE.nPE4 : 0.) );
      
      denom = 0.;
      denom  = ( (pmtQuality[0] ? alpha[0] * eta[0] : 0.) +
		 (pmtQuality[1] ? alpha[1] * eta[1] : 0.) +
		 (pmtQuality[2] ? alpha[2] * eta[2] : 0.) + 
		 (pmtQuality[3] ? alpha[3] * eta[3] : 0.) ); 
      
      t->ScintE.energy = t->EvisE = (denom!=0. ? numer/denom : 0.);
      t->ScintE.denergy = (denom!=0. ? sqrt(t->ScintE.energy/denom) : 0.);
      
      //WEST
      numer = denom = 0.;
      
      numer = ( (pmtQuality[4] ? t->ScintW.nPE1 : 0.) +
		(pmtQuality[5] ? t->ScintW.nPE2 : 0.) + 
		(pmtQuality[6] ? t->ScintW.nPE3 : 0.) + 
		(pmtQuality[7] ? t->ScintW.nPE4 : 0.) );
      
      denom  = ( (pmtQuality[4] ? alpha[4] * eta[4] : 0.) +
		 (pmtQuality[5] ? alpha[5] * eta[5] : 0.) +
		 (pmtQuality[6] ? alpha[6] * eta[6] : 0.) + 
		 (pmtQuality[7] ? alpha[7] * eta[7] : 0.) ); 
      
      
      t->ScintW.energy = t->EvisW = (denom!=0. ? numer/denom : 0.);
      t->ScintW.denergy = (denom!=0. ? sqrt(t->ScintW.energy/denom) : 0.);
      
      
      // Determine the reconstructed energy
      
      typeIndex = t->Type==0 ? 0:(t->Type==1 ? 1:2); //for retrieving the parameters from EQ2Etrue
      
      totalEvis=0.;
      
      if (t->Side==0) {
	totalEvis = t->Type==1 ? (t->EvisE+t->EvisW):t->EvisE;
	if (t->EvisE>0. && totalEvis>0.) {
	  t->Erecon = eRecon.getErecon(0,typeIndex,totalEvis);
	}
	else t->Erecon=-1.;
      }
      if (t->Side==1) {
	totalEvis = t->Type==1 ? (t->EvisE+t->EvisW):t->EvisW;
	if (t->EvisW>0. && totalEvis>0.) {
	  t->Erecon = eRecon.getErecon(1,typeIndex,totalEvis);
	}
	else t->Erecon=-1.;
      }

      // Last thing to do for electrons is position correct the anode signal and
      // apply the wirechamber energy calibration

      // Get the position response
      std::vector <Double_t> etaMWPC = anodeMap.getInterpolatedEta(t->xE.center,t->yE.center,
								   t->xW.center,t->yW.center);
      t->AnodeE = t->AnodeE / etaMWPC[0];
      t->AnodeW = t->AnodeW / etaMWPC[1];

      t->EMWPC_E = mwpcCal.applyCal( 0, t->AnodeE ) ;
      t->EMWPC_W = mwpcCal.applyCal( 1, t->AnodeW ) ;

    }

    // write out pedestal subtracted cathode values for all events
    
    for ( int ii = 0; ii<16; ++ii ) {
      t->Cathodes_Ex[ii] = t->Cathodes_Ex[ii] - pedPdc2[ii+16]; 
      t->Cathodes_Ey[ii] = t->Cathodes_Ey[ii] - pedPdc2[ii];
      t->Cathodes_Wx[ii] = t->Cathodes_Wx[ii] - pedPadc[ii+16];
      t->Cathodes_Wy[ii] = t->Cathodes_Wy[ii] - pedPadc[ii];
    }
    
    t->fillOutputTree();
    
  }
  // Write output ntuple
  t->writeOutputFile();
  
  delete t; //Closes files

  if ( checkIfReplayFileIsGood(std::string(tempOut)) != 1 ) {

    std::ofstream badRuns("badRuns.txt", std::fstream::app);
    badRuns << argv[1] << "\n";
    badRuns.close();

  }

  return 0;
}

