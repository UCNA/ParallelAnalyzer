//////////////////////////////////////////////////////////////////
// Code to apply cuts and do pedestal subtraction. 
// Ultimately produces first trees with branches of physical 
// meaning to be further processed in later replays.
//////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <vector>

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
#include "WireChamberResponse.h"
#include "DataTree.hh"
#include "MWPCPositionResponse.hh"

using namespace std;

bool OnlyReplayBadFiles = false;

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

int main(int argc, char *argv[])
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);
  TH1::AddDirectory(kFALSE);

  cout << "Run " << argv[1] << " ..." << endl;

  char tempOut[500];
  sprintf(tempOut, "%s/replay_pass1_%s.root",getenv("REPLAY_PASS1"), argv[1]);
  


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



  // Read cuts file
  char tempFileCuts[500];
  sprintf(tempFileCuts, "%s/cuts_%s.dat", getenv("CUTS"),argv[1]);
  cout << "... Reading: " << tempFileCuts << endl;

  ifstream fileCuts(tempFileCuts);
  fileCuts >> cutBeamBurstTime >> comment;
  fileCuts >> nCutsTimeWindows >> comment;
  if (nCutsTimeWindows > 0) {
    for (int i=0; i<nCutsTimeWindows; i++) {
      fileCuts >> cutTimeWindowLower[i] >> cutTimeWindowUpper[i];
    }
  }
  fileCuts >> cutEastAnode >> comment;
  fileCuts >> cutWestAnode >> comment;
  fileCuts >> cutEastTwoFold >> comment;
  fileCuts >> cutWestTwoFold >> comment;
  fileCuts >> cutEastTopVetoQADC >> comment;
  fileCuts >> cutEastTopVetoTDC >> comment;
  fileCuts >> cutEastDriftTubeTAC >> comment;
  fileCuts >> cutWestDriftTubeTAC >> comment;
  fileCuts >> cutEastBackingVetoQADC >> comment;
  fileCuts >> cutEastBackingVetoTDC >> comment;
  fileCuts >> cutWestBackingVetoQADC >> comment;
  fileCuts >> cutWestBackingVetoTDC >> comment;

  cout << "... Beam Burst T0 Cut: " << cutBeamBurstTime << endl;
  cout << "... Number of Time Windows Cuts: " << nCutsTimeWindows << endl;
  if (nCutsTimeWindows > 0) {
    for (int i=0; i<nCutsTimeWindows; i++) {
      cout << "        [" << cutTimeWindowLower[i] << ", " << cutTimeWindowUpper[i] << "]" << endl;
    }
  }
  cout << "... East MWPC Anode Cut: " << cutEastAnode << endl;
  cout << "... West MWPC Anode Cut: " << cutWestAnode << endl;
  cout << "... East Scintillator Two-Fold Trigger Cut: " << cutEastTwoFold << endl;
  cout << "... West Scintillator Two-Fold Trigger Cut: " << cutWestTwoFold << endl;
  cout << "... East Top Veto QADC Cut: " << cutEastTopVetoQADC << endl;
  cout << "... East Top Veto TDC Cut: " << cutEastTopVetoTDC << endl;
  cout << "... East Drift Tube TAC Cut: " << cutEastDriftTubeTAC << endl;
  cout << "... West Drift Tube TAC Cut: " << cutWestDriftTubeTAC << endl;
  cout << "... East Backing Veto QADC Cut: " << cutEastBackingVetoQADC << endl;
  cout << "... East Backing Veto TDC Cut: " << cutEastBackingVetoTDC << endl;
  cout << "... West Backing Veto QADC Cut: " << cutWestBackingVetoQADC << endl;
  cout << "... West Backing Veto TDC Cut: " << cutWestBackingVetoTDC << endl;

  //First load separate PMT pedestals as produced by the trigger thresholds

  std::vector < std::vector <Double_t> > pmtPedestals = loadPMTpedestals(atoi(argv[1]));
  

  // Read pedestals file
  char tempFilePed[500];
  int iRun;
  sprintf(tempFilePed, "%s/pedestals_%s.dat", getenv("PEDESTALS"), argv[1]);
  cout << "... Reading: " << tempFilePed << endl;

  ifstream filePed(tempFilePed);
  for (int i=0; i<8; i++) {
    filePed >> iRun >> pedQadc[i];   
    pedQadc[i] = pmtPedestals[i][0]; //replace pedQadc[i] with the pedestals separately loaded
  }
  for (int i=0; i<32; i++) {
    filePed >> iRun >> pedPdc2[i];
  }
  for (int i=0; i<32; i++) {
    filePed >> iRun >> pedPadc[i];
  }

  
  filePed >> iRun >> pedPdc30;
  filePed >> iRun >> pedPdc34;

  

  //cout << iRun << " " << pedPdc30 << endl;
  //cout << iRun << " " << pedPdc34 << endl;
  
  // Open output ntuple


  DataTree *t = new DataTree();
  t->makeOutputTree(std::string(tempOut),"pass1");  
  
  
  // Open input ntuple
  char tempIn[500];
  sprintf(tempIn, "%s/full%s.root",getenv("UCNA_RAW_DATA"), argv[1]);

  TFile *fileIn = new TFile(tempIn, "READ");
  TTree *Tin = (TTree*)(fileIn->Get("h1"));

  // Variables
  Tin->SetBranchAddress("Pdc30",  &Pdc30);
  Tin->SetBranchAddress("Pdc34",  &Pdc34);

  Tin->SetBranchAddress("Tdc016", &Tdc016);
  Tin->SetBranchAddress("Tdc017", &Tdc017);
  Tin->SetBranchAddress("Tdc00", &Tdc00);
  Tin->SetBranchAddress("Tdc01", &Tdc01);
  Tin->SetBranchAddress("Tdc02", &Tdc02);
  Tin->SetBranchAddress("Tdc03", &Tdc03);
  Tin->SetBranchAddress("Tdc08", &Tdc08);
  Tin->SetBranchAddress("Tdc09", &Tdc09);
  Tin->SetBranchAddress("Tdc014", &Tdc010); //Note that what used to be TDC10 is now TDC14
  Tin->SetBranchAddress("Tdc011", &Tdc011);

  Tin->SetBranchAddress("Sis00",  &Sis00);

  Tin->SetBranchAddress("Qadc0",  &Qadc[0]);
  Tin->SetBranchAddress("Qadc1",  &Qadc[1]);
  Tin->SetBranchAddress("Qadc2",  &Qadc[2]);
  Tin->SetBranchAddress("Qadc3",  &Qadc[3]);
  Tin->SetBranchAddress("Qadc4",  &Qadc[4]);
  Tin->SetBranchAddress("Qadc5",  &Qadc[5]);
  Tin->SetBranchAddress("Qadc6",  &Qadc[6]);
  Tin->SetBranchAddress("Qadc7",  &Qadc[7]);

  Tin->SetBranchAddress("Pdc20",  &Pdc2[0]);
  Tin->SetBranchAddress("Pdc21",  &Pdc2[1]);
  Tin->SetBranchAddress("Pdc22",  &Pdc2[2]);
  Tin->SetBranchAddress("Pdc23",  &Pdc2[3]);
  Tin->SetBranchAddress("Pdc24",  &Pdc2[4]);
  Tin->SetBranchAddress("Pdc25",  &Pdc2[5]);
  Tin->SetBranchAddress("Pdc26",  &Pdc2[6]);
  Tin->SetBranchAddress("Pdc27",  &Pdc2[7]);
  Tin->SetBranchAddress("Pdc28",  &Pdc2[8]);
  Tin->SetBranchAddress("Pdc29",  &Pdc2[9]);
  Tin->SetBranchAddress("Pdc210", &Pdc2[10]);
  Tin->SetBranchAddress("Pdc211", &Pdc2[11]);
  Tin->SetBranchAddress("Pdc212", &Pdc2[12]);
  Tin->SetBranchAddress("Pdc213", &Pdc2[13]);
  Tin->SetBranchAddress("Pdc214", &Pdc2[14]);
  Tin->SetBranchAddress("Pdc215", &Pdc2[15]);
  Tin->SetBranchAddress("Pdc216", &Pdc2[16]);
  Tin->SetBranchAddress("Pdc217", &Pdc2[17]);
  Tin->SetBranchAddress("Pdc218", &Pdc2[18]);
  Tin->SetBranchAddress("Pdc219", &Pdc2[19]);
  Tin->SetBranchAddress("Pdc220", &Pdc2[20]);
  Tin->SetBranchAddress("Pdc221", &Pdc2[21]);
  Tin->SetBranchAddress("Pdc222", &Pdc2[22]);
  Tin->SetBranchAddress("Pdc223", &Pdc2[23]);
  Tin->SetBranchAddress("Pdc224", &Pdc2[24]);
  Tin->SetBranchAddress("Pdc225", &Pdc2[25]);
  Tin->SetBranchAddress("Pdc226", &Pdc2[26]);
  Tin->SetBranchAddress("Pdc227", &Pdc2[27]);
  Tin->SetBranchAddress("Pdc228", &Pdc2[28]);
  Tin->SetBranchAddress("Pdc229", &Pdc2[29]);
  Tin->SetBranchAddress("Pdc230", &Pdc2[30]);
  Tin->SetBranchAddress("Pdc231", &Pdc2[31]);

  Tin->SetBranchAddress("Padc0",  &Padc[0]);
  Tin->SetBranchAddress("Padc1",  &Padc[1]);
  Tin->SetBranchAddress("Padc2",  &Padc[2]);
  Tin->SetBranchAddress("Padc3",  &Padc[3]);
  Tin->SetBranchAddress("Padc4",  &Padc[4]);
  Tin->SetBranchAddress("Padc5",  &Padc[5]);
  Tin->SetBranchAddress("Padc6",  &Padc[6]);
  Tin->SetBranchAddress("Padc7",  &Padc[7]);
  Tin->SetBranchAddress("Padc8",  &Padc[8]);
  Tin->SetBranchAddress("Padc9",  &Padc[9]);
  Tin->SetBranchAddress("Padc10", &Padc[10]);
  Tin->SetBranchAddress("Padc11", &Padc[11]);
  Tin->SetBranchAddress("Padc12", &Padc[12]);
  Tin->SetBranchAddress("Padc13", &Padc[13]);
  Tin->SetBranchAddress("Padc14", &Padc[14]);
  Tin->SetBranchAddress("Padc15", &Padc[15]);
  Tin->SetBranchAddress("Padc16", &Padc[16]);
  Tin->SetBranchAddress("Padc17", &Padc[17]);
  Tin->SetBranchAddress("Padc18", &Padc[18]);
  Tin->SetBranchAddress("Padc19", &Padc[19]);
  Tin->SetBranchAddress("Padc20", &Padc[20]);
  Tin->SetBranchAddress("Padc21", &Padc[21]);
  Tin->SetBranchAddress("Padc22", &Padc[22]);
  Tin->SetBranchAddress("Padc23", &Padc[23]);
  Tin->SetBranchAddress("Padc24", &Padc[24]);
  Tin->SetBranchAddress("Padc25", &Padc[25]);
  Tin->SetBranchAddress("Padc26", &Padc[26]);
  Tin->SetBranchAddress("Padc27", &Padc[27]);
  Tin->SetBranchAddress("Padc28", &Padc[28]);
  Tin->SetBranchAddress("Padc29", &Padc[29]);
  Tin->SetBranchAddress("Padc30", &Padc[30]);
  Tin->SetBranchAddress("Padc31", &Padc[31]);

  Tin->SetBranchAddress("S83028", &S83028);
  Tin->SetBranchAddress("S8200", &S8200);
  Tin->SetBranchAddress("Clk0", &Clk0);
  Tin->SetBranchAddress("Clk1", &Clk1);
  Tin->SetBranchAddress("Clk2", &Clk2);
  Tin->SetBranchAddress("Clk3", &Clk3);

  Tin->SetBranchAddress("Pdc38",  &Pdc38);
  Tin->SetBranchAddress("Pdc39",  &Pdc39);
  Tin->SetBranchAddress("Pdc310", &Pdc310);
  Tin->SetBranchAddress("Pdc311", &Pdc311);

  Tin->SetBranchAddress("Qadc9",  &Qadc9);
  Tin->SetBranchAddress("Tdc019", &Tdc019);
  Tin->SetBranchAddress("Pdc313", &Pdc313);
  Tin->SetBranchAddress("Pdc315", &Pdc315);
  Tin->SetBranchAddress("Qadc8",  &Qadc8);
  Tin->SetBranchAddress("Tdc018", &Tdc018);
  Tin->SetBranchAddress("Qadc10", &Qadc10);
  Tin->SetBranchAddress("Tdc020", &Tdc020);

  Tin->SetBranchAddress("Number", &Number);
  Tin->SetBranchAddress("Delt0", &Delt0);

  Tin->SetBranchAddress("Evnb0", &Evnb[0]);
  Tin->SetBranchAddress("Evnb1", &Evnb[1]);
  Tin->SetBranchAddress("Evnb2", &Evnb[2]);
  Tin->SetBranchAddress("Evnb3", &Evnb[3]);
  Tin->SetBranchAddress("Evnb4", &Evnb[4]);
  Tin->SetBranchAddress("Bkhf0", &Bkhf[0]);
  Tin->SetBranchAddress("Bkhf1", &Bkhf[1]);
  Tin->SetBranchAddress("Bkhf2", &Bkhf[2]);
  Tin->SetBranchAddress("Bkhf3", &Bkhf[3]);
  Tin->SetBranchAddress("Bkhf4", &Bkhf[4]);


  int nEvents = Tin->GetEntries();
  cout << "... Processing nEvents = " << nEvents << endl;

  //Get length of run from UNBLINDED TIME
  float runLengthBlind[2] = {0.};
  float runLengthTrue = 0.;
  Tin->GetEvent(nEvents-1);
  runLengthTrue = S83028*scalerCountsToTime;
  runLengthBlind[0] = Clk0*scalerCountsToTime;
  runLengthBlind[1] = Clk1*scalerCountsToTime;
  
  //Histograms to hold UCNMonRate
  int binWidth = 10;
  int nbins = (int)(runLengthTrue+20.)/binWidth;
  t->UCN_Mon_1_Rate = new TH1F("UCN_Mon_1_Rate","UCN Mon 1 Rate",nbins, 0., (float)nbins*binWidth);
  t->UCN_Mon_2_Rate = new TH1F("UCN_Mon_2_Rate","UCN Mon 2 Rate",nbins, 0., (float)nbins*binWidth);
  t->UCN_Mon_3_Rate = new TH1F("UCN_Mon_3_Rate","UCN Mon 3 Rate",nbins, 0., (float)nbins*binWidth);
  t->UCN_Mon_4_Rate = new TH1F("UCN_Mon_4_Rate","UCN Mon 4 Rate",nbins, 0., (float)nbins*binWidth);


  // Loop over events
  for (int i=0; i<nEvents; i++) {

    Tin->GetEvent(i);

    if ( i%10000==0 ) std::cout << i << std::endl;
    
    t->xE.nClipped = t->yE.nClipped = t->xW.nClipped = t->yW.nClipped = 0;

    Int_t iSis00 = (int) Sis00;

    // Calculate pedestal-subtracted PMT QADC values
    for (int j=0; j<8; j++) {
      pmt[j] = ((double) Qadc[j]) - pedQadc[j];
    }
    
    // Set cathode values to non-ped subtracted values
    for (int j=0; j<32; j++) {
      
      if (j<16) {
	t->Cathodes_Ey[j] = (double)Pdc2[j];
	t->Cathodes_Wy[j] = (double)Padc[j];
      }
      else {
	t->Cathodes_Ex[j-16] = (double)Pdc2[j];
	t->Cathodes_Wx[j-16] = (double)Padc[j];
      }
    }

    //////////////////////////////////////////////////////////
    // For making a max Cathode signal cut

    double CathMaxCut = 300.; // This is pedestal subtracted
    
    MWPCCathodeHandler cathResp(t->Cathodes_Ex,t->Cathodes_Ey,t->Cathodes_Wx,t->Cathodes_Wy,&pedPdc2[16],&pedPdc2[0],&pedPadc[16],&pedPadc[0]);
    
    //Get the max signal in each plane (ped subtracted)
    double cathMaxEX = cathResp.getMaxSignalEX();
    double cathMaxEY = cathResp.getMaxSignalEY();
    double cathMaxWX = cathResp.getMaxSignalWX();
    double cathMaxWY = cathResp.getMaxSignalWY();
    
    double maxCathSumE = cathMaxEX + cathMaxEY;
    double maxCathSumW = cathMaxWX + cathMaxWY;

    //Also saving one set of positions for making position maps
    cathResp.findAllPositions(true,false);

    
    std::vector<double> posex = cathResp.getPosEX();
    std::vector<double> posey = cathResp.getPosEY();
    std::vector<double> poswx = cathResp.getPosWX();
    std::vector<double> poswy = cathResp.getPosWY();
    
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

      

    ///////////////////////////////////////////////////////////
    
    
    // Calculate pedestal-subtracted MWPC Anode PADC values
    AnodeE = ((double) Pdc30) - pedPdc30;
    AnodeW = ((double) Pdc34) - pedPdc34;

    // UCN monitor events
    bool UCNMonitorTrigger = false;
    
    float time = S83028*scalerCountsToTime;
    if (iSis00==260) {t->UCN_Mon_1_Rate->Fill(time,1./binWidth); UCNMonitorTrigger = true;}
    else if (iSis00==516) {t->UCN_Mon_2_Rate->Fill(time,1./binWidth); UCNMonitorTrigger = true;}
    else if (iSis00==1028) {t->UCN_Mon_3_Rate->Fill(time,1./binWidth); UCNMonitorTrigger = true;}
    else if (iSis00==2052) {t->UCN_Mon_4_Rate->Fill(time,1./binWidth); UCNMonitorTrigger = true;}
    

    // Events with a muon hit
    bool muonHitEast = false;
    bool muonHitWest = false;
    if ( (((double) Qadc8)  > cutEastBackingVetoQADC) ||
         (((double) Tdc018) > cutEastBackingVetoTDC)  ||
         (((double) Qadc9)  > cutEastTopVetoQADC)     ||
         (((double) Tdc019) > cutEastTopVetoTDC)      ||
         (((double) Pdc313) > cutEastDriftTubeTAC) ) {
      muonHitEast = true;
    }
    if ( (((double) Qadc10) > cutWestBackingVetoQADC) ||
         (((double) Tdc020) > cutWestBackingVetoTDC)  ||
         (((double) Pdc315) > cutWestDriftTubeTAC) ) {
      muonHitWest = true;
    }
    
    t->TaggedBackE = (((double) Qadc8)  > cutEastBackingVetoQADC || ((double) Tdc018) > cutEastBackingVetoTDC)?true:false;
    t->TaggedBackW = (((double) Qadc10)  > cutWestBackingVetoQADC || ((double) Tdc020) > cutEastBackingVetoTDC)?true:false;
    t->TaggedTopE = (((double) Qadc9)  > cutEastTopVetoQADC || ((double) Tdc019) > cutEastTopVetoTDC)?true:false;
    t->TaggedTopW = false;
    t->TaggedDriftE = (((double) Pdc313)  > cutEastDriftTubeTAC)?true:false;
    t->TaggedDriftW = (((double) Pdc315)  > cutWestDriftTubeTAC)?true:false;

    t->EastBackADC = (double) Qadc8;
    t->WestBackADC = (double) Qadc10;
    t->EastBackTDC = (double) Tdc018;
    t->WestBackTDC = (double) Tdc020;
    t->EastDriftVetoADC = (double) Pdc313;
    t->WestDriftVetoADC = (double) Pdc315;
    t->EastTopVetoADC = (double) Qadc9;
    t->EastTopVetoTDC = (double) Tdc019;
    

    // LED trigger events
    bool LEDTrigger = false;
    if ( (iSis00 == 128) || (iSis00 == 129) || (iSis00 == 130) || (iSis00 == 131) || (iSis00 == 163) ) {
      LEDTrigger = true;
    }

    // Bi pulser trigger events
    bool bismuthPulser = false;
    if (iSis00 == 32) bismuthPulser = true;

    // Scintillator events
    bool mwpcHitEast = false;
    bool mwpcHitWest = false;
    bool scintillatorHitEast = false;
    bool scintillatorHitWest = false;
    bool scintillatorHitBoth = false;
    bool coincidenceEast = false;
    bool coincidenceWest = false;
    bool scintillatorHitFirstEast = false;
    bool scintillatorHitFirstWest = false;
    bool scintillatorHitBothBad = false;
    bool triggerEast = false;
    bool triggerWest = false;

    t->PassedAnoE = t->PassedAnoW = false;
    t->PassedCathE = t->PassedCathW = false;
    
    if ( ((double) Pdc30)  > cutEastAnode) { t->PassedAnoE=true; }
    if ( ((double) Pdc34)  > cutWestAnode) { t->PassedAnoW=true; }


    // Using the cathode sum of max wires to determine MWPC trigger
    if ( maxCathSumE > CathMaxCut ) { mwpcHitEast = true; t->PassedCathE=true; }
    if ( maxCathSumW > CathMaxCut ) { mwpcHitWest = true; t->PassedCathW=true; }

    double timeEastTwoFold = ((double) Tdc016)*tdcChannelToTime;
    double timeWestTwoFold = ((double) Tdc017)*tdcChannelToTime;
    if (timeEastTwoFold > 0.*tdcChannelToTime) scintillatorHitEast = true;
    if (timeWestTwoFold > 0.*tdcChannelToTime) scintillatorHitWest = true;
    if (scintillatorHitEast && scintillatorHitWest) {
      scintillatorHitBoth = true;
      if ( (timeEastTwoFold > cutEastTwoFold) && (timeWestTwoFold < cutWestTwoFold) ) {
        scintillatorHitFirstEast = true;
      }
      if ( (timeEastTwoFold < cutEastTwoFold) && (timeWestTwoFold > cutWestTwoFold) ) {
        scintillatorHitFirstWest = true;
      }
      if ( (timeEastTwoFold > cutEastTwoFold) && (timeWestTwoFold > cutWestTwoFold) ) {
        scintillatorHitBothBad = true;
      }
      if ( (timeEastTwoFold < cutEastTwoFold) && (timeWestTwoFold < cutWestTwoFold) ) {
        scintillatorHitBothBad = true;
      }
    }

    if (mwpcHitEast && scintillatorHitEast) coincidenceEast = true;
    if (mwpcHitWest && scintillatorHitWest) coincidenceWest = true;


    // Event PID logic

    type = 4; //not an electron event
    side = 2; //No scintillator triggers

    if (UCNMonitorTrigger) PID = 5;
    else if (LEDTrigger) PID = 3;
    else if (muonHitEast || muonHitWest) PID = 2;
    else if (bismuthPulser) PID = 4;
    else if ( (scintillatorHitEast && !mwpcHitEast && !scintillatorHitWest && !mwpcHitWest) ||
              (!scintillatorHitEast && !mwpcHitEast && scintillatorHitWest && !mwpcHitWest) ) { PID = 0; side = scintillatorHitEast ? 0 : 1; }
    else if ( (coincidenceEast && !scintillatorHitWest && !mwpcHitWest) ||
              (!scintillatorHitEast && !mwpcHitEast && coincidenceWest) ||
              (coincidenceEast && coincidenceWest && !scintillatorHitBothBad) ||
              (coincidenceEast && !scintillatorHitWest && mwpcHitWest) ||
              (!scintillatorHitEast && mwpcHitEast && coincidenceWest) ) PID = 1;
    else PID = 6;

    
    if (PID == 1) {
      if (coincidenceEast && !scintillatorHitWest && !mwpcHitWest) {
        type = 0;
        side = 0;
      }
      if (!scintillatorHitEast && !mwpcHitEast && coincidenceWest) {
        type = 0;
        side = 1;
      }
      if (coincidenceEast && coincidenceWest && !scintillatorHitBothBad) {
        type = 1;
        if (scintillatorHitFirstEast) side = 0;
        if (scintillatorHitFirstWest) side = 1;
      }
      if (coincidenceEast && !scintillatorHitWest && mwpcHitWest) {
        type = 2;
        side = 0;
      }
      if (!scintillatorHitEast && mwpcHitEast && coincidenceWest) {
        type = 2;
        side = 1;
      }
    }

    

    
    //////////////////////////////////////////////////////////////////////////////////
    /*  Swank's addition to the UK PA.  (2 of 2)
	this seems like a good spot since its around the position calculation. 
    */

    int xeRC,yeRC,xwRC,ywRC;
    WireChamberResponse * WCR = new WireChamberResponse();

    //changing the length 32 double array to length 16  float array.  using variables defined in WCR. 
    for(int j = 0; j<16; j++)
      {
	WCR->cathex[j]=Pdc2[16+j];
	WCR->cathey[j]=Pdc2[j];
	WCR->cathwx[j]=Padc[16+j];
	WCR->cathwy[j]=Padc[j];
      } 
    
    xeRC=WCR->ResponseType(WCR->cathex);   //response class x co-ordinate, East. 
    yeRC=WCR->ResponseType(WCR->cathey);   //response class y co-ordinate, East. 
    xwRC=WCR->ResponseType(WCR->cathwx);   //response class x co-ordinate, West. 
    ywRC=WCR->ResponseType(WCR->cathwy);   //response class y co-ordinate, West.

    delete WCR;
    /*
      End of swank's addtion (2 of 2)
    */
    ////////////////////////////////////////////////////////////////////////////////////
    
    // Calculate MWPC positions for electron events
    
    
    // Pass Everything to output tree
    t->TriggerNum = (int) Number;
    t->EvtN = i;
    t->Sis00 = iSis00;
    t->DeltaT = ((double)Delt0)*scalerCountsToTime;
    t->Tof = (double) S8200;
    t->TimeE = Clk0*scalerCountsToTime;
    t->TimeW = Clk1*scalerCountsToTime;
    t->Time = S83028*scalerCountsToTime;
    t->badTimeFlag = 0;
    t->oldTimeE = Clk0*scalerCountsToTime;
    t->oldTimeW = Clk1*scalerCountsToTime;
    t->oldTime = S83028*scalerCountsToTime;
    t->TDCE = (double) Tdc016;
    t->TDCW = (double) Tdc017;
    t->TDCE1 = (double) Tdc00;
    t->TDCE2 = (double) Tdc01;
    t->TDCE3 = (double) Tdc02;
    t->TDCE4 = (double) Tdc03;
    t->TDCW1 = (double) Tdc08;
    t->TDCW2 = (double) Tdc09;
    t->TDCW3 = (double) Tdc010;
    t->TDCW4 = (double) Tdc011;
    
    //Wirechambers
    
    t->xE.err = 0.;
    t->xE.rawCenter = 0.;
      
    t->yE.err = 0.;
    t->yE.rawCenter = 0.;
   
    t->xW.err = 0.;
    t->xW.rawCenter = 0.;
   
    t->yW.err = 0.;
    t->yW.rawCenter = 0.;
    
    t->ScintE.q1 = pmt[0];
    t->ScintE.q2 = pmt[1];
    t->ScintE.q3 = pmt[2];
    t->ScintE.q4 = pmt[3]; 
    t->ScintE.e1=t->ScintE.de1=t->ScintE.e2=t->ScintE.de2=t->ScintE.e3=t->ScintE.de3=t->ScintE.e4=t->ScintE.de4=0.;
    t->ScintE.energy=t->ScintE.denergy=0.;
    t->ScintE.nPE1=t->ScintE.nPE2=t->ScintE.nPE3=t->ScintE.nPE4=0.;
    
    t->ScintW.q1 = pmt[4];
    t->ScintW.q2 = pmt[5];
    t->ScintW.q3 = pmt[6];
    t->ScintW.q4 = pmt[7]; 
    t->ScintW.e1=t->ScintW.de1=t->ScintW.e2=t->ScintW.de2=t->ScintW.e3=t->ScintW.de3=t->ScintW.e4=t->ScintW.de4=0.;
    t->ScintW.energy=t->ScintW.denergy=0.;
    t->ScintW.nPE1=t->ScintW.nPE2=t->ScintW.nPE3=t->ScintW.nPE4=0.;
    
    t->EvisE = t->EvisW = 0.;
    
    t->CathSumE =  t->xE.cathSum+t->yE.cathSum;
    t->CathSumW =  t->xW.cathSum+t->yW.cathSum;
    
    t->CathMaxE = (t->xE.maxValue>t->yE.maxValue)?t->yE.maxValue:t->xE.maxValue;
    t->CathMaxW = (t->xW.maxValue>t->yW.maxValue)?t->yW.maxValue:t->xW.maxValue;
    
    t->EMWPC_E = t->EMWPC_W = 0.;
    
    t->AnodeE = AnodeE; // Pedestal subtracted
    t->AnodeW = AnodeW;
    
    t->EvnbGood = t->BkhfGood = true;
    for (Int_t i = 0; i<5; i++) {
      if ((int)Evnb[i]-t->TriggerNum) t->EvnbGood = false;
      if ((int)Bkhf[i]!=17) t->BkhfGood = false;
    }
    
    t->xeRC = xeRC;
    t->yeRC = yeRC;
    t->xwRC = xwRC;
    t->ywRC = ywRC;
    
    t->PID = PID;
    t->Type = type;
    t->Side = side;
    t->ProbIII = 0.;
    t->Erecon = 0.;
    t->old_Erecon = 0.;
    t->gaus_Erecon = 0.;


      /*timeE_BB = Clk2*scalerCountsToTime;
      timeW_BB = Clk3*scalerCountsToTime;
      UBtime = S83028*scalerCountsToTime;
      UBtime_BB = S8200*scalerCountsToTime;
      twoFoldE = Tdc016;
      twoFoldW = Tdc017;*/
  
    t->fillOutputTree();
    
  }
  
  fileIn->Close();

  
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

