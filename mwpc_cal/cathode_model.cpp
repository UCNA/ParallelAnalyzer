#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <sstream>

// ROOT libraries
#include <TF1.h>
#include <TH1D.h>
#include <TMath.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TCanvas.h>
#include <TStyle.h>

#include "runInfo.h"
#include "DataTree.hh"
#include "pedestals.h"
#include "MWPCPositionResponse.hh"
#include "peaks.hh"
#include "positionMapHandler.hh"


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

std::vector <int>  readOctetFile(int octet) {

  std::vector <int> runs;
  
  char fileName[200];
  sprintf(fileName,"%s/All_Octets/octet_list_%i.dat",getenv("OCTET_LIST"),octet);
  std::ifstream infile(fileName);
  std::string runTypeHold;
  int runNumberHold;
  int numRuns = 0;
  // Populate the map with the runType and runNumber from this octet
  while (infile >> runTypeHold >> runNumberHold) {

    if (runTypeHold=="A2" || runTypeHold=="A5" || runTypeHold=="A7" || runTypeHold=="A10" || 
	runTypeHold=="B2" || runTypeHold=="B5" || runTypeHold=="B7" || runTypeHold=="B10" )  {
     
      runs.push_back(runNumberHold);

    }
    numRuns++;
  }

  infile.close();
 
  std::cout << "Read in octet file for octet " << octet << "\n";
  return runs;

};


int main(int argc, char *argv[])
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  int runNumber = atoi(argv[1]);

  // Reading in pedestals file to get cathode pedestals
  // Read pedestals file
  char tempFilePed[500];
  int iRun;
  sprintf(tempFilePed, "%s/pedestals_%i.dat", getenv("PEDESTALS"), runNumber);
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
  

  std::vector <Double_t> cathEX_trig(16,0.); // Vectors to hold final model values
  std::vector <Double_t> cathEX_clip(16,0.);
  
  std::vector <Double_t> cathEY_trig(16,0.);
  std::vector <Double_t> cathEY_clip(16,0.);

  std::vector <Double_t> cathWX_trig(16,0.);
  std::vector <Double_t> cathWX_clip(16,0.);

  std::vector <Double_t> cathWY_trig(16,0.);
  std::vector <Double_t> cathWY_clip(16,0.);

  // Load the files
    
  TString infile = TString::Format("%s/replay_pass3_%i.root",getenv("REPLAY_PASS3"),runNumber);
  TFile *data = new TFile(infile, "READ");
  TTree *tData = (TTree*)data->Get("pass3");
  
  infile = TString::Format("%s/beta/revCalSim_%i_Beta.root",getenv("REVCALSIM"),runNumber);
  TFile *sim = new TFile(infile, "READ");
  TTree *tSim = (TTree*)sim->Get("revCalSim");

  Double_t triggerInc = 0.01;
  Double_t clippInc = 0.1;
  
  for ( int c=0; c<16; ++c ) { // loop over all cathodes

    std::cout << "//////////// Cathode "<< c << " /////////////\n";

    Double_t dataTriggerCut = 100.;
    Double_t dataClippingCut = 4090.;
    
    //First we do the threshold determination
    
    //East X
    Double_t dataTotalN = tData->GetEntries("PID==1 && ( Side==0 || ( Type==1 || Type==2 ) )");
    Double_t dataTriggeredN = tData->GetEntries(TString::Format("Cathodes_Ex[%i]>%f && PID==1 && ( Side==0 || ( Type==1 || Type==2 ) )",c,dataTriggerCut).Data());
    Double_t triggerFrac = dataTriggeredN / dataTotalN;
    
    Double_t simTotalN = tSim->GetEntries("PID==1 && ( side==0 || ( type==1 || type==2 ) )");
    Double_t simTriggeredN = simTotalN;
    Double_t simTriggerFrac = 1.;
    
    Double_t simTrigger = 0.;
    
    std::cout << "EX Threshold \n";
    while ( simTriggerFrac > triggerFrac ) {
      simTrigger += triggerInc;
      simTriggeredN = tSim->GetEntries(TString::Format("Cath_EX[%i]>%f && PID==1 && ( side==0 || ( type==1 || type==2 ) )",c,simTrigger).Data());
      simTriggerFrac = simTriggeredN / simTotalN ;
      cathEX_trig[c] = simTrigger;
    }
    
    //East Y
    dataTriggeredN = tData->GetEntries(TString::Format("Cathodes_Ey[%i]>%f && PID==1 && ( Side==0 || ( Type==1 || Type==2 ) )",c,dataTriggerCut).Data());
    triggerFrac = dataTriggeredN / dataTotalN;
    
    simTriggeredN = simTotalN;
    simTriggerFrac = 1.;
    
    simTrigger = 0.;

    std::cout << "EY Threshold \n";
    while ( simTriggerFrac > triggerFrac ) {
      simTrigger += triggerInc;
      simTriggeredN = tSim->GetEntries(TString::Format("Cath_EY[%i]>%f && PID==1 && ( side==0 || ( type==1 || type==2 ) )",c,simTrigger).Data());
      simTriggerFrac = simTriggeredN / simTotalN ;
      cathEY_trig[c] = simTrigger;
    }
    
    //West X
    dataTotalN = tData->GetEntries("PID==1 && ( Side==1 || ( Type==1 || Type==2 ) )");
    dataTriggeredN = tData->GetEntries(TString::Format("Cathodes_Wx[%i]>%f && PID==1 && ( Side==1 || ( Type==1 || Type==2 ) )",c,dataTriggerCut).Data());
    triggerFrac = dataTriggeredN / dataTotalN;
    
    simTotalN = tSim->GetEntries("PID==1 && ( side==1 || ( type==1 || type==2 ) )");
    simTriggeredN = simTotalN;
    simTriggerFrac = 1.;
    
    simTrigger = 0.;
    std::cout << "WX Threshold \n";
    while ( simTriggerFrac > triggerFrac ) {
      simTrigger += triggerInc;
      simTriggeredN = tSim->GetEntries(TString::Format("Cath_WX[%i]>%f && PID==1 && ( side==1 || ( type==1 || type==2 ) )",c,simTrigger).Data());
      simTriggerFrac = simTriggeredN / simTotalN ;
      cathWX_trig[c] = simTrigger;
    }
    
    //West Y
    dataTriggeredN = tData->GetEntries(TString::Format("Cathodes_Wy[%i]>%f && PID==1 && ( Side==1 || ( Type==1 || Type==2 ) )",c,dataTriggerCut).Data());
    triggerFrac = dataTriggeredN / dataTotalN;
    
    simTriggeredN = simTotalN;
    simTriggerFrac = 1.;
    
    simTrigger = 0.;
    std::cout << "WY Threshold\n";
    while ( simTriggerFrac > triggerFrac ) {
      simTrigger += triggerInc;
      simTriggeredN = tSim->GetEntries(TString::Format("Cath_WY[%i]>%f && PID==1 && ( side==1 || ( type==1 || type==2 ) )",c,simTrigger).Data());
      simTriggerFrac = simTriggeredN / simTotalN ;
      cathWY_trig[c] = simTrigger;
    }
    
    
    //Now we do the clipping determination
    
    //East X
    dataClippingCut = 4090. - pedPdc2[16+c]; // pedestal subtracted clipping
    
    dataTotalN = tData->GetEntries(TString::Format("PID==1 && ( Side==0 || ( Type==1 || Type==2 ) ) && Cathodes_Ex[%i]>%f",c,dataTriggerCut).Data());
    Double_t dataClippedN = tData->GetEntries(TString::Format("Cathodes_Ex[%i]>%f && PID==1 && ( Side==0 || ( Type==1 || Type==2 ) )",c,dataClippingCut).Data());
    Double_t ClippedFrac = dataClippedN / dataTotalN;
    
    simTotalN = tSim->GetEntries(TString::Format("PID==1 && ( side==0 || ( type==1 || type==2 ) ) && Cath_EX[%i]>%f",c,cathEX_trig[c]).Data());
    Double_t simClippedN = 0.;
    Double_t simClippedFrac = 0.;
    
    Double_t simClip = 9.;
    std::cout << "EX Clipping\n";// clipped Frac above " << dataClippingCut << ": "
    //	      << dataClippedN << "/" << dataTotalN << " = " << ClippedFrac << "\n";   
    while ( simClippedFrac < ClippedFrac ) {
      simClip -= clippInc;
      simClippedN = tSim->GetEntries(TString::Format("Cath_EX[%i]>%f && PID==1 && ( side==0 || ( type==1 || type==2 ) )",c,simClip).Data());
      simClippedFrac = simClippedN / simTotalN ;
    }
    cathEX_clip[c] = simClip;
    
    // East Y
    dataClippingCut = 4090. - pedPdc2[c]; // pedestal subtracted clipping
    
    dataTotalN = tData->GetEntries(TString::Format("PID==1 && ( Side==0 || ( Type==1 || Type==2 ) ) && Cathodes_Ey[%i]>%f",c,dataTriggerCut).Data());
    dataClippedN = tData->GetEntries(TString::Format("Cathodes_Ey[%i]>%f && PID==1 && ( Side==0 || ( Type==1 || Type==2 ) )",c,dataClippingCut).Data());
    ClippedFrac = dataClippedN / dataTotalN;
    
    simTotalN = tSim->GetEntries(TString::Format("PID==1 && ( side==0 || ( type==1 || type==2 ) ) && Cath_EY[%i]>%f",c,cathEX_trig[c]).Data());
    simClippedN = 0.;
    simClippedFrac = 0.;
    
    simClip = 9.;
    std::cout << "EY Clipping\n";    
    while ( simClippedFrac < ClippedFrac ) {
      simClip -= clippInc;
      simClippedN = tSim->GetEntries(TString::Format("Cath_EY[%i]>%f && PID==1 && ( side==0 || ( type==1 || type==2 ) )",c,simClip).Data());
      simClippedFrac = simClippedN / simTotalN ;
    }
    cathEY_clip[c] = simClip;
    
    // West X
    dataClippingCut = 4090. - pedPadc[16+c]; // pedestal subtracted clipping
    
    dataTotalN = tData->GetEntries(TString::Format("PID==1 && ( Side==1 || ( Type==1 || Type==2 ) ) && Cathodes_Wx[%i]>%f",c,dataTriggerCut).Data());
    dataClippedN = tData->GetEntries(TString::Format("Cathodes_Wx[%i]>%f && PID==1 && ( Side==1 || ( Type==1 || Type==2 ) )",c,dataClippingCut).Data());
    ClippedFrac = dataClippedN / dataTotalN;
    
    simTotalN = tSim->GetEntries(TString::Format("PID==1 && ( side==1 || ( type==1 || type==2 ) ) && Cath_WX[%i]>%f",c,cathEX_trig[c]).Data());
    simClippedN = 0.;
    simClippedFrac = 0.;
    
    simClip = 9.;
    std::cout << "WX Clipping\n";
    while ( simClippedFrac < ClippedFrac ) {
      simClip -= clippInc;
      simClippedN = tSim->GetEntries(TString::Format("Cath_WX[%i]>%f && PID==1 && ( side==1 || ( type==1 || type==2 ) )",c,simClip).Data());
      simClippedFrac = simClippedN / simTotalN ;
    }
    cathWX_clip[c] = simClip;
    
    // West Y
    dataClippingCut = 4090. - pedPadc[c]; // pedestal subtracted clipping
    
    dataTotalN = tData->GetEntries(TString::Format("PID==1 && ( Side==1 || ( Type==1 || Type==2 ) ) && Cathodes_Wy[%i]>%f",c,dataTriggerCut).Data());
    dataClippedN = tData->GetEntries(TString::Format("Cathodes_Wy[%i]>%f && PID==1 && ( Side==1 || ( Type==1 || Type==2 ) )",c,dataClippingCut).Data());
    ClippedFrac = dataClippedN / dataTotalN;
    
    simTotalN = tSim->GetEntries(TString::Format("PID==1 && ( side==1 || ( type==1 || type==2 ) ) && Cath_WY[%i]>%f",c,cathEX_trig[c]).Data());
    simClippedN = 0.;
    simClippedFrac = 0.;
    
    simClip = 9.;
    std::cout << "WY Clipping\n";
    while ( simClippedFrac < ClippedFrac ) {
      simClip -= clippInc;
      simClippedN = tSim->GetEntries(TString::Format("Cath_WY[%i]>%f && PID==1 && ( side==1 || ( type==1 || type==2 ) )",c,simClip).Data());
      simClippedFrac = simClippedN / simTotalN ;
    }
    cathWY_clip[c] = simClip;
    
  }

  
  delete sim;
  delete data;

  std::ofstream ofile(TString::Format("%s/cathode_model/cathode_model_%i.dat",getenv("MWPC_CALIBRATION"),runNumber).Data());

  ofile << "#cathEX_trigg\tcathEX_clip\tcathEY_trigg\tcathEY_clip\tcathWX_trigg\tcathWX_clip\tcathWY_trigg\tcathWY_clip\n";

  for ( int c=0; c<16; ++c ) {

    ofile << cathEX_trig[c] << "\t\t"
	  << cathEX_clip[c] << "\t\t"
	  << cathEY_trig[c] << "\t\t"
	  << cathEY_clip[c] << "\t\t"
	  << cathWX_trig[c] << "\t\t"
	  << cathWX_clip[c] << "\t\t"
	  << cathWY_trig[c] << "\t\t"
	  << cathWY_clip[c];
    if ( c<15 ) ofile << "\n";
  }

  ofile.close();
  
 
  return 0;
}
