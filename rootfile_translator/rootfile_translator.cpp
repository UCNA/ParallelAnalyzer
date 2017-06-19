//////////////////////////////////////////////////////////
// Takes replay_pass3_#.root files and creates a spec_#.root
// file to be used in other analyzers. The spec file doesn't
// include all the values as of yet... 
//////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <string>

// ROOT libraries
#include "TRandom3.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TH1D.h>

#include "DataTree.hh"
#include "DataTreeFLOAT.hh"

using namespace std;

int main(int argc, char *argv[])
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  cout << "Run " << argv[1] << " ..." << endl;
  cout << "... Creating spec_" << argv[1] << ".root from replay_pass3 data ..." << endl;

  // Open output ntuple
  char tempOut[500];
  sprintf(tempOut, "%s/spec_%s.root",getenv("UK_SPEC_REPLAY"), argv[1]);
  //sprintf(tempOut, "replay_pass3_%s.root", argv[1]);

  DataTreeFLOAT *spec = new DataTreeFLOAT();
  spec->makeOutputTree(std::string(tempOut),"phys");

  // Input ntuple
  char tempIn[500];
  sprintf(tempIn, "%s/replay_pass3_%s.root", getenv("REPLAY_PASS3"),argv[1]);
  //sprintf(tempIn, "../replay_pass2/replay_pass2_%s.root",argv[1]);

  DataTree *UK = new DataTree();
  UK->setupInputTree(std::string(tempIn),"pass3");

  int nEvents = UK->getEntries();
  cout << "... Processing nEvents = " << nEvents << endl;

  // Loop over events
  for (int i=0; i<nEvents; i++) {
    UK->getEvent(i);

    if (i%10000==0) std::cout << i << std::endl;

    spec->TriggerNum = UK->TriggerNum;
    spec->Sis00 = UK->Sis00;
    spec->DeltaT = (float)(UK->DeltaT);
    spec->Tof = (float)UK->Tof;
    spec->TimeE = (float)UK->TimeE;
    spec->TimeW = (float)UK->TimeW;
    spec->TDCE = (float)UK->TDCE;
    spec->TDCW = (float)UK->TDCW;
    spec->TDCE1 = (float)UK->TDCE1;
    spec->TDCE2 = (float)UK->TDCE2;
    spec->TDCE3 = (float)UK->TDCE3;
    spec->TDCE4 = (float)UK->TDCE4;
    spec->TDCW1 = (float)UK->TDCW1;
    spec->TDCW2 = (float)UK->TDCW2;
    spec->TDCW3 = (float)UK->TDCW3;
    spec->TDCW4 = (float)UK->TDCW4;
    spec->EvisE = (float)UK->EvisE;
    spec->EvisW = (float)UK->EvisW;
    spec->CathSumE = (float)UK->CathSumE;
    spec->CathSumW = (float)UK->CathSumW;
    spec->CathMaxE = (float)UK->CathMaxE;
    spec->CathMaxW = (float)UK->CathMaxW;
    spec->EMWPC_E = (float)UK->EMWPC_E;
    spec->EMWPC_W = (float)UK->EMWPC_W;
    spec->AnodeE = (float)UK->AnodeE;
    spec->AnodeW = (float)UK->AnodeW;
    spec->PassedAnoE = (int)UK->PassedAnoE;
    spec->PassedAnoW = (int)UK->PassedAnoW;
    spec->PassedCathE = (int)UK->PassedCathE;
    spec->PassedCathW = (int)UK->PassedCathW;
    spec->TaggedBackE = (int)UK->TaggedBackE;
    spec->TaggedBackW = (int)UK->TaggedBackW;
    spec->TaggedTopE = (int)UK->TaggedTopE;
    spec->TaggedTopW = (int)UK->TaggedTopW;
    spec->TaggedDriftE = (int)UK->TaggedDriftE;
    spec->TaggedDriftW = (int)UK->TaggedDriftW;
    spec->EastBackADC = (float)UK->EastBackADC;
    spec->WestBackADC = (float)UK->WestBackADC;
    spec->EastBackTDC = (float)UK->EastBackTDC;
    spec->WestBackTDC = (float)UK->WestBackTDC;
    spec->EastDriftVetoADC = (float)UK->EastDriftVetoADC;
    spec->WestDriftVetoADC = (float)UK->WestDriftVetoADC;
    spec->EastTopVetoADC = (float)UK->EastTopVetoADC;
    spec->EastTopVetoTDC = (float)UK->EastTopVetoTDC;
    spec->EvnbGood = (int)UK->EvnbGood;
    spec->BkhfGood = (int)UK->BkhfGood;
    spec->xeRC = UK->xeRC; 
    spec->yeRC = UK->yeRC;
    spec->xwRC = UK->xwRC;
    spec->ywRC = UK->ywRC;
    spec->PID = UK->PID; 
    spec->Type = UK->Type;
    spec->Side = UK->Side;
    spec->ProbIII = (float)UK->ProbIII;
    spec->Erecon = (float)UK->Erecon;

    spec->badTimeFlag = UK->badTimeFlag; 
    spec->oldTimeE = (float)UK->oldTimeE;
    spec->oldTimeW = (float)UK->oldTimeW;
    spec->oldTime = (float)UK->oldTime;

    copy(UK->Cathodes_Ex,UK->Cathodes_Ex+16,spec->Cathodes_Ex);
    copy(UK->Cathodes_Ey,UK->Cathodes_Ey+16,spec->Cathodes_Ey);
    copy(UK->Cathodes_Wx,UK->Cathodes_Wx+16,spec->Cathodes_Wx);
    copy(UK->Cathodes_Wy,UK->Cathodes_Wy+16,spec->Cathodes_Wy);

    spec->xE.center = (float)UK->xE.center;
    spec->xE.width = (float)UK->xE.width;
    spec->xE.cathSum = (float)UK->xE.cathSum;
    spec->xE.maxValue = (float)UK->xE.maxValue;
    spec->xE.maxWire = UK->xE.maxWire;
    spec->xE.mult = UK->xE.mult;
    spec->xE.nClipped = UK->xE.nClipped;
    spec->xE.err = UK->xE.err;
    spec->xE.rawCenter = (float)UK->xE.rawCenter;
    spec->xE.height = (float)UK->xE.height;

    spec->yE.center = (float)UK->yE.center;
    spec->yE.width = (float)UK->yE.width;
    spec->yE.cathSum = (float)UK->yE.cathSum;
    spec->yE.maxValue = (float)UK->yE.maxValue;
    spec->yE.maxWire = UK->yE.maxWire;
    spec->yE.mult = UK->yE.mult;
    spec->yE.nClipped = UK->yE.nClipped;
    spec->yE.err = UK->yE.err;
    spec->yE.rawCenter = (float)UK->yE.rawCenter;
    spec->yE.height = (float)UK->yE.height;

    spec->xW.center = (float)UK->xW.center;
    spec->xW.width = (float)UK->xW.width;
    spec->xW.cathSum = (float)UK->xW.cathSum;
    spec->xW.maxValue = (float)UK->xW.maxValue;
    spec->xW.maxWire = UK->xW.maxWire;
    spec->xW.mult = UK->xW.mult;
    spec->xW.nClipped = UK->xW.nClipped;
    spec->xW.err = UK->xW.err;
    spec->xW.rawCenter = (float)UK->xW.rawCenter;
    spec->xW.height = (float)UK->xW.height;

    spec->yW.center = (float)UK->yW.center;
    spec->yW.width = (float)UK->yW.width;
    spec->yW.cathSum = (float)UK->yW.cathSum;
    spec->yW.maxValue = (float)UK->yW.maxValue;
    spec->yW.maxWire = UK->yW.maxWire;
    spec->yW.mult = UK->yW.mult;
    spec->yW.nClipped = UK->yW.nClipped;
    spec->yW.err = UK->yW.err;
    spec->yW.rawCenter = (float)UK->yW.rawCenter;
    spec->yW.height = (float)UK->yW.height;

    spec->ScintE.q1 = (float)UK->ScintE.q1;
    spec->ScintE.q2 = (float)UK->ScintE.q2;
    spec->ScintE.q3 = (float)UK->ScintE.q3;
    spec->ScintE.q4 = (float)UK->ScintE.q4;
    spec->ScintE.e1 = (float)UK->ScintE.e1;
    spec->ScintE.de1 = (float)UK->ScintE.de1;
    spec->ScintE.e2 = (float)UK->ScintE.e2;
    spec->ScintE.de2 = (float)UK->ScintE.de2;
    spec->ScintE.e3 = (float)UK->ScintE.e3;
    spec->ScintE.de3 = (float)UK->ScintE.de3;
    spec->ScintE.e4 = (float)UK->ScintE.e4;
    spec->ScintE.de4 = (float)UK->ScintE.de4;
    spec->ScintE.energy = (float)UK->ScintE.energy;
    spec->ScintE.denergy = (float)UK->ScintE.denergy;
    spec->ScintE.nPE1 = (float)UK->ScintE.nPE1;
    spec->ScintE.nPE2 = (float)UK->ScintE.nPE2;
    spec->ScintE.nPE3 = (float)UK->ScintE.nPE3;
    spec->ScintE.nPE4 = (float)UK->ScintE.nPE4;

    spec->ScintW.q1 = (float)UK->ScintW.q1;
    spec->ScintW.q2 = (float)UK->ScintW.q2;
    spec->ScintW.q3 = (float)UK->ScintW.q3;
    spec->ScintW.q4 = (float)UK->ScintW.q4;
    spec->ScintW.e1 = (float)UK->ScintW.e1;
    spec->ScintW.de1 = (float)UK->ScintW.de1;
    spec->ScintW.e2 = (float)UK->ScintW.e2;
    spec->ScintW.de2 = (float)UK->ScintW.de2;
    spec->ScintW.e3 = (float)UK->ScintW.e3;
    spec->ScintW.de3 = (float)UK->ScintW.de3;
    spec->ScintW.e4 = (float)UK->ScintW.e4;
    spec->ScintW.de4 = (float)UK->ScintW.de4;
    spec->ScintW.energy = (float)UK->ScintW.energy;
    spec->ScintW.denergy = (float)UK->ScintW.denergy;
    spec->ScintW.nPE1 = (float)UK->ScintW.nPE1;
    spec->ScintW.nPE2 = (float)UK->ScintW.nPE2;
    spec->ScintW.nPE3 = (float)UK->ScintW.nPE3;
    spec->ScintW.nPE4 = (float)UK->ScintW.nPE4;

    spec->fillOutputTree();
  }

  spec->UCN_Mon_1_Rate = UK->UCN_Mon_1_Rate;
  spec->UCN_Mon_2_Rate = UK->UCN_Mon_2_Rate;
  spec->UCN_Mon_3_Rate = UK->UCN_Mon_3_Rate;
  spec->UCN_Mon_4_Rate = UK->UCN_Mon_4_Rate;


  // Write output ntuple
  spec->writeOutputFile();
  
  delete spec; //Closes files
  delete UK;

  return 0;
}


  
