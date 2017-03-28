#ifndef DATATREEFLOAT_H
#define DATATREEFLOAT_H

#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <string>
#include <iostream>

struct MWPC_f {
  Float_t center;
  Float_t width;
  Float_t cathSum;
  Float_t maxValue;
  Int_t maxWire;
  Int_t mult;
  Int_t nClipped;
  Int_t err;
  Float_t rawCenter;
  Float_t height;
};

struct Scint_f {
  Float_t q1, q2, q3, q4;
  Float_t e1, de1, e2, de2, e3, de3, e4, de4;
  Float_t energy, denergy; 
  Float_t nPE1, nPE2, nPE3, nPE4;
};

class DataTreeFLOAT {
private:
  TFile *inputFile, *outputFile; 
  TTree *inputTree, *outputTree;

public:
  DataTreeFLOAT(); //Constructor
  ~DataTreeFLOAT();

  void makeOutputTree(std::string outputFile, std::string outputTree); 
  void fillOutputTree() {if (outputTree) outputTree->Fill();};
  void writeOutputFile();
  void setupInputTree(std::string inputFile, std::string inputTree);
  void getEvent(UInt_t N) {inputTree->GetEvent(N);}
  Int_t getEntries() {return inputTree->GetEntriesFast();}

  //List of final variables
  Int_t TriggerNum; //trigger number from full.root file
  Int_t EvtN;       // Event number when read in 
  Int_t Sis00;      // Trigger source
  Float_t DeltaT;  // Time since last event
  Float_t Tof;     // time since last beam burst
  Float_t TimeE;   // blinded time East
  Float_t TimeW;   // blinded time West
  Float_t TDCE, TDCW; // Trigger TDC by side
  Float_t TDCE1, TDCE2, TDCE3, TDCE4; //East ind PMT trigger TDC
  Float_t TDCW1, TDCW2, TDCW3, TDCW4; //West ind PMT trigger TDC
  MWPC_f xE, xW, yE, yW; //MWPC structures
  Float_t Cathodes_Ex[16], Cathodes_Ey[16]; //array of cathode values
  Float_t Cathodes_Wx[16], Cathodes_Wy[16]; //array of cathode values
  Scint_f ScintE, ScintW; //Scintillator variable information
  Float_t EvisE, EvisW; // Visible (quenched) energy on each side
  Float_t CathSumE, CathSumW; //sum of both x-y cathode planes
  Float_t CathMaxE, CathMaxW; //smaller max-value of each plane
  Float_t EMWPC_E, EMWPC_W; // estimated wirechamber energy deposition
  Float_t AnodeE, AnodeW; // Wirechamber anode signals
  Int_t PassedAnoE, PassedAnoW; //Whether or not passed anode threshold
  Int_t PassedCathE, PassedCathW; //Whether or not passed cathode threshold
  Int_t TaggedBackE, TaggedBackW; // Whether or not tagged by backing veto
  Int_t TaggedTopE, TaggedTopW;  // Whether tagged as muon by top veto
                                // which was only on East side so west is 0
  Int_t TaggedDriftE, TaggedDriftW; // Whether tagged by muon veto drift tubes
  Float_t EastBackADC, WestBackADC; // muon backing veto ADC
  Float_t EastBackTDC, WestBackTDC; // muon backing veto TDC
  Float_t EastDriftVetoADC, WestDriftVetoADC; // muon drift veto ADC
  Float_t EastTopVetoADC; //East Top veto ADC (No west)
  Float_t EastTopVetoTDC; // East top veto TDC (No west)
  Int_t EvnbGood, BkhfGood; // DAQ header and footer quality

  Int_t xeRC, yeRC, xwRC, ywRC; //Swank's wirechamber response class variables
  
  Int_t PID; //Particle ID
  Int_t Type; //Event type
  Int_t Side; //Earlier Trigger side
  Float_t ProbIII; //Probability of type 3 event
  Float_t Erecon; //Final reconstructed energy of an event

  // New variabled for beam cuts
  Int_t badTimeFlag; //This is 0 for good events, 1 for bad events
  Float_t oldTimeE;   // blinded time East
  Float_t oldTimeW;   // blinded time West
  Float_t oldTime;    // UNBLINDED time

  TH1F *UCN_Mon_1_Rate, *UCN_Mon_2_Rate, *UCN_Mon_3_Rate, *UCN_Mon_4_Rate;

  //std::string inputTreeName, outputTreeName; 
};
  



#endif


