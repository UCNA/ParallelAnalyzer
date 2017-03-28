#include <iostream>
#include <iomanip>
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

#include "fullTreeVariables.h"
#include "pedestals.h"

using namespace std;

struct cuts {
  double cutBeamBurstTime; // Beam Burst T0                                                                                                    
  int nCutsTimeWindows; // Number of Time Window Cuts                                                                                          
  std::vector <Double_t> cutTimeWindowLower;
  std::vector <Double_t> cutTimeWindowUpper;
  Double_t cutEastAnode;
  Double_t cutWestAnode;
  Double_t cutEastTwoFold;
  Double_t cutWestTwoFold;
  Double_t cutEastTopVetoQADC;
  Double_t cutEastTopVetoTDC;
  Double_t cutEastDriftTubeTAC;
  Double_t cutWestDriftTubeTAC;
  Double_t cutEastBackingVetoQADC;
  Double_t cutEastBackingVetoTDC;
  Double_t cutWestBackingVetoQADC;
  Double_t cutWestBackingVetoTDC;
};

void loadCuts(Int_t runNumber, cuts* Cuts);


int main(int argc, char *argv[])
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  // Open output ntuple
  char tempOut[500];
  sprintf(tempOut, "%s/pedestals_%s.root",getenv("PEDESTALS"),argv[1]);
  TFile *fileOut = new TFile(tempOut,"RECREATE");

  Int_t run = atoi(argv[1]);

  // Define output histograms


  // NOTE: These are not the PMT pedestals used in replay pass1!
  //PMT pedestals are calculated when doing the trigger function
  TH1F *his11 = new TH1F("his11", "", 4000,0.5,4000.5); // East PMT #1
  TH1F *his12 = new TH1F("his12", "", 4000,0.5,4000.5); // East PMT #2
  TH1F *his13 = new TH1F("his13", "", 4000,0.5,4000.5); // East PMT #3
  TH1F *his14 = new TH1F("his14", "", 4000,0.5,4000.5); // East PMT #4
  TH1F *his15 = new TH1F("his15", "", 4000,0.5,4000.5); // West PMT #1
  TH1F *his16 = new TH1F("his16", "", 4000,0.5,4000.5); // West PMT #2
  TH1F *his17 = new TH1F("his17", "", 4000,0.5,4000.5); // West PMT #3
  TH1F *his18 = new TH1F("his18", "", 4000,0.5,4000.5); // West PMT #4

  TH1F *his101 = new TH1F("his101", "", 4000,0.5,4000.5); // Pdc2[0]
  TH1F *his102 = new TH1F("his102", "", 4000,0.5,4000.5); // Pdc2[1]
  TH1F *his103 = new TH1F("his103", "", 4000,0.5,4000.5); // Pdc2[2]
  TH1F *his104 = new TH1F("his104", "", 4000,0.5,4000.5); // Pdc2[3]
  TH1F *his105 = new TH1F("his105", "", 4000,0.5,4000.5); // Pdc2[4]
  TH1F *his106 = new TH1F("his106", "", 4000,0.5,4000.5); // Pdc2[5]
  TH1F *his107 = new TH1F("his107", "", 4000,0.5,4000.5); // Pdc2[6]
  TH1F *his108 = new TH1F("his108", "", 4000,0.5,4000.5); // Pdc2[7]
  TH1F *his109 = new TH1F("his109", "", 4000,0.5,4000.5); // Pdc2[8]
  TH1F *his110 = new TH1F("his110", "", 4000,0.5,4000.5); // Pdc2[9]
  TH1F *his111 = new TH1F("his111", "", 4000,0.5,4000.5); // Pdc2[10]
  TH1F *his112 = new TH1F("his112", "", 4000,0.5,4000.5); // Pdc2[11]
  TH1F *his113 = new TH1F("his113", "", 4000,0.5,4000.5); // Pdc2[12]
  TH1F *his114 = new TH1F("his114", "", 4000,0.5,4000.5); // Pdc2[13]
  TH1F *his115 = new TH1F("his115", "", 4000,0.5,4000.5); // Pdc2[14]
  TH1F *his116 = new TH1F("his116", "", 4000,0.5,4000.5); // Pdc2[15]
  TH1F *his117 = new TH1F("his117", "", 4000,0.5,4000.5); // Pdc2[16]
  TH1F *his118 = new TH1F("his118", "", 4000,0.5,4000.5); // Pdc2[17]
  TH1F *his119 = new TH1F("his119", "", 4000,0.5,4000.5); // Pdc2[18]
  TH1F *his120 = new TH1F("his120", "", 4000,0.5,4000.5); // Pdc2[19]
  TH1F *his121 = new TH1F("his121", "", 4000,0.5,4000.5); // Pdc2[20]
  TH1F *his122 = new TH1F("his122", "", 4000,0.5,4000.5); // Pdc2[21]
  TH1F *his123 = new TH1F("his123", "", 4000,0.5,4000.5); // Pdc2[22]
  TH1F *his124 = new TH1F("his124", "", 4000,0.5,4000.5); // Pdc2[23]
  TH1F *his125 = new TH1F("his125", "", 4000,0.5,4000.5); // Pdc2[24]
  TH1F *his126 = new TH1F("his126", "", 4000,0.5,4000.5); // Pdc2[25]
  TH1F *his127 = new TH1F("his127", "", 4000,0.5,4000.5); // Pdc2[26]
  TH1F *his128 = new TH1F("his128", "", 4000,0.5,4000.5); // Pdc2[27]
  TH1F *his129 = new TH1F("his129", "", 4000,0.5,4000.5); // Pdc2[28]
  TH1F *his130 = new TH1F("his130", "", 4000,0.5,4000.5); // Pdc2[29]
  TH1F *his131 = new TH1F("his131", "", 4000,0.5,4000.5); // Pdc2[30]
  TH1F *his132 = new TH1F("his132", "", 4000,0.5,4000.5); // Pdc2[31]

  TH1F *his201 = new TH1F("his201", "", 4000,0.5,4000.5); // Padc2[0]
  TH1F *his202 = new TH1F("his202", "", 4000,0.5,4000.5); // Padc2[1]
  TH1F *his203 = new TH1F("his203", "", 4000,0.5,4000.5); // Padc2[2]
  TH1F *his204 = new TH1F("his204", "", 4000,0.5,4000.5); // Padc2[3]
  TH1F *his205 = new TH1F("his205", "", 4000,0.5,4000.5); // Padc2[4]
  TH1F *his206 = new TH1F("his206", "", 4000,0.5,4000.5); // Padc2[5]
  TH1F *his207 = new TH1F("his207", "", 4000,0.5,4000.5); // Padc2[6]
  TH1F *his208 = new TH1F("his208", "", 4000,0.5,4000.5); // Padc2[7]
  TH1F *his209 = new TH1F("his209", "", 4000,0.5,4000.5); // Padc2[8]
  TH1F *his210 = new TH1F("his210", "", 4000,0.5,4000.5); // Padc2[9]
  TH1F *his211 = new TH1F("his211", "", 4000,0.5,4000.5); // Padc2[10]
  TH1F *his212 = new TH1F("his212", "", 4000,0.5,4000.5); // Padc2[11]
  TH1F *his213 = new TH1F("his213", "", 4000,0.5,4000.5); // Padc2[12]
  TH1F *his214 = new TH1F("his214", "", 4000,0.5,4000.5); // Padc2[13]
  TH1F *his215 = new TH1F("his215", "", 4000,0.5,4000.5); // Padc2[14]
  TH1F *his216 = new TH1F("his216", "", 4000,0.5,4000.5); // Padc2[15]
  TH1F *his217 = new TH1F("his217", "", 4000,0.5,4000.5); // Padc2[16]
  TH1F *his218 = new TH1F("his218", "", 4000,0.5,4000.5); // Padc2[17]
  TH1F *his219 = new TH1F("his219", "", 4000,0.5,4000.5); // Padc2[18]
  TH1F *his220 = new TH1F("his220", "", 4000,0.5,4000.5); // Padc2[19]
  TH1F *his221 = new TH1F("his221", "", 4000,0.5,4000.5); // Padc2[20]
  TH1F *his222 = new TH1F("his222", "", 4000,0.5,4000.5); // Padc2[21]
  TH1F *his223 = new TH1F("his223", "", 4000,0.5,4000.5); // Padc2[22]
  TH1F *his224 = new TH1F("his224", "", 4000,0.5,4000.5); // Padc2[23]
  TH1F *his225 = new TH1F("his225", "", 4000,0.5,4000.5); // Padc2[24]
  TH1F *his226 = new TH1F("his226", "", 4000,0.5,4000.5); // Padc2[25]
  TH1F *his227 = new TH1F("his227", "", 4000,0.5,4000.5); // Padc2[26]
  TH1F *his228 = new TH1F("his228", "", 4000,0.5,4000.5); // Padc2[27]
  TH1F *his229 = new TH1F("his229", "", 4000,0.5,4000.5); // Padc2[28]
  TH1F *his230 = new TH1F("his230", "", 4000,0.5,4000.5); // Padc2[29]
  TH1F *his231 = new TH1F("his231", "", 4000,0.5,4000.5); // Padc2[30]
  TH1F *his232 = new TH1F("his232", "", 4000,0.5,4000.5); // Padc2[31]
  
  TH1F *his300 = new TH1F("his300", "", 4000,0.5,4000.5); // Pdc30
  TH1F *his301 = new TH1F("his301", "", 4000,0.5,4000.5); // Pdc34

  // Open input ntuple
  char tempIn[500];
  sprintf(tempIn, "/extern/mabrow05/ucna/rawdata/full%s.root", argv[1]);

  //ADDING IN TDC VALUES
  Float_t TdcE1, TdcE2, TdcE3, TdcE4, TdcW1, TdcW2, TdcW3, TdcW4;
  Float_t TdcEast2fold, TdcWest2fold;

  TFile *fileIn = new TFile(tempIn, "READ");
  TTree *Tin = (TTree*)(fileIn->Get("h1"));

  // Variables
  Tin->SetBranchAddress("Sis00",  &Sis00);

  Tin->SetBranchAddress("Tdc016", &TdcEast2fold);
  Tin->SetBranchAddress("Tdc017", &TdcWest2fold);
  Tin->SetBranchAddress("Tdc00", &TdcE1);
  Tin->SetBranchAddress("Tdc01", &TdcE2);
  Tin->SetBranchAddress("Tdc02", &TdcE3);
  Tin->SetBranchAddress("Tdc03", &TdcE4);
  Tin->SetBranchAddress("Tdc08", &TdcW1);
  Tin->SetBranchAddress("Tdc09", &TdcW2);
  Tin->SetBranchAddress("Tdc014", &TdcW3);
  Tin->SetBranchAddress("Tdc011", &TdcW4);

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

  Tin->SetBranchAddress("Pdc30", &Pdc30);
  Tin->SetBranchAddress("Pdc34", &Pdc34);

  //Load cuts for this run
  cuts Cuts;

  loadCuts(atoi(argv[1]),&Cuts);

  
  int nEvents = Tin->GetEntries();
  cout << "Run " << argv[1] << " ..." << endl;
  cout << "... Processing nEvents = " << nEvents << endl;


  // Loop over events
  for (int i=0; i<nEvents; i++) {

    Tin->GetEvent(i);
    int iSis00 = (int) Sis00;

    bool triggerEast = false;
    bool triggerWest = false;
    bool triggerBoth = false;
    bool triggerUCN  = false;
    bool triggerBiPulser = false;

    if (iSis00 == 1) triggerEast = true;
    else if (iSis00 == 2) triggerWest = true;
    else if (iSis00 == 3) triggerBoth = true;
    else if (iSis00 == 32) triggerBiPulser = true;
    else if (iSis00 == 260 || iSis00 == 516 || iSis00 == 1028 || iSis00 == 2052) triggerUCN = true;

    //Computing  pedestals
    if ( triggerBiPulser ) { // Choosing Bi Pulser events

      //East
      if (TdcE1<0.0001 && TdcE2<0.0001 && TdcE3<0.0001 && TdcE4<0.0001) {

	// East PMTs
	//his11->Fill(Qadc[0]);
	//his12->Fill(Qadc[1]);
	//his13->Fill(Qadc[2]);
	//his14->Fill(Qadc[3]);

	//East Wirechambers
	his101->Fill(Pdc2[0]);
	his102->Fill(Pdc2[1]);
	his103->Fill(Pdc2[2]);
	his104->Fill(Pdc2[3]);
	his105->Fill(Pdc2[4]);
	his106->Fill(Pdc2[5]);
	his107->Fill(Pdc2[6]);
	his108->Fill(Pdc2[7]);
	his109->Fill(Pdc2[8]);
	his110->Fill(Pdc2[9]);
	his111->Fill(Pdc2[10]);
	his112->Fill(Pdc2[11]);
	his113->Fill(Pdc2[12]);
	his114->Fill(Pdc2[13]);
	his115->Fill(Pdc2[14]);
	his116->Fill(Pdc2[15]);
	his117->Fill(Pdc2[16]);
	his118->Fill(Pdc2[17]);
	his119->Fill(Pdc2[18]);
	his120->Fill(Pdc2[19]);
	his121->Fill(Pdc2[20]);
	his122->Fill(Pdc2[21]);
	his123->Fill(Pdc2[22]);
	his124->Fill(Pdc2[23]);
	his125->Fill(Pdc2[24]);
	his126->Fill(Pdc2[25]);
	his127->Fill(Pdc2[26]);
	his128->Fill(Pdc2[27]);
	his129->Fill(Pdc2[28]);
	his130->Fill(Pdc2[29]);
	his131->Fill(Pdc2[30]);
	his132->Fill(Pdc2[31]);

	his300->Fill(Pdc30);



      }

      //West 
      if (TdcW1<0.0001 && TdcW2<0.0001 && TdcW3<0.0001 && TdcW4<0.0001) {

	//West PMTs
	//his15->Fill(Qadc[4]);
	//his16->Fill(Qadc[5]);
	//his17->Fill(Qadc[6]);
	//his18->Fill(Qadc[7]);

	//West Wirechamber
	his201->Fill(Padc[0]);
	his202->Fill(Padc[1]);
	his203->Fill(Padc[2]);
	his204->Fill(Padc[3]);
	his205->Fill(Padc[4]);
	his206->Fill(Padc[5]);
	his207->Fill(Padc[6]);
	his208->Fill(Padc[7]);
	his209->Fill(Padc[8]);
	his210->Fill(Padc[9]);
	his211->Fill(Padc[10]);
	his212->Fill(Padc[11]);
	his213->Fill(Padc[12]);
	his214->Fill(Padc[13]);
	his215->Fill(Padc[14]);
	his216->Fill(Padc[15]);
	his217->Fill(Padc[16]);
	his218->Fill(Padc[17]);
	his219->Fill(Padc[18]);
	his220->Fill(Padc[19]);
	his221->Fill(Padc[20]);
	his222->Fill(Padc[21]);
	his223->Fill(Padc[22]);
	his224->Fill(Padc[23]);
	his225->Fill(Padc[24]);
	his226->Fill(Padc[25]);
	his227->Fill(Padc[26]);
	his228->Fill(Padc[27]);
	his229->Fill(Padc[28]);
	his230->Fill(Padc[29]);
	his231->Fill(Padc[30]);
	his232->Fill(Padc[31]);

	his301->Fill(Pdc34);
      }
      
    }

    // Now handle PMT pedestals using electron events
    else if ( triggerEast || triggerWest || triggerBoth ) {

      if ( TdcE1<0.000001 ) his11->Fill(Qadc[0]);
      if ( TdcE2<0.000001 ) his12->Fill(Qadc[1]);
      if ( TdcE3<0.000001 ) his13->Fill(Qadc[2]);
      if ( TdcE4<0.000001 ) his14->Fill(Qadc[3]);

      if ( TdcW1<0.000001 ) his15->Fill(Qadc[4]);
      if ( TdcW2<0.000001 ) his16->Fill(Qadc[5]);
      if ( TdcW3<0.000001 ) his17->Fill(Qadc[6]);
      if ( TdcW4<0.000001 ) his18->Fill(Qadc[7]);	

    }

  }


  pedQadc[0] = his11->GetXaxis()->GetBinCenter(his11->GetMaximumBin()); his11->GetXaxis()->SetRangeUser(pedQadc[0]-15., pedQadc[0]+15.);
  pedQadc[1] = his12->GetXaxis()->GetBinCenter(his12->GetMaximumBin()); his12->GetXaxis()->SetRangeUser(pedQadc[1]-15., pedQadc[1]+15.);
  pedQadc[2] = his13->GetXaxis()->GetBinCenter(his13->GetMaximumBin()); his13->GetXaxis()->SetRangeUser(pedQadc[2]-15., pedQadc[2]+15.);
  pedQadc[3] = his14->GetXaxis()->GetBinCenter(his14->GetMaximumBin()); his14->GetXaxis()->SetRangeUser(pedQadc[3]-15., pedQadc[3]+15.);
  pedQadc[4] = his15->GetXaxis()->GetBinCenter(his15->GetMaximumBin()); his15->GetXaxis()->SetRangeUser(pedQadc[4]-15., pedQadc[4]+15.);
  pedQadc[5] = his16->GetXaxis()->GetBinCenter(his16->GetMaximumBin()); his16->GetXaxis()->SetRangeUser(pedQadc[5]-15., pedQadc[5]+15.);
  pedQadc[6] = his17->GetXaxis()->GetBinCenter(his17->GetMaximumBin()); his17->GetXaxis()->SetRangeUser(pedQadc[6]-15., pedQadc[6]+15.);
  pedQadc[7] = his18->GetXaxis()->GetBinCenter(his18->GetMaximumBin()); his18->GetXaxis()->SetRangeUser(pedQadc[7]-15., pedQadc[7]+15.);

  pedPdc2[0]  = his101->GetXaxis()->GetBinCenter(his101->GetMaximumBin()); his101->GetXaxis()->SetRangeUser(pedPdc2[0]-50., pedPdc2[0]+50.);
  pedPdc2[1]  = his102->GetXaxis()->GetBinCenter(his102->GetMaximumBin()); his102->GetXaxis()->SetRangeUser(pedPdc2[1]-50., pedPdc2[1]+50.);
  pedPdc2[2]  = his103->GetXaxis()->GetBinCenter(his103->GetMaximumBin()); his103->GetXaxis()->SetRangeUser(pedPdc2[2]-50., pedPdc2[2]+50.);
  pedPdc2[3]  = his104->GetXaxis()->GetBinCenter(his104->GetMaximumBin()); his104->GetXaxis()->SetRangeUser(pedPdc2[3]-50., pedPdc2[3]+50.);
  pedPdc2[4]  = his105->GetXaxis()->GetBinCenter(his105->GetMaximumBin()); his105->GetXaxis()->SetRangeUser(pedPdc2[4]-50., pedPdc2[4]+50.);
  pedPdc2[5]  = his106->GetXaxis()->GetBinCenter(his106->GetMaximumBin()); his106->GetXaxis()->SetRangeUser(pedPdc2[5]-50., pedPdc2[5]+50.);
  pedPdc2[6]  = his107->GetXaxis()->GetBinCenter(his107->GetMaximumBin()); his107->GetXaxis()->SetRangeUser(pedPdc2[6]-50., pedPdc2[6]+50.);
  pedPdc2[7]  = his108->GetXaxis()->GetBinCenter(his108->GetMaximumBin()); his108->GetXaxis()->SetRangeUser(pedPdc2[7]-50., pedPdc2[7]+50.);
  pedPdc2[8]  = his109->GetXaxis()->GetBinCenter(his109->GetMaximumBin()); his109->GetXaxis()->SetRangeUser(pedPdc2[8]-50., pedPdc2[8]+50.);
  pedPdc2[9]  = his110->GetXaxis()->GetBinCenter(his110->GetMaximumBin()); his110->GetXaxis()->SetRangeUser(pedPdc2[9]-50., pedPdc2[9]+50.);
  pedPdc2[10] = his111->GetXaxis()->GetBinCenter(his111->GetMaximumBin()); his111->GetXaxis()->SetRangeUser(pedPdc2[10]-50., pedPdc2[10]+50.);
  pedPdc2[11] = his112->GetXaxis()->GetBinCenter(his112->GetMaximumBin()); his112->GetXaxis()->SetRangeUser(pedPdc2[11]-50., pedPdc2[11]+50.);
  pedPdc2[12] = his113->GetXaxis()->GetBinCenter(his113->GetMaximumBin()); his113->GetXaxis()->SetRangeUser(pedPdc2[12]-50., pedPdc2[12]+50.);
  pedPdc2[13] = his114->GetXaxis()->GetBinCenter(his114->GetMaximumBin()); his114->GetXaxis()->SetRangeUser(pedPdc2[13]-50., pedPdc2[13]+50.);
  pedPdc2[14] = his115->GetXaxis()->GetBinCenter(his115->GetMaximumBin()); his115->GetXaxis()->SetRangeUser(pedPdc2[14]-50., pedPdc2[14]+50.);
  pedPdc2[15] = his116->GetXaxis()->GetBinCenter(his116->GetMaximumBin()); his116->GetXaxis()->SetRangeUser(pedPdc2[15]-50., pedPdc2[15]+50.);
  pedPdc2[16] = his117->GetXaxis()->GetBinCenter(his117->GetMaximumBin()); his117->GetXaxis()->SetRangeUser(pedPdc2[16]-50., pedPdc2[16]+50.);
  pedPdc2[17] = his118->GetXaxis()->GetBinCenter(his118->GetMaximumBin()); his118->GetXaxis()->SetRangeUser(pedPdc2[17]-50., pedPdc2[17]+50.);
  pedPdc2[18] = his119->GetXaxis()->GetBinCenter(his119->GetMaximumBin()); his119->GetXaxis()->SetRangeUser(pedPdc2[18]-50., pedPdc2[18]+50.);
  pedPdc2[19] = his120->GetXaxis()->GetBinCenter(his120->GetMaximumBin()); his120->GetXaxis()->SetRangeUser(pedPdc2[19]-50., pedPdc2[19]+50.);
  pedPdc2[20] = his121->GetXaxis()->GetBinCenter(his121->GetMaximumBin()); his121->GetXaxis()->SetRangeUser(pedPdc2[20]-50., pedPdc2[20]+50.);
  pedPdc2[21] = his122->GetXaxis()->GetBinCenter(his122->GetMaximumBin()); his122->GetXaxis()->SetRangeUser(pedPdc2[21]-50., pedPdc2[21]+50.);
  pedPdc2[22] = his123->GetXaxis()->GetBinCenter(his123->GetMaximumBin()); his123->GetXaxis()->SetRangeUser(pedPdc2[22]-50., pedPdc2[22]+50.);
  pedPdc2[23] = his124->GetXaxis()->GetBinCenter(his124->GetMaximumBin()); his124->GetXaxis()->SetRangeUser(pedPdc2[23]-50., pedPdc2[23]+50.);
  pedPdc2[24] = his125->GetXaxis()->GetBinCenter(his125->GetMaximumBin()); his125->GetXaxis()->SetRangeUser(pedPdc2[24]-50., pedPdc2[24]+50.);
  pedPdc2[25] = his126->GetXaxis()->GetBinCenter(his126->GetMaximumBin()); his126->GetXaxis()->SetRangeUser(pedPdc2[25]-50., pedPdc2[25]+50.);
  pedPdc2[26] = his127->GetXaxis()->GetBinCenter(his127->GetMaximumBin()); his127->GetXaxis()->SetRangeUser(pedPdc2[26]-50., pedPdc2[26]+50.);
  pedPdc2[27] = his128->GetXaxis()->GetBinCenter(his128->GetMaximumBin()); his128->GetXaxis()->SetRangeUser(pedPdc2[27]-50., pedPdc2[27]+50.);
  pedPdc2[28] = his129->GetXaxis()->GetBinCenter(his129->GetMaximumBin()); his129->GetXaxis()->SetRangeUser(pedPdc2[28]-50., pedPdc2[28]+50.);
  pedPdc2[29] = his130->GetXaxis()->GetBinCenter(his130->GetMaximumBin()); his130->GetXaxis()->SetRangeUser(pedPdc2[29]-50., pedPdc2[29]+50.);
  pedPdc2[30] = his131->GetXaxis()->GetBinCenter(his131->GetMaximumBin()); his131->GetXaxis()->SetRangeUser(pedPdc2[30]-50., pedPdc2[30]+50.);
  pedPdc2[31] = his132->GetXaxis()->GetBinCenter(his132->GetMaximumBin()); his132->GetXaxis()->SetRangeUser(pedPdc2[31]-50., pedPdc2[31]+50.);

  pedPadc[0]  = his201->GetXaxis()->GetBinCenter(his201->GetMaximumBin()); his201->GetXaxis()->SetRangeUser(pedPadc[0]-50., pedPadc[0]+50.);
  pedPadc[1]  = his202->GetXaxis()->GetBinCenter(his202->GetMaximumBin()); his202->GetXaxis()->SetRangeUser(pedPadc[1]-50., pedPadc[1]+50.);
  pedPadc[2]  = his203->GetXaxis()->GetBinCenter(his203->GetMaximumBin()); his203->GetXaxis()->SetRangeUser(pedPadc[2]-50., pedPadc[2]+50.);
  pedPadc[3]  = his204->GetXaxis()->GetBinCenter(his204->GetMaximumBin()); his204->GetXaxis()->SetRangeUser(pedPadc[3]-50., pedPadc[3]+50.);
  pedPadc[4]  = his205->GetXaxis()->GetBinCenter(his205->GetMaximumBin()); his205->GetXaxis()->SetRangeUser(pedPadc[4]-50., pedPadc[4]+50.);
  pedPadc[5]  = his206->GetXaxis()->GetBinCenter(his206->GetMaximumBin()); his206->GetXaxis()->SetRangeUser(pedPadc[5]-50., pedPadc[5]+50.);
  pedPadc[6]  = his207->GetXaxis()->GetBinCenter(his207->GetMaximumBin()); his207->GetXaxis()->SetRangeUser(pedPadc[6]-50., pedPadc[6]+50.);
  pedPadc[7]  = his208->GetXaxis()->GetBinCenter(his208->GetMaximumBin()); his208->GetXaxis()->SetRangeUser(pedPadc[7]-50., pedPadc[7]+50.);
  pedPadc[8]  = his209->GetXaxis()->GetBinCenter(his209->GetMaximumBin()); his209->GetXaxis()->SetRangeUser(pedPadc[8]-50., pedPadc[8]+50.);
  pedPadc[9]  = his210->GetXaxis()->GetBinCenter(his210->GetMaximumBin()); his210->GetXaxis()->SetRangeUser(pedPadc[9]-50., pedPadc[9]+50.);
  pedPadc[10] = his211->GetXaxis()->GetBinCenter(his211->GetMaximumBin()); his211->GetXaxis()->SetRangeUser(pedPadc[10]-50., pedPadc[10]+50.);
  pedPadc[11] = his212->GetXaxis()->GetBinCenter(his212->GetMaximumBin()); his212->GetXaxis()->SetRangeUser(pedPadc[11]-50., pedPadc[11]+50.);
  pedPadc[12] = his213->GetXaxis()->GetBinCenter(his213->GetMaximumBin()); his213->GetXaxis()->SetRangeUser(pedPadc[12]-50., pedPadc[12]+50.);
  pedPadc[13] = his214->GetXaxis()->GetBinCenter(his214->GetMaximumBin()); his214->GetXaxis()->SetRangeUser(pedPadc[13]-50., pedPadc[13]+50.);
  pedPadc[14] = his215->GetXaxis()->GetBinCenter(his215->GetMaximumBin()); his215->GetXaxis()->SetRangeUser(pedPadc[14]-50., pedPadc[14]+50.);
  pedPadc[15] = his216->GetXaxis()->GetBinCenter(his216->GetMaximumBin()); his216->GetXaxis()->SetRangeUser(pedPadc[15]-50., pedPadc[15]+50.);
  pedPadc[16] = his217->GetXaxis()->GetBinCenter(his217->GetMaximumBin()); his217->GetXaxis()->SetRangeUser(pedPadc[16]-50., pedPadc[16]+50.);
  pedPadc[17] = his218->GetXaxis()->GetBinCenter(his218->GetMaximumBin()); his218->GetXaxis()->SetRangeUser(pedPadc[17]-50., pedPadc[17]+50.);
  pedPadc[18] = his219->GetXaxis()->GetBinCenter(his219->GetMaximumBin()); his219->GetXaxis()->SetRangeUser(pedPadc[18]-50., pedPadc[18]+50.);
  pedPadc[19] = his220->GetXaxis()->GetBinCenter(his220->GetMaximumBin()); his220->GetXaxis()->SetRangeUser(pedPadc[19]-50., pedPadc[19]+50.);
  pedPadc[20] = his221->GetXaxis()->GetBinCenter(his221->GetMaximumBin()); his221->GetXaxis()->SetRangeUser(pedPadc[20]-50., pedPadc[20]+50.);
  pedPadc[21] = his222->GetXaxis()->GetBinCenter(his222->GetMaximumBin()); his222->GetXaxis()->SetRangeUser(pedPadc[21]-50., pedPadc[21]+50.);
  pedPadc[22] = his223->GetXaxis()->GetBinCenter(his223->GetMaximumBin()); his223->GetXaxis()->SetRangeUser(pedPadc[22]-50., pedPadc[22]+50.);
  pedPadc[23] = his224->GetXaxis()->GetBinCenter(his224->GetMaximumBin()); his224->GetXaxis()->SetRangeUser(pedPadc[23]-50., pedPadc[23]+50.);
  pedPadc[24] = his225->GetXaxis()->GetBinCenter(his225->GetMaximumBin()); his225->GetXaxis()->SetRangeUser(pedPadc[24]-50., pedPadc[24]+50.);
  pedPadc[25] = his226->GetXaxis()->GetBinCenter(his226->GetMaximumBin()); his226->GetXaxis()->SetRangeUser(pedPadc[25]-50., pedPadc[25]+50.);
  pedPadc[26] = his227->GetXaxis()->GetBinCenter(his227->GetMaximumBin()); his227->GetXaxis()->SetRangeUser(pedPadc[26]-50., pedPadc[26]+50.);
  pedPadc[27] = his228->GetXaxis()->GetBinCenter(his228->GetMaximumBin()); his228->GetXaxis()->SetRangeUser(pedPadc[27]-50., pedPadc[27]+50.);
  pedPadc[28] = his229->GetXaxis()->GetBinCenter(his229->GetMaximumBin()); his229->GetXaxis()->SetRangeUser(pedPadc[28]-50., pedPadc[28]+50.);
  pedPadc[29] = his230->GetXaxis()->GetBinCenter(his230->GetMaximumBin()); his230->GetXaxis()->SetRangeUser(pedPadc[29]-50., pedPadc[29]+50.);
  pedPadc[30] = his231->GetXaxis()->GetBinCenter(his231->GetMaximumBin()); his231->GetXaxis()->SetRangeUser(pedPadc[30]-50., pedPadc[30]+50.);
  pedPadc[31] = his232->GetXaxis()->GetBinCenter(his232->GetMaximumBin()); his232->GetXaxis()->SetRangeUser(pedPadc[31]-50., pedPadc[31]+50.);

  pedPdc30 = his300->GetXaxis()->GetBinCenter(his300->GetMaximumBin());his300->GetXaxis()->SetRangeUser(pedPdc30-50., pedPdc30+50.);
  pedPdc34 = his301->GetXaxis()->GetBinCenter(his301->GetMaximumBin());his301->GetXaxis()->SetRangeUser(pedPdc34-50., pedPdc34+50.);

  //calc mean again after Adjusting Range

  pedQadc[0] = his11->GetMean();
  pedQadc[1] = his12->GetMean();
  pedQadc[2] = his13->GetMean();
  pedQadc[3] = his14->GetMean();
  pedQadc[4] = his15->GetMean();
  pedQadc[5] = his16->GetMean();
  pedQadc[6] = his17->GetMean();
  pedQadc[7] = his18->GetMean();

  pedPdc2[0]  = his101->GetMean();
  pedPdc2[1]  = his102->GetMean();
  pedPdc2[2]  = his103->GetMean();
  pedPdc2[3]  = his104->GetMean();
  pedPdc2[4]  = his105->GetMean();
  pedPdc2[5]  = his106->GetMean();
  pedPdc2[6]  = his107->GetMean();
  pedPdc2[7]  = his108->GetMean();
  pedPdc2[8]  = his109->GetMean();
  pedPdc2[9]  = his110->GetMean();
  pedPdc2[10] = his111->GetMean();
  pedPdc2[11] = his112->GetMean();
  pedPdc2[12] = his113->GetMean();
  pedPdc2[13] = his114->GetMean();
  pedPdc2[14] = his115->GetMean();
  pedPdc2[15] = his116->GetMean();
  pedPdc2[16] = his117->GetMean();
  pedPdc2[17] = his118->GetMean();
  pedPdc2[18] = his119->GetMean();
  pedPdc2[19] = his120->GetMean();
  pedPdc2[20] = his121->GetMean();
  pedPdc2[21] = his122->GetMean();
  pedPdc2[22] = his123->GetMean();
  pedPdc2[23] = his124->GetMean();
  pedPdc2[24] = his125->GetMean();
  pedPdc2[25] = his126->GetMean();
  pedPdc2[26] = his127->GetMean();
  pedPdc2[27] = his128->GetMean();
  pedPdc2[28] = his129->GetMean();
  pedPdc2[29] = his130->GetMean();
  pedPdc2[30] = his131->GetMean();
  pedPdc2[31] = his132->GetMean();

  pedPadc[0]  = his201->GetMean();
  pedPadc[1]  = his202->GetMean();
  pedPadc[2]  = his203->GetMean();
  pedPadc[3]  = his204->GetMean();
  pedPadc[4]  = his205->GetMean();
  pedPadc[5]  = his206->GetMean();
  pedPadc[6]  = his207->GetMean();
  pedPadc[7]  = his208->GetMean();
  pedPadc[8]  = his209->GetMean();
  pedPadc[9]  = his210->GetMean();
  pedPadc[10] = his211->GetMean();
  pedPadc[11] = his212->GetMean();
  pedPadc[12] = his213->GetMean();
  pedPadc[13] = his214->GetMean();
  pedPadc[14] = his215->GetMean();
  pedPadc[15] = his216->GetMean();
  pedPadc[16] = his217->GetMean();
  pedPadc[17] = his218->GetMean();
  pedPadc[18] = his219->GetMean();
  pedPadc[19] = his220->GetMean();
  pedPadc[20] = his221->GetMean();
  pedPadc[21] = his222->GetMean();
  pedPadc[22] = his223->GetMean();
  pedPadc[23] = his224->GetMean();
  pedPadc[24] = his225->GetMean();
  pedPadc[25] = his226->GetMean();
  pedPadc[26] = his227->GetMean();
  pedPadc[27] = his228->GetMean();
  pedPadc[28] = his229->GetMean();
  pedPadc[29] = his230->GetMean();
  pedPadc[30] = his231->GetMean();
  pedPadc[31] = his232->GetMean();

  pedPdc30 = his300->GetMean();
  pedPdc34 = his301->GetMean();


  // Write pedestals to file
  char tempPedFile[500];
  sprintf(tempPedFile, "%s/pedestals_%s.dat", getenv("PEDESTALS"),argv[1]);
  ofstream outPedFile(tempPedFile);

  outPedFile << std::fixed << std::setprecision(7);

  for (int i=0; i<8; i++) {
    outPedFile <<  argv[1] << " " << pedQadc[i] << endl;
  }
  for (int i=0; i<32; i++) {
    outPedFile << argv[1] << " " << pedPdc2[i] << endl;
  }
  for (int i=0; i<32; i++) {
    outPedFile << argv[1] << " " << pedPadc[i] << endl;
  }
  outPedFile << argv[1] << " " << pedPdc30 << endl;
  outPedFile << argv[1] << " " << pedPdc34 << endl;

  outPedFile.close();


  ofstream outWidthFile(TString::Format("%s/PMT_pedestals_%i.dat", getenv("PEDESTALS"),run).Data());

  outWidthFile << std::fixed << std::setprecision(7);


  outWidthFile << run << " " << his11->GetMean() << " " << his11->GetRMS() << std::endl;
  outWidthFile << run << " " << his12->GetMean() << " " << his12->GetRMS() << std::endl;
  outWidthFile << run << " " << his13->GetMean() << " " << his13->GetRMS() << std::endl;
  outWidthFile << run << " " << his14->GetMean() << " " << his14->GetRMS() << std::endl;
  outWidthFile << run << " " << his15->GetMean() << " " << his15->GetRMS() << std::endl;
  outWidthFile << run << " " << his16->GetMean() << " " << his16->GetRMS() << std::endl;
  outWidthFile << run << " " << his17->GetMean() << " " << his17->GetRMS() << std::endl;
  outWidthFile << run << " " << his18->GetMean() << " " << his18->GetRMS();
  
  outWidthFile.close();



  // Write output ntuple
  fileOut->Write();
  fileOut->Close();

  return 0;
}

void loadCuts(Int_t runNumber, cuts* Cuts) {

  // Read cuts file                                                                                                                            
  char tempFileCuts[500];
  sprintf(tempFileCuts, "%s/cuts_%i.dat", getenv("CUTS"),runNumber);
  std::cout << "... Reading: " << tempFileCuts << std::endl;

  ifstream fileCuts(tempFileCuts);
  Char_t comment[500];
  fileCuts >> Cuts->cutBeamBurstTime >> comment;
  fileCuts >> Cuts->nCutsTimeWindows >> comment;

  Double_t cutLowerHold, cutUpperHold;

  if (Cuts->nCutsTimeWindows > 0) {
    for (int i=0; i<Cuts->nCutsTimeWindows; i++) {
      fileCuts >> cutLowerHold >> cutUpperHold;
      Cuts->cutTimeWindowLower.push_back(cutLowerHold);
      Cuts->cutTimeWindowUpper.push_back(cutUpperHold);
    }
  }
  fileCuts >> Cuts->cutEastAnode >> comment;
  fileCuts >> Cuts->cutWestAnode >> comment;
  fileCuts >> Cuts->cutEastTwoFold >> comment;
  fileCuts >> Cuts->cutWestTwoFold >> comment;
  fileCuts >> Cuts->cutEastTopVetoQADC >> comment;
  fileCuts >> Cuts->cutEastTopVetoTDC >> comment;
  fileCuts >> Cuts->cutEastDriftTubeTAC >> comment;
  fileCuts >> Cuts->cutWestDriftTubeTAC >> comment;
  fileCuts >> Cuts->cutEastBackingVetoQADC >> comment;
  fileCuts >> Cuts->cutEastBackingVetoTDC >> comment;
  fileCuts >> Cuts->cutWestBackingVetoQADC >> comment;
  fileCuts >> Cuts->cutWestBackingVetoTDC >> comment;

  std::cout << "... Beam Burst T0 Cut: " << Cuts->cutBeamBurstTime << std::endl;
  std::cout << "... Number of Time Windows Cuts: " << Cuts->nCutsTimeWindows << std::endl;
  if (Cuts->nCutsTimeWindows > 0) {
    for (int i=0; i<Cuts->nCutsTimeWindows; i++) {
      std::cout << "        [" << Cuts->cutTimeWindowLower[i] << ", " << Cuts->cutTimeWindowUpper[i] << "]" << std::endl;
    }
  }
  std::cout << "... East MWPC Anode Cut: " << Cuts->cutEastAnode << std::endl;
  std::cout << "... West MWPC Anode Cut: " << Cuts->cutWestAnode << std::endl;
  std::cout << "... East Scintillator Two-Fold Trigger Cut: " << Cuts->cutEastTwoFold << std::endl;
  std::cout << "... West Scintillator Two-Fold Trigger Cut: " << Cuts->cutWestTwoFold << std::endl;
  std::cout << "... East Top Veto QADC Cut: " << Cuts->cutEastTopVetoQADC << std::endl;
  std::cout << "... East Top Veto TDC Cut: " << Cuts->cutEastTopVetoTDC << std::endl;
  std::cout << "... East Drift Tube TAC Cut: " << Cuts->cutEastDriftTubeTAC << std::endl;
  std::cout << "... West Drift Tube TAC Cut: " << Cuts->cutWestDriftTubeTAC << std::endl;
  std::cout << "... East Backing Veto QADC Cut: " << Cuts->cutEastBackingVetoQADC << std::endl;
  std::cout << "... East Backing Veto TDC Cut: " << Cuts->cutEastBackingVetoTDC << std::endl;
  std::cout << "... West Backing Veto QADC Cut: " << Cuts->cutWestBackingVetoQADC << std::endl;
  std::cout << "... West Backing Veto TDC Cut: " << Cuts->cutWestBackingVetoTDC << std::endl;

}


