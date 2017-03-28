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

#include "fullTreeVariables.h"
#include "MWPCGeometry.h"

using namespace std;

int main(int argc, char *argv[])
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  // Open output ntuple
  char tempOut[500];
  sprintf(tempOut, "%s/basic_histograms_%s.root", getenv("BASIC_HISTOGRAMS"),argv[1]);
  TFile *fileOut = new TFile(tempOut,"RECREATE");

  // Define output histograms
  TH1F *his1  = new TH1F("his1",  "", 1000,0.0,1000.0); // East MWPC Anode PADC
  TH1F *his2  = new TH1F("his2",  "", 1000,0.0,1000.0); // West MWPC Anode PADC
  TH1F *his3  = new TH1F("his3",  "",  180,0.0, 180.0); // East Two-Fold Timing for Scintillator Triggers
  TH1F *his4  = new TH1F("his4",  "",  180,0.0, 180.0); // West Two-Fold Timing for Scintillator Triggers

  TH1F *his11 = new TH1F("his11", "", 4000,0.0,4000.0); // East PMT #1 QADC Scintillator Triggers
  TH1F *his12 = new TH1F("his12", "", 4000,0.0,4000.0); // East PMT #2 QADC Scintillator Triggers
  TH1F *his13 = new TH1F("his13", "", 4000,0.0,4000.0); // East PMT #3 QADC Scintillator Triggers
  TH1F *his14 = new TH1F("his14", "", 4000,0.0,4000.0); // East PMT #4 QADC Scintillator Triggers
  TH1F *his15 = new TH1F("his15", "", 4000,0.0,4000.0); // West PMT #1 QADC Scintillator Triggers
  TH1F *his16 = new TH1F("his16", "", 4000,0.0,4000.0); // West PMT #2 QADC Scintillator Triggers
  TH1F *his17 = new TH1F("his17", "", 4000,0.0,4000.0); // West PMT #3 QADC Scintillator Triggers
  TH1F *his18 = new TH1F("his18", "", 4000,0.0,4000.0); // West PMT #4 QADC Scintillator Triggers

  TH1F *his21 = new TH1F("his21", "", 4000,0.0,4000.0); // East PMT #1 QADC Bi Pulser Triggers
  TH1F *his22 = new TH1F("his22", "", 4000,0.0,4000.0); // East PMT #2 QADC Bi Pulser Triggers
  TH1F *his23 = new TH1F("his23", "", 4000,0.0,4000.0); // East PMT #3 QADC Bi Pulser Triggers
  TH1F *his24 = new TH1F("his24", "", 4000,0.0,4000.0); // East PMT #4 QADC Bi Pulser Triggers
  TH1F *his25 = new TH1F("his25", "", 4000,0.0,4000.0); // West PMT #1 QADC Bi Pulser Triggers
  TH1F *his26 = new TH1F("his26", "", 4000,0.0,4000.0); // West PMT #2 QADC Bi Pulser Triggers
  TH1F *his27 = new TH1F("his27", "", 4000,0.0,4000.0); // West PMT #3 QADC Bi Pulser Triggers
  TH1F *his28 = new TH1F("his28", "", 4000,0.0,4000.0); // West PMT #4 QADC Bi Pulser Triggers

  TH1F *his31 = new TH1F("his31", "", 4000,0.0,4000.0); // East PMT #1 QADC LED Pulser Triggers
  TH1F *his32 = new TH1F("his32", "", 4000,0.0,4000.0); // East PMT #2 QADC LED Pulser Triggers
  TH1F *his33 = new TH1F("his33", "", 4000,0.0,4000.0); // East PMT #3 QADC LED Pulser Triggers
  TH1F *his34 = new TH1F("his34", "", 4000,0.0,4000.0); // East PMT #4 QADC LED Pulser Triggers
  TH1F *his35 = new TH1F("his35", "", 4000,0.0,4000.0); // West PMT #1 QADC LED Pulser Triggers
  TH1F *his36 = new TH1F("his36", "", 4000,0.0,4000.0); // West PMT #2 QADC LED Pulser Triggers
  TH1F *his37 = new TH1F("his37", "", 4000,0.0,4000.0); // West PMT #3 QADC LED Pulser Triggers
  TH1F *his38 = new TH1F("his38", "", 4000,0.0,4000.0); // West PMT #4 QADC LED Pulser Triggers

  TH1F *his101 = new TH1F("his101", "", 4000,0.0,4000.0); // Pdc2[0]
  TH1F *his102 = new TH1F("his102", "", 4000,0.0,4000.0); // Pdc2[1]
  TH1F *his103 = new TH1F("his103", "", 4000,0.0,4000.0); // Pdc2[2]
  TH1F *his104 = new TH1F("his104", "", 4000,0.0,4000.0); // Pdc2[3]
  TH1F *his105 = new TH1F("his105", "", 4000,0.0,4000.0); // Pdc2[4]
  TH1F *his106 = new TH1F("his106", "", 4000,0.0,4000.0); // Pdc2[5]
  TH1F *his107 = new TH1F("his107", "", 4000,0.0,4000.0); // Pdc2[6]
  TH1F *his108 = new TH1F("his108", "", 4000,0.0,4000.0); // Pdc2[7]
  TH1F *his109 = new TH1F("his109", "", 4000,0.0,4000.0); // Pdc2[8]
  TH1F *his110 = new TH1F("his110", "", 4000,0.0,4000.0); // Pdc2[9]
  TH1F *his111 = new TH1F("his111", "", 4000,0.0,4000.0); // Pdc2[10]
  TH1F *his112 = new TH1F("his112", "", 4000,0.0,4000.0); // Pdc2[11]
  TH1F *his113 = new TH1F("his113", "", 4000,0.0,4000.0); // Pdc2[12]
  TH1F *his114 = new TH1F("his114", "", 4000,0.0,4000.0); // Pdc2[13]
  TH1F *his115 = new TH1F("his115", "", 4000,0.0,4000.0); // Pdc2[14]
  TH1F *his116 = new TH1F("his116", "", 4000,0.0,4000.0); // Pdc2[15]
  TH1F *his117 = new TH1F("his117", "", 4000,0.0,4000.0); // Pdc2[16]
  TH1F *his118 = new TH1F("his118", "", 4000,0.0,4000.0); // Pdc2[17]
  TH1F *his119 = new TH1F("his119", "", 4000,0.0,4000.0); // Pdc2[18]
  TH1F *his120 = new TH1F("his120", "", 4000,0.0,4000.0); // Pdc2[19]
  TH1F *his121 = new TH1F("his121", "", 4000,0.0,4000.0); // Pdc2[20]
  TH1F *his122 = new TH1F("his122", "", 4000,0.0,4000.0); // Pdc2[21]
  TH1F *his123 = new TH1F("his123", "", 4000,0.0,4000.0); // Pdc2[22]
  TH1F *his124 = new TH1F("his124", "", 4000,0.0,4000.0); // Pdc2[23]
  TH1F *his125 = new TH1F("his125", "", 4000,0.0,4000.0); // Pdc2[24]
  TH1F *his126 = new TH1F("his126", "", 4000,0.0,4000.0); // Pdc2[25]
  TH1F *his127 = new TH1F("his127", "", 4000,0.0,4000.0); // Pdc2[26]
  TH1F *his128 = new TH1F("his128", "", 4000,0.0,4000.0); // Pdc2[27]
  TH1F *his129 = new TH1F("his129", "", 4000,0.0,4000.0); // Pdc2[28]
  TH1F *his130 = new TH1F("his130", "", 4000,0.0,4000.0); // Pdc2[29]
  TH1F *his131 = new TH1F("his131", "", 4000,0.0,4000.0); // Pdc2[30]
  TH1F *his132 = new TH1F("his132", "", 4000,0.0,4000.0); // Pdc2[31]

  TH1F *his201 = new TH1F("his201", "", 4000,0.0,4000.0); // Padc2[0]
  TH1F *his202 = new TH1F("his202", "", 4000,0.0,4000.0); // Padc2[1]
  TH1F *his203 = new TH1F("his203", "", 4000,0.0,4000.0); // Padc2[2]
  TH1F *his204 = new TH1F("his204", "", 4000,0.0,4000.0); // Padc2[3]
  TH1F *his205 = new TH1F("his205", "", 4000,0.0,4000.0); // Padc2[4]
  TH1F *his206 = new TH1F("his206", "", 4000,0.0,4000.0); // Padc2[5]
  TH1F *his207 = new TH1F("his207", "", 4000,0.0,4000.0); // Padc2[6]
  TH1F *his208 = new TH1F("his208", "", 4000,0.0,4000.0); // Padc2[7]
  TH1F *his209 = new TH1F("his209", "", 4000,0.0,4000.0); // Padc2[8]
  TH1F *his210 = new TH1F("his210", "", 4000,0.0,4000.0); // Padc2[9]
  TH1F *his211 = new TH1F("his211", "", 4000,0.0,4000.0); // Padc2[10]
  TH1F *his212 = new TH1F("his212", "", 4000,0.0,4000.0); // Padc2[11]
  TH1F *his213 = new TH1F("his213", "", 4000,0.0,4000.0); // Padc2[12]
  TH1F *his214 = new TH1F("his214", "", 4000,0.0,4000.0); // Padc2[13]
  TH1F *his215 = new TH1F("his215", "", 4000,0.0,4000.0); // Padc2[14]
  TH1F *his216 = new TH1F("his216", "", 4000,0.0,4000.0); // Padc2[15]
  TH1F *his217 = new TH1F("his217", "", 4000,0.0,4000.0); // Padc2[16]
  TH1F *his218 = new TH1F("his218", "", 4000,0.0,4000.0); // Padc2[17]
  TH1F *his219 = new TH1F("his219", "", 4000,0.0,4000.0); // Padc2[18]
  TH1F *his220 = new TH1F("his220", "", 4000,0.0,4000.0); // Padc2[19]
  TH1F *his221 = new TH1F("his221", "", 4000,0.0,4000.0); // Padc2[20]
  TH1F *his222 = new TH1F("his222", "", 4000,0.0,4000.0); // Padc2[21]
  TH1F *his223 = new TH1F("his223", "", 4000,0.0,4000.0); // Padc2[22]
  TH1F *his224 = new TH1F("his224", "", 4000,0.0,4000.0); // Padc2[23]
  TH1F *his225 = new TH1F("his225", "", 4000,0.0,4000.0); // Padc2[24]
  TH1F *his226 = new TH1F("his226", "", 4000,0.0,4000.0); // Padc2[25]
  TH1F *his227 = new TH1F("his227", "", 4000,0.0,4000.0); // Padc2[26]
  TH1F *his228 = new TH1F("his228", "", 4000,0.0,4000.0); // Padc2[27]
  TH1F *his229 = new TH1F("his229", "", 4000,0.0,4000.0); // Padc2[28]
  TH1F *his230 = new TH1F("his230", "", 4000,0.0,4000.0); // Padc2[29]
  TH1F *his231 = new TH1F("his231", "", 4000,0.0,4000.0); // Padc2[30]
  TH1F *his232 = new TH1F("his232", "", 4000,0.0,4000.0); // Padc2[31]

  TH1F *his301 = new TH1F("his301", "", 5000,0.0,5000.0); // Timing Spectrum for Scintillator Triggers
  TH1F *his302 = new TH1F("his302", "", 5000,0.0,5000.0); // Timing Spectrum for Gate Valve UCN Monitor
  TH1F *his303 = new TH1F("his303", "", 5000,0.0,5000.0); // Timing Spectrum for Switcher UCN Monitor
  TH1F *his304 = new TH1F("his304", "", 5000,0.0,5000.0); // Timing Spectrum for AFP Fe Foil UCN Monitor
  TH1F *his305 = new TH1F("his305", "", 5000,0.0,5000.0); // Timing Spectrum for SCS UCN Monitor

  TH1F *his401 = new TH1F("his401", "", 4000,0.0,4000.0); // East Top Veto QADC
  TH1F *his402 = new TH1F("his402", "", 4000,0.0,4000.0); // East Top Veto TDC
  TH1F *his403 = new TH1F("his403", "", 4000,0.0,4000.0); // East Drift Tube Veto TAC
  TH1F *his404 = new TH1F("his404", "", 4000,0.0,4000.0); // West Drift Tube Veto TAC
  TH1F *his405 = new TH1F("his405", "", 4000,0.0,4000.0); // East Backing Veto QADC
  TH1F *his406 = new TH1F("his406", "", 4000,0.0,4000.0); // East Backing Veto TDC
  TH1F *his407 = new TH1F("his407", "", 4000,0.0,4000.0); // West Backing Veto QADC
  TH1F *his408 = new TH1F("his408", "", 4000,0.0,4000.0); // West Backing Veto TDC

  TH2F *his501 = new TH2F("his501", "", 120,-60.0,60.0, 120,-60.0,60.0); // Coarse East (x,y) Position Map
  TH2F *his502 = new TH2F("his502", "", 120,-60.0,60.0, 120,-60.0,60.0); // Coarse West (x,y) Position Map

  // East MWPC y                                                                  
  posPdc2[0]  =  76.20;
  posPdc2[1]  =  66.04;
  posPdc2[2]  =  55.88;
  posPdc2[3]  =  45.72;
  posPdc2[4]  =  35.56;
  posPdc2[5]  =  25.40;
  posPdc2[6]  =  15.24;
  posPdc2[7]  =   5.08;
  posPdc2[8]  =  -5.08;
  posPdc2[9]  = -15.24;
  posPdc2[10] = -25.40;
  posPdc2[11] = -35.56;
  posPdc2[12] = -45.72;
  posPdc2[13] = -55.88;
  posPdc2[14] = -66.04;
  posPdc2[15] = -76.20;

  // East MWPC x                                                                  
  posPdc2[16] =  76.20;
  posPdc2[17] =  66.04;
  posPdc2[18] =  55.88;
  posPdc2[19] =  45.72;
  posPdc2[20] =  35.56;
  posPdc2[21] =  25.40;
  posPdc2[22] =  15.24;
  posPdc2[23] =   5.08;
  posPdc2[24] =  -5.08;
  posPdc2[25] = -15.24;
  posPdc2[26] = -25.40;
  posPdc2[27] = -35.56;
  posPdc2[28] = -45.72;
  posPdc2[29] = -55.88;
  posPdc2[30] = -66.04;
  posPdc2[31] = -76.20;

  // West MWPC y                                                                  
  posPadc[0]  =  76.20;
  posPadc[1]  =  66.04;
  posPadc[2]  =  55.88;
  posPadc[3]  =  45.72;
  posPadc[4]  =  35.56;
  posPadc[5]  =  25.40;
  posPadc[6]  =  15.24;
  posPadc[7]  =   5.08;
  posPadc[8]  =  -5.08;
  posPadc[9]  = -15.24;
  posPadc[10] = -25.40;
  posPadc[11] = -35.56;
  posPadc[12] = -45.72;
  posPadc[13] = -55.88;
  posPadc[14] = -66.04;
  posPadc[15] = -76.20;

  // West MWPC x                                                                  
  posPadc[16] = -76.20;
  posPadc[17] = -66.04;
  posPadc[18] = -55.88;
  posPadc[19] = -45.72;
  posPadc[20] = -35.56;
  posPadc[21] = -25.40;
  posPadc[22] = -15.24;
  posPadc[23] =  -5.08;
  posPadc[24] =   5.08;
  posPadc[25] =  15.24;
  posPadc[26] =  25.40;
  posPadc[27] =  35.56;
  posPadc[28] =  45.72;
  posPadc[29] =  55.88;
  posPadc[30] =  66.04;
  posPadc[31] =  76.20;

  // Open input ntuple
  char tempIn[500];
  sprintf(tempIn, "/extern/mabrow05/ucna/rawdata/full%s.root", argv[1]);

  TFile *fileIn = new TFile(tempIn, "READ");
  TTree *Tin = (TTree*)(fileIn->Get("h1"));

  // Variables
  Tin->SetBranchAddress("Pdc30",  &Pdc30);
  Tin->SetBranchAddress("Pdc34",  &Pdc34);

  Tin->SetBranchAddress("Tdc016", &Tdc016);
  Tin->SetBranchAddress("Tdc017", &Tdc017);

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

  int nEvents = Tin->GetEntries();
  cout << "Run " << argv[1] << " ..." << endl;
  cout << "... Processing nEvents = " << nEvents << endl;

  // Loop over events
  for (int i=0; i<nEvents; i++) {
    Tin->GetEvent(i);
    int iSis00 = (int) Sis00;

    his1->Fill(Pdc30);
    his2->Fill(Pdc34);

    if (iSis00 < 4) {
      his3->Fill(Tdc016*tdcChannelToTime);
      his4->Fill(Tdc017*tdcChannelToTime);
    }

    if (iSis00 < 4) {
      his11->Fill(Qadc[0]);
      his12->Fill(Qadc[1]);
      his13->Fill(Qadc[2]);
      his14->Fill(Qadc[3]);
      his15->Fill(Qadc[4]);
      his16->Fill(Qadc[5]);
      his17->Fill(Qadc[6]);
      his18->Fill(Qadc[7]);
    }

    if (iSis00 == 32) {
      his21->Fill(Qadc[0]);
      his22->Fill(Qadc[1]);
      his23->Fill(Qadc[2]);
      his24->Fill(Qadc[3]);
      his25->Fill(Qadc[4]);
      his26->Fill(Qadc[5]);
      his27->Fill(Qadc[6]);
      his28->Fill(Qadc[7]);
    }

    if (iSis00 == 163) {
      his31->Fill(Qadc[0]);
      his32->Fill(Qadc[1]);
      his33->Fill(Qadc[2]);
      his34->Fill(Qadc[3]);
      his35->Fill(Qadc[4]);
      his36->Fill(Qadc[5]);
      his37->Fill(Qadc[6]);
      his38->Fill(Qadc[7]);
    }

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

    if (iSis00 < 4) {
      his301->Fill(S83028*scalerCountsToTime);
    }

    if (iSis00 == 260)  his302->Fill(S83028*scalerCountsToTime);
    if (iSis00 == 516)  his303->Fill(S83028*scalerCountsToTime);
    if (iSis00 == 1028) his304->Fill(S83028*scalerCountsToTime);
    if (iSis00 == 2052) his305->Fill(S83028*scalerCountsToTime);

    his401->Fill(Qadc9);
    his402->Fill(Tdc019);
    his403->Fill(Pdc313);
    his404->Fill(Pdc315);
    his405->Fill(Qadc8);
    his406->Fill(Tdc018);
    his407->Fill(Qadc10);
    his408->Fill(Tdc020);

    double xMWPCEast, yMWPCEast, xMWPCWest, yMWPCWest;
    double maxEastxPdc2, maxEastyPdc2, maxWestxPadc, maxWestyPadc;
    maxEastxPdc2 = -1.0e99;
    maxEastyPdc2 = -1.0e99;
    maxWestxPadc = -1.0e99;
    maxWestyPadc = -1.0e99;
    xMWPCEast = -100.;
    yMWPCEast = -100.;
    xMWPCWest = -100.;
    yMWPCWest = -100.;
    for (int j=0; j<16; j++) {
      if (Pdc2[j+16] > maxEastxPdc2) {
        xMWPCEast    = posPdc2[j+16];
        maxEastxPdc2 = Pdc2[j+16];
      }
      if (Pdc2[j]    > maxEastyPdc2) {
        yMWPCEast    = posPdc2[j];
        maxEastyPdc2 = Pdc2[j];
      }
      if (Padc[j+16] > maxWestxPadc) {
        xMWPCWest    = posPadc[j+16];
        maxWestxPadc = Padc[j+16];
      }
      if (Padc[j]    > maxWestyPadc) {
        yMWPCWest    = posPadc[j];
        maxWestyPadc = Padc[j];
      }
    }
    xMWPCEast *= positionProjection;
    yMWPCEast *= positionProjection;
    xMWPCWest *= positionProjection;
    yMWPCWest *= positionProjection;

    if (iSis00 < 4) {
      his501->Fill(xMWPCEast, yMWPCEast);
      his502->Fill(xMWPCWest, yMWPCWest);
    }

  }

  // Write output ntuple
  fileOut->Write();
  fileOut->Close();

  return 0;
}
