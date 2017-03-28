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
#include "pedestals.h"

using namespace std;

int main(int argc, char *argv[])
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);
  
  // Define output histograms
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

  TFile *fileIn = new TFile(tempIn, "READ");
  TTree *Tin = (TTree*)(fileIn->Get("h1"));

  // Variables
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

  Tin->SetBranchAddress("Pdc30", &Pdc30);
  Tin->SetBranchAddress("Pdc34", &Pdc34);
  
  int nEvents = Tin->GetEntries();
  cout << "Run " << argv[1] << " ..." << endl;
  cout << "... Processing nEvents = " << nEvents << endl;

  int nEastPed = 0;
  int nWestPed = 0;

  // Loop over events
  for (int i=0; i<nEvents; i++) {
    Tin->GetEvent(i);
    int iSis00 = (int) Sis00;

    bool triggerEast = false;
    bool triggerWest = false;
    bool triggerUCN  = false;

    if (iSis00 == 1) triggerEast = true;
    if (iSis00 == 2) triggerWest = true;
    if (iSis00 == 260 || iSis00 == 516 || iSis00 == 1028 || iSis00 == 2052) triggerUCN = true;

    if (triggerWest) { // || triggerUCN) {
      his11->Fill(Qadc[0]);
      his12->Fill(Qadc[1]);
      his13->Fill(Qadc[2]);
      his14->Fill(Qadc[3]);

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

      nEastPed++;
    }

    if (triggerEast) { // || triggerUCN) {
      his15->Fill(Qadc[4]);
      his16->Fill(Qadc[5]);
      his17->Fill(Qadc[6]);
      his18->Fill(Qadc[7]);

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

      nWestPed++;
    }

  }

  // Extract mean values of pedestals
  for (int i=0; i<8; i++) {
    pedQadc[i] = 0.;
  }
  

  
  // Write pedestals to file
  char tempFile[500];
  sprintf(tempFile, "%s/pedestal_widths_%s.dat", getenv("PEDESTALS"),argv[1]);
  ofstream outWidthFile(tempFile);


  outWidthFile << argv[1] << " " << his11->GetMean() << " " << his11->GetRMS() << endl;
  outWidthFile << argv[1] << " " << his12->GetMean() << " " << his12->GetRMS() << endl;
  outWidthFile << argv[1] << " " << his13->GetMean() << " " << his13->GetRMS() << endl;
  outWidthFile << argv[1] << " " << his14->GetMean() << " " << his14->GetRMS() << endl;
  outWidthFile << argv[1] << " " << his15->GetMean() << " " << his15->GetRMS() << endl;
  outWidthFile << argv[1] << " " << his16->GetMean() << " " << his16->GetRMS() << endl;
  outWidthFile << argv[1] << " " << his17->GetMean() << " " << his17->GetRMS() << endl;
  outWidthFile << argv[1] << " " << his18->GetMean() << " " << his18->GetRMS();
  
  outWidthFile.close();
 

  return 0;
}
