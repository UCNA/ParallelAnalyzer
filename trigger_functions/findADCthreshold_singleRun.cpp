#include <vector>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <TString.h>
#include <TH2D.h>
#include <TMinuit.h>
#include <TRandom3.h>
#include <TH1.h>
#include <TF1.h>
#include <TH1F.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TMath.h>

#include "DataTree.hh"




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

std::vector < std::vector <Double_t> > loadPMTpedestals(Int_t runNumber);

void loadCuts(Int_t runNumber, cuts* Cuts);

std::vector <Double_t> loadGainFactors(Int_t runNumber);

void findDiscriminatorThresh(Int_t rn); 

std::vector<Double_t> fitThreshold(TGraphAsymmErrors *trigg, std::vector <Double_t> paramGuesses = std::vector <Double_t>(0));


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {
  
  std::vector <Int_t> evtTypes;

  if (argc!=2) { std::cout << "USAGE: ./findADCthreshold.exe [run]\n"; exit(0); }

  findDiscriminatorThresh(atoi(argv[1]));

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::vector < std::vector <Double_t> > loadPMTpedestals(Int_t runNumber) {

  Char_t temp[500];
  std::vector < std::vector < Double_t > > peds (8,std::vector <Double_t> (2,0.));
  sprintf(temp,"%s/PMT_pedestals_%i.dat",getenv("PEDESTALS"),runNumber);
  std::ifstream infile;
  infile.open(temp);

  Int_t i = 0;
  Int_t run;

  while (infile >> run >> peds[i][0] >> peds[i][1]) { std::cout << "Pedestal " << i << ": " << peds[i][0] << " " << peds[i][1] << std::endl; i++; }
  return peds;

};


void findDiscriminatorThresh(Int_t rn) {

  

  Int_t run = rn;
  
  //TFile *f = NULL;
  //TTree *t = NULL;

  Double_t lower_limit = -60.1;
  Int_t hisBinWidth = 2;
  Double_t upper_limit = 159.9;
  Int_t nBins = (int)((upper_limit-lower_limit)/hisBinWidth);
  //nBins = 100;

  TH1F *totalE1 = new TH1F("totalE1", "total events E1", nBins, lower_limit, upper_limit);
  TH1F *triggerE1 = new TH1F("triggerE1", "trigger events E1", nBins, lower_limit, upper_limit);
  TGraphAsymmErrors *probE1 = new TGraphAsymmErrors();
  
  TH1F *totalE2 = new TH1F("totalE2", "total events E2", nBins, lower_limit, upper_limit);
  TH1F *triggerE2 = new TH1F("triggerE2", "trigger events E2", nBins, lower_limit, upper_limit);
  TGraphAsymmErrors *probE2 = new TGraphAsymmErrors();

  TH1F *totalE3 = new TH1F("totalE3", "total events E3", nBins, lower_limit, upper_limit);
  TH1F *triggerE3 = new TH1F("triggerE3", "trigger events E3", nBins, lower_limit, upper_limit);
  TGraphAsymmErrors *probE3 = new TGraphAsymmErrors();

  TH1F *totalE4 = new TH1F("totalE4", "total events E4", nBins, lower_limit, upper_limit);
  TH1F *triggerE4 = new TH1F("triggerE4", "trigger events E4", nBins, lower_limit, upper_limit);
  TGraphAsymmErrors *probE4 = new TGraphAsymmErrors();

  TH1F *totalW1 = new TH1F("totalW1", "total events W1", nBins, lower_limit, upper_limit);
  TH1F *triggerW1 = new TH1F("triggerW1", "trigger events W1", nBins, lower_limit, upper_limit);
  TGraphAsymmErrors *probW1 = new TGraphAsymmErrors();
  
  TH1F *totalW2 = new TH1F("totalW2", "total events W2", nBins, lower_limit, upper_limit);
  TH1F *triggerW2 = new TH1F("triggerW2", "trigger events W2", nBins, lower_limit, upper_limit);
  TGraphAsymmErrors *probW2 = new TGraphAsymmErrors();

  TH1F *totalW3 = new TH1F("totalW3", "total events W3", nBins, lower_limit, upper_limit);
  TH1F *triggerW3 = new TH1F("triggerW3", "trigger events W3", nBins, lower_limit, upper_limit);
  TGraphAsymmErrors *probW3 = new TGraphAsymmErrors();

  TH1F *totalW4 = new TH1F("totalW4", "total events W4", nBins, lower_limit, upper_limit);
  TH1F *triggerW4 = new TH1F("triggerW4", "trigger events W4", nBins, lower_limit, upper_limit);
  TGraphAsymmErrors *probW4 = new TGraphAsymmErrors();

  //PMT Pedestal Histograms
  TH1F *pedE1 = new TH1F("pedE1", "Pedestal E1", nBins, lower_limit, upper_limit);
  TH1F *pedE2 = new TH1F("pedE2", "Pedestal E2", nBins, lower_limit, upper_limit);
  TH1F *pedE3 = new TH1F("pedE3", "Pedestal E3", nBins, lower_limit, upper_limit);
  TH1F *pedE4 = new TH1F("pedE4", "Pedestal E4", nBins, lower_limit, upper_limit);
  TH1F *pedW1 = new TH1F("pedW1", "Pedestal W1", nBins, lower_limit, upper_limit);
  TH1F *pedW2 = new TH1F("pedW2", "Pedestal W2", nBins, lower_limit, upper_limit);
  TH1F *pedW3 = new TH1F("pedW3", "Pedestal W3", nBins, lower_limit, upper_limit);
  TH1F *pedW4 = new TH1F("pedW4", "Pedestal W4", nBins, lower_limit, upper_limit);
    

  // All pertinent variables to be read in from data files
  /*  Float_t QadcE1, QadcE2, QadcE3, QadcE4, QadcW1, QadcW2, QadcW3, QadcW4;
  Float_t PdcWest, PdcEast; //Wirechamber Anode values
  Float_t TdcE1, TdcE2, TdcE3, TdcE4, TdcW1, TdcW2, TdcW3, TdcW4; 
  Float_t TdcEast2fold, TdcWest2fold;
  Float_t Sis00;
  */
  //Load cuts for this run
  cuts Cuts;
  loadCuts(run,&Cuts);

  std::vector < std::vector <Double_t> > peds = loadPMTpedestals(run);
  std::vector <Double_t> gain = loadGainFactors(run);
    
	 
  DataTree *t = new DataTree();

  char tempIn[500];
  sprintf(tempIn, "%s/replay_pass1_%i.root", getenv("REPLAY_PASS1"),rn);
  //sprintf(tempIn, "../replay_pass1/replay_pass1_%s.root", argv[1]);
  t->setupInputTree(std::string(tempIn),"pass1");

  
  
    
  Int_t nevents = t->getEntries();
    
  std::cout << "nEvents = " << nevents << std::endl;

  for (int evt=0; evt<nevents; evt++) {
      
    t->getEvent(evt);
    Int_t iSis00 = t->Sis00;

    if (evt%100000==0) std::cout << evt << " in run " << run << std::endl;

      
    //if ((int)Sis00>3) /* && !((int)Sis00>=32 && (int)Sis00<=35))*/ continue; //Not a beta trigger or Bi pulser event

    //Int_t side = (int)Sis00 - 1;

    //if ( ((int)Sis00==1 || (int)Sis00==3) && PdcEast<Cuts.cutEastAnode) continue; //not an electron event
    //if ( ((int)Sis00==2  || (int)Sis00==3) && PdcWest<Cuts.cutWestAnode) continue; //not an electron event
      

    // East Side Triggers (including backscatters)
    if ( iSis00==1 ) { //&& t->xeRC>0 && t->xeRC<4 && t->yeRC>0 && t->yeRC<4 ) {  //TdcEast2fold > 0.00001 ) {

      //Fill West side total histograms
      //if (PdcWest>Cuts.cutWestAnode) {
      //std::cout <<  (QadcW1-peds[4][0])*gain[4] << std::endl;

      totalW1->Fill((t->ScintW.q1)*gain[4]);
      totalW2->Fill((t->ScintW.q2)*gain[5]);
      totalW3->Fill((t->ScintW.q3)*gain[6]);
      totalW4->Fill((t->ScintW.q4)*gain[7]);
      
      if ( t->TDCW1>0.00001 ) triggerW1->Fill((t->ScintW.q1)*gain[4]);//-pedW1)*gainW1);
      else pedW1->Fill((t->ScintW.q1)*gain[4]);
      if ( t->TDCW2>0.00001 ) triggerW2->Fill((t->ScintW.q2)*gain[5]);//-pedW2)*gainW2);
      else pedW2->Fill((t->ScintW.q2)*gain[5]);
      if ( t->TDCW3>0.00001 ) triggerW3->Fill((t->ScintW.q3)*gain[6]);//-pedW3)*gainW3);
      else pedW3->Fill((t->ScintW.q3)*gain[6]);
      if ( t->TDCW4>0.00001 ) triggerW4->Fill((t->ScintW.q4)*gain[7]);//-pedW4)*gainW4);
      else pedW4->Fill((t->ScintW.q4)*gain[7]);
      //}
      
      
      
      
      
      //3 PMTs trigger
      if ( (t->TDCE1>0.00001 && t->TDCE2>0.00001 && t->TDCE3>0.00001) || (t->TDCE2>0.00001 && t->TDCE3>0.00001 && t->TDCE4>0.00001) || (t->TDCE1>0.00001 && t->TDCE3>0.00001 && t->TDCE4>0.00001) || (t->TDCE1>0.00001 && t->TDCE2>0.00001 && t->TDCE4>0.00001) ) { 
	
	//East side 3-fold trigger events
	totalE1->Fill((t->ScintE.q1)*gain[0]);
	totalE2->Fill((t->ScintE.q2)*gain[1]);
	totalE3->Fill((t->ScintE.q3)*gain[2]);
	totalE4->Fill((t->ScintE.q4)*gain[3]);
	
	if ( t->TDCE1>0.00001 ) triggerE1->Fill((t->ScintE.q1)*gain[0]);//-pedE1)*gainE1);
	else pedE1->Fill((t->ScintE.q1)*gain[0]);
	if ( t->TDCE2>0.00001 ) triggerE2->Fill((t->ScintE.q2)*gain[1]);//-pedE2)*gainE2);
	else pedE2->Fill((t->ScintE.q2)*gain[1]);
	if ( t->TDCE3>0.00001 ) triggerE3->Fill((t->ScintE.q3)*gain[2]);//-pedE3)*gainE3);
	else pedE3->Fill((t->ScintE.q3)*gain[2]);
	if ( t->TDCE4>0.00001 ) triggerE4->Fill((t->ScintE.q4)*gain[3]);//-pedE4)*gainE4);
	else pedE4->Fill((t->ScintE.q4)*gain[3]);

      }
      
    }
    
    // West Side Triggers (including backscatters)
    if ( iSis00==2 ) { //&& t->xwRC>0 && t->xwRC<4 && t->ywRC>0 && t->ywRC<4 )  { //t->TDCWest2fold > 0.00001 ) {
      
      //Fill East side total histograms
      //if (PdcEast>Cuts.cutEastAnode) {
      totalE1->Fill((t->ScintE.q1)*gain[0]);
      totalE2->Fill((t->ScintE.q2)*gain[1]);
      totalE3->Fill((t->ScintE.q3)*gain[2]);
      totalE4->Fill((t->ScintE.q4)*gain[3]);
      
      if ( t->TDCE1>0.00001 ) triggerE1->Fill((t->ScintE.q1)*gain[0]);//-pedE1)*gainE1);
      else pedE1->Fill((t->ScintE.q1)*gain[0]);
      if ( t->TDCE2>0.00001 ) triggerE2->Fill((t->ScintE.q2)*gain[1]);//-pedE2)*gainE2);
      else pedE2->Fill((t->ScintE.q2)*gain[1]);
      if ( t->TDCE3>0.00001 ) triggerE3->Fill((t->ScintE.q3)*gain[2]);//-pedE3)*gainE3);
      else pedE3->Fill((t->ScintE.q3)*gain[2]);
      if ( t->TDCE4>0.00001 ) triggerE4->Fill((t->ScintE.q4)*gain[3]);//-pedE4)*gainE4);
      else pedE4->Fill((t->ScintE.q4)*gain[3]);
      //}
      
      
      if ( (t->TDCW1>0.00001 && t->TDCW2>0.00001 && t->TDCW3>0.00001) || (t->TDCW2>0.00001 && t->TDCW3>0.00001 && t->TDCW4>0.00001) || (t->TDCW1>0.00001 && t->TDCW3>0.00001 && t->TDCW4>0.00001) || (t->TDCW1>0.00001 && t->TDCW2>0.00001 && t->TDCW4>0.00001) ) { 
	
	//West side events
	totalW1->Fill((t->ScintW.q1)*gain[4]);
	totalW2->Fill((t->ScintW.q2)*gain[5]);
	totalW3->Fill((t->ScintW.q3)*gain[6]);
	totalW4->Fill((t->ScintW.q4)*gain[7]);
      
	if ( t->TDCW1>0.00001 ) triggerW1->Fill((t->ScintW.q1)*gain[4]);//-pedW1)*gainW1);
	else pedW1->Fill((t->ScintW.q1)*gain[4]);
	if ( t->TDCW2>0.00001 ) triggerW2->Fill((t->ScintW.q2)*gain[5]);//-pedW2)*gainW2);
	else pedW2->Fill((t->ScintW.q2)*gain[5]);
	if ( t->TDCW3>0.00001 ) triggerW3->Fill((t->ScintW.q3)*gain[6]);//-pedW3)*gainW3);
	else pedW3->Fill((t->ScintW.q3)*gain[6]);
	if ( t->TDCW4>0.00001 ) triggerW4->Fill((t->ScintW.q4)*gain[7]);//-pedW4)*gainW4);
	else pedW4->Fill((t->ScintW.q4)*gain[7]);
	  
      }
    }

    if (iSis00==3) { // && t->xeRC>0 && t->xeRC<4 && t->yeRC>0 && t->yeRC<4 && t->xwRC>0 && t->xwRC<4 && t->ywRC>0 && t->ywRC<4) {
      
      if ( (t->TDCE1>0.00001 && t->TDCE2>0.00001 && t->TDCE3>0.00001) || (t->TDCE2>0.00001 && t->TDCE3>0.00001 && t->TDCE4>0.00001) || (t->TDCE1>0.00001 && t->TDCE3>0.00001 && t->TDCE4>0.00001) || (t->TDCE1>0.00001 && t->TDCE2>0.00001 && t->TDCE4>0.00001) ) { 
	
	//East side 3-fold trigger events
	totalE1->Fill((t->ScintE.q1)*gain[0]);
	totalE2->Fill((t->ScintE.q2)*gain[1]);
	totalE3->Fill((t->ScintE.q3)*gain[2]);
	totalE4->Fill((t->ScintE.q4)*gain[3]);
	
	if ( t->TDCE1>0.00001 ) triggerE1->Fill((t->ScintE.q1)*gain[0]);//-pedE1)*gainE1);
	else pedE1->Fill((t->ScintE.q1)*gain[0]);
	if ( t->TDCE2>0.00001 ) triggerE2->Fill((t->ScintE.q2)*gain[1]);//-pedE2)*gainE2);
	else pedE2->Fill((t->ScintE.q2)*gain[1]);
	if ( t->TDCE3>0.00001 ) triggerE3->Fill((t->ScintE.q3)*gain[2]);//-pedE3)*gainE3);
	else pedE3->Fill((t->ScintE.q3)*gain[2]);
	if ( t->TDCE4>0.00001 ) triggerE4->Fill((t->ScintE.q4)*gain[3]);//-pedE4)*gainE4);
	else pedE4->Fill((t->ScintE.q4)*gain[3]);
	
      }

      
      if ( (t->TDCW1>0.00001 && t->TDCW2>0.00001 && t->TDCW3>0.00001) || (t->TDCW2>0.00001 && t->TDCW3>0.00001 && t->TDCW4>0.00001) || (t->TDCW1>0.00001 && t->TDCW3>0.00001 && t->TDCW4>0.00001) || (t->TDCW1>0.00001 && t->TDCW2>0.00001 && t->TDCW4>0.00001) ) { 
	
	//West side events
	totalW1->Fill((t->ScintW.q1)*gain[4]);
	totalW2->Fill((t->ScintW.q2)*gain[5]);
	totalW3->Fill((t->ScintW.q3)*gain[6]);
	totalW4->Fill((t->ScintW.q4)*gain[7]);
	
	if ( t->TDCW1>0.00001 ) triggerW1->Fill((t->ScintW.q1)*gain[4]);//-pedW1)*gainW1);
	else pedW1->Fill((t->ScintW.q1)*gain[4]);
	if ( t->TDCW2>0.00001 ) triggerW2->Fill((t->ScintW.q2)*gain[5]);//-pedW2)*gainW2);
	else pedW2->Fill((t->ScintW.q2)*gain[5]);
	if ( t->TDCW3>0.00001 ) triggerW3->Fill((t->ScintW.q3)*gain[6]);//-pedW3)*gainW3);
	else pedW3->Fill((t->ScintW.q3)*gain[6]);
	if ( t->TDCW4>0.00001 ) triggerW4->Fill((t->ScintW.q4)*gain[7]);//-pedW4)*gainW4);
	else pedW4->Fill((t->ScintW.q4)*gain[7]);
      }
    }
    
  }
  delete t;



  triggerE1->GetXaxis()->SetTitle("ADC channels above pedestal");
  triggerE2->GetXaxis()->SetTitle("ADC channels above pedestal");
  triggerE3->GetXaxis()->SetTitle("ADC channels above pedestal");
  triggerE4->GetXaxis()->SetTitle("ADC channels above pedestal");
  triggerW1->GetXaxis()->SetTitle("ADC channels above pedestal");
  triggerW2->GetXaxis()->SetTitle("ADC channels above pedestal");
  triggerW3->GetXaxis()->SetTitle("ADC channels above pedestal");
  triggerW4->GetXaxis()->SetTitle("ADC channels above pedestal");

  totalE1->GetXaxis()->SetTitle("ADC channels above pedestal");
  totalE2->GetXaxis()->SetTitle("ADC channels above pedestal");
  totalE3->GetXaxis()->SetTitle("ADC channels above pedestal");
  totalE4->GetXaxis()->SetTitle("ADC channels above pedestal");
  totalW1->GetXaxis()->SetTitle("ADC channels above pedestal");
  totalW2->GetXaxis()->SetTitle("ADC channels above pedestal");
  totalW3->GetXaxis()->SetTitle("ADC channels above pedestal");
  totalW4->GetXaxis()->SetTitle("ADC channels above pedestal");

  
  triggerE1->SetLineColor(kRed);
  triggerE1->SetLineStyle(7);
  triggerE2->SetLineColor(kRed);
  triggerE2->SetLineStyle(7);
  triggerE3->SetLineColor(kRed);
  triggerE3->SetLineStyle(7);
  triggerE4->SetLineColor(kRed);
  triggerE4->SetLineStyle(7);
  triggerW1->SetLineColor(kRed);
  triggerW1->SetLineStyle(7);
  triggerW2->SetLineColor(kRed);
  triggerW2->SetLineStyle(7);
  triggerW3->SetLineColor(kRed);
  triggerW3->SetLineStyle(7);
  triggerW4->SetLineColor(kRed);
  triggerW4->SetLineStyle(7);

  pedE1->SetLineColor(kGreen);
  pedE1->SetLineStyle(7);
  pedE2->SetLineColor(kGreen);
  pedE2->SetLineStyle(7);
  pedE3->SetLineColor(kGreen);
  pedE3->SetLineStyle(7);
  pedE4->SetLineColor(kGreen);
  pedE4->SetLineStyle(7);
  pedW1->SetLineColor(kGreen);
  pedW1->SetLineStyle(7);
  pedW2->SetLineColor(kGreen);
  pedW2->SetLineStyle(7);
  pedW3->SetLineColor(kGreen);
  pedW3->SetLineStyle(7);
  pedW4->SetLineColor(kGreen);
  pedW4->SetLineStyle(7);

  probE1->Divide(triggerE1, totalE1);
  probE2->Divide(triggerE2, totalE2);
  probE3->Divide(triggerE3, totalE3);
  probE4->Divide(triggerE4, totalE4);

  probW1->Divide(triggerW1, totalW1);
  probW2->Divide(triggerW2, totalW2);
  probW3->Divide(triggerW3, totalW3);
  probW4->Divide(triggerW4, totalW4);


  probE1->SetMarkerStyle(20);
  probE2->SetMarkerStyle(20);
  probE3->SetMarkerStyle(20);
  probE4->SetMarkerStyle(20);
  probW1->SetMarkerStyle(20);
  probW2->SetMarkerStyle(20);
  probW3->SetMarkerStyle(20);
  probW4->SetMarkerStyle(20);

  probE1->SetTitle("Trigger Probabability PMT East 1");
  probE2->SetTitle("Trigger Probabability PMT East 2");
  probE3->SetTitle("Trigger Probabability PMT East 3");
  probE4->SetTitle("Trigger Probabability PMT East 4");
  probW1->SetTitle("Trigger Probabability PMT West 1");
  probW2->SetTitle("Trigger Probabability PMT West 2");
  probW3->SetTitle("Trigger Probabability PMT West 3");
  probW4->SetTitle("Trigger Probabability PMT West 4");

  probE1->GetXaxis()->SetTitle("ADC channels above pedestal");
  probE2->GetXaxis()->SetTitle("ADC channels above pedestal");
  probE3->GetXaxis()->SetTitle("ADC channels above pedestal");
  probE4->GetXaxis()->SetTitle("ADC channels above pedestal");
  probW1->GetXaxis()->SetTitle("ADC channels above pedestal");
  probW2->GetXaxis()->SetTitle("ADC channels above pedestal");
  probW3->GetXaxis()->SetTitle("ADC channels above pedestal");
  probW4->GetXaxis()->SetTitle("ADC channels above pedestal");

  

  
  //Fit the probability histograms and write to file
  TString filename = TString::Format("%s/runs/%i_thresholds.dat",getenv("TRIGGER_FUNC"),run);
  std::ofstream funcfile(filename.Data());
  
  std::vector <Double_t> params(8,0);
  
  params = fitThreshold(probE1);
  funcfile << params[0] << " " << params[1] << " " << params[2] << " " << params[3] 
  	   << " " << params[4] << " " << params[5] << " " << params[6] << " " << params[7] << std::endl;
  params = fitThreshold(probE2);
  funcfile << params[0] << " " << params[1] << " " << params[2] << " " << params[3] 
  	   << " " << params[4] << " " << params[5] << " " << params[6] << " " << params[7] << std::endl;
  params = fitThreshold(probE3);
  funcfile << params[0] << " " << params[1] << " " << params[2] << " " << params[3] 
  	   << " " << params[4] << " " << params[5] << " " << params[6] << " " << params[7] << std::endl;
  params =   fitThreshold(probE4);
  funcfile << params[0] << " " << params[1] << " " << params[2] << " " << params[3] 
  	   << " " << params[4] << " " << params[5] << " " << params[6] << " " << params[7] << std::endl;
  params = fitThreshold(probW1);
  funcfile << params[0] << " " << params[1] << " " << params[2] << " " << params[3] 
  	   << " " << params[4] << " " << params[5] << " " << params[6] << " " << params[7] << std::endl;
  params = fitThreshold(probW2);
  funcfile << params[0] << " " << params[1] << " " << params[2] << " " << params[3] 
  	   << " " << params[4] << " " << params[5] << " " << params[6] << " " << params[7] << std::endl;
  params = fitThreshold(probW3);
  funcfile << params[0] << " " << params[1] << " " << params[2] << " " << params[3] 
  	   << " " << params[4] << " " << params[5] << " " << params[6] << " " << params[7] << std::endl;
  params = fitThreshold(probW4);
  funcfile << params[0] << " " << params[1] << " " << params[2] << " " << params[3] 
  	   << " " << params[4] << " " << params[5] << " " << params[6] << " " << params[7] << std::endl;

  funcfile.close();

  probE1->SetMinimum(-0.1);
  probE2->SetMinimum(-0.1);
  probE3->SetMinimum(-0.1);
  probE4->SetMinimum(-0.1);
  probW1->SetMinimum(-0.1);
  probW2->SetMinimum(-0.1);
  probW3->SetMinimum(-0.1);
  probW4->SetMinimum(-0.1);


  TString pdf_file = TString::Format("%s/runs/%i_Thresholds.pdf",getenv("TRIGGER_FUNC"),run);

  TCanvas *cEtrigg = new TCanvas("cEtrigg"," ",2000,800);
  cEtrigg->Divide(4,1);

  cEtrigg->cd(1);
  gPad->SetLogy();
  totalE1->Draw();
  triggerE1->Draw("SAME");
  pedE1->Draw("SAME");
  cEtrigg->cd(2);
  gPad->SetLogy();
  totalE2->Draw();
  triggerE2->Draw("SAME");
  pedE2->Draw("SAME");
  cEtrigg->cd(3);
  gPad->SetLogy();
  totalE3->Draw();
  triggerE3->Draw("SAME");
  pedE3->Draw("SAME");
  cEtrigg->cd(4);
  gPad->SetLogy();
  totalE4->Draw();
  triggerE4->Draw("SAME");
  pedE4->Draw("SAME");

  cEtrigg->Print(pdf_file+"(");

  TCanvas *cWtrigg = new TCanvas("cWtrigg"," ",2000,800);
  cWtrigg->Divide(4,1);

  cWtrigg->cd(1);
  gPad->SetLogy();
  totalW1->Draw();
  triggerW1->Draw("SAME");
  pedW1->Draw("SAME");
  cWtrigg->cd(2);
  gPad->SetLogy();
  totalW2->Draw();
  triggerW2->Draw("SAME");
  pedW2->Draw("SAME");
  cWtrigg->cd(3);
  gPad->SetLogy();
  totalW3->Draw();
  triggerW3->Draw("SAME");
  pedW3->Draw("SAME");
  cWtrigg->cd(4);
  gPad->SetLogy();
  totalW4->Draw();
  triggerW4->Draw("SAME");
  pedW4->Draw("SAME");

  cWtrigg->Print(pdf_file);

  TCanvas *cEprob = new TCanvas("cEprob"," ",2000,800);
  cEprob->Divide(4,1);

  cEprob->cd(1);
  probE1->Draw("AP");
  cEprob->cd(2);
  probE2->Draw("AP");
  cEprob->cd(3);
  probE3->Draw("AP");
  cEprob->cd(4);
  probE4->Draw("AP");

  cEprob->Print(pdf_file);

  TCanvas *cWprob = new TCanvas("cWprob"," ",2000,800);
  cWprob->Divide(4,1);

  cWprob->cd(1);
  probW1->Draw("AP");
  cWprob->cd(2);
  probW2->Draw("AP");
  cWprob->cd(3);
  probW3->Draw("AP");
  cWprob->cd(4);
  probW4->Draw("AP");

  cWprob->Print(pdf_file+")");

  //Clean-up
  delete cEtrigg; delete cWtrigg; delete cEprob; delete cWprob;

  delete triggerE1; delete triggerE2; delete triggerE3; delete triggerE4;
  delete triggerW1; delete triggerW2; delete triggerW3; delete triggerW4;

  delete totalE1; delete totalE2; delete totalE3; delete totalE4;
  delete totalW1; delete totalW2; delete totalW3; delete totalW4;

  delete probE1; delete probE2; delete probE3; delete probE4;
  delete probW1; delete probW2; delete probW3; delete probW4;


  /*//Process Pedestal Histograms

  Double_t pedQadc[8] = {0.}; //PMT pedestals
  
  pedQadc[0] = pedE1->GetXaxis()->GetBinCenter(pedE1->GetMaximumBin()); pedE1->GetXaxis()->SetRangeUser(pedQadc[0]-50., pedQadc[0]+50.);
  pedQadc[1] = pedE2->GetXaxis()->GetBinCenter(pedE2->GetMaximumBin()); pedE2->GetXaxis()->SetRangeUser(pedQadc[1]-50., pedQadc[1]+50.);
  pedQadc[2] = pedE3->GetXaxis()->GetBinCenter(pedE3->GetMaximumBin()); pedE3->GetXaxis()->SetRangeUser(pedQadc[2]-50., pedQadc[2]+50.);
  pedQadc[3] = pedE4->GetXaxis()->GetBinCenter(pedE4->GetMaximumBin()); pedE4->GetXaxis()->SetRangeUser(pedQadc[3]-50., pedQadc[3]+50.);
  pedQadc[4] = pedW1->GetXaxis()->GetBinCenter(pedW1->GetMaximumBin()); pedW1->GetXaxis()->SetRangeUser(pedQadc[4]-50., pedQadc[4]+50.);
  pedQadc[5] = pedW2->GetXaxis()->GetBinCenter(pedW2->GetMaximumBin()); pedW2->GetXaxis()->SetRangeUser(pedQadc[5]-50., pedQadc[5]+50.);
  pedQadc[6] = pedW3->GetXaxis()->GetBinCenter(pedW3->GetMaximumBin()); pedW3->GetXaxis()->SetRangeUser(pedQadc[6]-50., pedQadc[6]+50.);
  pedQadc[7] = pedW4->GetXaxis()->GetBinCenter(pedW4->GetMaximumBin()); pedW4->GetXaxis()->SetRangeUser(pedQadc[7]-50., pedQadc[7]+50.);

  //calc mean again after Adjusting Range

  std::ofstream outWidthFile(TString::Format("%s/PMT_pedestals_%i.dat", getenv("PEDESTALS"),run).Data());

  outWidthFile << std::fixed << std::setprecision(7);


  outWidthFile << run << " " << pedE1->GetMean() << " " << pedE1->GetRMS() << std::endl;
  outWidthFile << run << " " << pedE2->GetMean() << " " << pedE2->GetRMS() << std::endl;
  outWidthFile << run << " " << pedE3->GetMean() << " " << pedE3->GetRMS() << std::endl;
  outWidthFile << run << " " << pedE4->GetMean() << " " << pedE4->GetRMS() << std::endl;
  outWidthFile << run << " " << pedW1->GetMean() << " " << pedW1->GetRMS() << std::endl;
  outWidthFile << run << " " << pedW2->GetMean() << " " << pedW2->GetRMS() << std::endl;
  outWidthFile << run << " " << pedW3->GetMean() << " " << pedW3->GetRMS() << std::endl;
  outWidthFile << run << " " << pedW4->GetMean() << " " << pedW4->GetRMS();
  
  outWidthFile.close();
  */

  delete pedE1; delete pedE2; delete pedE3; delete pedE4;
  delete pedW1; delete pedW2; delete pedW3; delete pedW4;

}


void loadCuts(Int_t runNumber, cuts* Cuts) {

  // Read cuts file
  char tempFileCuts[500];
  sprintf(tempFileCuts, "%s/cuts_%i.dat", getenv("CUTS"),runNumber);
  std::cout << "... Reading: " << tempFileCuts << std::endl;

  std::ifstream fileCuts(tempFileCuts);
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

std::vector <Double_t> fitThreshold(TGraphAsymmErrors *trigg, std::vector <Double_t> params) {

  bool defaultParams = true;
  if (params.size()==8) defaultParams = false;

  TF1 *erf = new TF1("erf","([0]+[1]*TMath::Erf((x-[2])/[3]))*(0.5-.5*TMath::TanH((x-[2])/[4]))+(0.5+.5*TMath::TanH((x-[2])/[4]))*([5]+[6]*TMath::TanH((x-[2])/[7]))",-50.,150.);
  //TF1 *erf = new TF1("erf","([0]+[1]*TMath::Erf((x-[2])/[3]))",-30.,150.);

  if (defaultParams) {
    params.resize(8);
    erf->SetParameter(0,0.5); //Constant offset of erf                                                                                       
    erf->SetParameter(1,0.5); //Scaling of erf                                                                                               
    erf->SetParameter(2,30.); //Mean of gaussian integrated for erf                                                                          
    erf->SetParameter(3,7.); //std. dev. of gaussian integrated for erf                                                                      
    erf->SetParameter(4,9.); //severity of transition function "turn on"                                                                     
    erf->SetParameter(5,0.5); //constant offset of second tanh                                                                               
    erf->SetParameter(6,0.5); //Scaling of second tanh                                                                                        
    erf->SetParameter(7,12.); // stretching factor of second tanh                                                                            
    
  }
  else {
    erf->SetParameter(0,params[0]); //Constant offset of erf                                                                            
    erf->SetParameter(1,params[1]); //Scaling of erf                                                                                   
    erf->SetParameter(2,params[2]); //Mean of gaussian integrated for erf                                                               
    erf->SetParameter(3,params[3]); //std. dev. of gaussian integrated for erf                                                         
    erf->SetParameter(4,params[4]); //severity of transition function "turn on"                                                        
    erf->SetParameter(5,params[5]); //constant offset of second tanh                                                                   
    erf->SetParameter(6,params[6]); //Scaling of second tanh                                                                           
    erf->SetParameter(7,params[7]); // stretching factor of second tanh
  }
  
  erf->FixParameter(0,0.5);
  erf->FixParameter(1,0.5);
  erf->SetParLimits(3,1.,100.);
  erf->SetParLimits(4,1.,100.);
  erf->FixParameter(5,0.5);                                                                                                              
  erf->FixParameter(6,0.5);
  //erf->SetParLimits(7,0.,50.); //Constant offset of erf     
  
    
  trigg->Fit("erf","R");

  Int_t attempt = 0;

  while ( gMinuit->fCstatu!=TString("CONVERGED ") && attempt<20 ) {

    TRandom3 rand(0);
                                                                                            
    erf->SetParameter(2, rand.Gaus(30., 15.)); //Mean of gaussian integrated for erf                                                                         
    erf->SetParameter(3, rand.Gaus(10.,2.)); //std. dev. of gaussian integrated for erf    
                                                                  
    erf->SetParameter(4, rand.Gaus(15.,7.)); //severity of transition function "turn on"                                                                     
                              
    erf->SetParameter(7,rand.Gaus(12.,3.)); // stretching factor of second tanh 
    
    erf->SetParLimits(3,3.,100.);
    
    trigg->Fit("erf","R");

    attempt+=1;
  }

  for (UInt_t i=0; i<params.size();i++) params[i] = erf->GetParameter(i);

  delete erf;

  return params;
}

std::vector <Double_t> loadGainFactors(Int_t runNumber) {

  // Read gain corrections file                                                                                                                
  char tempFileGain[500];
  sprintf(tempFileGain, "%s/gain_bismuth_%i.dat",getenv("GAIN_BISMUTH"), runNumber);
  std::cout << "... Reading: " << tempFileGain << std::endl;

  std::vector <Double_t> gainCorrection(8,1.);
  std::vector <Double_t> fitMean(8,1.);
  std::ifstream fileGain(tempFileGain);
  for (int i=0; i<8; i++) {
    fileGain >> fitMean[i] >> gainCorrection[i];
  }
  std::cout << "...   PMT E1: " << gainCorrection[0] << std::endl;
  std::cout << "...   PMT E2: " << gainCorrection[1] << std::endl;
  std::cout << "...   PMT E3: " << gainCorrection[2] << std::endl;
  std::cout << "...   PMT E4: " << gainCorrection[3] << std::endl;
  std::cout << "...   PMT W1: " << gainCorrection[4] << std::endl;
  std::cout << "...   PMT W2: " << gainCorrection[5] << std::endl;
  std::cout << "...   PMT W3: " << gainCorrection[6] << std::endl;
  std::cout << "...   PMT W4: " << gainCorrection[7] << std::endl;

  return gainCorrection;
}
