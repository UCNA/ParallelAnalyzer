#include <iostream>
#include <vector>
#include <cstdlib>

void trigger_threshold(Int_t XeRunPeriod, bool mpmData=false) {
  
  Char_t temp[500];
  //Create file to write parameters of fit.
  sprintf(temp,"%s/trigger_functions_XePeriod_%i.dat",getenv("TRIGGER_FUNC"),XeRunPeriod);
  ofstream triggFunc(temp);

  // Read in the Xe runs in this runPeriod
  vector <Int_t> XeRuns;
  sprintf(temp,"../run_lists/Xenon_Calibration_Run_Period_%i.dat",XeRunPeriod);
  ifstream XeRunList(temp);

  Int_t run, numRuns=0;
  while (XeRunList >> run) {
    XeRuns.push_back(run);
    numRuns++;
    //if (numRuns>30) break;
  }

  if (numRuns!=XeRuns.size()) {
    cout << "Runs!=XeRuns.size()\n";
    exit(0);
  }
  XeRunList.close();
  cout << "There are " << numRuns << " in this Xe run period\n";

  TChain *chain;
  if (mpmData) chain = new TChain("phys");
  else chain = new TChain("pass3");

  if (mpmData) {
    for (Int_t i=0; i<numRuns; i++) {
      sprintf(temp,"%s/hists/spec_%i.root",getenv("UCNAOUTPUTDIR"),XeRuns[i]);
      chain->AddFile(temp);
    }
  }
  else {
    for (Int_t i=0; i<numRuns; i++) {
      sprintf(temp,"%s/replay_pass3_%i.root",getenv("REPLAY_PASS3"),XeRuns[i]);
      chain->AddFile(temp);
    }
  }
 
  //sprintf(temp,"/extern/UCNA/replay_pass4_MB/replay_pass4_%i.root",runNumber);
  //TFile *file = new TFile(temp,"READ");

  //TTree *t = (TTree*)file->Get("pass4");

  //Int_t nbins = 75;
  Double_t lower_limit = -20.;
  Int_t binWidth = 1;
  Double_t East_upper_limit = 150.;
  Double_t West_upper_limit = 150.;
  Int_t nbinsE = (int)(East_upper_limit/binWidth);
  Int_t nbinsW = (int)(West_upper_limit/binWidth);
  TH1F *Etrigg = new TH1F("Etrigg","West Type 1: EvisE",nbinsE,lower_limit,East_upper_limit);
  TH1F *EnoTrigg = new TH1F("EnoTrigg","West Type 2/3: EvisE",nbinsE,lower_limit,East_upper_limit);
  TH1F *Etotal = new TH1F("Etotal","West Type 1,2/3: EvisE",nbinsE,lower_limit,East_upper_limit);
  TH1F *EtriggFunc = new TH1F("EtriggFunc","East Trigger Probability",nbinsE,lower_limit,East_upper_limit);
  EtriggFunc->SetMarkerStyle(20);

  TH1F *Wtrigg = new TH1F("Wtrigg","East Type 1: EvisW",nbinsW,lower_limit,West_upper_limit);
  TH1F *WnoTrigg = new TH1F("WnoTrigg","East Type 2/3: EvisW",nbinsW,lower_limit,West_upper_limit);
  TH1F *Wtotal = new TH1F("Wtotal","East Type 1,2/3: EvisW",nbinsW,lower_limit,West_upper_limit);
  TH1F *WtriggFunc = new TH1F("WtriggFunc","West Trigger Probability",nbinsW,lower_limit,West_upper_limit);
  WtriggFunc->SetMarkerStyle(20);

  //TF1 *erf = new TF1("erf","([5]*TMath::Erf((x-[0])/[1])+0.5)+[2]*TMath::Gaus(x,[3],[4])",0.,150.);
  //TF1 *erf = new TF1("erf","([0]+[1]*TMath::Erf((x-[2])/[3]))+[4]*TMath::Gaus(x,[5],[6])",0.,150.);
  //TF1 *erf = new TF1("erf","([0]+[1]*TMath::Erf((x-[2])/[3])*TMath::TanH([7]*x)+[4]*TMath::Gaus(x,[5],[6]))",0.,150.);

  // Best fit function. This is a shifted erf which goes into a shifted tanh, where the smooth transition is done via application of another shifted tanh
  TF1 *erf = new TF1("erf","([0]+[1]*TMath::Erf((x-[2])/[3]))*(0.5-.5*TMath::TanH((x-[2])/[4]))+(0.5+.5*TMath::TanH((x-[2])/[4]))*([5]+[6]*TMath::TanH((x-[2])/[7]))",-20.,150.);
  erf->SetParameter(0,0.5); //Constant offset of erf
  erf->FixParameter(0,0.5); 
  erf->SetParameter(1,0.5); //Scaling of erf
  erf->FixParameter(1,0.5); 
  erf->SetParameter(2,30.); //Mean of gaussian integrated for erf
  erf->SetParameter(3,7..); //std. dev. of gaussian integrated for erf
  erf->SetParameter(4,15.); //severity of transition function "turn on"
  erf->SetParLimits(4,1.,25.);
  erf->SetParameter(5,0.5); //constant offset of second tanh
  erf->SetParameter(6,0.5); //Scaling of second tanh
  //erf->FixParameter(5,0.5); 
  //erf->FixParameter(6,0.5); 
  erf->SetParameter(7,20.); // stretching factor of second tanh

  //Working Parameters of the east trigger function for erf*tanh+gaus
  /*erf->SetParameter(0,0.5); //Constant offset of erf
  erf->SetParameter(1,0.5); //Scaling of erf
  erf->SetParameter(2,30.); //Mean of gaussian integrated for erf
  erf->SetParameter(3,7..); //std. dev. of gaussian integrated for erf
  erf->SetParameter(4,6.); //severity of transition function "turn on"
  erf->SetParLimits(4,1.,25.);
  erf->SetParameter(5,0.5); //constant offset of second tanh
  erf->SetParameter(6,0.5); //Scaling of second tanh
  erf->SetParameter(7,13.); // stretching factor of second tanh
  */
  TCanvas *c1 = new TCanvas("c1"," ",1200.,1600.);
  c1->Divide(2,2);
  c1->cd(1);
  chain->Draw("EvisE>>Etrigg","(Type==1 && Side==1) || (Type==0 && Side==0) && PID==1");
  c1->cd(2);  
  chain->Draw("EvisE>>EnoTrigg","Type==2 && Side==1 && PID==1");
  c1->cd(3);
  Etotal->Add(Etrigg,EnoTrigg);
  Etotal->Draw();
  c1->cd(4);
  EtriggFunc->Divide(Etrigg,Etotal);
  EtriggFunc->SetStats(0);
  //Etrigg->Draw("P");

  EtriggFunc->Fit("erf","","",0.,East_upper_limit);
  
  EtriggFunc->Draw("P");
  triggFunc << erf->GetParameter(0) << " " << erf->GetParameter(1) << " " << erf->GetParameter(2) << " " << erf->GetParameter(3) << " " << erf->GetParameter(4) << " " << erf->GetParameter(5) << " " << erf->GetParameter(6) << " " << erf->GetParameter(7) << endl;


  TCanvas *c2 = new TCanvas("c2"," ",1200.,1600.);
  c2->Divide(2,2);
  c2->cd(1);
  chain->Draw("EvisW>>Wtrigg","(Type==1 && Side==0) || (Type==0 && Side==1) && PID==1");
  c2->cd(2);  
  chain->Draw("EvisW>>WnoTrigg","Type==2 && Side==0 && PID==1");
  c2->cd(3);
  Wtotal->Add(Wtrigg,WnoTrigg);
  Wtotal->Draw();
  c2->cd(4);
  WtriggFunc->Divide(Wtrigg,Wtotal);
  WtriggFunc->SetStats(0);
  //Wtrigg->Draw("P");

  //erf->SetParameter(0,0.5); //Constant offset of erf
  //erf->SetParameter(1,0.5); //Scaling of erf
  erf->SetParameter(2,23.); //Mean of gaussian integrated for erf
  erf->SetParameter(3,8.); //std. dev. of gaussian integrated for erf
  erf->SetParameter(4,10.); //severity of transition function "turn on"
  erf->SetParameter(5,0.5); //constant offset of second tanh
  erf->SetParameter(6,0.5); //Scaling of second tanh
  erf->SetParameter(7,8.); // stretching factor of second tanh

  //Old Parameters
  //erf->SetParameter(0,0.5); //Constant offset of erf
  //erf->SetParameter(1,0.5); //Scaling of erf
  //erf->SetParameter(2,10.); //Mean of gaussian integrated for erf
  //erf->SetParameter(3,5.); //std. dev. of gaussian integrated for erf
  //erf->SetParameter(4,-0.1); //Scaling of additional gaussian for cancellation of erf overshooting distribution
  //erf->SetParameter(5,14.); //Mean of additional gaussian
  //erf->SetParameter(6,6.); //Std. dev of additional gaussian
  //erf->SetParameter(7,25.);
  //erf->SetParameter(7,0.036);
  //erf->SetParameter(8,25.);


  WtriggFunc->Fit("erf","","",-20.,West_upper_limit);
  //WtriggFunc->Fit("erf","R");

  WtriggFunc->Draw("P");

  triggFunc << erf->GetParameter(0) << " " << erf->GetParameter(1) << " " << erf->GetParameter(2) << " " << erf->GetParameter(3) << " " << erf->GetParameter(4) << " " << erf->GetParameter(5) << " " << erf->GetParameter(6) << " " << erf->GetParameter(7) << endl;

  triggFunc.close();
  
  //Write the canvases to file
  sprintf(temp,"%i",XeRunPeriod);
  TString canvas1 = TString(getenv("TRIGGER_FUNC")) + TString("trigger_functions_XePeriod_")+TString(temp)+TString(".pdf(");
  TString canvas2 = TString(getenv("TRIGGER_FUNC")) + TString("trigger_functions_XePeriod_")+TString(temp)+TString(".pdf)");
  c1->Print(canvas1);
  c2->Print(canvas2);
  
}
