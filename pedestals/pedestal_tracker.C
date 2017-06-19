void pedestal_tracker(Int_t run) {
  Int_t runNumber = run;//22822;//22312;//19381;//20519;//20574;
  
  TFile *fileIn = new TFile(TString::Format("/extern/mabrow05/ucna/rawdata/full%i.root",runNumber), "READ");
  TTree *Tin = (TTree*)(fileIn->Get("h1"));
  
  cout << "Here\n";

  TH1F *hisOppBiE1 = new TH1F("hisOppBiE1","Opposite Side Bi Pulser",301, -0.5, 300.5);
  TH1F *hisOpp2foldE1 = new TH1F("hisOpp2foldE1","Opposite Side 2-fold",301, -0.5, 300.5);
  TH1F *hisMonE1 = new TH1F("hisMonE1","UCN Monitor Peds",301, -0.5, 300.5);
  TH1F *hisMPME1 = new TH1F("hisMPME1","Opposite 2-fold + UCNMon",301, -0.5, 300.5);

  hisOppBiE1->SetLineColor(1);
  hisOpp2foldE1->SetLineColor(3);
  hisMonE1->SetLineColor(4);
  hisMPME1->SetLineColor(2);

  TH1F *hisOppBiE2 = new TH1F("hisOppBiE2","Opposite Side Bi Pulser",301, -0.5, 300.5);
  TH1F *hisOpp2foldE2 = new TH1F("hisOpp2foldE2","Opposite Side 2-fold",301, -0.5, 300.5);
  TH1F *hisMonE2 = new TH1F("hisMonE2","UCN Monitor Peds",301, -0.5, 300.5);
  TH1F *hisMPME2 = new TH1F("hisMPME2","Opposite 2-fold + UCNMon",301, -0.5, 300.5);

  hisOppBiE2->SetLineColor(1);
  hisOpp2foldE2->SetLineColor(3);
  hisMonE2->SetLineColor(4);
  hisMPME2->SetLineColor(2);

  TH1F *hisOppBiE3 = new TH1F("hisOppBiE3","Opposite Side Bi Pulser",301, -0.5, 300.5);
  TH1F *hisOpp2foldE3 = new TH1F("hisOpp2foldE3","Opposite Side 2-fold",301, -0.5, 300.5);
  TH1F *hisMonE3 = new TH1F("hisMonE3","UCN Monitor Peds",301, -0.5, 300.5);
  TH1F *hisMPME3 = new TH1F("hisMPME3","Opposite 2-fold + UCNMon",301, -0.5, 300.5);

  hisOppBiE3->SetLineColor(1);
  hisOpp2foldE3->SetLineColor(3);
  hisMonE3->SetLineColor(4);
  hisMPME3->SetLineColor(2);

  TH1F *hisOppBiE4 = new TH1F("hisOppBiE4","Opposite Side Bi Pulser",301, -0.5, 300.5);
  TH1F *hisOpp2foldE4 = new TH1F("hisOpp2foldE4","Opposite Side 2-fold",301, -0.5, 300.5);
  TH1F *hisMonE4 = new TH1F("hisMonE4","UCN Monitor Peds",301, -0.5, 300.5);
  TH1F *hisMPME4 = new TH1F("hisMPME4","Opposite 2-fold + UCNMon",301, -0.5, 300.5);

  hisOppBiE4->SetLineColor(1);
  hisOpp2foldE4->SetLineColor(3);
  hisMonE4->SetLineColor(4);
  hisMPME4->SetLineColor(2);

  TH1F *hisOppBiW1 = new TH1F("hisOppBiW1","Opposite Side Bi Pulser",301, -0.5, 300.5);
  TH1F *hisOpp2foldW1 = new TH1F("hisOpp2foldW1","Opposite Side 2-fold",301, -0.5, 300.5);
  TH1F *hisMonW1 = new TH1F("hisMonW1","UCN Monitor Peds",301, -0.5, 300.5);
  TH1F *hisMPMW1 = new TH1F("hisMPMW1","Opposite 2-fold + UCNMon",301, -0.5, 300.5);

  hisOppBiW1->SetLineColor(1);
  hisOpp2foldW1->SetLineColor(3);
  hisMonW1->SetLineColor(4);
  hisMPMW1->SetLineColor(2);

  TH1F *hisOppBiW2 = new TH1F("hisOppBiW2","Opposite Side Bi Pulser",301, -0.5, 300.5);
  TH1F *hisOpp2foldW2 = new TH1F("hisOpp2foldW2","Opposite Side 2-fold",301, -0.5, 300.5);
  TH1F *hisMonW2 = new TH1F("hisMonW2","UCN Monitor Peds",301, -0.5, 300.5);
  TH1F *hisMPMW2 = new TH1F("hisMPMW2","Opposite 2-fold + UCNMon",301, -0.5, 300.5);

  hisOppBiW2->SetLineColor(1);
  hisOpp2foldW2->SetLineColor(3);
  hisMonW2->SetLineColor(4);
  hisMPMW2->SetLineColor(2);

  TH1F *hisOppBiW3 = new TH1F("hisOppBiW3","Opposite Side Bi Pulser",301, -0.5, 300.5);
  TH1F *hisOpp2foldW3 = new TH1F("hisOpp2foldW3","Opposite Side 2-fold",301, -0.5, 300.5);
  TH1F *hisMonW3 = new TH1F("hisMonW3","UCN Monitor Peds",301, -0.5, 300.5);
  TH1F *hisMPMW3 = new TH1F("hisMPMW3","Opposite 2-fold + UCNMon",301, -0.5, 300.5);

  hisOppBiW3->SetLineColor(1);
  hisOpp2foldW3->SetLineColor(3);
  hisMonW3->SetLineColor(4);
  hisMPMW3->SetLineColor(2);

  TH1F *hisOppBiW4 = new TH1F("hisOppBiW4","Opposite Side Bi Pulser",301, -0.5, 300.5);
  TH1F *hisOpp2foldW4 = new TH1F("hisOpp2foldW4","Opposite Side 2-fold",301, -0.5, 300.5);
  TH1F *hisMonW4 = new TH1F("hisMonW4","UCN Monitor Peds",301, -0.5, 300.5);
  TH1F *hisMPMW4 = new TH1F("hisMPMW4","Opposite 2-fold + UCNMon",301, -0.5, 300.5);

  hisOppBiW4->SetLineColor(1);
  hisOpp2foldW4->SetLineColor(3);
  hisMonW4->SetLineColor(4);
  hisMPMW4->SetLineColor(2);
  

  
 
    

  TCanvas *cE1 = new TCanvas("cE1","cE1",900, 900);
  cE1->Divide(2,2);

  
  cE1->cd(1);
  Tin->Draw("Qadc0>>hisOppBiE1","(int(Sis00)==32) && Tdc00<0.0001 && Tdc01<0.0001 && Tdc02<0.0001 && Tdc03<0.0001");
  cE1->cd(2);
  Tin->Draw("Qadc0>>hisOpp2foldE1","(int(Sis00)==2 || int(Sis00)==1 || int(Sis00)==3) && Tdc00<0.0001 ");
  cE1->cd(3);
  Tin->Draw("Qadc0>>hisMonE1","Tdc00<0.0001 && (int(Sis00) == 260 || int(Sis00) == 516 || int(Sis00) == 1028 || int(Sis00) == 2052 || int(Sis00) == (2052+260) || int(Sis00) == (2052+516) || int(Sis00) == (2052+1028) || int(Sis00) == (260+516) || int(Sis00) == (260+1028) || int(Sis00) == (516+1028) )");
  cE1->cd(4);
  Tin->Draw("Qadc0>>hisMPME1","Tdc00<0.0001 && (int(Sis00) == 260 || int(Sis00) == 516 || int(Sis00) == 1028 || int(Sis00) == 2052 || int(Sis00) == (2052+260) || int(Sis00) == (2052+516) || int(Sis00) == (2052+1028) || int(Sis00) == (260+516) || int(Sis00) == (260+1028) || int(Sis00) == (516+1028) || int(Sis00)==2)");


  TCanvas *cE2 = new TCanvas("cE2","cE2",900, 900);
  cE2->Divide(2,2);

  
  cE2->cd(1);
  Tin->Draw("Qadc1>>hisOppBiE2","(int(Sis00)==32) && Tdc00<0.0001 && Tdc01<0.0001 && Tdc02<0.0001 && Tdc03<0.0001");
  cE2->cd(2);
  Tin->Draw("Qadc1>>hisOpp2foldE2","(int(Sis00)==2 || int(Sis00)==1 || int(Sis00)==3) && Tdc01<0.0001 ");
  cE2->cd(3);
  Tin->Draw("Qadc1>>hisMonE2","Tdc01<0.0001 && (int(Sis00) == 260 || int(Sis00) == 516 || int(Sis00) == 1028 || int(Sis00) == 2052 || int(Sis00) == (2052+260) || int(Sis00) == (2052+516) || int(Sis00) == (2052+1028) || int(Sis00) == (260+516) || int(Sis00) == (260+1028) || int(Sis00) == (516+1028) )");
  cE2->cd(4);
  Tin->Draw("Qadc1>>hisMPME2","Tdc01<0.0001 && (int(Sis00) == 260 || int(Sis00) == 516 || int(Sis00) == 1028 || int(Sis00) == 2052 || int(Sis00) == (2052+260) || int(Sis00) == (2052+516) || int(Sis00) == (2052+1028) || int(Sis00) == (260+516) || int(Sis00) == (260+1028) || int(Sis00) == (516+1028) || int(Sis00)==2)");

  TCanvas *cE3 = new TCanvas("cE3","cE3",900, 900);
  cE3->Divide(2,2);

  
  cE3->cd(1);
  Tin->Draw("Qadc2>>hisOppBiE3","(int(Sis00)==32) && Tdc00<0.0001 && Tdc01<0.0001 && Tdc02<0.0001 && Tdc03<0.0001");
  cE3->cd(2);
  Tin->Draw("Qadc2>>hisOpp2foldE3","(int(Sis00)==2 || int(Sis00)==1 || int(Sis00)==3) && Tdc02<0.0001 ");
  cE3->cd(3);
  Tin->Draw("Qadc2>>hisMonE3","Tdc02<0.0001 && (int(Sis00) == 260 || int(Sis00) == 516 || int(Sis00) == 1028 || int(Sis00) == 2052 || int(Sis00) == (2052+260) || int(Sis00) == (2052+516) || int(Sis00) == (2052+1028) || int(Sis00) == (260+516) || int(Sis00) == (260+1028) || int(Sis00) == (516+1028) )");
  cE3->cd(4);
  Tin->Draw("Qadc2>>hisMPME3","Tdc02<0.0001 && (int(Sis00) == 260 || int(Sis00) == 516 || int(Sis00) == 1028 || int(Sis00) == 2052 || int(Sis00) == (2052+260) || int(Sis00) == (2052+516) || int(Sis00) == (2052+1028) || int(Sis00) == (260+516) || int(Sis00) == (260+1028) || int(Sis00) == (516+1028) || int(Sis00)==2)");

  TCanvas *cE4 = new TCanvas("cE4","cE4",900, 900);
  cE4->Divide(2,2);

  
  cE4->cd(1);
  Tin->Draw("Qadc3>>hisOppBiE4","(int(Sis00)==32) && Tdc00<0.0001 && Tdc01<0.0001 && Tdc02<0.0001 && Tdc03<0.0001");
  cE4->cd(2);
  Tin->Draw("Qadc3>>hisOpp2foldE4","(int(Sis00)==2 || int(Sis00)==1 || int(Sis00)==3) && Tdc03<0.0001 ");
  cE4->cd(3);
  Tin->Draw("Qadc3>>hisMonE4","Tdc03<0.0001 && (int(Sis00) == 260 || int(Sis00) == 516 || int(Sis00) == 1028 || int(Sis00) == 2052 || int(Sis00) == (2052+260) || int(Sis00) == (2052+516) || int(Sis00) == (2052+1028) || int(Sis00) == (260+516) || int(Sis00) == (260+1028) || int(Sis00) == (516+1028) )");
  cE4->cd(4);
  Tin->Draw("Qadc3>>hisMPME4","Tdc03<0.0001 && (int(Sis00) == 260 || int(Sis00) == 516 || int(Sis00) == 1028 || int(Sis00) == 2052 || int(Sis00) == (2052+260) || int(Sis00) == (2052+516) || int(Sis00) == (2052+1028) || int(Sis00) == (260+516) || int(Sis00) == (260+1028) || int(Sis00) == (516+1028) || int(Sis00)==2)");
  

  //West side


  TCanvas *cW1 = new TCanvas("cW1","cW1",900, 900);
  cW1->Divide(2,2);

  
  cW1->cd(1);
  Tin->Draw("Qadc4>>hisOppBiW1","(int(Sis00)==32) && Tdc08<0.0001 && Tdc09<0.0001 && Tdc014<0.0001 && Tdc011<0.0001");
  cW1->cd(2);
  Tin->Draw("Qadc4>>hisOpp2foldW1","(int(Sis00)==2 || int(Sis00)==1 || int(Sis00)==3) && Tdc08<0.0001 ");
  cW1->cd(3);
  Tin->Draw("Qadc4>>hisMonW1","Tdc08<0.0001 && (int(Sis00) == 260 || int(Sis00) == 516 || int(Sis00) == 1028 || int(Sis00) == 2052 || int(Sis00) == (2052+260) || int(Sis00) == (2052+516) || int(Sis00) == (2052+1028) || int(Sis00) == (260+516) || int(Sis00) == (260+1028) || int(Sis00) == (516+1028) )");
  cW1->cd(4);
  Tin->Draw("Qadc4>>hisMPMW1","Tdc08<0.0001 && (int(Sis00) == 260 || int(Sis00) == 516 || int(Sis00) == 1028 || int(Sis00) == 2052 || int(Sis00) == (2052+260) || int(Sis00) == (2052+516) || int(Sis00) == (2052+1028) || int(Sis00) == (260+516) || int(Sis00) == (260+1028) || int(Sis00) == (516+1028) || int(Sis00)==1)");


  TCanvas *cW2 = new TCanvas("cW2","cW2",900, 900);
  cW2->Divide(2,2);

  
  cW2->cd(1);
  Tin->Draw("Qadc5>>hisOppBiW2","(int(Sis00)==32) && Tdc08<0.0001 && Tdc09<0.0001 && Tdc014<0.0001 && Tdc011<0.0001");
  cW2->cd(2);
  Tin->Draw("Qadc5>>hisOpp2foldW2","(int(Sis00)==2 || int(Sis00)==1 || int(Sis00)==3) && Tdc09<0.0001 ");
  cW2->cd(3);
  Tin->Draw("Qadc5>>hisMonW2","Tdc09<0.0001 && (int(Sis00) == 260 || int(Sis00) == 516 || int(Sis00) == 1028 || int(Sis00) == 2052 || int(Sis00) == (2052+260) || int(Sis00) == (2052+516) || int(Sis00) == (2052+1028) || int(Sis00) == (260+516) || int(Sis00) == (260+1028) || int(Sis00) == (516+1028) )");
  cW2->cd(4);
  Tin->Draw("Qadc5>>hisMPMW2","Tdc09<0.0001 && (int(Sis00) == 260 || int(Sis00) == 516 || int(Sis00) == 1028 || int(Sis00) == 2052 || int(Sis00) == (2052+260) || int(Sis00) == (2052+516) || int(Sis00) == (2052+1028) || int(Sis00) == (260+516) || int(Sis00) == (260+1028) || int(Sis00) == (516+1028) || int(Sis00)==1)");

  TCanvas *cW3 = new TCanvas("cW3","cW3",900, 900);
  cW3->Divide(2,2);

  
  cW3->cd(1);
  Tin->Draw("Qadc6>>hisOppBiW3","(int(Sis00)==32) && Tdc08<0.0001 && Tdc09<0.0001 && Tdc014<0.0001 && Tdc011<0.0001");
  cW3->cd(2);
  Tin->Draw("Qadc6>>hisOpp2foldW3","(int(Sis00)==2 || int(Sis00)==1 || int(Sis00)==3) && Tdc014<0.0001 ");
  cW3->cd(3);
  Tin->Draw("Qadc6>>hisMonW3","Tdc014<0.0001 && (int(Sis00) == 260 || int(Sis00) == 516 || int(Sis00) == 1028 || int(Sis00) == 2052 || int(Sis00) == (2052+260) || int(Sis00) == (2052+516) || int(Sis00) == (2052+1028) || int(Sis00) == (260+516) || int(Sis00) == (260+1028) || int(Sis00) == (516+1028) )");
  cW3->cd(4);
  Tin->Draw("Qadc6>>hisMPMW3","Tdc014<0.0001 && (int(Sis00) == 260 || int(Sis00) == 516 || int(Sis00) == 1028 || int(Sis00) == 2052 || int(Sis00) == (2052+260) || int(Sis00) == (2052+516) || int(Sis00) == (2052+1028) || int(Sis00) == (260+516) || int(Sis00) == (260+1028) || int(Sis00) == (516+1028) || int(Sis00)==1)");

  TCanvas *cW4 = new TCanvas("cW4","cW4",900, 900);
  cW4->Divide(2,2);

  
  cW4->cd(1);
  Tin->Draw("Qadc7>>hisOppBiW4","(int(Sis00)==32) && Tdc08<0.0001 && Tdc09<0.0001 && Tdc014<0.0001 && Tdc011<0.0001");
  cW4->cd(2);
  Tin->Draw("Qadc7>>hisOpp2foldW4","(int(Sis00)==2 || int(Sis00)==1 || int(Sis00)==3) && Tdc011<0.0001 ");
  cW4->cd(3);
  Tin->Draw("Qadc7>>hisMonW4","Tdc011<0.0001 && (int(Sis00) == 260 || int(Sis00) == 516 || int(Sis00) == 1028 || int(Sis00) == 2052 || int(Sis00) == (2052+260) || int(Sis00) == (2052+516) || int(Sis00) == (2052+1028) || int(Sis00) == (260+516) || int(Sis00) == (260+1028) || int(Sis00) == (516+1028) )");
  cW4->cd(4);
  Tin->Draw("Qadc7>>hisMPMW4","Tdc011<0.0001 && (int(Sis00) == 260 || int(Sis00) == 516 || int(Sis00) == 1028 || int(Sis00) == 2052 || int(Sis00) == (2052+260) || int(Sis00) == (2052+516) || int(Sis00) == (2052+1028) || int(Sis00) == (260+516) || int(Sis00) == (260+1028) || int(Sis00) == (516+1028) || int(Sis00)==1)");
  

 


  
}
    
