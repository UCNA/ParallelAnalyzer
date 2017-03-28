#include <iostream>
#include <vector>
#include <cstdlib>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1D.h>
#include "MBUtils.hh"

#include "TriggerMap.hh"

void trigger_threshold(Int_t XeRunPeriod, Double_t binWidth, bool mpmData=false);

int main(int argc,char *argv[]) {
  if (argc<3) {
    std::cout << "Usage: ./triggerThreshold.exe [XeRunPeriod] [BinWidth in mm]\n";
    exit(0);
  }

  trigger_threshold(atoi(argv[1]), atof(argv[2]));

}


void trigger_threshold(Int_t XeRunPeriod, Double_t binWidth, bool mpmData) {
  
  //Create TriggerMap object
  Int_t nParams = 8;
  
  TriggerMap tmap(binWidth, nParams);

  Int_t nBinsXY = tmap.getNbinsXY();
  

  // Read in the Xe runs in this runPeriod
  std::vector <Int_t> XeRuns;
  Char_t temp[200];
  sprintf(temp,"../run_lists/Xenon_Calibration_Run_Period_%i.dat",XeRunPeriod);
  ifstream XeRunList(temp);

  Int_t run;
  UInt_t numRuns=0;
  while (XeRunList >> run) {
    XeRuns.push_back(run);
    numRuns++;
    //if (numRuns>30) break;
  }

  if (numRuns!=XeRuns.size()) {
    std::cout << "Runs!=XeRuns.size()\n";
    exit(0);
  }
  XeRunList.close();
  std::cout << "There are " << numRuns << " in this Xe run period\n";

  TChain *chain;
  if (mpmData) chain = new TChain("phys");
  else chain = new TChain("pass3");

  if (mpmData) {
    for (UInt_t i=0; i<numRuns; i++) {
      sprintf(temp,"%s/hists/spec_%i.root",getenv("UCNAOUTPUTDIR"),XeRuns[i]);
      chain->AddFile(temp);
    }
  }
  else {
    for (UInt_t i=0; i<numRuns; i++) {
      sprintf(temp,"%s/replay_pass3_%i.root",getenv("REPLAY_PASS3"),XeRuns[i]);
      chain->AddFile(temp);
    }
  }
 
  
  //Int_t nbins = 75;
  Double_t lower_limit = -50.;
  Int_t hisBinWidth = 1;
  Double_t East_upper_limit = 150.;
  Double_t West_upper_limit = 150.;
  Int_t nbinsE = (int)(East_upper_limit/hisBinWidth);
  Int_t nbinsW = (int)(West_upper_limit/hisBinWidth);

  std::vector < std::vector <TH1D*> > Etrigg(nBinsXY, std::vector <TH1D*> (nBinsXY, NULL));
  std::vector < std::vector <TH1D*> > EnoTrigg(nBinsXY, std::vector <TH1D*> (nBinsXY, NULL));
  std::vector < std::vector <TH1D*> > Etotal(nBinsXY, std::vector <TH1D*> (nBinsXY, NULL));
  std::vector < std::vector <TH1D*> > EtriggFunc(nBinsXY, std::vector <TH1D*> (nBinsXY, NULL));

  std::vector < std::vector <TH1D*> > Wtrigg(nBinsXY, std::vector <TH1D*> (nBinsXY, NULL));
  std::vector < std::vector <TH1D*> > WnoTrigg(nBinsXY, std::vector <TH1D*> (nBinsXY, NULL));
  std::vector < std::vector <TH1D*> > Wtotal(nBinsXY, std::vector <TH1D*> (nBinsXY, NULL));
  std::vector < std::vector <TH1D*> > WtriggFunc(nBinsXY, std::vector <TH1D*> (nBinsXY, NULL));
  
  for (Int_t xb=0; xb<nBinsXY; xb++) {
    for (Int_t yb=0; yb<nBinsXY; yb++) {

      Double_t xpos = tmap.getBinCenter(xb);
      Double_t ypos = tmap.getBinCenter(yb);

      Etrigg[xb][yb] = new TH1D(TString::Format("Etrigg_%0.1f_%0.1f",xpos, ypos),TString::Format("West Type 1 & East Type 0,1,2/3: EvisE (%0.1f, %0.1f)",xpos,ypos),nbinsE,lower_limit,East_upper_limit);
      EnoTrigg[xb][yb] = new TH1D(TString::Format("EnoTrigg_%0.1f_%0.1f",xpos, ypos),TString::Format("West Type 2/3: EvisE (%0.1f, %0.1f)",xpos,ypos),nbinsE,lower_limit,East_upper_limit);
      Etotal[xb][yb] = new TH1D(TString::Format("Etotal_%0.1f_%0.1f",xpos, ypos),TString::Format("West Type 1,2/3 & East Type 0,1,2/3: EvisE (%0.1f, %0.1f)",xpos,ypos),nbinsE,lower_limit,East_upper_limit);
      EtriggFunc[xb][yb] = new TH1D(TString::Format("EtriggFunc_%0.1f_%0.1f",xpos, ypos),TString::Format("East Trigger Probability (%0.1f, %0.1f)",xpos,ypos),nbinsE,lower_limit,East_upper_limit);
      EtriggFunc[xb][yb]->SetMarkerStyle(20);

      Wtrigg[xb][yb] = new TH1D(TString::Format("Wtrigg_%0.1f_%0.1f",xpos, ypos),TString::Format("East Type 1 & West Type 0,1,2/3: EvisW (%0.1f, %0.1f)",xpos,ypos),nbinsW,lower_limit,West_upper_limit);
      WnoTrigg[xb][yb] = new TH1D(TString::Format("WnoTrigg_%0.1f_%0.1f",xpos, ypos),TString::Format("East Type 2/3: EvisW (%0.1f, %0.1f)",xpos,ypos),nbinsW,lower_limit,West_upper_limit);
      Wtotal[xb][yb] = new TH1D(TString::Format("Wtotal_%0.1f_%0.1f",xpos, ypos),TString::Format("East Type 1,2/3 & West Type 0,1,2/3: EvisW (%0.1f, %0.1f)",xpos,ypos),nbinsW,lower_limit,West_upper_limit);
      WtriggFunc[xb][yb] = new TH1D(TString::Format("WtriggFunc_%0.1f_%0.1f",xpos, ypos),TString::Format("West Trigger Probability (%0.1f, %0.1f)",xpos,ypos),nbinsW,lower_limit,West_upper_limit);
      WtriggFunc[xb][yb]->SetMarkerStyle(20);

    }
  }
      
      
  // Best fit function. This is a shifted erf which goes into a shifted tanh, where the smooth transition is done via application of another shifted tanh
  TF1 *erf = new TF1("erf","([0]+[1]*TMath::Erf((x-[2])/[3]))*(0.5-.5*TMath::TanH((x-[2])/[4]))+(0.5+.5*TMath::TanH((x-[2])/[4]))*([5]+[6]*TMath::TanH((x-[2])/[7]))",-20.,150.);
  
 

  TString pdfFileBase = TString(getenv("TRIGGER_FUNC")) + TString("trigger_functions_XePeriod_")+itos(XeRunPeriod)+"_"+ftos(binWidth)+TString("mm.pdf");
  

  for (Int_t xb=0; xb<nBinsXY; xb++) {
    for (Int_t yb=0; yb<nBinsXY; yb++) {

      Double_t xpos = tmap.getBinCenter(xb);
      Double_t ypos = tmap.getBinCenter(yb);
      
      TString posCutE = TString::Format("xE.center<%f && xE.center>%f && yE.center<%f && yE.center>%f && (xE.center*xE.center+yE.center*yE.center)<2500.",tmap.getBinUpper(xb), tmap.getBinLower(xb),tmap.getBinUpper(yb), tmap.getBinLower(yb));
      TString posCutW = TString::Format("xW.center<%f && xW.center>%f && yW.center<%f && yW.center>%f && (xW.center*xW.center+yW.center*yW.center)<2500.",tmap.getBinUpper(xb), tmap.getBinLower(xb),tmap.getBinUpper(yb), tmap.getBinLower(yb));

      std::vector <Double_t> fitParams(nParams,0.);
  
      TCanvas *c1 = new TCanvas("c1"," ",2000,1600);
      c1->Divide(4,2);


      c1->cd(1);
      chain->Draw(TString::Format("EvisE>>Etrigg_%0.1f_%0.1f",xpos, ypos),TString::Format("(Type==1 && Side==1) || (Side==0) && PID==1 && %s",posCutE.Data()));
      c1->cd(2);  
      chain->Draw(TString::Format("EvisE>>EnoTrigg_%0.1f_%0.1f",xpos, ypos),TString::Format("Type==2 && Side==1 && PID==1 && %s",posCutE.Data()));
      c1->cd(3);
      Etotal[xb][yb]->Add(Etrigg[xb][yb],EnoTrigg[xb][yb]);
      Etotal[xb][yb]->Draw();
      c1->cd(4);
      EtriggFunc[xb][yb]->Divide(Etrigg[xb][yb],Etotal[xb][yb]);
      EtriggFunc[xb][yb]->SetStats(0);
      
      erf->SetParameter(0,0.5); //Constant offset of erf
      erf->FixParameter(0,0.5); 
      erf->SetParameter(1,0.5); //Scaling of erf
      erf->FixParameter(1,0.5); 
      erf->SetParameter(2,25.); //Mean of gaussian integrated for erf
      erf->SetParameter(3,7.); //std. dev. of gaussian integrated for erf
      erf->SetParameter(4,9.); //severity of transition function "turn on"
      erf->SetParLimits(4,1.,25.);
      erf->SetParameter(5,0.5); //constant offset of second tanh
      erf->SetParameter(6,0.5); //Scaling of second tanh
      erf->SetParLimits(4,0.,2.);
      //erf->FixParameter(5,0.5); 
      //erf->FixParameter(6,0.5); 
      erf->SetParameter(7,12.); // stretching factor of second tanh
      erf->SetParLimits(7,0.,50.);


      EtriggFunc[xb][yb]->Fit("erf","","",-20.,East_upper_limit);
      
      EtriggFunc[xb][yb]->Draw("P");
      
      for (Int_t p=0; p<nParams; p++) fitParams[p] = erf->GetParameter(p);
      
      tmap.setTriggerMapPointEast(xb,yb,fitParams);


      c1->cd(5);
      chain->Draw(TString::Format("EvisW>>Wtrigg_%0.1f_%0.1f",xpos, ypos),TString::Format("(Type==1 && Side==0) || (Side==1) && PID==1 && %s",posCutW.Data()));
      c1->cd(6);  
      chain->Draw(TString::Format("EvisW>>WnoTrigg_%0.1f_%0.1f",xpos, ypos),TString::Format("Type==2 && Side==0 && PID==1 && %s",posCutW.Data()));
      c1->cd(7);
      Wtotal[xb][yb]->Add(Wtrigg[xb][yb],WnoTrigg[xb][yb]);
      Wtotal[xb][yb]->Draw();
      c1->cd(8);
      WtriggFunc[xb][yb]->Divide(Wtrigg[xb][yb],Wtotal[xb][yb]);
      WtriggFunc[xb][yb]->SetStats(0);

      erf->SetParameter(2,15.); //Mean of gaussian integrated for erf
      erf->SetParameter(3,5.); //std. dev. of gaussian integrated for erf
      erf->SetParameter(4,7.); //severity of transition function "turn on"
      erf->SetParameter(5,0.5); //constant offset of second tanh
      erf->SetParameter(6,0.5); //Scaling of second tanh
      erf->SetParameter(7,8.); // stretching factor of second tanh
      
      WtriggFunc[xb][yb]->Fit("erf","","",-20.,West_upper_limit);
      
      WtriggFunc[xb][yb]->Draw("P");
      
      for (Int_t p=0; p<nParams; p++) fitParams[p] = erf->GetParameter(p);
      
      tmap.setTriggerMapPointWest(xb,yb,fitParams);
      
      if (xb==0 && yb==0 && xb!=(nBinsXY-1) && yb!=(nBinsXY-1)) c1->Print(pdfFileBase+TString("("));
      else if (xb==(nBinsXY-1) && yb==(nBinsXY-1) && xb!=0 && yb!=0) c1->Print(pdfFileBase+TString(")"));
      else c1->Print(pdfFileBase);
      
      delete c1;
    }
  }
  
   for (Int_t xb=0; xb<nBinsXY; xb++) {
    for (Int_t yb=0; yb<nBinsXY; yb++) {
      delete Etrigg[xb][yb]; delete EnoTrigg[xb][yb]; delete Etotal[xb][yb]; delete EtriggFunc[xb][yb]; 
      delete Wtrigg[xb][yb]; delete WnoTrigg[xb][yb]; delete Wtotal[xb][yb]; delete WtriggFunc[xb][yb];
    }
   }  

  tmap.writeTriggerMap(XeRunPeriod);     
  
}
