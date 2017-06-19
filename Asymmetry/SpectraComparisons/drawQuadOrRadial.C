
#include <vector>

void drawQuadOrRadial(bool quad=true, TString year="2011-2012", TString anaCh="C") {

  TFile *f1 = new TFile(TString::Format("Octets_%s_DATA_anaCh%s_%s.root",
					year==TString("2011-2012")?"0-59":"60-121",
					anaCh.Data(),quad?"Quadrants":"Radial"),"READ");
  TFile *f2 = new TFile(TString::Format("Octets_%s_SIM_anaCh%s_%s.root",
					year==TString("2011-2012")?"0-59":"60-121",
					anaCh.Data(),quad?"Quadrants":"Radial"),"READ");

  int numAsymms = quad?4:6;

  std::vector <TH1D*> data(numAsymms,NULL);
  std::vector <TH1D*> sim(numAsymms,NULL);
  std::vector <TCanvas*> c(numAsymms,NULL);

  TCanvas *c1 = new TCanvas("c1","c1",1000, 800);
  c1->Divide(numAsymms/2,2);

  for (int i=0; i<numAsymms; ++i) {

    c1->cd(i+1);
    
    data[i] = (TH1D*)f1->Get(TString::Format("FG_%s_%i",quad?"Quadrants":"Radial",i));
    sim[i] = (TH1D*)f2->Get(TString::Format("FG_%s_%i",quad?"Quadrants":"Radial",i));
    sim[i]->SetLineColor(kRed);
    sim[i]->Scale(data[i]->Integral()/sim[i]->Integral());

    data[i]->Draw();
    sim[i]->Draw("SAME");
  }

}

  



