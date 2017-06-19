#include <fstream>

void betaAsymmetryFitter(TString year="2011-2012") {

  gStyle->SetTitleSize(0.06,"t");
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetTitleYSize(0.05);
  gStyle->SetTitleYOffset(1.0);
  gStyle->SetTitleXSize(0.06);
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetLabelSize(0.035,"xyz");
  gStyle->SetOptFit(0);//(1111);
  gStyle->SetStatY(0.85);
  gStyle->SetStatX(0.975);
  gStyle->SetStatW(.09);

  gStyle->SetFillStyle(0000); 
  //gStyle->SetStatStyle(0); 
  //gStyle->SetTitleStyle(0); 
  //gStyle->SetCanvasBorderSize(0); 
  //gStyle->SetFrameBorderSize(0); 
  gStyle->SetLegendBorderSize(0); 
  //gStyle->SetStatBorderSize(0); 
  //gStyle->SetTitleBorderSize(0);
  
  TString octets = year==TString("2011-2012") ? "0-59" : "60-121";
  
  std::vector <Double_t> energy;
  std::vector <Double_t> AsymWithEnDep;
  std::vector <Double_t> AsymWithEnDepErr;
  std::vector <Double_t> Asym;
  std::vector <Double_t> AsymErr;

  Int_t en = 0;
  Double_t a = 0.;
  Double_t aErr = 0.;

  // Load the beta/2 corrected data first
  TString fname = TString::Format("%s/Asymmetries/UnCorr_OctetAsymmetries_AnaChC_Octets_%s_BinByBin.txt",getenv("SIM_ANALYSIS_RESULTS"),octets.Data());
  //UNBLINDED_

  std::ifstream infile(fname.Data());
  
  while ( infile >> en >> a >> aErr ) {
    //cout << en << "\t" << a << "\t" << aErr << "\t\n";
    energy.push_back(en);
    Asym.push_back(a);
    AsymErr.push_back(aErr);
  }
  
  infile.close();

  
  
  fname = TString::Format("%s/Asymmetries/UnCorr_OctetAsymmetries_AnaChC_Octets_%s_BinByBin_withEnergyDependence.txt",getenv("SIM_ANALYSIS_RESULTS"),octets.Data());

  infile.open(fname.Data());
  
  while ( infile >> en >> a >> aErr ) {
    AsymWithEnDep.push_back(a);
    AsymWithEnDepErr.push_back(aErr);
  }
  
  infile.close();
  
  TCanvas *c1 = new TCanvas("c1");

  TGraphErrors *gr1 = new TGraphErrors(80,&energy[4],&AsymWithEnDep[4],NULL,&AsymWithEnDepErr[4]);
  gr1->SetTitle(TString::Format("A(E) vs. Energy for %s",year.Data()).Data());
  gr1->GetXaxis()->SetTitle("Energy (keV)");
  gr1->GetYaxis()->SetTitle("A(E)");
  gr1->GetXaxis()->CenterTitle();
  gr1->GetYaxis()->CenterTitle();
  gr1->SetMarkerStyle(21);
  gr1->SetMarkerSize(0.6);
  gr1->SetMinimum(-0.1);
  gr1->SetMaximum(0.);
  
  Double_t A0 = -0.1184;
  TF1 *f1 = new TF1("f1","(-0.1184/2.*sqrt(x*x+2.*x*510.99892)/(x+510.99892))*[0]*(1.+[1]*x)*(1+[2]/(1.+exp((x-[3])/[4])))",50.,780.);
  f1->SetParameters(1.,0.,2.,50.,10.);
  //f1->SetParLimits(4,1.,20.);

  gr1->Fit("f1","LR");

  gr1->Draw("AP");
  
  TF1 *f2 = new TF1("f2","(-0.1184/2.*sqrt(x*x+2.*x*510.99892)/(x+510.99892))",0., 1000.);
  f2->SetLineColor(kBlack);
  f2->SetLineStyle(2);

  f2->Draw("SAME");

  TCanvas *c2 = new TCanvas("c2");

  TGraphErrors *gr2 = new TGraphErrors(80,&energy[4],&Asym[4],NULL,&AsymErr[4]);
  gr2->SetTitle(TString::Format("A_{0} vs. Energy for %s",year.Data()).Data());
  gr2->GetXaxis()->SetTitle("Energy (keV)");
  gr2->GetYaxis()->SetTitle("A_{0}");
  gr2->GetXaxis()->CenterTitle();
  gr2->GetYaxis()->CenterTitle();
  gr2->SetMarkerStyle(21);
  gr2->SetMarkerSize(0.6);
  gr2->SetMinimum(-0.14);
  gr2->SetMaximum(-0.10);
  
  TF1 *f3 = new TF1("f3","[0]",50.,750.);
  f3->SetParameter(0,-0.1184);

  gr2->Fit("f3","LR");

  gr2->Draw("AP");
  
  

}

  

 
