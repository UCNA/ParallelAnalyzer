

#include <vector>

void EnergyDependentCorrections() {

  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.08,"t");
  gStyle->SetPadLeftMargin(0.15);
  //gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleYSize(0.05);
  gStyle->SetTitleYOffset(1.4);
  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleSize(0.04,"Z");
  gStyle->SetTitleOffset(1.6,"Z");
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetFillStyle(0);

  // read in the corrections and plot their deltaA/A for each bin
  // This can be done for both theory and MC

  //Int_t octmin = year==2011?0:60;
  //Int_t octmax = year==2011?59:121;

  std::vector <Double_t> energy;
  std::vector <Double_t> UnCorr2011;
  std::vector <Double_t> TheoryCorr2011;
  std::vector <Double_t> AllCorr2011;

  std::vector <Double_t> UnCorr2012;
  std::vector <Double_t> TheoryCorr2012;
  std::vector <Double_t> AllCorr2012;

  Double_t en, val, err;
  
  std::ifstream infile(TString::Format("%s/Asymmetries/UnCorr_OctetAsymmetries_AnaChC_Octets_%i-%i_BinByBin.txt",getenv("ANALYSIS_RESULTS"),0,59));
  
  while (infile >> en >> val >> err ) {
    energy.push_back(en);
    UnCorr2011.push_back(val);
  }

  infile.close();

  infile.open(TString::Format("%s/Asymmetries/DeltaTheoryOnly_OctetAsymmetries_AnaChC_Octets_%i-%i_BinByBin.txt",getenv("ANALYSIS_RESULTS"),0,59));
  
  while (infile >> en >> val >> err) {
    TheoryCorr2011.push_back(val);
  }

  infile.close();

  infile.open(TString::Format("%s/Asymmetries/AllCorr_OctetAsymmetries_AnaChC_Octets_%i-%i_BinByBin.txt",getenv("ANALYSIS_RESULTS"),0,59));
  
  while (infile >> en >> val >> err) {
    AllCorr2011.push_back(val);
  }

  infile.close();
    
  infile.open(TString::Format("%s/Asymmetries/UnCorr_OctetAsymmetries_AnaChC_Octets_%i-%i_BinByBin.txt",getenv("ANALYSIS_RESULTS"),60,121));
  
  while (infile >> en >> val >> err ) {
    energy.push_back(en);
    UnCorr2012.push_back(val);
  }

  infile.close();

  infile.open(TString::Format("%s/Asymmetries/DeltaTheoryOnly_OctetAsymmetries_AnaChC_Octets_%i-%i_BinByBin.txt",getenv("ANALYSIS_RESULTS"),60,121));
  
  while (infile >> en >> val >> err) {
    TheoryCorr2012.push_back(val);
  }

  infile.close();

  infile.open(TString::Format("%s/Asymmetries/AllCorr_OctetAsymmetries_AnaChC_Octets_%i-%i_BinByBin.txt",getenv("ANALYSIS_RESULTS"),60,121));
  
  while (infile >> en >> val >> err) {
    AllCorr2012.push_back(val);
  }

  infile.close();


  
  
  
  

  std::vector <Double_t> en_Th;
  std::vector <Double_t> corr_Th;
  std::vector <Double_t> en_MC2011;
  std::vector <Double_t> corr_MC2011;
  std::vector <Double_t> en_MC2012;
  std::vector <Double_t> corr_MC2012;
  
  
  for (int i=2;i<80;++i) {

    if (UnCorr2011[i]!=0.) en_Th.push_back(energy[i]), corr_Th.push_back(100.*(TheoryCorr2011[i]/UnCorr2011[i]-1.));
    if (TheoryCorr2011[i]!=0.) en_MC2011.push_back(energy[i]), corr_MC2011.push_back(100.*(AllCorr2011[i]/TheoryCorr2011[i]-1.));
    if (TheoryCorr2012[i]!=0.) en_MC2012.push_back(energy[i]), corr_MC2012.push_back(100.*(AllCorr2012[i]/TheoryCorr2012[i]-1.));
  }

  

  TCanvas *c2 = new TCanvas("c2","c2",1200,800);
  c2->Divide(2,1);

  c2->cd(1);

  TGraph *gTheory = new TGraph(en_Th.size(),&en_Th[0],&corr_Th[0]);
  gTheory->SetTitle("Theory Corrections vs. Energy");
  gTheory->SetMarkerStyle(kFullCircle);
  gTheory->SetMinimum(-5.);
  gTheory->SetMaximum(0.);
  gTheory->SetLineWidth(3);
  gTheory->SetLineColor(kBlack);
  gTheory->GetYaxis()->SetTitle("#DeltaA/A (%)");
  gTheory->GetYaxis()->CenterTitle();
  gTheory->GetXaxis()->SetTitle("Energy (keV)");
  gTheory->GetXaxis()->CenterTitle();
  gTheory->Draw("AL");
  
  c2->cd(2);

  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("MC Corrections vs. Energy");
  
  TGraph *gMC2011 = new TGraph(en_MC2011.size(),&en_MC2011[0],&corr_MC2011[0]);
  gMC2011->SetMarkerStyle(kFullCircle);
  //gMC2011->SetMinimum(-10.);
  //gMC2011->SetMaximum(6.);
  gMC2011->SetLineWidth(3);
  gMC2011->SetLineColor(kBlue);
  gMC2011->GetYaxis()->SetTitle("#DeltaA/A (%)");
  gMC2011->GetYaxis()->CenterTitle();
  gMC2011->GetXaxis()->SetTitle("Energy (keV)");
  gMC2011->GetXaxis()->CenterTitle();
  //gMC2011->Draw("AC");

  TGraph *gMC2012 = new TGraph(en_MC2012.size(),&en_MC2012[0],&corr_MC2012[0]);
  gMC2012->SetMarkerStyle(kFullCircle);
  //gMC2012->SetMinimum(-10.);
  //gMC2012->SetMaximum(6.);
  gMC2012->SetLineWidth(3);
  gMC2012->SetLineColor(kGreen);
  //gMC2012->GetYaxis()->SetTitle("#DeltaA/A (%)");
  //gMC2012->GetYaxis()->CenterTitle();
  //gMC2012->GetXaxis()->SetTitle("Energy (keV)");
  //gMC2012->GetXaxis()->CenterTitle();
  //gMC2012->Draw("AC");
  
  mg->Add(gMC2011);
  mg->Add(gMC2012);
  mg->SetMinimum(-8.);
  mg->SetMaximum(6.);
  
  mg->Draw("AC");
  mg->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mg->GetYaxis()->CenterTitle();
  mg->GetXaxis()->SetTitle("Energy (keV)");
  mg->GetXaxis()->CenterTitle();
  gPad->Modified();

  TLegend *leg = new TLegend(0.55,0.75,0.85,0.85);
  leg->AddEntry(gMC2011,"2011-2012","l");
  leg->AddEntry(gMC2012,"2012-2013","l");
  leg->SetTextSize(0.05);
  leg->Draw("SAME");
}
