

void foilThicknessTest(TString src) {


  gStyle->SetOptLogy();

  TChain *geom2011 = new TChain("anaTree");
  TChain *geom2012 = new TChain("anaTree");
  TChain *geom2012_iso = new TChain("anaTree");

  for (UInt_t i=0; i<10; ++i ) {

    geom2011->AddFile(TString::Format("%s/%s/analyzed_%i.root",getenv("SIM_2011_2012"),src.Data()));
    geom2012->AddFile(TString::Format("%s/%s/analyzed_%i.root",getenv("SIM_2012_2013"),src.Data()));
    geom2012_iso->AddFile(TString::Format("%s/%s/analyzed_%i.root",getenv("SIM_2012_2013_ISOBUTANE"),src.Data()));
	
  }

  Double_t histMin = 0.;
  Double_t histMax = 5.;

  // histograms for saving energy deposited in foils
  TH1D *edepFoil2011E = new TH1D("edepFoil2011E","Energy Deposited in East Endcap 2011",100, histMin, histMax);
  edepFoil2011E->GetXaxis()->SetTitle("Energy (keV)");
  //edepFoil2011E->SetLineColor(kRed);
  TH1D *edepFoil2011W = new TH1D("edepFoil2011W","Energy Deposited in West Endcap 2011",100, histMin, histMax);
  edepFoil2011W->GetXaxis()->SetTitle("Energy (keV)");

  TH1D *edepFoil2012E = new TH1D("edepFoil2012E","Energy Deposited in East Endcap 2012",100, histMin, histMax);
  edepFoil2012E->GetXaxis()->SetTitle("Energy (keV)");
  //edepFoil2012E->SetLineColor(kRed);
  TH1D *edepFoil2012W = new TH1D("edepFoil2012W","Energy Deposited in West Endcap 2012",100, histMin, histMax);
  edepFoil2012W->GetXaxis()->SetTitle("Energy (keV)");

  TH1D *edepFoil2012E_iso = new TH1D("edepFoil2012E_iso","Energy Deposited in East Endcap 2012_iso",100, histMin, histMax);
  edepFoil2012E_iso->GetXaxis()->SetTitle("Energy (keV)");
  //edepFoil2012E_iso->SetLineColor(kRed);
  TH1D *edepFoil2012W_iso = new TH1D("edepFoil2012W_iso","Energy Deposited in West Endcap 2012_iso",100, histMin, histMax);
  edepFoil2012W_iso->GetXaxis()->SetTitle("Energy (keV)");


  TString cutE = TString::Format("primTheta>TMath::Pi()/2. && EdepE>0. && hitCountSD[5]==1");
  TString cutW = TString::Format("primTheta<TMath::Pi()/2. && EdepW>0. && hitCountSD[15]==1");


  TCanvas *c1 = new TCanvas("c1","c1",1600, 1200);
  c1->Divide(3,2);
  
  c1->cd(1);
  geom2011->Draw("EdepSD[5]>>edepFoil2011E",cutE.Data());
  c1->cd(4);
  geom2011->Draw("EdepSD[15]>>edepFoil2011W",cutW.Data());

  c1->cd(2);
  geom2012->Draw("EdepSD[5]>>edepFoil2012E",cutE.Data());
  c1->cd(5);
  geom2012->Draw("EdepSD[15]>>edepFoil2012W",cutW.Data());

  c1->cd(3);
  geom2012_iso->Draw("EdepSD[5]>>edepFoil2012E_iso",cutE.Data());
  c1->cd(6);
  geom2012_iso->Draw("EdepSD[15]>>edepFoil2012W_iso",cutW.Data());


  std::cout << "2011-2012 East Mean Edep = " << edepFoil2011E->GetMean() << "\n";
  std::cout << "2011-2012 West Mean Edep = " << edepFoil2011W->GetMean() << "\n\n";

  std::cout << "2012-2013 East Mean Edep = " << edepFoil2012E->GetMean() << "\n";
  std::cout << "2012-2013 West Mean Edep = " << edepFoil2012W->GetMean() << "\n\n";

  std::cout << "2012-2013_iso East Mean Edep = " << edepFoil2012E_iso->GetMean() << "\n";
  std::cout << "2012-2013_iso West Mean Edep = " << edepFoil2012W_iso->GetMean() << "\n\n";


};
