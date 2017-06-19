

void type4Fraction(TString geom) {


  TChain *chain = new TChain("anaTree");

  for (int i=0;i<200; ++i) {
    std::cout << "File " << i << std::endl;
    chain->AddFile(TString::Format("%s/Beta_polE/analyzed_%i.root",
				   geom==TString("2011-2012")?getenv("SIM_2011_2012"):
				   (geom==TString("2012-2013")?getenv("SIM_2012_2013"):
				    getenv("SIM_2012_2013_ISOBUTANE")),i));
    chain->AddFile(TString::Format("%s/Beta_polW/analyzed_%i.root",
				   geom==TString("2011-2012")?getenv("SIM_2011_2012"):
				   (geom==TString("2012-2013")?getenv("SIM_2012_2013"):
				    getenv("SIM_2012_2013_ISOBUTANE")),i));
  }

  TFile *file = new TFile(TString::Format("lost_events_%s.root",geom.Data()),"RECREATE");
  
  TH1D *lost = new TH1D("lost","Lost Events",80,0.,800.);
  TH1D *type0 = new TH1D("type0","Type 0 Events",80,0.,800.);
  TH1D *percLost = new TH1D("percLost","Lost/Type0",80,0.,800.);

  chain->Draw("primKE>>lost","EdepQE<0.001 && EdepQW<0.001 && primKE<800.","goff");
  chain->Draw("primKE>>type0","((EdepQE>0.001 && EdepQW==0. && MWPCEnergyW==0.) || (EdepQW>0.001 && EdepQE==0. && MWPCEnergyE==0.)) && primKE<800.","goff");

  percLost->Divide(lost,type0);

  //  percLost->SetMarkerStyle(kFullCircle);
  //percLost->Draw("P");

  file->Write();
  delete file;
}
