

void AsymmByAnach(TString year, TString corrections, bool withPOL, Int_t ebinLow=220, Int_t ebinHigh=680) {

  bool readInAsymms = true;
  
  const Int_t numAnaChT0 = 4;
  const Int_t numAnaChBacksc = 5;
  TString anChT0[numAnaChT0] = {"A","B","C","D"};//
  TString anChBacksc[numAnaChBacksc] = {"F","G","H","J","K"};
  
  std::vector <TString> anaChoicesT0;
  std::vector <Double_t> indicesT0;
  std::vector <Double_t> AsymmDataT0;
  std::vector <Double_t> AsymmDataErrT0; //Statistical Errors
  std::vector <Double_t> AsymmSimT0;
  std::vector <Double_t> AsymmSimErrT0; //Statistical Errors

  std::vector <TString> anaChoicesBacksc;
  std::vector <Double_t> indicesBacksc;
  std::vector <Double_t> AsymmDataBacksc;
  std::vector <Double_t> AsymmDataErrBacksc; //Statistical Errors
  std::vector <Double_t> AsymmSimBacksc;
  std::vector <Double_t> AsymmSimErrBacksc; //Statistical Errors

  //Fill anaChoices
  for (Int_t i=0; i< sizeof(anChT0)/sizeof(TString); i++) 
    { anaChoicesT0.push_back(anChT0[i]); indicesT0.push_back(i+1); }

  for (Int_t i=0; i< sizeof(anChBacksc)/sizeof(TString); i++) 
    { anaChoicesBacksc.push_back(anChBacksc[i]); indicesBacksc.push_back(i+1); }
  
  
  //Read in Asymmetries
  ifstream infile;
  std::string strHold;
  
  if ( readInAsymms ) {
    
    for (Int_t i=0; i<anaChoicesT0.size(); i++) {
      infile.open(TString::Format("%s/Asymmetries/%s%s_OctetAsymmetries_AnaCh%s_%i-%i_Octets_%s.txt",
				  getenv("ANALYSIS_RESULTS"),corrections.Data(),
				  (withPOL?"_withPOL":""),anaChoicesT0[i].Data(),
				  ebinLow,ebinHigh,year==TString("2011-2012")?"0-59":"60-121").Data());
      
      Double_t AsymHold = 0.;
      Double_t AsymErrHold = 0.;
      if (infile.is_open()) {
	for (Int_t j=0; j<3; j++) infile >> strHold >> AsymHold >> AsymErrHold;
      }
      AsymmDataT0.push_back(AsymHold);
      AsymmDataErrT0.push_back(AsymErrHold);
      
      infile.close();
      
      infile.open(TString::Format("%s/Asymmetries/%s_OctetAsymmetries_AnaCh%s_%i-%i_Octets_%s.txt",
				  getenv("SIM_ANALYSIS_RESULTS"),corrections.Data(),
				  anaChoicesT0[i].Data(),ebinLow,ebinHigh,
				  year==TString("2011-2012")?"0-59":"60-121").Data() );
      
      AsymHold = 0.;
      AsymErrHold = 0.;
      if (infile.is_open()) {
	for (Int_t j=0; j<3; j++) infile >> strHold >> AsymHold >> AsymErrHold;
      }
      AsymmSimT0.push_back(AsymHold);
      AsymmSimErrT0.push_back(AsymErrHold);
      
      infile.close();
    }

    for (Int_t i=0; i<anaChoicesBacksc.size(); i++) {
      infile.open(TString::Format("%s/Asymmetries/%s%s_OctetAsymmetries_AnaCh%s_%i-%i_Octets_%s.txt",
				  getenv("ANALYSIS_RESULTS"),corrections.Data(),
				  (withPOL?"_withPOL":""),anaChoicesBacksc[i].Data(),
				  ebinLow,ebinHigh,year==TString("2011-2012")?"0-59":"60-121").Data());
      
      Double_t AsymHold = 0.;
      Double_t AsymErrHold = 0.;
      if (infile.is_open()) {
	for (Int_t j=0; j<3; j++) infile >> strHold >> AsymHold >> AsymErrHold;
      }
      AsymmDataBacksc.push_back(AsymHold);
      AsymmDataErrBacksc.push_back(AsymErrHold);
      
      infile.close();
      
      infile.open(TString::Format("%s/Asymmetries/%s_OctetAsymmetries_AnaCh%s_%i-%i_Octets_%s.txt",
				  getenv("SIM_ANALYSIS_RESULTS"),corrections.Data(),
				  anaChoicesBacksc[i].Data(),ebinLow,ebinHigh,
				  year==TString("2011-2012")?"0-59":"60-121").Data() );
      
      AsymHold = 0.;
      AsymErrHold = 0.;
      if (infile.is_open()) {
	for (Int_t j=0; j<3; j++) infile >> strHold >> AsymHold >> AsymErrHold;
      }
      AsymmSimBacksc.push_back(AsymHold);
      AsymmSimErrBacksc.push_back(AsymErrHold);
      
      infile.close();
    }
  }    
  else {
    
    AsymmDataT0.push_back(-0.12404); AsymmDataErrT0.push_back(0.00079);
    AsymmDataT0.push_back(-0.12755); AsymmDataErrT0.push_back(0.00079);
    AsymmDataT0.push_back(-0.0925); AsymmDataErrT0.push_back(0.0056);
    AsymmDataT0.push_back(-0.0083); AsymmDataErrT0.push_back(0.0079);
    
    AsymmSimT0.push_back(-0.12291); AsymmSimErrT0.push_back(0.00069);
    AsymmSimT0.push_back(-0.12551); AsymmSimErrT0.push_back(0.00069);
    AsymmSimT0.push_back(-0.0998); AsymmSimErrT0.push_back(0.0046);
    AsymmSimT0.push_back(-0.0251); AsymmSimErrT0.push_back(0.0074);
  }
  

  

  Double_t xerr[10] = {0.};

  //cout << anaChoices.size() << " " << anaChoices[0] << " " << AsymmData[0] << " " << xerr[0] << " " << AsymmDataErr[0] << endl;
    

  gStyle->SetLegendBorderSize(0);
  gStyle->SetFillStyle(0);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetTitleXOffset(1.2);

  TCanvas *c1 = new TCanvas("c1","demo bin labels",10,10,600,600);

  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);

  
  TH1F * data = new TH1F("data","Analysis Choices with Type 0",AsymmDataT0.size(),0.5,AsymmDataT0.size()+0.5);
  data->SetMarkerColor(kBlue);
  data->SetLineColor(kBlue);
  data->SetLineWidth(3);
  data->SetMarkerStyle(kFullCircle);
  //data->GetXaxis()->SetNdivisions();
  //data->SetCanExtend(TH1::kAllAxes);
  data->SetStats(0);

  for (Int_t i=0;i<AsymmDataT0.size();++i) {
    data->SetBinContent(i+1,AsymmDataT0[i]);
    data->SetBinError(i+1,AsymmDataErrT0[i]);
  }
  //data->LabelsDeflate("X");
  //data->LabelsDeflate("Y");
  //data->LabelsOption("v");
  data->GetXaxis()->SetBinLabel(1,"A");
  data->GetXaxis()->SetBinLabel(2,"B");
  data->GetXaxis()->SetBinLabel(3,"C");
  data->GetXaxis()->SetBinLabel(4,"D");
  data->GetXaxis()->SetLabelSize(0.07);
  data->GetXaxis()->SetTitle("Analysis Choice");
  data->GetYaxis()->SetTitle("Asymmetry");
  data->GetYaxis()->SetTitleOffset(2.);
  data->GetXaxis()->SetTitleSize(0.05);
  data->GetXaxis()->CenterTitle();
  data->GetYaxis()->CenterTitle();
  data->SetMinimum(-0.13);
  data->SetMaximum(-0.118);

  data->Draw("EX0");

  TH1F * sim = new TH1F("","Sim vs MC Blinded Asymmetries",AsymmSimT0.size(),0.5,AsymmSimT0.size()+0.5);
  sim->SetMarkerColor(kRed);
  sim->SetLineColor(kRed);
  sim->SetLineWidth(3);
  sim->SetMarkerStyle(kFullSquare);
  //sim->GetXaxis()->SetNdivisions();
  //sim->SetCanExtend(TH1::kAllAxes);
  sim->SetStats(0);

  for (Int_t i=0;i<AsymmSimT0.size();++i) {
    sim->SetBinContent(i+1,AsymmSimT0[i]);
    sim->SetBinError(i+1,AsymmSimErrT0[i]);
  }
  
  sim->Draw("SAME EX0");
  
  
  TLegend *leg = new TLegend(0.55,0.75,0.9,0.9);
  //leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  leg->AddEntry(data,"Data");
  leg->AddEntry(sim,"Monte Carlo");
   //leg->AddEntry("gr","Graph with error bars","lep");
  leg->Draw();


  ///////////////// backscatter anachoices /////////////////////////

  TCanvas *c2 = new TCanvas("c2","demo bin labels",10,10,600,600);

  c2->SetLeftMargin(0.15);
  c2->SetBottomMargin(0.15);

  
  TH1F * dataBacksc = new TH1F("dataBacksc","Backscattering Analysis Choices",AsymmDataBacksc.size(),0.5,AsymmDataBacksc.size()+0.5);
  dataBacksc->SetMarkerColor(kBlue);
  dataBacksc->SetLineColor(kBlue);
  dataBacksc->SetLineWidth(3);
  dataBacksc->SetMarkerStyle(kFullCircle);
  //dataBacksc->GetXaxis()->SetNdivisions();
  //dataBacksc->SetCanExtend(TH1::kAllAxes);
  dataBacksc->SetStats(0);

  for (Int_t i=0;i<AsymmDataBacksc.size();++i) {
    dataBacksc->SetBinContent(i+1,AsymmDataBacksc[i]);
    dataBacksc->SetBinError(i+1,AsymmDataErrBacksc[i]);
  }
  //dataBacksc->LabelsDeflate("X");
  //dataBacksc->LabelsDeflate("Y");
  //dataBacksc->LabelsOption("v");
  dataBacksc->GetXaxis()->SetBinLabel(1,"F");
  dataBacksc->GetXaxis()->SetBinLabel(2,"G");
  dataBacksc->GetXaxis()->SetBinLabel(3,"H");
  dataBacksc->GetXaxis()->SetBinLabel(4,"J");
  dataBacksc->GetXaxis()->SetBinLabel(5,"K");
  dataBacksc->GetXaxis()->SetLabelSize(0.07);
  dataBacksc->GetXaxis()->SetTitle("Analysis Choice");
  dataBacksc->GetYaxis()->SetTitle("Asymmetry");
  dataBacksc->GetYaxis()->SetTitleOffset(2.);
  dataBacksc->GetXaxis()->SetTitleSize(0.05);
  dataBacksc->GetXaxis()->CenterTitle();
  dataBacksc->GetYaxis()->CenterTitle();
  dataBacksc->SetMinimum(-0.13);
  dataBacksc->SetMaximum(0.02);

  dataBacksc->Draw("EX0");

  TH1F * simBacksc = new TH1F("","Sim vs MC Blinded Asymmetries",AsymmSimBacksc.size(),0.5,AsymmSimBacksc.size()+0.5);
  simBacksc->SetMarkerColor(kRed);
  simBacksc->SetLineColor(kRed);
  simBacksc->SetLineWidth(3);
  simBacksc->SetMarkerStyle(kFullSquare);
  //simBacksc->GetXaxis()->SetNdivisions();
  //simBacksc->SetCanExtend(TH1::kAllAxes);
  simBacksc->SetStats(0);

  for (Int_t i=0;i<AsymmSimBacksc.size();++i) {
    simBacksc->SetBinContent(i+1,AsymmSimBacksc[i]);
    simBacksc->SetBinError(i+1,AsymmSimErrBacksc[i]);
  }
  
  simBacksc->Draw("SAME EX0");
 
  TLegend *leg = new TLegend(0.65,0.65,0.95,0.8);
  //leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  leg->AddEntry(dataBacksc,"Data");
  leg->AddEntry(simBacksc,"Monte Carlo");
   //leg->AddEntry("gr","Graph with error bars","lep");
  //leg->Draw();

  
  
}
