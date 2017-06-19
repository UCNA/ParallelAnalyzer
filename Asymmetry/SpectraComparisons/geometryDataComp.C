

void geometryDataComp(TString simORdata) {
  gStyle->SetOptStat(0);

  Int_t col2011 = 1;
  Int_t col2012 = 4;

  TString normType = "ALL";

  Double_t normLow = 200.;
  Double_t normHigh = 780.;
  
  Double_t xAxisMax = 1200.;

  //Storing event fractions for data, E/W and sfON/OFF
  Double_t t0E_sfON, t0E_sfOFF, t1E_sfON, t1E_sfOFF, t2E_sfON, t2E_sfOFF, t3E_sfON, t3E_sfOFF;
  Double_t t0W_sfON, t0W_sfOFF, t1W_sfON, t1W_sfOFF, t2W_sfON, t2W_sfOFF, t3W_sfON, t3W_sfOFF;
  t0E_sfON = t0E_sfOFF = t1E_sfON = t1E_sfOFF = t23E_sfON = t2E_sfOFF = t3E_sfOFF = 0;
  t0W_sfON = t0W_sfOFF = t1W_sfON = t1W_sfOFF = t23W_sfON = t2W_sfOFF = t3E_sfOFF = 0.;
  
  TFile *geom2011_file = new TFile(TString::Format("Octets_0-59_%s.root",simORdata.Data()),"READ");
  TFile *geom2012_file = new TFile(TString::Format("Octets_60-121_%s.root",simORdata.Data()),"READ");


  //Load the histograms

  TH1D *geom2011ALL = (TH1D*)geom2011_file->Get("EreconALL");
  TH1D *geom20110 = (TH1D*)geom2011_file->Get("Erecon0");
  TH1D *geom20111 = (TH1D*)geom2011_file->Get("Erecon1");
  TH1D *geom20112 = (TH1D*)geom2011_file->Get("Erecon2");
  TH1D *geom20113 = (TH1D*)geom2011_file->Get("Erecon3");
  
  TH1D *geom2012ALL = (TH1D*)geom2012_file->Get("EreconALL");
  TH1D *geom20120 = (TH1D*)geom2012_file->Get("Erecon0");
  TH1D *geom20121 = (TH1D*)geom2012_file->Get("Erecon1");
  TH1D *geom20122 = (TH1D*)geom2012_file->Get("Erecon2");
  TH1D *geom20123 = (TH1D*)geom2012_file->Get("Erecon3");
  
  TH1D *geom2011ALL_sfOFF_E = (TH1D*)geom2011_file->Get("EreconALL_sfOFF_E");
  TH1D *geom20110_sfOFF_E = (TH1D*)geom2011_file->Get("Erecon0_sfOFF_E");
  TH1D *geom20111_sfOFF_E = (TH1D*)geom2011_file->Get("Erecon1_sfOFF_E");
  TH1D *geom20112_sfOFF_E = (TH1D*)geom2011_file->Get("Erecon2_sfOFF_E");
  TH1D *geom20113_sfOFF_E = (TH1D*)geom2011_file->Get("Erecon3_sfOFF_E");
  
  TH1D *geom2012ALL_sfOFF_E = (TH1D*)geom2012_file->Get("EreconALL_sfOFF_E");
  TH1D *geom20120_sfOFF_E = (TH1D*)geom2012_file->Get("Erecon0_sfOFF_E");
  TH1D *geom20121_sfOFF_E = (TH1D*)geom2012_file->Get("Erecon1_sfOFF_E");
  TH1D *geom20122_sfOFF_E = (TH1D*)geom2012_file->Get("Erecon2_sfOFF_E");
  TH1D *geom20123_sfOFF_E = (TH1D*)geom2012_file->Get("Erecon3_sfOFF_E");

  TH1D *geom2011ALL_sfOFF_W = (TH1D*)geom2011_file->Get("EreconALL_sfOFF_W");
  TH1D *geom20110_sfOFF_W = (TH1D*)geom2011_file->Get("Erecon0_sfOFF_W");
  TH1D *geom20111_sfOFF_W = (TH1D*)geom2011_file->Get("Erecon1_sfOFF_W");
  TH1D *geom20112_sfOFF_W = (TH1D*)geom2011_file->Get("Erecon2_sfOFF_W");
  TH1D *geom20113_sfOFF_W = (TH1D*)geom2011_file->Get("Erecon3_sfOFF_W");
  
  TH1D *geom2012ALL_sfOFF_W = (TH1D*)geom2012_file->Get("EreconALL_sfOFF_W");
  TH1D *geom20120_sfOFF_W = (TH1D*)geom2012_file->Get("Erecon0_sfOFF_W");
  TH1D *geom20121_sfOFF_W = (TH1D*)geom2012_file->Get("Erecon1_sfOFF_W");
  TH1D *geom20122_sfOFF_W = (TH1D*)geom2012_file->Get("Erecon2_sfOFF_W");
  TH1D *geom20123_sfOFF_W = (TH1D*)geom2012_file->Get("Erecon3_sfOFF_W");
  
  TH1D *geom2011ALL_sfON_E = (TH1D*)geom2011_file->Get("EreconALL_sfON_E");
  TH1D *geom20110_sfON_E = (TH1D*)geom2011_file->Get("Erecon0_sfON_E");
  TH1D *geom20111_sfON_E = (TH1D*)geom2011_file->Get("Erecon1_sfON_E");
  TH1D *geom20112_sfON_E = (TH1D*)geom2011_file->Get("Erecon2_sfON_E");
  TH1D *geom20113_sfON_E = (TH1D*)geom2011_file->Get("Erecon3_sfON_E");
  
  TH1D *geom2012ALL_sfON_E = (TH1D*)geom2012_file->Get("EreconALL_sfON_E");
  TH1D *geom20120_sfON_E = (TH1D*)geom2012_file->Get("Erecon0_sfON_E");
  TH1D *geom20121_sfON_E = (TH1D*)geom2012_file->Get("Erecon1_sfON_E");
  TH1D *geom20122_sfON_E = (TH1D*)geom2012_file->Get("Erecon2_sfON_E");
  TH1D *geom20123_sfON_E = (TH1D*)geom2012_file->Get("Erecon3_sfON_E");

  TH1D *geom2011ALL_sfON_W = (TH1D*)geom2011_file->Get("EreconALL_sfON_W");
  TH1D *geom20110_sfON_W = (TH1D*)geom2011_file->Get("Erecon0_sfON_W");
  TH1D *geom20111_sfON_W = (TH1D*)geom2011_file->Get("Erecon1_sfON_W");
  TH1D *geom20112_sfON_W = (TH1D*)geom2011_file->Get("Erecon2_sfON_W");
  TH1D *geom20113_sfON_W = (TH1D*)geom2011_file->Get("Erecon3_sfON_W");
  
  TH1D *geom2012ALL_sfON_W = (TH1D*)geom2012_file->Get("EreconALL_sfON_W");
  TH1D *geom20120_sfON_W = (TH1D*)geom2012_file->Get("Erecon0_sfON_W");
  TH1D *geom20121_sfON_W = (TH1D*)geom2012_file->Get("Erecon1_sfON_W");
  TH1D *geom20122_sfON_W = (TH1D*)geom2012_file->Get("Erecon2_sfON_W");
  TH1D *geom20123_sfON_W = (TH1D*)geom2012_file->Get("Erecon3_sfON_W");
  

  t0E_sfOFF = geom20110_sfOFF_E->Integral(18,79);
  t0E_sfON = geom20110_sfON_E->Integral(18,79);
  t1E_sfOFF = geom20111_sfOFF_E->Integral(18,79);
  t1E_sfON = geom20111_sfON_E->Integral(18,79);
  t2E_sfOFF = geom20112_sfOFF_E->Integral(18,79);
  t2E_sfON = geom20112_sfON_E->Integral(18,79);
  t3E_sfOFF = geom20113_sfOFF_E->Integral(18,79);
  t3E_sfON = geom20113_sfON_E->Integral(18,79);
  t0W_sfOFF = geom20110_sfOFF_W->Integral(18,79);
  t0W_sfON = geom20110_sfON_W->Integral(18,79);
  t1W_sfOFF = geom20111_sfOFF_W->Integral(18,79);
  t1W_sfON = geom20111_sfON_W->Integral(18,79);
  t2W_sfOFF = geom20112_sfOFF_W->Integral(18,79);
  t2W_sfON = geom20112_sfON_W->Integral(18,79);
  t3W_sfOFF = geom20113_sfOFF_W->Integral(18,79);
  t3W_sfON = geom20113_sfON_W->Integral(18,79);
  

  double totalE_sfOFF = t0E_sfOFF + t1E_sfOFF + t2E_sfOFF + t3E_sfOFF;
  double totalE_sfON = t0E_sfON + t1E_sfON + t2E_sfON + t3E_sfON;
  double totalW_sfOFF = t0W_sfOFF + t1W_sfOFF + t2W_sfOFF + t3W_sfOFF;
  double totalW_sfON = t0W_sfON + t1W_sfON + t2W_sfON + t3W_sfON;

  
  std::cout << "********************************************\n"
	    << "          Event ratios (180-780 keV)\n\n"
	    << "\t\tsfON\t\t\tsfOFF\n"
	    << "Type 0E:\t" << t0E_sfON/totalE_sfON << "\t\t" << t0E_sfOFF/totalE_sfOFF << "\n"
	    << "Type 0W:\t" << t0W_sfON/totalW_sfON << "\t\t" << t0W_sfOFF/totalW_sfOFF << "\n\n"
	    << "Type 1E:\t" << t1E_sfON/totalE_sfON << "\t\t" << t1E_sfOFF/totalE_sfOFF << "\n"
	    << "Type 1W:\t" << t1W_sfON/totalW_sfON << "\t\t" << t1W_sfOFF/totalW_sfOFF << "\n\n"
	    << "Type 2E:\t" << t2E_sfON/totalE_sfON << "\t\t" << t2E_sfOFF/totalE_sfOFF << "\n"
	    << "Type 2W:\t" << t2W_sfON/totalW_sfON << "\t\t" << t2W_sfOFF/totalW_sfOFF << "\n\n"
	    << "Type 3E:\t" << t3E_sfON/totalE_sfON << "\t\t" << t3E_sfOFF/totalE_sfOFF << "\n"
	    << "Type 3W:\t" << t3W_sfON/totalW_sfON << "\t\t" << t3W_sfOFF/totalW_sfOFF << "\n\n"
	    << "********************************************\n";
    
 
    
  //Normalize
  Double_t normFactor = 1.;
  Double_t normFactor_sfOFF_E = 1.;
  Double_t normFactor_sfON_E = 1.;
  Double_t normFactor_sfOFF_W = 1.;
  Double_t normFactor_sfON_W = 1.;
  
  if (normType==TString("ALL")) {
    normFactor = 1./(geom2011ALL->Integral(geom2011ALL->GetXaxis()->FindFixBin(normLow),geom2011ALL->GetXaxis()->FindFixBin(normHigh))/geom2012ALL->Integral(geom2012ALL->GetXaxis()->FindFixBin(normLow),geom2012ALL->GetXaxis()->FindFixBin(normHigh)));
    normFactor_sfOFF_E = 1./(geom2011ALL_sfOFF_E->Integral(geom2011ALL_sfOFF_E->GetXaxis()->FindFixBin(normLow),geom2011ALL_sfOFF_E->GetXaxis()->FindFixBin(normHigh))/geom2012ALL_sfOFF_E->Integral(geom2012ALL_sfOFF_E->GetXaxis()->FindFixBin(normLow),geom2012ALL_sfOFF_E->GetXaxis()->FindFixBin(normHigh)));
    normFactor_sfON_E = 1./(geom2011ALL_sfON_E->Integral(geom2011ALL_sfON_E->GetXaxis()->FindFixBin(normLow),geom2011ALL_sfON_E->GetXaxis()->FindFixBin(normHigh))/geom2012ALL_sfON_E->Integral(geom2012ALL_sfON_E->GetXaxis()->FindFixBin(normLow),geom2012ALL_sfON_E->GetXaxis()->FindFixBin(normHigh)));
    normFactor_sfOFF_W = 1./(geom2011ALL_sfOFF_W->Integral(geom2011ALL_sfOFF_W->GetXaxis()->FindFixBin(normLow),geom2011ALL_sfOFF_W->GetXaxis()->FindFixBin(normHigh))/geom2012ALL_sfOFF_W->Integral(geom2012ALL_sfOFF_W->GetXaxis()->FindFixBin(normLow),geom2012ALL_sfOFF_W->GetXaxis()->FindFixBin(normHigh)));
    normFactor_sfON_W = 1./(geom2011ALL_sfON_W->Integral(geom2011ALL_sfON_W->GetXaxis()->FindFixBin(normLow),geom2011ALL_sfON_W->GetXaxis()->FindFixBin(normHigh))/geom2012ALL_sfON_W->Integral(geom2012ALL_sfON_W->GetXaxis()->FindFixBin(normLow),geom2012ALL_sfON_W->GetXaxis()->FindFixBin(normHigh)));
    
  }

  if (normType==TString("0")) {
    normFactor = geom20110->Integral(geom20110->GetXaxis()->FindFixBin(normLow),geom20110->GetXaxis()->FindFixBin(normHigh))/geom20120->Integral(geom20120->GetXaxis()->FindFixBin(normLow),geom20120->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfOFF_E = geom20110_sfOFF_E->Integral(geom20110_sfOFF_E->GetXaxis()->FindFixBin(normLow),geom20110_sfOFF_E->GetXaxis()->FindFixBin(normHigh))/geom20120_sfOFF_E->Integral(geom20120_sfOFF_E->GetXaxis()->FindFixBin(normLow),geom20120_sfOFF_E->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfON_E = geom20110_sfON_E->Integral(geom20110_sfON_E->GetXaxis()->FindFixBin(normLow),geom20110_sfON_E->GetXaxis()->FindFixBin(normHigh))/geom20120_sfON_E->Integral(geom20120_sfON_E->GetXaxis()->FindFixBin(normLow),geom20120_sfON_E->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfOFF_W = geom20110_sfOFF_W->Integral(geom20110_sfOFF_W->GetXaxis()->FindFixBin(normLow),geom20110_sfOFF_W->GetXaxis()->FindFixBin(normHigh))/geom20120_sfOFF_W->Integral(geom20120_sfOFF_W->GetXaxis()->FindFixBin(normLow),geom20120_sfOFF_W->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfON_W = geom20110_sfON_W->Integral(geom20110_sfON_W->GetXaxis()->FindFixBin(normLow),geom20110_sfON_W->GetXaxis()->FindFixBin(normHigh))/geom20120_sfON_W->Integral(geom20120_sfON_W->GetXaxis()->FindFixBin(normLow),geom20120_sfON_W->GetXaxis()->FindFixBin(normHigh));
  }
  
  if (normType==TString("1")) {
    normFactor = geom20111->Integral(geom20111->GetXaxis()->FindFixBin(normLow),geom20111->GetXaxis()->FindFixBin(normHigh))/geom20121->Integral(geom20121->GetXaxis()->FindFixBin(normLow),geom20121->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfOFF_E = geom20111_sfOFF_E->Integral(geom20111_sfOFF_E->GetXaxis()->FindFixBin(normLow),geom20111_sfOFF_E->GetXaxis()->FindFixBin(normHigh))/geom20121_sfOFF_E->Integral(geom20121_sfOFF_E->GetXaxis()->FindFixBin(normLow),geom20121_sfOFF_E->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfON_E = geom20111_sfON_E->Integral(geom20111_sfON_E->GetXaxis()->FindFixBin(normLow),geom20111_sfON_E->GetXaxis()->FindFixBin(normHigh))/geom20121_sfON_E->Integral(geom20121_sfON_E->GetXaxis()->FindFixBin(normLow),geom20121_sfON_E->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfOFF_W = geom20111_sfOFF_W->Integral(geom20111_sfOFF_W->GetXaxis()->FindFixBin(normLow),geom20111_sfOFF_W->GetXaxis()->FindFixBin(normHigh))/geom20121_sfOFF_W->Integral(geom20121_sfOFF_W->GetXaxis()->FindFixBin(normLow),geom20121_sfOFF_W->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfON_W = geom20111_sfON_W->Integral(geom20111_sfON_W->GetXaxis()->FindFixBin(normLow),geom20111_sfON_W->GetXaxis()->FindFixBin(normHigh))/geom20121_sfON_W->Integral(geom20121_sfON_W->GetXaxis()->FindFixBin(normLow),geom20121_sfON_W->GetXaxis()->FindFixBin(normHigh));
  }
  

  if (normType==TString("2")) {
    normFactor = geom20112->Integral(geom20112->GetXaxis()->FindFixBin(normLow),geom20112->GetXaxis()->FindFixBin(normHigh))/geom20122->Integral(geom20122->GetXaxis()->FindFixBin(normLow),geom20122->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfOFF_E = geom20112_sfOFF_E->Integral(geom20112_sfOFF_E->GetXaxis()->FindFixBin(normLow),geom20112_sfOFF_E->GetXaxis()->FindFixBin(normHigh))/geom20122_sfOFF_E->Integral(geom20122_sfOFF_E->GetXaxis()->FindFixBin(normLow),geom20122_sfOFF_E->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfON_E = geom20112_sfON_E->Integral(geom20112_sfON_E->GetXaxis()->FindFixBin(normLow),geom20112_sfON_E->GetXaxis()->FindFixBin(normHigh))/geom20122_sfON_E->Integral(geom20122_sfON_E->GetXaxis()->FindFixBin(normLow),geom20122_sfON_E->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfOFF_W = geom20112_sfOFF_W->Integral(geom20112_sfOFF_W->GetXaxis()->FindFixBin(normLow),geom20112_sfOFF_W->GetXaxis()->FindFixBin(normHigh))/geom20122_sfOFF_W->Integral(geom20122_sfOFF_W->GetXaxis()->FindFixBin(normLow),geom20122_sfOFF_W->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfON_W = geom20112_sfON_W->Integral(geom20112_sfON_W->GetXaxis()->FindFixBin(normLow),geom20112_sfON_W->GetXaxis()->FindFixBin(normHigh))/geom20122_sfON_W->Integral(geom20122_sfON_W->GetXaxis()->FindFixBin(normLow),geom20122_sfON_W->GetXaxis()->FindFixBin(normHigh));
  }

  if (normType==TString("3")) {
    normFactor = geom20113->Integral(geom20113->GetXaxis()->FindFixBin(normLow),geom20113->GetXaxis()->FindFixBin(normHigh))/geom20123->Integral(geom20123->GetXaxis()->FindFixBin(normLow),geom20123->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfOFF_E = geom20113_sfOFF_E->Integral(geom20113_sfOFF_E->GetXaxis()->FindFixBin(normLow),geom20113_sfOFF_E->GetXaxis()->FindFixBin(normHigh))/geom20123_sfOFF_E->Integral(geom20123_sfOFF_E->GetXaxis()->FindFixBin(normLow),geom20123_sfOFF_E->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfON_E = geom20113_sfON_E->Integral(geom20113_sfON_E->GetXaxis()->FindFixBin(normLow),geom20113_sfON_E->GetXaxis()->FindFixBin(normHigh))/geom20123_sfON_E->Integral(geom20123_sfON_E->GetXaxis()->FindFixBin(normLow),geom20123_sfON_E->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfOFF_W = geom20113_sfOFF_W->Integral(geom20113_sfOFF_W->GetXaxis()->FindFixBin(normLow),geom20113_sfOFF_W->GetXaxis()->FindFixBin(normHigh))/geom20123_sfOFF_W->Integral(geom20123_sfOFF_W->GetXaxis()->FindFixBin(normLow),geom20123_sfOFF_W->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfON_W = geom20113_sfON_W->Integral(geom20113_sfON_W->GetXaxis()->FindFixBin(normLow),geom20113_sfON_W->GetXaxis()->FindFixBin(normHigh))/geom20123_sfON_W->Integral(geom20123_sfON_W->GetXaxis()->FindFixBin(normLow),geom20123_sfON_W->GetXaxis()->FindFixBin(normHigh));
  }
  
  geom2011ALL->Scale(normFactor);
  geom20110->Scale(normFactor);
  geom20111->Scale(normFactor);
  geom20112->Scale(normFactor);
  geom20113->Scale(normFactor);
  geom2011ALL_sfOFF_E->Scale(normFactor_sfOFF_E);
  geom20110_sfOFF_E->Scale(normFactor_sfOFF_E);
  geom20111_sfOFF_E->Scale(normFactor_sfOFF_E);
  geom20112_sfOFF_E->Scale(normFactor_sfOFF_E);
  geom20113_sfOFF_E->Scale(normFactor_sfOFF_E);
  geom2011ALL_sfON_E->Scale(normFactor_sfON_E);
  geom20110_sfON_E->Scale(normFactor_sfON_E);
  geom20111_sfON_E->Scale(normFactor_sfON_E);
  geom20112_sfON_E->Scale(normFactor_sfON_E);
  geom20113_sfON_E->Scale(normFactor_sfON_E);
  geom2011ALL_sfOFF_W->Scale(normFactor_sfOFF_W);
  geom20110_sfOFF_W->Scale(normFactor_sfOFF_W);
  geom20111_sfOFF_W->Scale(normFactor_sfOFF_W);
  geom20112_sfOFF_W->Scale(normFactor_sfOFF_W);
  geom20113_sfOFF_W->Scale(normFactor_sfOFF_W);
  geom2011ALL_sfON_W->Scale(normFactor_sfON_W);
  geom20110_sfON_W->Scale(normFactor_sfON_W);
  geom20111_sfON_W->Scale(normFactor_sfON_W);
  geom20112_sfON_W->Scale(normFactor_sfON_W);
  geom20113_sfON_W->Scale(normFactor_sfON_W);
  

  ////////////////////// SUPER_SUMS ////////////////////////////////////////////

  //Residuals...
  TCanvas *c1 = new TCanvas("c1","c1", 1600., 600.);
  c1->Divide(2,1);
  c1->cd(2);
  Int_t nBins = geom2011ALL->GetNbinsX();
  Double_t Min = geom2011ALL->GetXaxis()->GetBinLowEdge(geom2011ALL->GetXaxis()->GetFirst());
  Double_t Max = geom2011ALL->GetXaxis()->GetBinUpEdge(geom2011ALL->GetXaxis()->GetLast());

  //residual
  TH1D *resid = new TH1D("resid","Residuals: geom2011 - geom2012", nBins, Min, Max);
  resid->Add(geom2012ALL,geom2011ALL,-1,1);
  resid->GetXaxis()->SetRangeUser(0., xAxisMax);
  resid->GetYaxis()->SetTitle("event rate (mHz/keV)");
  //resid->SetMaximum(2.);
  //resid->SetMinimum(-2.);
  resid->SetLineWidth(2);
  resid->SetMarkerColor(col2011);
  resid->SetMarkerStyle(34);
  resid->SetMarkerSize(1.);

  TH1D *perc_resid = new TH1D("perc_resid","East Fractional Residuals: (geom2012/geom2013 - geom2011/geom2012)/MC", nBins, Min, Max);
  perc_resid->Add(geom2012ALL,geom2011ALL,1,-1);
  perc_resid->Divide(perc_resid,geom2012ALL,1.,1.);
  perc_resid->GetXaxis()->SetRangeUser(0., xAxisMax);
  //perc_resid->GetYaxis()->SetTitle("event rate (mHz/keV)");
  perc_resid->SetMaximum(.1);
  perc_resid->SetMinimum(-.1);
  perc_resid->SetLineWidth(2);
  


  //perc_resid->Draw();
  resid->Draw("E0");
  c1->Update();
  TLine *l = new TLine(Min, 0., xAxisMax, 0.);
  l->SetLineStyle(8);
  l->Draw();

  c1->cd(1);
  
  geom2011ALL->SetTitle("Super-Sum All Event Types");
  geom2011ALL->SetMarkerColor(col2011);
  geom2011ALL->SetMarkerStyle(22);
  geom2011ALL->SetMarkerSize(0.75);
  geom2011ALL->SetLineColor(col2011);
  //geom2011ALL->SetFillStyle(3002);
  //geom2011ALL->SetFillColor(col2011);
  geom2011ALL->SetLineWidth(3);
  geom2012ALL->SetMarkerColor(col2012);
  geom2012ALL->SetLineColor(col2012);
  geom2012ALL->SetMarkerSize(0.75);
  geom2012ALL->SetMarkerStyle(20);
  geom2011ALL->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom2011ALL->SetMaximum(geom2011ALL->GetMaximum()*1.2);
  geom2011ALL->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom2012ALL->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom2011ALL->Draw("HISTE0");
  geom2012ALL->Draw("SAMEE0");


  TCanvas *c2 = new TCanvas("c2", "c2", 1800., 600.);
  c2->Divide(4,1);
  c2->cd(1);
  
  

  geom20110->SetMarkerColor(col2011);
  geom20110->SetMarkerStyle(22);
  geom20110->SetMarkerSize(0.75);
  geom20110->SetLineColor(col2011);
  //geom20110->SetFillStyle(3002);
  //geom20110->SetFillColor(col2011);
  geom20110->SetLineWidth(3);
  geom20120->SetMarkerColor(col2012);
  geom20120->SetLineColor(col2012);
  geom20120->SetMarkerSize(0.75);
  geom20120->SetMarkerStyle(20);
  geom20110->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20110->SetMaximum(geom20110->GetMaximum()*1.2);
  geom20110->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom20120->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20110->Draw("HISTE0");
  geom20120->Draw("SAMEE0");
  
  c2->cd(2);
  

  geom20111->SetMarkerColor(col2011);
  geom20111->SetMarkerStyle(22);
  geom20111->SetMarkerSize(0.75);
  geom20111->SetLineColor(col2011);
  //geom20111->SetFillStyle(3002);
  //geom20111->SetFillColor(col2011);
  geom20111->SetLineWidth(3);
  geom20111->SetMinimum(0.);
  geom20121->SetMarkerColor(col2012);
  geom20121->SetLineColor(col2012);
  geom20121->SetMarkerSize(0.75);
  geom20121->SetMarkerStyle(20);
  geom20111->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20111->SetMaximum(geom20111->GetMaximum()*1.2);
  geom20111->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom20121->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20111->Draw("HISTE0");
  geom20121->Draw("SAMEE0");
  

  c2->cd(3);

  
  geom20112->SetMarkerColor(col2011);
  geom20112->SetMarkerStyle(22);
  geom20112->SetMarkerSize(0.75);
  geom20112->SetLineColor(col2011);
  //geom20112->SetFillStyle(3002);
  //geom20112->SetFillColor(col2011);
  geom20112->SetLineWidth(3);
  geom20112->SetMinimum(0.);
  geom20122->SetMarkerColor(col2012);
  geom20122->SetLineColor(col2012);
  geom20122->SetMarkerSize(0.75);
  geom20122->SetMarkerStyle(20);
  geom20112->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20112->SetMaximum(geom20112->GetMaximum()*1.2);
  geom20112->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom20122->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20112->Draw("HISTE0");
  geom20122->Draw("SAMEE0");

   c2->cd(4);

  
  geom20113->SetMarkerColor(col2011);
  geom20113->SetMarkerStyle(22);
  geom20113->SetMarkerSize(0.75);
  geom20113->SetLineColor(col2011);
  //geom20113->SetFillStyle(3002);
  //geom20113->SetFillColor(col2011);
  geom20113->SetLineWidth(3);
  geom20113->SetMinimum(0.);
  geom20123->SetMarkerColor(col2012);
  geom20123->SetLineColor(col2012);
  geom20123->SetMarkerSize(0.75);
  geom20123->SetMarkerStyle(20);
  geom20113->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20113->SetMaximum(geom20113->GetMaximum()*1.2);
  geom20113->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom20123->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20113->Draw("HISTE0");
  geom20123->Draw("SAMEE0");


  TString pdfFile = TString::Format("spectraComp_%s_Type%s_%0.0f-%0.0f.pdf",simORdata.Data(),normType.Data(), normLow, normHigh);
  //TString pdfFileStart = pdfFile + TString("(");
  //TString pdfFileEnd = pdfFile + TString
  c1->Print(TString::Format("%s(",pdfFile.Data()));
  c2->Print(TString::Format("%s)",pdfFile.Data()));


  ////////////////////// flipper OFF East ////////////////////////////////////////////

  //Residuals...
  TCanvas *c1_sfOFF_E = new TCanvas("c1_sfOFF_E","c1_sfOFF_E", 1600., 600.);
  c1_sfOFF_E->Divide(2,1);
  c1_sfOFF_E->cd(2);
  nBins = geom2011ALL_sfOFF_E->GetNbinsX();
  Min = geom2011ALL_sfOFF_E->GetXaxis()->GetBinLowEdge(geom2011ALL->GetXaxis()->GetFirst());
  Max = geom2011ALL_sfOFF_E->GetXaxis()->GetBinUpEdge(geom2011ALL->GetXaxis()->GetLast());

  //residual
  TH1D *resid_sfOFF_E = new TH1D("resid_sfOFF_E","Residuals: geom2012-geom2011", nBins, Min, Max);
  resid_sfOFF_E->Add(geom2012ALL_sfOFF_E,geom2011ALL_sfOFF_E,1,-1);
  resid_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  resid_sfOFF_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  //resid_sfOFF_E->SetMaximum(2.);
  //resid_sfOFF_E->SetMinimum(-2.);
  resid_sfOFF_E->SetLineWidth(2);
  resid_sfOFF_E->SetMarkerColor(col2011);
  resid_sfOFF_E->SetMarkerStyle(34);
  resid_sfOFF_E->SetMarkerSize(1.);

  TH1D *perc_resid_sfOFF_E = new TH1D("perc_resid_sfOFF_E","East Fractional Residuals: (MC-geom2011)/MC", nBins, Min, Max);
  perc_resid_sfOFF_E->Add(geom2012ALL_sfOFF_E,geom2011ALL_sfOFF_E,1,-1);
  perc_resid_sfOFF_E->Divide(perc_resid,geom2012ALL_sfOFF_E,1.,1.);
  perc_resid_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  //perc_resid_sfOFF_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  perc_resid_sfOFF_E->SetMaximum(.1);
  perc_resid_sfOFF_E->SetMinimum(-.1);
  perc_resid_sfOFF_E->SetLineWidth(2);
  


  //perc_resid_sfOFF_E->Draw();
  resid_sfOFF_E->Draw("E0");
  c1_sfOFF_E->Update();
  l = new TLine(Min, 0., xAxisMax, 0.);
  l->SetLineStyle(8);
  l->Draw();

  c1_sfOFF_E->cd(1);
  
  geom2011ALL_sfOFF_E->SetTitle("East Spin Flipper OFF: All Event Types");
  geom2011ALL_sfOFF_E->SetMarkerColor(col2011);
  geom2011ALL_sfOFF_E->SetMarkerStyle(22);
  geom2011ALL_sfOFF_E->SetMarkerSize(0.75);
  geom2011ALL_sfOFF_E->SetLineColor(col2011);
  //geom2011ALL_sfOFF_E->SetFillStyle(3002);
  //geom2011ALL_sfOFF_E->SetFillColor(col2011);
  geom2011ALL_sfOFF_E->SetLineWidth(3);
  geom2012ALL_sfOFF_E->SetMarkerColor(col2012);
  geom2012ALL_sfOFF_E->SetLineColor(col2012);
  geom2012ALL_sfOFF_E->SetMarkerSize(0.75);
  geom2012ALL_sfOFF_E->SetMarkerStyle(20);
  geom2011ALL_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom2011ALL_sfOFF_E->SetMaximum(geom2011ALL_sfOFF_E->GetMaximum()*1.2);
  geom2011ALL_sfOFF_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom2012ALL_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom2011ALL_sfOFF_E->Draw("HISTE0");
  geom2012ALL_sfOFF_E->Draw("SAMEE0");


  TCanvas *c2_sfOFF_E = new TCanvas("c2_sfOFF_E", "c2_sfOFF_E", 1600., 600.);
  c2_sfOFF_E->Divide(4,1);
  c2_sfOFF_E->cd(1);
  
  

  geom20110_sfOFF_E->SetMarkerColor(col2011);
  geom20110_sfOFF_E->SetMarkerStyle(22);
  geom20110_sfOFF_E->SetMarkerSize(0.75);
  geom20110_sfOFF_E->SetLineColor(col2011);
  //geom20110_sfOFF_E->SetFillStyle(3002);
  //geom20110_sfOFF_E->SetFillColor(col2011);
  geom20110_sfOFF_E->SetLineWidth(3);
  geom20120_sfOFF_E->SetMarkerColor(col2012);
  geom20120_sfOFF_E->SetLineColor(col2012);
  geom20120_sfOFF_E->SetMarkerSize(0.75);
  geom20120_sfOFF_E->SetMarkerStyle(20);
  geom20110_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20110_sfOFF_E->SetMaximum(geom20110_sfOFF_E->GetMaximum()*1.2);
  geom20110_sfOFF_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom20120_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20110_sfOFF_E->Draw("HISTE0");
  geom20120_sfOFF_E->Draw("SAMEE0");
  
  c2_sfOFF_E->cd(2);
  

  geom20111_sfOFF_E->SetMarkerColor(col2011);
  geom20111_sfOFF_E->SetMarkerStyle(22);
  geom20111_sfOFF_E->SetMarkerSize(0.75);
  geom20111_sfOFF_E->SetLineColor(col2011);
  //geom20111_sfOFF_E->SetFillStyle(3002);
  //geom20111_sfOFF_E->SetFillColor(col2011);
  geom20111_sfOFF_E->SetLineWidth(3);
  geom20111_sfOFF_E->SetMinimum(0.);
  geom20121_sfOFF_E->SetMarkerColor(col2012);
  geom20121_sfOFF_E->SetLineColor(col2012);
  geom20121_sfOFF_E->SetMarkerSize(0.75);
  geom20121_sfOFF_E->SetMarkerStyle(20);
  geom20111_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20111_sfOFF_E->SetMaximum(geom20111_sfOFF_E->GetMaximum()*1.2);
  geom20111_sfOFF_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom20121_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20111_sfOFF_E->Draw("HISTE0");
  geom20121_sfOFF_E->Draw("SAMEE0");
  

  c2_sfOFF_E->cd(3);

  
  geom20112_sfOFF_E->SetMarkerColor(col2011);
  geom20112_sfOFF_E->SetMarkerStyle(22);
  geom20112_sfOFF_E->SetMarkerSize(0.75);
  geom20112_sfOFF_E->SetLineColor(col2011);
  //geom20112_sfOFF_E->SetFillStyle(3002);
  //geom20112_sfOFF_E->SetFillColor(col2011);
  geom20112_sfOFF_E->SetLineWidth(3);
  geom20112_sfOFF_E->SetMinimum(0.);
  geom20122_sfOFF_E->SetMarkerColor(col2012);
  geom20122_sfOFF_E->SetLineColor(col2012);
  geom20122_sfOFF_E->SetMarkerSize(0.75);
  geom20122_sfOFF_E->SetMarkerStyle(20);
  geom20112_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20112_sfOFF_E->SetMaximum(geom20112_sfOFF_E->GetMaximum()*1.2);
  geom20112_sfOFF_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom20122_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20112_sfOFF_E->Draw("HISTE0");
  geom20122_sfOFF_E->Draw("SAMEE0");

  c2_sfOFF_E->cd(4);

  
  geom20113_sfOFF_E->SetMarkerColor(col2011);
  geom20113_sfOFF_E->SetMarkerStyle(22);
  geom20113_sfOFF_E->SetMarkerSize(0.75);
  geom20113_sfOFF_E->SetLineColor(col2011);
  //geom20113_sfOFF_E->SetFillStyle(3002);
  //geom20113_sfOFF_E->SetFillColor(col2011);
  geom20113_sfOFF_E->SetLineWidth(3);
  geom20113_sfOFF_E->SetMinimum(0.);
  geom20123_sfOFF_E->SetMarkerColor(col2012);
  geom20123_sfOFF_E->SetLineColor(col2012);
  geom20123_sfOFF_E->SetMarkerSize(0.75);
  geom20123_sfOFF_E->SetMarkerStyle(20);
  geom20113_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20113_sfOFF_E->SetMaximum(geom20113_sfOFF_E->GetMaximum()*1.2);
  geom20113_sfOFF_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom20123_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20113_sfOFF_E->Draw("HISTE0");
  geom20123_sfOFF_E->Draw("SAMEE0");


  pdfFile = TString::Format("spectraComp_%s_Type%s_%0.0f-%0.0f_sfOFF_E.pdf",simORdata.Data(),normType.Data(), normLow, normHigh);
  //TString pdfFileStart = pdfFile + TString("(");
  //TString pdfFileEnd = pdfFile + TString
  c1_sfOFF_E->Print(TString::Format("%s(",pdfFile.Data()));
  c2_sfOFF_E->Print(TString::Format("%s)",pdfFile.Data()));


  ////////////////////// flipper ON East ////////////////////////////////////////////

  //Residuals...
  TCanvas *c1_sfON_E = new TCanvas("c1_sfON_E","c1_sfON_E", 1600., 600.);
  c1_sfON_E->Divide(2,1);
  c1_sfON_E->cd(2);
  nBins = geom2011ALL_sfON_E->GetNbinsX();
  Min = geom2011ALL_sfON_E->GetXaxis()->GetBinLowEdge(geom2011ALL->GetXaxis()->GetFirst());
  Max = geom2011ALL_sfON_E->GetXaxis()->GetBinUpEdge(geom2011ALL->GetXaxis()->GetLast());

  //residual
  TH1D *resid_sfON_E = new TH1D("resid_sfON_E","Residuals: geom2012-geom2011", nBins, Min, Max);
  resid_sfON_E->Add(geom2012ALL_sfON_E,geom2011ALL_sfON_E,1,-1);
  resid_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  resid_sfON_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  //resid_sfON_E->SetMaximum(2.);
  //resid_sfON_E->SetMinimum(-2.);
  resid_sfON_E->SetLineWidth(2);
  resid_sfON_E->SetMarkerColor(col2011);
  resid_sfON_E->SetMarkerStyle(34);
  resid_sfON_E->SetMarkerSize(1.);

  TH1D *perc_resid_sfON_E = new TH1D("perc_resid_sfON_E","East Fractional Residuals: (MC-geom2011)/MC", nBins, Min, Max);
  perc_resid_sfON_E->Add(geom2012ALL_sfON_E,geom2011ALL_sfON_E,1,-1);
  perc_resid_sfON_E->Divide(perc_resid,geom2012ALL_sfON_E,1.,1.);
  perc_resid_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  //perc_resid_sfON_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  perc_resid_sfON_E->SetMaximum(.1);
  perc_resid_sfON_E->SetMinimum(-.1);
  perc_resid_sfON_E->SetLineWidth(2);
  


  //perc_resid_sfON_E->Draw();
  resid_sfON_E->Draw("E0");
  c1_sfON_E->Update();
  l = new TLine(Min, 0., xAxisMax, 0.);
  l->SetLineStyle(8);
  l->Draw();

  c1_sfON_E->cd(1);
  
  geom2011ALL_sfON_E->SetTitle("East Spin Flipper ON: All Event Types");
  geom2011ALL_sfON_E->SetMarkerColor(col2011);
  geom2011ALL_sfON_E->SetMarkerStyle(22);
  geom2011ALL_sfON_E->SetMarkerSize(0.75);
  geom2011ALL_sfON_E->SetLineColor(col2011);
  //geom2011ALL_sfON_E->SetFillStyle(3002);
  //geom2011ALL_sfON_E->SetFillColor(col2011);
  geom2011ALL_sfON_E->SetLineWidth(3);
  geom2012ALL_sfON_E->SetMarkerColor(col2012);
  geom2012ALL_sfON_E->SetLineColor(col2012);
  geom2012ALL_sfON_E->SetMarkerSize(0.75);
  geom2012ALL_sfON_E->SetMarkerStyle(20);
  geom2011ALL_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom2011ALL_sfON_E->SetMaximum(geom2011ALL_sfON_E->GetMaximum()*1.2);
  geom2011ALL_sfON_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom2012ALL_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom2011ALL_sfON_E->Draw("HISTE0");
  geom2012ALL_sfON_E->Draw("SAMEE0");


  TCanvas *c2_sfON_E = new TCanvas("c2_sfON_E", "c2_sfON_E", 1600., 600.);
  c2_sfON_E->Divide(4,1);
  c2_sfON_E->cd(1);
  
  

  geom20110_sfON_E->SetMarkerColor(col2011);
  geom20110_sfON_E->SetMarkerStyle(22);
  geom20110_sfON_E->SetMarkerSize(0.75);
  geom20110_sfON_E->SetLineColor(col2011);
  //geom20110_sfON_E->SetFillStyle(3002);
  //geom20110_sfON_E->SetFillColor(col2011);
  geom20110_sfON_E->SetLineWidth(3);
  geom20120_sfON_E->SetMarkerColor(col2012);
  geom20120_sfON_E->SetLineColor(col2012);
  geom20120_sfON_E->SetMarkerSize(0.75);
  geom20120_sfON_E->SetMarkerStyle(20);
  geom20110_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20110_sfON_E->SetMaximum(geom20110_sfON_E->GetMaximum()*1.2);
  geom20110_sfON_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom20120_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20110_sfON_E->Draw("HISTE0");
  geom20120_sfON_E->Draw("SAMEE0");
  
  c2_sfON_E->cd(2);
  

  geom20111_sfON_E->SetMarkerColor(col2011);
  geom20111_sfON_E->SetMarkerStyle(22);
  geom20111_sfON_E->SetMarkerSize(0.75);
  geom20111_sfON_E->SetLineColor(col2011);
  //geom20111_sfON_E->SetFillStyle(3002);
  //geom20111_sfON_E->SetFillColor(col2011);
  geom20111_sfON_E->SetLineWidth(3);
  geom20111_sfON_E->SetMinimum(0.);
  geom20121_sfON_E->SetMarkerColor(col2012);
  geom20121_sfON_E->SetLineColor(col2012);
  geom20121_sfON_E->SetMarkerSize(0.75);
  geom20121_sfON_E->SetMarkerStyle(20);
  geom20111_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20111_sfON_E->SetMaximum(geom20111_sfON_E->GetMaximum()*1.2);
  geom20111_sfON_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom20121_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20111_sfON_E->Draw("HISTE0");
  geom20121_sfON_E->Draw("SAMEE0");
  

  c2_sfON_E->cd(3);

  
  geom20112_sfON_E->SetMarkerColor(col2011);
  geom20112_sfON_E->SetMarkerStyle(22);
  geom20112_sfON_E->SetMarkerSize(0.75);
  geom20112_sfON_E->SetLineColor(col2011);
  //geom20112_sfON_E->SetFillStyle(3002);
  //geom20112_sfON_E->SetFillColor(col2011);
  geom20112_sfON_E->SetLineWidth(3);
  geom20112_sfON_E->SetMinimum(0.);
  geom20122_sfON_E->SetMarkerColor(col2012);
  geom20122_sfON_E->SetLineColor(col2012);
  geom20122_sfON_E->SetMarkerSize(0.75);
  geom20122_sfON_E->SetMarkerStyle(20);
  geom20112_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20112_sfON_E->SetMaximum(geom20112_sfON_E->GetMaximum()*1.2);
  geom20112_sfON_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom20122_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20112_sfON_E->Draw("HISTE0");
  geom20122_sfON_E->Draw("SAMEE0");


  c2_sfON_E->cd(4);

  
  geom20113_sfON_E->SetMarkerColor(col2011);
  geom20113_sfON_E->SetMarkerStyle(22);
  geom20113_sfON_E->SetMarkerSize(0.75);
  geom20113_sfON_E->SetLineColor(col2011);
  //geom20113_sfON_E->SetFillStyle(3002);
  //geom20113_sfON_E->SetFillColor(col2011);
  geom20113_sfON_E->SetLineWidth(3);
  geom20113_sfON_E->SetMinimum(0.);
  geom20123_sfON_E->SetMarkerColor(col2012);
  geom20123_sfON_E->SetLineColor(col2012);
  geom20123_sfON_E->SetMarkerSize(0.75);
  geom20123_sfON_E->SetMarkerStyle(20);
  geom20113_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20113_sfON_E->SetMaximum(geom20113_sfON_E->GetMaximum()*1.2);
  geom20113_sfON_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom20123_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20113_sfON_E->Draw("HISTE0");
  geom20123_sfON_E->Draw("SAMEE0");


  pdfFile = TString::Format("spectraComp_%s_Type%s_%0.0f-%0.0f_sfON_E.pdf",simORdata.Data(),normType.Data(), normLow, normHigh);
  //TString pdfFileStart = pdfFile + TString("(");
  //TString pdfFileEnd = pdfFile + TString
  c1_sfON_E->Print(TString::Format("%s(",pdfFile.Data()));
  c2_sfON_E->Print(TString::Format("%s)",pdfFile.Data()));


  ////////////////////// flipper OFF West ////////////////////////////////////////////

  //Residuals...
  TCanvas *c1_sfOFF_W = new TCanvas("c1_sfOFF_W","c1_sfOFF_W", 1600., 600.);
  c1_sfOFF_W->Divide(2,1);
  c1_sfOFF_W->cd(2);
  nBins = geom2011ALL_sfOFF_W->GetNbinsX();
  Min = geom2011ALL_sfOFF_W->GetXaxis()->GetBinLowEdge(geom2011ALL->GetXaxis()->GetFirst());
  Max = geom2011ALL_sfOFF_W->GetXaxis()->GetBinUpEdge(geom2011ALL->GetXaxis()->GetLast());

  //residual
  TH1D *resid_sfOFF_W = new TH1D("resid_sfOFF_W","Residuals: geom2012-geom2011", nBins, Min, Max);
  resid_sfOFF_W->Add(geom2012ALL_sfOFF_W,geom2011ALL_sfOFF_W,1,-1);
  resid_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  resid_sfOFF_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  //resid_sfOFF_W->SetMaximum(2.);
  //resid_sfOFF_W->SetMinimum(-2.);
  resid_sfOFF_W->SetLineWidth(2);
  resid_sfOFF_W->SetMarkerColor(col2011);
  resid_sfOFF_W->SetMarkerStyle(34);
  resid_sfOFF_W->SetMarkerSize(1.);

  TH1D *perc_resid_sfOFF_W = new TH1D("perc_resid_sfOFF_W","East Fractional Residuals: (MC-geom2011)/MC", nBins, Min, Max);
  perc_resid_sfOFF_W->Add(geom2012ALL_sfOFF_W,geom2011ALL_sfOFF_W,1,-1);
  perc_resid_sfOFF_W->Divide(perc_resid,geom2012ALL_sfOFF_W,1.,1.);
  perc_resid_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  //perc_resid_sfOFF_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  perc_resid_sfOFF_W->SetMaximum(.1);
  perc_resid_sfOFF_W->SetMinimum(-.1);
  perc_resid_sfOFF_W->SetLineWidth(2);
  


  //perc_resid_sfOFF_W->Draw();
  resid_sfOFF_W->Draw("E0");
  c1_sfOFF_W->Update();
  l = new TLine(Min, 0., xAxisMax, 0.);
  l->SetLineStyle(8);
  l->Draw();

  c1_sfOFF_W->cd(1);
  
  geom2011ALL_sfOFF_W->SetTitle("West Spin Flipper OFF: All Event Types");
  geom2011ALL_sfOFF_W->SetMarkerColor(col2011);
  geom2011ALL_sfOFF_W->SetMarkerStyle(22);
  geom2011ALL_sfOFF_W->SetMarkerSize(0.75);
  geom2011ALL_sfOFF_W->SetLineColor(col2011);
  //geom2011ALL_sfOFF_W->SetFillStyle(3002);
  //geom2011ALL_sfOFF_W->SetFillColor(col2011);
  geom2011ALL_sfOFF_W->SetLineWidth(3);
  geom2012ALL_sfOFF_W->SetMarkerColor(col2012);
  geom2012ALL_sfOFF_W->SetLineColor(col2012);
  geom2012ALL_sfOFF_W->SetMarkerSize(0.75);
  geom2012ALL_sfOFF_W->SetMarkerStyle(20);
  geom2011ALL_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom2011ALL_sfOFF_W->SetMaximum(geom2011ALL_sfOFF_W->GetMaximum()*1.2);
  geom2011ALL_sfOFF_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom2012ALL_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom2011ALL_sfOFF_W->Draw("HISTE0");
  geom2012ALL_sfOFF_W->Draw("SAMEE0");


  TCanvas *c2_sfOFF_W = new TCanvas("c2_sfOFF_W", "c2_sfOFF_W", 1600., 600.);
  c2_sfOFF_W->Divide(4,1);
  c2_sfOFF_W->cd(1);
  
  

  geom20110_sfOFF_W->SetMarkerColor(col2011);
  geom20110_sfOFF_W->SetMarkerStyle(22);
  geom20110_sfOFF_W->SetMarkerSize(0.75);
  geom20110_sfOFF_W->SetLineColor(col2011);
  //geom20110_sfOFF_W->SetFillStyle(3002);
  //geom20110_sfOFF_W->SetFillColor(col2011);
  geom20110_sfOFF_W->SetLineWidth(3);
  geom20120_sfOFF_W->SetMarkerColor(col2012);
  geom20120_sfOFF_W->SetLineColor(col2012);
  geom20120_sfOFF_W->SetMarkerSize(0.75);
  geom20120_sfOFF_W->SetMarkerStyle(20);
  geom20110_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20110_sfOFF_W->SetMaximum(geom20110_sfOFF_W->GetMaximum()*1.2);
  geom20110_sfOFF_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom20120_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20110_sfOFF_W->Draw("HISTE0");
  geom20120_sfOFF_W->Draw("SAMEE0");
  
  c2_sfOFF_W->cd(2);
  

  geom20111_sfOFF_W->SetMarkerColor(col2011);
  geom20111_sfOFF_W->SetMarkerStyle(22);
  geom20111_sfOFF_W->SetMarkerSize(0.75);
  geom20111_sfOFF_W->SetLineColor(col2011);
  //geom20111_sfOFF_W->SetFillStyle(3002);
  //geom20111_sfOFF_W->SetFillColor(col2011);
  geom20111_sfOFF_W->SetLineWidth(3);
  geom20111_sfOFF_W->SetMinimum(0.);
  geom20121_sfOFF_W->SetMarkerColor(col2012);
  geom20121_sfOFF_W->SetLineColor(col2012);
  geom20121_sfOFF_W->SetMarkerSize(0.75);
  geom20121_sfOFF_W->SetMarkerStyle(20);
  geom20111_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20111_sfOFF_W->SetMaximum(geom20111_sfOFF_W->GetMaximum()*1.2);
  geom20111_sfOFF_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom20121_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20111_sfOFF_W->Draw("HISTE0");
  geom20121_sfOFF_W->Draw("SAMEE0");
  

  c2_sfOFF_W->cd(3);

  
  geom20112_sfOFF_W->SetMarkerColor(col2011);
  geom20112_sfOFF_W->SetMarkerStyle(22);
  geom20112_sfOFF_W->SetMarkerSize(0.75);
  geom20112_sfOFF_W->SetLineColor(col2011);
  //geom20112_sfOFF_W->SetFillStyle(3002);
  //geom20112_sfOFF_W->SetFillColor(col2011);
  geom20112_sfOFF_W->SetLineWidth(3);
  geom20112_sfOFF_W->SetMinimum(0.);
  geom20122_sfOFF_W->SetMarkerColor(col2012);
  geom20122_sfOFF_W->SetLineColor(col2012);
  geom20122_sfOFF_W->SetMarkerSize(0.75);
  geom20122_sfOFF_W->SetMarkerStyle(20);
  geom20112_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20112_sfOFF_W->SetMaximum(geom20112_sfOFF_W->GetMaximum()*1.2);
  geom20112_sfOFF_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom20122_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20112_sfOFF_W->Draw("HISTE0");
  geom20122_sfOFF_W->Draw("SAMEE0");


  c2_sfOFF_W->cd(4);

  
  geom20113_sfOFF_W->SetMarkerColor(col2011);
  geom20113_sfOFF_W->SetMarkerStyle(22);
  geom20113_sfOFF_W->SetMarkerSize(0.75);
  geom20113_sfOFF_W->SetLineColor(col2011);
  //geom20113_sfOFF_W->SetFillStyle(3002);
  //geom20113_sfOFF_W->SetFillColor(col2011);
  geom20113_sfOFF_W->SetLineWidth(3);
  geom20113_sfOFF_W->SetMinimum(0.);
  geom20123_sfOFF_W->SetMarkerColor(col2012);
  geom20123_sfOFF_W->SetLineColor(col2012);
  geom20123_sfOFF_W->SetMarkerSize(0.75);
  geom20123_sfOFF_W->SetMarkerStyle(20);
  geom20113_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20113_sfOFF_W->SetMaximum(geom20113_sfOFF_W->GetMaximum()*1.2);
  geom20113_sfOFF_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom20123_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20113_sfOFF_W->Draw("HISTE0");
  geom20123_sfOFF_W->Draw("SAMEE0");


  pdfFile = TString::Format("spectraComp_%s_Type%s_%0.0f-%0.0f_sfOFF_W.pdf",simORdata.Data(),normType.Data(), normLow, normHigh);
  //TString pdfFileStart = pdfFile + TString("(");
  //TString pdfFileEnd = pdfFile + TString
  c1_sfOFF_W->Print(TString::Format("%s(",pdfFile.Data()));
  c2_sfOFF_W->Print(TString::Format("%s)",pdfFile.Data()));


  ////////////////////// flipper ON West ////////////////////////////////////////////

  //Residuals...
  TCanvas *c1_sfON_W = new TCanvas("c1_sfON_W","c1_sfON_W", 1600., 600.);
  c1_sfON_W->Divide(2,1);
  c1_sfON_W->cd(2);
  nBins = geom2011ALL_sfON_W->GetNbinsX();
  Min = geom2011ALL_sfON_W->GetXaxis()->GetBinLowEdge(geom2011ALL->GetXaxis()->GetFirst());
  Max = geom2011ALL_sfON_W->GetXaxis()->GetBinUpEdge(geom2011ALL->GetXaxis()->GetLast());

  //residual
  TH1D *resid_sfON_W = new TH1D("resid_sfON_W","Residuals: geom2012-geom2011", nBins, Min, Max);
  resid_sfON_W->Add(geom2012ALL_sfON_W,geom2011ALL_sfON_W,1,-1);
  resid_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  resid_sfON_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  //resid_sfON_W->SetMaximum(2.);
  //resid_sfON_W->SetMinimum(-2.);
  resid_sfON_W->SetLineWidth(2);
  resid_sfON_W->SetMarkerColor(col2011);
  resid_sfON_W->SetMarkerStyle(34);
  resid_sfON_W->SetMarkerSize(1.);

  TH1D *perc_resid_sfON_W = new TH1D("perc_resid_sfON_W","East Fractional Residuals: (MC-geom2011)/MC", nBins, Min, Max);
  perc_resid_sfON_W->Add(geom2012ALL_sfON_W,geom2011ALL_sfON_W,1,-1);
  perc_resid_sfON_W->Divide(perc_resid,geom2012ALL_sfON_W,1.,1.);
  perc_resid_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  //perc_resid_sfON_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  perc_resid_sfON_W->SetMaximum(.1);
  perc_resid_sfON_W->SetMinimum(-.1);
  perc_resid_sfON_W->SetLineWidth(2);
  


  //perc_resid_sfON_W->Draw();
  resid_sfON_W->Draw("E0");
  c1_sfON_W->Update();
  l = new TLine(Min, 0., xAxisMax, 0.);
  l->SetLineStyle(8);
  l->Draw();

  c1_sfON_W->cd(1);
  
  geom2011ALL_sfON_W->SetTitle("West Spin Flipper ON: All Event Types");
  geom2011ALL_sfON_W->SetMarkerColor(col2011);
  geom2011ALL_sfON_W->SetMarkerStyle(22);
  geom2011ALL_sfON_W->SetMarkerSize(0.75);
  geom2011ALL_sfON_W->SetLineColor(col2011);
  //geom2011ALL_sfON_W->SetFillStyle(3002);
  //geom2011ALL_sfON_W->SetFillColor(col2011);
  geom2011ALL_sfON_W->SetLineWidth(3);
  geom2012ALL_sfON_W->SetMarkerColor(col2012);
  geom2012ALL_sfON_W->SetLineColor(col2012);
  geom2012ALL_sfON_W->SetMarkerSize(0.75);
  geom2012ALL_sfON_W->SetMarkerStyle(20);
  geom2011ALL_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom2011ALL_sfON_W->SetMaximum(geom2011ALL_sfON_W->GetMaximum()*1.2);
  geom2011ALL_sfON_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom2012ALL_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom2011ALL_sfON_W->Draw("HISTE0");
  geom2012ALL_sfON_W->Draw("SAMEE0");


  TCanvas *c2_sfON_W = new TCanvas("c2_sfON_W", "c2_sfON_W", 1600., 600.);
  c2_sfON_W->Divide(4,1);
  c2_sfON_W->cd(1);
  
  

  geom20110_sfON_W->SetMarkerColor(col2011);
  geom20110_sfON_W->SetMarkerStyle(22);
  geom20110_sfON_W->SetMarkerSize(0.75);
  geom20110_sfON_W->SetLineColor(col2011);
  //geom20110_sfON_W->SetFillStyle(3002);
  //geom20110_sfON_W->SetFillColor(col2011);
  geom20110_sfON_W->SetLineWidth(3);
  geom20120_sfON_W->SetMarkerColor(col2012);
  geom20120_sfON_W->SetLineColor(col2012);
  geom20120_sfON_W->SetMarkerSize(0.75);
  geom20120_sfON_W->SetMarkerStyle(20);
  geom20110_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20110_sfON_W->SetMaximum(geom20110_sfON_W->GetMaximum()*1.2);
  geom20110_sfON_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom20120_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20110_sfON_W->Draw("HISTE0");
  geom20120_sfON_W->Draw("SAMEE0");
  
  c2_sfON_W->cd(2);
  

  geom20111_sfON_W->SetMarkerColor(col2011);
  geom20111_sfON_W->SetMarkerStyle(22);
  geom20111_sfON_W->SetMarkerSize(0.75);
  geom20111_sfON_W->SetLineColor(col2011);
  //geom20111_sfON_W->SetFillStyle(3002);
  //geom20111_sfON_W->SetFillColor(col2011);
  geom20111_sfON_W->SetLineWidth(3);
  geom20111_sfON_W->SetMinimum(0.);
  geom20121_sfON_W->SetMarkerColor(col2012);
  geom20121_sfON_W->SetLineColor(col2012);
  geom20121_sfON_W->SetMarkerSize(0.75);
  geom20121_sfON_W->SetMarkerStyle(20);
  geom20111_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20111_sfON_W->SetMaximum(geom20111_sfON_W->GetMaximum()*1.2);
  geom20111_sfON_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom20121_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20111_sfON_W->Draw("HISTE0");
  geom20121_sfON_W->Draw("SAMEE0");
  

  c2_sfON_W->cd(3);

  
  geom20112_sfON_W->SetMarkerColor(col2011);
  geom20112_sfON_W->SetMarkerStyle(22);
  geom20112_sfON_W->SetMarkerSize(0.75);
  geom20112_sfON_W->SetLineColor(col2011);
  //geom20112_sfON_W->SetFillStyle(3002);
  //geom20112_sfON_W->SetFillColor(col2011);
  geom20112_sfON_W->SetLineWidth(3);
  geom20112_sfON_W->SetMinimum(0.);
  geom20122_sfON_W->SetMarkerColor(col2012);
  geom20122_sfON_W->SetLineColor(col2012);
  geom20122_sfON_W->SetMarkerSize(0.75);
  geom20122_sfON_W->SetMarkerStyle(20);
  geom20112_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20112_sfON_W->SetMaximum(geom20112_sfON_W->GetMaximum()*1.2);
  geom20112_sfON_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom20122_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20112_sfON_W->Draw("HISTE0");
  geom20122_sfON_W->Draw("SAMEE0");


  c2_sfON_W->cd(4);


  geom20113_sfON_W->SetMarkerColor(col2011);
  geom20113_sfON_W->SetMarkerStyle(22);
  geom20113_sfON_W->SetMarkerSize(0.75);
  geom20113_sfON_W->SetLineColor(col2011);
  //geom20113_sfON_W->SetFillStyle(3002);
  //geom20113_sfON_W->SetFillColor(col2011);
  geom20113_sfON_W->SetLineWidth(3);
  geom20113_sfON_W->SetMinimum(0.);
  geom20123_sfON_W->SetMarkerColor(col2012);
  geom20123_sfON_W->SetLineColor(col2012);
  geom20123_sfON_W->SetMarkerSize(0.75);
  geom20123_sfON_W->SetMarkerStyle(20);
  geom20113_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20113_sfON_W->SetMaximum(geom20113_sfON_W->GetMaximum()*1.2);
  geom20113_sfON_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  geom20123_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  geom20113_sfON_W->Draw("HISTE0");
  geom20123_sfON_W->Draw("SAMEE0");


  pdfFile = TString::Format("spectraComp_%s_Type%s_%0.0f-%0.0f_sfON_W.pdf",simORdata.Data(),normType.Data(), normLow, normHigh);
  //TString pdfFileStart = pdfFile + TString("(");
  //TString pdfFileEnd = pdfFile + TString
  c1_sfON_W->Print(TString::Format("%s(",pdfFile.Data()));
  c2_sfON_W->Print(TString::Format("%s)",pdfFile.Data()));
}
  
