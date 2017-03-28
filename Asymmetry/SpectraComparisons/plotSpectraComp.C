

{
  gStyle->SetOptStat(0);

  int octetStart = 0;
  int octetEnd = 59;

  TString normType = "ALL";

  Double_t normLow = 200.;
  Double_t normHigh = 780.;
  
  Double_t xAxisMax = 1200.;

  //Storing event fractions for data, E/W and sfON/OFF
  Double_t t0E_sfON, t0E_sfOFF, t1E_sfON, t1E_sfOFF, t2E_sfON, t2E_sfOFF, t3E_sfON, t3E_sfOFF;
  Double_t t0W_sfON, t0W_sfOFF, t1W_sfON, t1W_sfOFF, t2W_sfON, t2W_sfOFF, t3W_sfON, t3W_sfOFF;
  t0E_sfON = t0E_sfOFF = t1E_sfON = t1E_sfOFF = t23E_sfON = t2E_sfOFF = t3E_sfOFF = 0;
  t0W_sfON = t0W_sfOFF = t1W_sfON = t1W_sfOFF = t23W_sfON = t2W_sfOFF = t3E_sfOFF = 0.;
  
  TFile *data_file = new TFile(TString::Format("Octets_%i-%i_DATA.root",octetStart,octetEnd),"READ");
  TFile *sim_file = new TFile(TString::Format("Octets_%i-%i_SIM.root",octetStart,octetEnd),"READ");


  //Load the histograms

  TH1D *dataALL = (TH1D*)data_file->Get("EreconALL");
  TH1D *data0 = (TH1D*)data_file->Get("Erecon0");
  TH1D *data1 = (TH1D*)data_file->Get("Erecon1");
  TH1D *data2 = (TH1D*)data_file->Get("Erecon2");
  TH1D *data3 = (TH1D*)data_file->Get("Erecon3");
  
  TH1D *simALL = (TH1D*)sim_file->Get("EreconALL");
  TH1D *sim0 = (TH1D*)sim_file->Get("Erecon0");
  TH1D *sim1 = (TH1D*)sim_file->Get("Erecon1");
  TH1D *sim2 = (TH1D*)sim_file->Get("Erecon2");
  TH1D *sim3 = (TH1D*)sim_file->Get("Erecon3");
  
  TH1D *dataALL_sfOFF_E = (TH1D*)data_file->Get("EreconALL_sfOFF_E");
  TH1D *data0_sfOFF_E = (TH1D*)data_file->Get("Erecon0_sfOFF_E");
  TH1D *data1_sfOFF_E = (TH1D*)data_file->Get("Erecon1_sfOFF_E");
  TH1D *data2_sfOFF_E = (TH1D*)data_file->Get("Erecon2_sfOFF_E");
  TH1D *data3_sfOFF_E = (TH1D*)data_file->Get("Erecon3_sfOFF_E");
  
  TH1D *simALL_sfOFF_E = (TH1D*)sim_file->Get("EreconALL_sfOFF_E");
  TH1D *sim0_sfOFF_E = (TH1D*)sim_file->Get("Erecon0_sfOFF_E");
  TH1D *sim1_sfOFF_E = (TH1D*)sim_file->Get("Erecon1_sfOFF_E");
  TH1D *sim2_sfOFF_E = (TH1D*)sim_file->Get("Erecon2_sfOFF_E");
  TH1D *sim3_sfOFF_E = (TH1D*)sim_file->Get("Erecon3_sfOFF_E");

  TH1D *dataALL_sfOFF_W = (TH1D*)data_file->Get("EreconALL_sfOFF_W");
  TH1D *data0_sfOFF_W = (TH1D*)data_file->Get("Erecon0_sfOFF_W");
  TH1D *data1_sfOFF_W = (TH1D*)data_file->Get("Erecon1_sfOFF_W");
  TH1D *data2_sfOFF_W = (TH1D*)data_file->Get("Erecon2_sfOFF_W");
  TH1D *data3_sfOFF_W = (TH1D*)data_file->Get("Erecon3_sfOFF_W");
  
  TH1D *simALL_sfOFF_W = (TH1D*)sim_file->Get("EreconALL_sfOFF_W");
  TH1D *sim0_sfOFF_W = (TH1D*)sim_file->Get("Erecon0_sfOFF_W");
  TH1D *sim1_sfOFF_W = (TH1D*)sim_file->Get("Erecon1_sfOFF_W");
  TH1D *sim2_sfOFF_W = (TH1D*)sim_file->Get("Erecon2_sfOFF_W");
  TH1D *sim3_sfOFF_W = (TH1D*)sim_file->Get("Erecon3_sfOFF_W");
  
  TH1D *dataALL_sfON_E = (TH1D*)data_file->Get("EreconALL_sfON_E");
  TH1D *data0_sfON_E = (TH1D*)data_file->Get("Erecon0_sfON_E");
  TH1D *data1_sfON_E = (TH1D*)data_file->Get("Erecon1_sfON_E");
  TH1D *data2_sfON_E = (TH1D*)data_file->Get("Erecon2_sfON_E");
  TH1D *data3_sfON_E = (TH1D*)data_file->Get("Erecon3_sfON_E");
  
  TH1D *simALL_sfON_E = (TH1D*)sim_file->Get("EreconALL_sfON_E");
  TH1D *sim0_sfON_E = (TH1D*)sim_file->Get("Erecon0_sfON_E");
  TH1D *sim1_sfON_E = (TH1D*)sim_file->Get("Erecon1_sfON_E");
  TH1D *sim2_sfON_E = (TH1D*)sim_file->Get("Erecon2_sfON_E");
  TH1D *sim3_sfON_E = (TH1D*)sim_file->Get("Erecon3_sfON_E");

  TH1D *dataALL_sfON_W = (TH1D*)data_file->Get("EreconALL_sfON_W");
  TH1D *data0_sfON_W = (TH1D*)data_file->Get("Erecon0_sfON_W");
  TH1D *data1_sfON_W = (TH1D*)data_file->Get("Erecon1_sfON_W");
  TH1D *data2_sfON_W = (TH1D*)data_file->Get("Erecon2_sfON_W");
  TH1D *data3_sfON_W = (TH1D*)data_file->Get("Erecon3_sfON_W");
  
  TH1D *simALL_sfON_W = (TH1D*)sim_file->Get("EreconALL_sfON_W");
  TH1D *sim0_sfON_W = (TH1D*)sim_file->Get("Erecon0_sfON_W");
  TH1D *sim1_sfON_W = (TH1D*)sim_file->Get("Erecon1_sfON_W");
  TH1D *sim2_sfON_W = (TH1D*)sim_file->Get("Erecon2_sfON_W");
  TH1D *sim3_sfON_W = (TH1D*)sim_file->Get("Erecon3_sfON_W");
  

  t0E_sfOFF = data0_sfOFF_E->Integral(18,79);
  t0E_sfON = data0_sfON_E->Integral(18,79);
  t1E_sfOFF = data1_sfOFF_E->Integral(18,79);
  t1E_sfON = data1_sfON_E->Integral(18,79);
  t2E_sfOFF = data2_sfOFF_E->Integral(18,79);
  t2E_sfON = data2_sfON_E->Integral(18,79);
  t3E_sfOFF = data3_sfOFF_E->Integral(18,79);
  t3E_sfON = data3_sfON_E->Integral(18,79);
  t0W_sfOFF = data0_sfOFF_W->Integral(18,79);
  t0W_sfON = data0_sfON_W->Integral(18,79);
  t1W_sfOFF = data1_sfOFF_W->Integral(18,79);
  t1W_sfON = data1_sfON_W->Integral(18,79);
  t2W_sfOFF = data2_sfOFF_W->Integral(18,79);
  t2W_sfON = data2_sfON_W->Integral(18,79);
  t3W_sfOFF = data3_sfOFF_W->Integral(18,79);
  t3W_sfON = data3_sfON_W->Integral(18,79);
  

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
    normFactor = dataALL->Integral(dataALL->GetXaxis()->FindFixBin(normLow),dataALL->GetXaxis()->FindFixBin(normHigh))/simALL->Integral(simALL->GetXaxis()->FindFixBin(normLow),simALL->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfOFF_E = dataALL_sfOFF_E->Integral(dataALL_sfOFF_E->GetXaxis()->FindFixBin(normLow),dataALL_sfOFF_E->GetXaxis()->FindFixBin(normHigh))/simALL_sfOFF_E->Integral(simALL_sfOFF_E->GetXaxis()->FindFixBin(normLow),simALL_sfOFF_E->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfON_E = dataALL_sfON_E->Integral(dataALL_sfON_E->GetXaxis()->FindFixBin(normLow),dataALL_sfON_E->GetXaxis()->FindFixBin(normHigh))/simALL_sfON_E->Integral(simALL_sfON_E->GetXaxis()->FindFixBin(normLow),simALL_sfON_E->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfOFF_W = dataALL_sfOFF_W->Integral(dataALL_sfOFF_W->GetXaxis()->FindFixBin(normLow),dataALL_sfOFF_W->GetXaxis()->FindFixBin(normHigh))/simALL_sfOFF_W->Integral(simALL_sfOFF_W->GetXaxis()->FindFixBin(normLow),simALL_sfOFF_W->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfON_W = dataALL_sfON_W->Integral(dataALL_sfON_W->GetXaxis()->FindFixBin(normLow),dataALL_sfON_W->GetXaxis()->FindFixBin(normHigh))/simALL_sfON_W->Integral(simALL_sfON_W->GetXaxis()->FindFixBin(normLow),simALL_sfON_W->GetXaxis()->FindFixBin(normHigh));
    
  }

  if (normType==TString("0")) {
    normFactor = data0->Integral(data0->GetXaxis()->FindFixBin(normLow),data0->GetXaxis()->FindFixBin(normHigh))/sim0->Integral(sim0->GetXaxis()->FindFixBin(normLow),sim0->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfOFF_E = data0_sfOFF_E->Integral(data0_sfOFF_E->GetXaxis()->FindFixBin(normLow),data0_sfOFF_E->GetXaxis()->FindFixBin(normHigh))/sim0_sfOFF_E->Integral(sim0_sfOFF_E->GetXaxis()->FindFixBin(normLow),sim0_sfOFF_E->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfON_E = data0_sfON_E->Integral(data0_sfON_E->GetXaxis()->FindFixBin(normLow),data0_sfON_E->GetXaxis()->FindFixBin(normHigh))/sim0_sfON_E->Integral(sim0_sfON_E->GetXaxis()->FindFixBin(normLow),sim0_sfON_E->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfOFF_W = data0_sfOFF_W->Integral(data0_sfOFF_W->GetXaxis()->FindFixBin(normLow),data0_sfOFF_W->GetXaxis()->FindFixBin(normHigh))/sim0_sfOFF_W->Integral(sim0_sfOFF_W->GetXaxis()->FindFixBin(normLow),sim0_sfOFF_W->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfON_W = data0_sfON_W->Integral(data0_sfON_W->GetXaxis()->FindFixBin(normLow),data0_sfON_W->GetXaxis()->FindFixBin(normHigh))/sim0_sfON_W->Integral(sim0_sfON_W->GetXaxis()->FindFixBin(normLow),sim0_sfON_W->GetXaxis()->FindFixBin(normHigh));
  }
  
  if (normType==TString("1")) {
    normFactor = data1->Integral(data1->GetXaxis()->FindFixBin(normLow),data1->GetXaxis()->FindFixBin(normHigh))/sim1->Integral(sim1->GetXaxis()->FindFixBin(normLow),sim1->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfOFF_E = data1_sfOFF_E->Integral(data1_sfOFF_E->GetXaxis()->FindFixBin(normLow),data1_sfOFF_E->GetXaxis()->FindFixBin(normHigh))/sim1_sfOFF_E->Integral(sim1_sfOFF_E->GetXaxis()->FindFixBin(normLow),sim1_sfOFF_E->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfON_E = data1_sfON_E->Integral(data1_sfON_E->GetXaxis()->FindFixBin(normLow),data1_sfON_E->GetXaxis()->FindFixBin(normHigh))/sim1_sfON_E->Integral(sim1_sfON_E->GetXaxis()->FindFixBin(normLow),sim1_sfON_E->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfOFF_W = data1_sfOFF_W->Integral(data1_sfOFF_W->GetXaxis()->FindFixBin(normLow),data1_sfOFF_W->GetXaxis()->FindFixBin(normHigh))/sim1_sfOFF_W->Integral(sim1_sfOFF_W->GetXaxis()->FindFixBin(normLow),sim1_sfOFF_W->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfON_W = data1_sfON_W->Integral(data1_sfON_W->GetXaxis()->FindFixBin(normLow),data1_sfON_W->GetXaxis()->FindFixBin(normHigh))/sim1_sfON_W->Integral(sim1_sfON_W->GetXaxis()->FindFixBin(normLow),sim1_sfON_W->GetXaxis()->FindFixBin(normHigh));
  }
  

  if (normType==TString("2")) {
    normFactor = data2->Integral(data2->GetXaxis()->FindFixBin(normLow),data2->GetXaxis()->FindFixBin(normHigh))/sim2->Integral(sim2->GetXaxis()->FindFixBin(normLow),sim2->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfOFF_E = data2_sfOFF_E->Integral(data2_sfOFF_E->GetXaxis()->FindFixBin(normLow),data2_sfOFF_E->GetXaxis()->FindFixBin(normHigh))/sim2_sfOFF_E->Integral(sim2_sfOFF_E->GetXaxis()->FindFixBin(normLow),sim2_sfOFF_E->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfON_E = data2_sfON_E->Integral(data2_sfON_E->GetXaxis()->FindFixBin(normLow),data2_sfON_E->GetXaxis()->FindFixBin(normHigh))/sim2_sfON_E->Integral(sim2_sfON_E->GetXaxis()->FindFixBin(normLow),sim2_sfON_E->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfOFF_W = data2_sfOFF_W->Integral(data2_sfOFF_W->GetXaxis()->FindFixBin(normLow),data2_sfOFF_W->GetXaxis()->FindFixBin(normHigh))/sim2_sfOFF_W->Integral(sim2_sfOFF_W->GetXaxis()->FindFixBin(normLow),sim2_sfOFF_W->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfON_W = data2_sfON_W->Integral(data2_sfON_W->GetXaxis()->FindFixBin(normLow),data2_sfON_W->GetXaxis()->FindFixBin(normHigh))/sim2_sfON_W->Integral(sim2_sfON_W->GetXaxis()->FindFixBin(normLow),sim2_sfON_W->GetXaxis()->FindFixBin(normHigh));
  }

  if (normType==TString("3")) {
    normFactor = data3->Integral(data3->GetXaxis()->FindFixBin(normLow),data3->GetXaxis()->FindFixBin(normHigh))/sim3->Integral(sim3->GetXaxis()->FindFixBin(normLow),sim3->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfOFF_E = data3_sfOFF_E->Integral(data3_sfOFF_E->GetXaxis()->FindFixBin(normLow),data3_sfOFF_E->GetXaxis()->FindFixBin(normHigh))/sim3_sfOFF_E->Integral(sim3_sfOFF_E->GetXaxis()->FindFixBin(normLow),sim3_sfOFF_E->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfON_E = data3_sfON_E->Integral(data3_sfON_E->GetXaxis()->FindFixBin(normLow),data3_sfON_E->GetXaxis()->FindFixBin(normHigh))/sim3_sfON_E->Integral(sim3_sfON_E->GetXaxis()->FindFixBin(normLow),sim3_sfON_E->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfOFF_W = data3_sfOFF_W->Integral(data3_sfOFF_W->GetXaxis()->FindFixBin(normLow),data3_sfOFF_W->GetXaxis()->FindFixBin(normHigh))/sim3_sfOFF_W->Integral(sim3_sfOFF_W->GetXaxis()->FindFixBin(normLow),sim3_sfOFF_W->GetXaxis()->FindFixBin(normHigh));
    normFactor_sfON_W = data3_sfON_W->Integral(data3_sfON_W->GetXaxis()->FindFixBin(normLow),data3_sfON_W->GetXaxis()->FindFixBin(normHigh))/sim3_sfON_W->Integral(sim3_sfON_W->GetXaxis()->FindFixBin(normLow),sim3_sfON_W->GetXaxis()->FindFixBin(normHigh));
  }
  
  simALL->Scale(normFactor);
  sim0->Scale(normFactor);
  sim1->Scale(normFactor);
  sim2->Scale(normFactor);
  sim3->Scale(normFactor);
  simALL_sfOFF_E->Scale(normFactor_sfOFF_E);
  sim0_sfOFF_E->Scale(normFactor_sfOFF_E);
  sim1_sfOFF_E->Scale(normFactor_sfOFF_E);
  sim2_sfOFF_E->Scale(normFactor_sfOFF_E);
  sim3_sfOFF_E->Scale(normFactor_sfOFF_E);
  simALL_sfON_E->Scale(normFactor_sfON_E);
  sim0_sfON_E->Scale(normFactor_sfON_E);
  sim1_sfON_E->Scale(normFactor_sfON_E);
  sim2_sfON_E->Scale(normFactor_sfON_E);
  sim3_sfON_E->Scale(normFactor_sfON_E);
  simALL_sfOFF_W->Scale(normFactor_sfOFF_W);
  sim0_sfOFF_W->Scale(normFactor_sfOFF_W);
  sim1_sfOFF_W->Scale(normFactor_sfOFF_W);
  sim2_sfOFF_W->Scale(normFactor_sfOFF_W);
  sim3_sfOFF_W->Scale(normFactor_sfOFF_W);
  simALL_sfON_W->Scale(normFactor_sfON_W);
  sim0_sfON_W->Scale(normFactor_sfON_W);
  sim1_sfON_W->Scale(normFactor_sfON_W);
  sim2_sfON_W->Scale(normFactor_sfON_W);
  sim3_sfON_W->Scale(normFactor_sfON_W);
  

  ////////////////////// SUPER_SUMS ////////////////////////////////////////////

  //Residuals...
  TCanvas *c1 = new TCanvas("c1","c1", 1600., 600.);
  c1->Divide(2,1);
  c1->cd(2);
  Int_t nBins = dataALL->GetNbinsX();
  Double_t Min = dataALL->GetXaxis()->GetBinLowEdge(dataALL->GetXaxis()->GetFirst());
  Double_t Max = dataALL->GetXaxis()->GetBinUpEdge(dataALL->GetXaxis()->GetLast());

  //residual
  TH1D *resid = new TH1D("resid","Residuals: MC-Data", nBins, Min, Max);
  resid->Add(simALL,dataALL,1,-1);
  resid->GetXaxis()->SetRangeUser(0., xAxisMax);
  resid->GetYaxis()->SetTitle("event rate (mHz/keV)");
  //resid->SetMaximum(2.);
  //resid->SetMinimum(-2.);
  resid->SetLineWidth(2);
  resid->SetMarkerColor(kBlue);
  resid->SetMarkerStyle(34);
  resid->SetMarkerSize(1.);

  TH1D *perc_resid = new TH1D("perc_resid","East Fractional Residuals: (MC-Data)/MC", nBins, Min, Max);
  perc_resid->Add(simALL,dataALL,1,-1);
  perc_resid->Divide(perc_resid,simALL,1.,1.);
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
  
  dataALL->SetTitle("Super-Sum All Event Types");
  dataALL->SetMarkerColor(kBlue);
  dataALL->SetMarkerStyle(22);
  dataALL->SetMarkerSize(0.75);
  dataALL->SetLineColor(kBlue);
  //dataALL->SetFillStyle(3002);
  //dataALL->SetFillColor(kBlue);
  dataALL->SetLineWidth(3);
  simALL->SetMarkerColor(kRed);
  simALL->SetLineColor(kRed);
  simALL->SetMarkerSize(0.75);
  simALL->SetMarkerStyle(20);
  dataALL->GetXaxis()->SetRangeUser(0., xAxisMax);
  dataALL->SetMaximum(dataALL->GetMaximum()*1.2);
  dataALL->GetYaxis()->SetTitle("event rate (mHz/keV)");
  simALL->GetXaxis()->SetRangeUser(0., xAxisMax);
  dataALL->Draw("HISTE0");
  simALL->Draw("SAMEE0");


  TCanvas *c2 = new TCanvas("c2", "c2", 1800., 600.);
  c2->Divide(4,1);
  c2->cd(1);
  
  

  data0->SetMarkerColor(kBlue);
  data0->SetMarkerStyle(22);
  data0->SetMarkerSize(0.75);
  data0->SetLineColor(kBlue);
  //data0->SetFillStyle(3002);
  //data0->SetFillColor(kBlue);
  data0->SetLineWidth(3);
  sim0->SetMarkerColor(kRed);
  sim0->SetLineColor(kRed);
  sim0->SetMarkerSize(0.75);
  sim0->SetMarkerStyle(20);
  data0->GetXaxis()->SetRangeUser(0., xAxisMax);
  data0->SetMaximum(data0->GetMaximum()*1.2);
  data0->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim0->GetXaxis()->SetRangeUser(0., xAxisMax);
  data0->Draw("HISTE0");
  sim0->Draw("SAMEE0");
  
  c2->cd(2);
  

  data1->SetMarkerColor(kBlue);
  data1->SetMarkerStyle(22);
  data1->SetMarkerSize(0.75);
  data1->SetLineColor(kBlue);
  //data1->SetFillStyle(3002);
  //data1->SetFillColor(kBlue);
  data1->SetLineWidth(3);
  data1->SetMinimum(0.);
  sim1->SetMarkerColor(kRed);
  sim1->SetLineColor(kRed);
  sim1->SetMarkerSize(0.75);
  sim1->SetMarkerStyle(20);
  data1->GetXaxis()->SetRangeUser(0., xAxisMax);
  data1->SetMaximum(data1->GetMaximum()*1.2);
  data1->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim1->GetXaxis()->SetRangeUser(0., xAxisMax);
  data1->Draw("HISTE0");
  sim1->Draw("SAMEE0");
  

  c2->cd(3);

  
  data2->SetMarkerColor(kBlue);
  data2->SetMarkerStyle(22);
  data2->SetMarkerSize(0.75);
  data2->SetLineColor(kBlue);
  //data2->SetFillStyle(3002);
  //data2->SetFillColor(kBlue);
  data2->SetLineWidth(3);
  data2->SetMinimum(0.);
  sim2->SetMarkerColor(kRed);
  sim2->SetLineColor(kRed);
  sim2->SetMarkerSize(0.75);
  sim2->SetMarkerStyle(20);
  data2->GetXaxis()->SetRangeUser(0., xAxisMax);
  data2->SetMaximum(data2->GetMaximum()*1.2);
  data2->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim2->GetXaxis()->SetRangeUser(0., xAxisMax);
  data2->Draw("HISTE0");
  sim2->Draw("SAMEE0");

   c2->cd(4);

  
  data3->SetMarkerColor(kBlue);
  data3->SetMarkerStyle(22);
  data3->SetMarkerSize(0.75);
  data3->SetLineColor(kBlue);
  //data3->SetFillStyle(3002);
  //data3->SetFillColor(kBlue);
  data3->SetLineWidth(3);
  data3->SetMinimum(0.);
  sim3->SetMarkerColor(kRed);
  sim3->SetLineColor(kRed);
  sim3->SetMarkerSize(0.75);
  sim3->SetMarkerStyle(20);
  data3->GetXaxis()->SetRangeUser(0., xAxisMax);
  data3->SetMaximum(data3->GetMaximum()*1.2);
  data3->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim3->GetXaxis()->SetRangeUser(0., xAxisMax);
  data3->Draw("HISTE0");
  sim3->Draw("SAMEE0");


  TString pdfFile = TString::Format("spectraComp_%i-%i_Type%s_%0.0f-%0.0f.pdf",octetStart,octetEnd,normType.Data(), normLow, normHigh);
  //TString pdfFileStart = pdfFile + TString("(");
  //TString pdfFileEnd = pdfFile + TString
  c1->Print(TString::Format("%s(",pdfFile.Data()));
  c2->Print(TString::Format("%s)",pdfFile.Data()));


  ////////////////////// flipper OFF East ////////////////////////////////////////////

  //Residuals...
  TCanvas *c1_sfOFF_E = new TCanvas("c1_sfOFF_E","c1_sfOFF_E", 1600., 600.);
  c1_sfOFF_E->Divide(2,1);
  c1_sfOFF_E->cd(2);
  nBins = dataALL_sfOFF_E->GetNbinsX();
  Min = dataALL_sfOFF_E->GetXaxis()->GetBinLowEdge(dataALL->GetXaxis()->GetFirst());
  Max = dataALL_sfOFF_E->GetXaxis()->GetBinUpEdge(dataALL->GetXaxis()->GetLast());

  //residual
  TH1D *resid_sfOFF_E = new TH1D("resid_sfOFF_E","Residuals: MC-Data", nBins, Min, Max);
  resid_sfOFF_E->Add(simALL_sfOFF_E,dataALL_sfOFF_E,1,-1);
  resid_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  resid_sfOFF_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  //resid_sfOFF_E->SetMaximum(2.);
  //resid_sfOFF_E->SetMinimum(-2.);
  resid_sfOFF_E->SetLineWidth(2);
  resid_sfOFF_E->SetMarkerColor(kBlue);
  resid_sfOFF_E->SetMarkerStyle(34);
  resid_sfOFF_E->SetMarkerSize(1.);

  TH1D *perc_resid_sfOFF_E = new TH1D("perc_resid_sfOFF_E","East Fractional Residuals: (MC-Data)/MC", nBins, Min, Max);
  perc_resid_sfOFF_E->Add(simALL_sfOFF_E,dataALL_sfOFF_E,1,-1);
  perc_resid_sfOFF_E->Divide(perc_resid,simALL_sfOFF_E,1.,1.);
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
  
  dataALL_sfOFF_E->SetTitle("East Spin Flipper OFF: All Event Types");
  dataALL_sfOFF_E->SetMarkerColor(kBlue);
  dataALL_sfOFF_E->SetMarkerStyle(22);
  dataALL_sfOFF_E->SetMarkerSize(0.75);
  dataALL_sfOFF_E->SetLineColor(kBlue);
  //dataALL_sfOFF_E->SetFillStyle(3002);
  //dataALL_sfOFF_E->SetFillColor(kBlue);
  dataALL_sfOFF_E->SetLineWidth(3);
  simALL_sfOFF_E->SetMarkerColor(kRed);
  simALL_sfOFF_E->SetLineColor(kRed);
  simALL_sfOFF_E->SetMarkerSize(0.75);
  simALL_sfOFF_E->SetMarkerStyle(20);
  dataALL_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  dataALL_sfOFF_E->SetMaximum(dataALL_sfOFF_E->GetMaximum()*1.2);
  dataALL_sfOFF_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  simALL_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  dataALL_sfOFF_E->Draw("HISTE0");
  simALL_sfOFF_E->Draw("SAMEE0");


  TCanvas *c2_sfOFF_E = new TCanvas("c2_sfOFF_E", "c2_sfOFF_E", 1600., 600.);
  c2_sfOFF_E->Divide(4,1);
  c2_sfOFF_E->cd(1);
  
  

  data0_sfOFF_E->SetMarkerColor(kBlue);
  data0_sfOFF_E->SetMarkerStyle(22);
  data0_sfOFF_E->SetMarkerSize(0.75);
  data0_sfOFF_E->SetLineColor(kBlue);
  //data0_sfOFF_E->SetFillStyle(3002);
  //data0_sfOFF_E->SetFillColor(kBlue);
  data0_sfOFF_E->SetLineWidth(3);
  sim0_sfOFF_E->SetMarkerColor(kRed);
  sim0_sfOFF_E->SetLineColor(kRed);
  sim0_sfOFF_E->SetMarkerSize(0.75);
  sim0_sfOFF_E->SetMarkerStyle(20);
  data0_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  data0_sfOFF_E->SetMaximum(data0_sfOFF_E->GetMaximum()*1.2);
  data0_sfOFF_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim0_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  data0_sfOFF_E->Draw("HISTE0");
  sim0_sfOFF_E->Draw("SAMEE0");
  
  c2_sfOFF_E->cd(2);
  

  data1_sfOFF_E->SetMarkerColor(kBlue);
  data1_sfOFF_E->SetMarkerStyle(22);
  data1_sfOFF_E->SetMarkerSize(0.75);
  data1_sfOFF_E->SetLineColor(kBlue);
  //data1_sfOFF_E->SetFillStyle(3002);
  //data1_sfOFF_E->SetFillColor(kBlue);
  data1_sfOFF_E->SetLineWidth(3);
  data1_sfOFF_E->SetMinimum(0.);
  sim1_sfOFF_E->SetMarkerColor(kRed);
  sim1_sfOFF_E->SetLineColor(kRed);
  sim1_sfOFF_E->SetMarkerSize(0.75);
  sim1_sfOFF_E->SetMarkerStyle(20);
  data1_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  data1_sfOFF_E->SetMaximum(data1_sfOFF_E->GetMaximum()*1.2);
  data1_sfOFF_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim1_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  data1_sfOFF_E->Draw("HISTE0");
  sim1_sfOFF_E->Draw("SAMEE0");
  

  c2_sfOFF_E->cd(3);

  
  data2_sfOFF_E->SetMarkerColor(kBlue);
  data2_sfOFF_E->SetMarkerStyle(22);
  data2_sfOFF_E->SetMarkerSize(0.75);
  data2_sfOFF_E->SetLineColor(kBlue);
  //data2_sfOFF_E->SetFillStyle(3002);
  //data2_sfOFF_E->SetFillColor(kBlue);
  data2_sfOFF_E->SetLineWidth(3);
  data2_sfOFF_E->SetMinimum(0.);
  sim2_sfOFF_E->SetMarkerColor(kRed);
  sim2_sfOFF_E->SetLineColor(kRed);
  sim2_sfOFF_E->SetMarkerSize(0.75);
  sim2_sfOFF_E->SetMarkerStyle(20);
  data2_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  data2_sfOFF_E->SetMaximum(data2_sfOFF_E->GetMaximum()*1.2);
  data2_sfOFF_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim2_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  data2_sfOFF_E->Draw("HISTE0");
  sim2_sfOFF_E->Draw("SAMEE0");

  c2_sfOFF_E->cd(4);

  
  data3_sfOFF_E->SetMarkerColor(kBlue);
  data3_sfOFF_E->SetMarkerStyle(22);
  data3_sfOFF_E->SetMarkerSize(0.75);
  data3_sfOFF_E->SetLineColor(kBlue);
  //data3_sfOFF_E->SetFillStyle(3002);
  //data3_sfOFF_E->SetFillColor(kBlue);
  data3_sfOFF_E->SetLineWidth(3);
  data3_sfOFF_E->SetMinimum(0.);
  sim3_sfOFF_E->SetMarkerColor(kRed);
  sim3_sfOFF_E->SetLineColor(kRed);
  sim3_sfOFF_E->SetMarkerSize(0.75);
  sim3_sfOFF_E->SetMarkerStyle(20);
  data3_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  data3_sfOFF_E->SetMaximum(data3_sfOFF_E->GetMaximum()*1.2);
  data3_sfOFF_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim3_sfOFF_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  data3_sfOFF_E->Draw("HISTE0");
  sim3_sfOFF_E->Draw("SAMEE0");


  pdfFile = TString::Format("spectraComp_%i-%i_Type%s_%0.0f-%0.0f_sfOFF_E.pdf",octetStart,octetEnd,normType.Data(), normLow, normHigh);
  //TString pdfFileStart = pdfFile + TString("(");
  //TString pdfFileEnd = pdfFile + TString
  c1_sfOFF_E->Print(TString::Format("%s(",pdfFile.Data()));
  c2_sfOFF_E->Print(TString::Format("%s)",pdfFile.Data()));


  ////////////////////// flipper ON East ////////////////////////////////////////////

  //Residuals...
  TCanvas *c1_sfON_E = new TCanvas("c1_sfON_E","c1_sfON_E", 1600., 600.);
  c1_sfON_E->Divide(2,1);
  c1_sfON_E->cd(2);
  nBins = dataALL_sfON_E->GetNbinsX();
  Min = dataALL_sfON_E->GetXaxis()->GetBinLowEdge(dataALL->GetXaxis()->GetFirst());
  Max = dataALL_sfON_E->GetXaxis()->GetBinUpEdge(dataALL->GetXaxis()->GetLast());

  //residual
  TH1D *resid_sfON_E = new TH1D("resid_sfON_E","Residuals: MC-Data", nBins, Min, Max);
  resid_sfON_E->Add(simALL_sfON_E,dataALL_sfON_E,1,-1);
  resid_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  resid_sfON_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  //resid_sfON_E->SetMaximum(2.);
  //resid_sfON_E->SetMinimum(-2.);
  resid_sfON_E->SetLineWidth(2);
  resid_sfON_E->SetMarkerColor(kBlue);
  resid_sfON_E->SetMarkerStyle(34);
  resid_sfON_E->SetMarkerSize(1.);

  TH1D *perc_resid_sfON_E = new TH1D("perc_resid_sfON_E","East Fractional Residuals: (MC-Data)/MC", nBins, Min, Max);
  perc_resid_sfON_E->Add(simALL_sfON_E,dataALL_sfON_E,1,-1);
  perc_resid_sfON_E->Divide(perc_resid,simALL_sfON_E,1.,1.);
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
  
  dataALL_sfON_E->SetTitle("East Spin Flipper ON: All Event Types");
  dataALL_sfON_E->SetMarkerColor(kBlue);
  dataALL_sfON_E->SetMarkerStyle(22);
  dataALL_sfON_E->SetMarkerSize(0.75);
  dataALL_sfON_E->SetLineColor(kBlue);
  //dataALL_sfON_E->SetFillStyle(3002);
  //dataALL_sfON_E->SetFillColor(kBlue);
  dataALL_sfON_E->SetLineWidth(3);
  simALL_sfON_E->SetMarkerColor(kRed);
  simALL_sfON_E->SetLineColor(kRed);
  simALL_sfON_E->SetMarkerSize(0.75);
  simALL_sfON_E->SetMarkerStyle(20);
  dataALL_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  dataALL_sfON_E->SetMaximum(dataALL_sfON_E->GetMaximum()*1.2);
  dataALL_sfON_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  simALL_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  dataALL_sfON_E->Draw("HISTE0");
  simALL_sfON_E->Draw("SAMEE0");


  TCanvas *c2_sfON_E = new TCanvas("c2_sfON_E", "c2_sfON_E", 1600., 600.);
  c2_sfON_E->Divide(4,1);
  c2_sfON_E->cd(1);
  
  

  data0_sfON_E->SetMarkerColor(kBlue);
  data0_sfON_E->SetMarkerStyle(22);
  data0_sfON_E->SetMarkerSize(0.75);
  data0_sfON_E->SetLineColor(kBlue);
  //data0_sfON_E->SetFillStyle(3002);
  //data0_sfON_E->SetFillColor(kBlue);
  data0_sfON_E->SetLineWidth(3);
  sim0_sfON_E->SetMarkerColor(kRed);
  sim0_sfON_E->SetLineColor(kRed);
  sim0_sfON_E->SetMarkerSize(0.75);
  sim0_sfON_E->SetMarkerStyle(20);
  data0_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  data0_sfON_E->SetMaximum(data0_sfON_E->GetMaximum()*1.2);
  data0_sfON_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim0_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  data0_sfON_E->Draw("HISTE0");
  sim0_sfON_E->Draw("SAMEE0");
  
  c2_sfON_E->cd(2);
  

  data1_sfON_E->SetMarkerColor(kBlue);
  data1_sfON_E->SetMarkerStyle(22);
  data1_sfON_E->SetMarkerSize(0.75);
  data1_sfON_E->SetLineColor(kBlue);
  //data1_sfON_E->SetFillStyle(3002);
  //data1_sfON_E->SetFillColor(kBlue);
  data1_sfON_E->SetLineWidth(3);
  data1_sfON_E->SetMinimum(0.);
  sim1_sfON_E->SetMarkerColor(kRed);
  sim1_sfON_E->SetLineColor(kRed);
  sim1_sfON_E->SetMarkerSize(0.75);
  sim1_sfON_E->SetMarkerStyle(20);
  data1_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  data1_sfON_E->SetMaximum(data1_sfON_E->GetMaximum()*1.2);
  data1_sfON_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim1_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  data1_sfON_E->Draw("HISTE0");
  sim1_sfON_E->Draw("SAMEE0");
  

  c2_sfON_E->cd(3);

  
  data2_sfON_E->SetMarkerColor(kBlue);
  data2_sfON_E->SetMarkerStyle(22);
  data2_sfON_E->SetMarkerSize(0.75);
  data2_sfON_E->SetLineColor(kBlue);
  //data2_sfON_E->SetFillStyle(3002);
  //data2_sfON_E->SetFillColor(kBlue);
  data2_sfON_E->SetLineWidth(3);
  data2_sfON_E->SetMinimum(0.);
  sim2_sfON_E->SetMarkerColor(kRed);
  sim2_sfON_E->SetLineColor(kRed);
  sim2_sfON_E->SetMarkerSize(0.75);
  sim2_sfON_E->SetMarkerStyle(20);
  data2_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  data2_sfON_E->SetMaximum(data2_sfON_E->GetMaximum()*1.2);
  data2_sfON_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim2_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  data2_sfON_E->Draw("HISTE0");
  sim2_sfON_E->Draw("SAMEE0");


  c2_sfON_E->cd(4);

  
  data3_sfON_E->SetMarkerColor(kBlue);
  data3_sfON_E->SetMarkerStyle(22);
  data3_sfON_E->SetMarkerSize(0.75);
  data3_sfON_E->SetLineColor(kBlue);
  //data3_sfON_E->SetFillStyle(3002);
  //data3_sfON_E->SetFillColor(kBlue);
  data3_sfON_E->SetLineWidth(3);
  data3_sfON_E->SetMinimum(0.);
  sim3_sfON_E->SetMarkerColor(kRed);
  sim3_sfON_E->SetLineColor(kRed);
  sim3_sfON_E->SetMarkerSize(0.75);
  sim3_sfON_E->SetMarkerStyle(20);
  data3_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  data3_sfON_E->SetMaximum(data3_sfON_E->GetMaximum()*1.2);
  data3_sfON_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim3_sfON_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  data3_sfON_E->Draw("HISTE0");
  sim3_sfON_E->Draw("SAMEE0");


  pdfFile = TString::Format("spectraComp_%i-%i_Type%s_%0.0f-%0.0f_sfON_E.pdf",octetStart,octetEnd,normType.Data(), normLow, normHigh);
  //TString pdfFileStart = pdfFile + TString("(");
  //TString pdfFileEnd = pdfFile + TString
  c1_sfON_E->Print(TString::Format("%s(",pdfFile.Data()));
  c2_sfON_E->Print(TString::Format("%s)",pdfFile.Data()));


  ////////////////////// flipper OFF West ////////////////////////////////////////////

  //Residuals...
  TCanvas *c1_sfOFF_W = new TCanvas("c1_sfOFF_W","c1_sfOFF_W", 1600., 600.);
  c1_sfOFF_W->Divide(2,1);
  c1_sfOFF_W->cd(2);
  nBins = dataALL_sfOFF_W->GetNbinsX();
  Min = dataALL_sfOFF_W->GetXaxis()->GetBinLowEdge(dataALL->GetXaxis()->GetFirst());
  Max = dataALL_sfOFF_W->GetXaxis()->GetBinUpEdge(dataALL->GetXaxis()->GetLast());

  //residual
  TH1D *resid_sfOFF_W = new TH1D("resid_sfOFF_W","Residuals: MC-Data", nBins, Min, Max);
  resid_sfOFF_W->Add(simALL_sfOFF_W,dataALL_sfOFF_W,1,-1);
  resid_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  resid_sfOFF_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  //resid_sfOFF_W->SetMaximum(2.);
  //resid_sfOFF_W->SetMinimum(-2.);
  resid_sfOFF_W->SetLineWidth(2);
  resid_sfOFF_W->SetMarkerColor(kBlue);
  resid_sfOFF_W->SetMarkerStyle(34);
  resid_sfOFF_W->SetMarkerSize(1.);

  TH1D *perc_resid_sfOFF_W = new TH1D("perc_resid_sfOFF_W","East Fractional Residuals: (MC-Data)/MC", nBins, Min, Max);
  perc_resid_sfOFF_W->Add(simALL_sfOFF_W,dataALL_sfOFF_W,1,-1);
  perc_resid_sfOFF_W->Divide(perc_resid,simALL_sfOFF_W,1.,1.);
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
  
  dataALL_sfOFF_W->SetTitle("West Spin Flipper OFF: All Event Types");
  dataALL_sfOFF_W->SetMarkerColor(kBlue);
  dataALL_sfOFF_W->SetMarkerStyle(22);
  dataALL_sfOFF_W->SetMarkerSize(0.75);
  dataALL_sfOFF_W->SetLineColor(kBlue);
  //dataALL_sfOFF_W->SetFillStyle(3002);
  //dataALL_sfOFF_W->SetFillColor(kBlue);
  dataALL_sfOFF_W->SetLineWidth(3);
  simALL_sfOFF_W->SetMarkerColor(kRed);
  simALL_sfOFF_W->SetLineColor(kRed);
  simALL_sfOFF_W->SetMarkerSize(0.75);
  simALL_sfOFF_W->SetMarkerStyle(20);
  dataALL_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  dataALL_sfOFF_W->SetMaximum(dataALL_sfOFF_W->GetMaximum()*1.2);
  dataALL_sfOFF_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  simALL_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  dataALL_sfOFF_W->Draw("HISTE0");
  simALL_sfOFF_W->Draw("SAMEE0");


  TCanvas *c2_sfOFF_W = new TCanvas("c2_sfOFF_W", "c2_sfOFF_W", 1600., 600.);
  c2_sfOFF_W->Divide(4,1);
  c2_sfOFF_W->cd(1);
  
  

  data0_sfOFF_W->SetMarkerColor(kBlue);
  data0_sfOFF_W->SetMarkerStyle(22);
  data0_sfOFF_W->SetMarkerSize(0.75);
  data0_sfOFF_W->SetLineColor(kBlue);
  //data0_sfOFF_W->SetFillStyle(3002);
  //data0_sfOFF_W->SetFillColor(kBlue);
  data0_sfOFF_W->SetLineWidth(3);
  sim0_sfOFF_W->SetMarkerColor(kRed);
  sim0_sfOFF_W->SetLineColor(kRed);
  sim0_sfOFF_W->SetMarkerSize(0.75);
  sim0_sfOFF_W->SetMarkerStyle(20);
  data0_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  data0_sfOFF_W->SetMaximum(data0_sfOFF_W->GetMaximum()*1.2);
  data0_sfOFF_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim0_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  data0_sfOFF_W->Draw("HISTE0");
  sim0_sfOFF_W->Draw("SAMEE0");
  
  c2_sfOFF_W->cd(2);
  

  data1_sfOFF_W->SetMarkerColor(kBlue);
  data1_sfOFF_W->SetMarkerStyle(22);
  data1_sfOFF_W->SetMarkerSize(0.75);
  data1_sfOFF_W->SetLineColor(kBlue);
  //data1_sfOFF_W->SetFillStyle(3002);
  //data1_sfOFF_W->SetFillColor(kBlue);
  data1_sfOFF_W->SetLineWidth(3);
  data1_sfOFF_W->SetMinimum(0.);
  sim1_sfOFF_W->SetMarkerColor(kRed);
  sim1_sfOFF_W->SetLineColor(kRed);
  sim1_sfOFF_W->SetMarkerSize(0.75);
  sim1_sfOFF_W->SetMarkerStyle(20);
  data1_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  data1_sfOFF_W->SetMaximum(data1_sfOFF_W->GetMaximum()*1.2);
  data1_sfOFF_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim1_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  data1_sfOFF_W->Draw("HISTE0");
  sim1_sfOFF_W->Draw("SAMEE0");
  

  c2_sfOFF_W->cd(3);

  
  data2_sfOFF_W->SetMarkerColor(kBlue);
  data2_sfOFF_W->SetMarkerStyle(22);
  data2_sfOFF_W->SetMarkerSize(0.75);
  data2_sfOFF_W->SetLineColor(kBlue);
  //data2_sfOFF_W->SetFillStyle(3002);
  //data2_sfOFF_W->SetFillColor(kBlue);
  data2_sfOFF_W->SetLineWidth(3);
  data2_sfOFF_W->SetMinimum(0.);
  sim2_sfOFF_W->SetMarkerColor(kRed);
  sim2_sfOFF_W->SetLineColor(kRed);
  sim2_sfOFF_W->SetMarkerSize(0.75);
  sim2_sfOFF_W->SetMarkerStyle(20);
  data2_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  data2_sfOFF_W->SetMaximum(data2_sfOFF_W->GetMaximum()*1.2);
  data2_sfOFF_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim2_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  data2_sfOFF_W->Draw("HISTE0");
  sim2_sfOFF_W->Draw("SAMEE0");


  c2_sfOFF_W->cd(4);

  
  data3_sfOFF_W->SetMarkerColor(kBlue);
  data3_sfOFF_W->SetMarkerStyle(22);
  data3_sfOFF_W->SetMarkerSize(0.75);
  data3_sfOFF_W->SetLineColor(kBlue);
  //data3_sfOFF_W->SetFillStyle(3002);
  //data3_sfOFF_W->SetFillColor(kBlue);
  data3_sfOFF_W->SetLineWidth(3);
  data3_sfOFF_W->SetMinimum(0.);
  sim3_sfOFF_W->SetMarkerColor(kRed);
  sim3_sfOFF_W->SetLineColor(kRed);
  sim3_sfOFF_W->SetMarkerSize(0.75);
  sim3_sfOFF_W->SetMarkerStyle(20);
  data3_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  data3_sfOFF_W->SetMaximum(data3_sfOFF_W->GetMaximum()*1.2);
  data3_sfOFF_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim3_sfOFF_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  data3_sfOFF_W->Draw("HISTE0");
  sim3_sfOFF_W->Draw("SAMEE0");


  pdfFile = TString::Format("spectraComp_%i-%i_Type%s_%0.0f-%0.0f_sfOFF_W.pdf",octetStart,octetEnd,normType.Data(), normLow, normHigh);
  //TString pdfFileStart = pdfFile + TString("(");
  //TString pdfFileEnd = pdfFile + TString
  c1_sfOFF_W->Print(TString::Format("%s(",pdfFile.Data()));
  c2_sfOFF_W->Print(TString::Format("%s)",pdfFile.Data()));


  ////////////////////// flipper ON West ////////////////////////////////////////////

  //Residuals...
  TCanvas *c1_sfON_W = new TCanvas("c1_sfON_W","c1_sfON_W", 1600., 600.);
  c1_sfON_W->Divide(2,1);
  c1_sfON_W->cd(2);
  nBins = dataALL_sfON_W->GetNbinsX();
  Min = dataALL_sfON_W->GetXaxis()->GetBinLowEdge(dataALL->GetXaxis()->GetFirst());
  Max = dataALL_sfON_W->GetXaxis()->GetBinUpEdge(dataALL->GetXaxis()->GetLast());

  //residual
  TH1D *resid_sfON_W = new TH1D("resid_sfON_W","Residuals: MC-Data", nBins, Min, Max);
  resid_sfON_W->Add(simALL_sfON_W,dataALL_sfON_W,1,-1);
  resid_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  resid_sfON_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  //resid_sfON_W->SetMaximum(2.);
  //resid_sfON_W->SetMinimum(-2.);
  resid_sfON_W->SetLineWidth(2);
  resid_sfON_W->SetMarkerColor(kBlue);
  resid_sfON_W->SetMarkerStyle(34);
  resid_sfON_W->SetMarkerSize(1.);

  TH1D *perc_resid_sfON_W = new TH1D("perc_resid_sfON_W","East Fractional Residuals: (MC-Data)/MC", nBins, Min, Max);
  perc_resid_sfON_W->Add(simALL_sfON_W,dataALL_sfON_W,1,-1);
  perc_resid_sfON_W->Divide(perc_resid,simALL_sfON_W,1.,1.);
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
  
  dataALL_sfON_W->SetTitle("West Spin Flipper ON: All Event Types");
  dataALL_sfON_W->SetMarkerColor(kBlue);
  dataALL_sfON_W->SetMarkerStyle(22);
  dataALL_sfON_W->SetMarkerSize(0.75);
  dataALL_sfON_W->SetLineColor(kBlue);
  //dataALL_sfON_W->SetFillStyle(3002);
  //dataALL_sfON_W->SetFillColor(kBlue);
  dataALL_sfON_W->SetLineWidth(3);
  simALL_sfON_W->SetMarkerColor(kRed);
  simALL_sfON_W->SetLineColor(kRed);
  simALL_sfON_W->SetMarkerSize(0.75);
  simALL_sfON_W->SetMarkerStyle(20);
  dataALL_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  dataALL_sfON_W->SetMaximum(dataALL_sfON_W->GetMaximum()*1.2);
  dataALL_sfON_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  simALL_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  dataALL_sfON_W->Draw("HISTE0");
  simALL_sfON_W->Draw("SAMEE0");


  TCanvas *c2_sfON_W = new TCanvas("c2_sfON_W", "c2_sfON_W", 1600., 600.);
  c2_sfON_W->Divide(4,1);
  c2_sfON_W->cd(1);
  
  

  data0_sfON_W->SetMarkerColor(kBlue);
  data0_sfON_W->SetMarkerStyle(22);
  data0_sfON_W->SetMarkerSize(0.75);
  data0_sfON_W->SetLineColor(kBlue);
  //data0_sfON_W->SetFillStyle(3002);
  //data0_sfON_W->SetFillColor(kBlue);
  data0_sfON_W->SetLineWidth(3);
  sim0_sfON_W->SetMarkerColor(kRed);
  sim0_sfON_W->SetLineColor(kRed);
  sim0_sfON_W->SetMarkerSize(0.75);
  sim0_sfON_W->SetMarkerStyle(20);
  data0_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  data0_sfON_W->SetMaximum(data0_sfON_W->GetMaximum()*1.2);
  data0_sfON_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim0_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  data0_sfON_W->Draw("HISTE0");
  sim0_sfON_W->Draw("SAMEE0");
  
  c2_sfON_W->cd(2);
  

  data1_sfON_W->SetMarkerColor(kBlue);
  data1_sfON_W->SetMarkerStyle(22);
  data1_sfON_W->SetMarkerSize(0.75);
  data1_sfON_W->SetLineColor(kBlue);
  //data1_sfON_W->SetFillStyle(3002);
  //data1_sfON_W->SetFillColor(kBlue);
  data1_sfON_W->SetLineWidth(3);
  data1_sfON_W->SetMinimum(0.);
  sim1_sfON_W->SetMarkerColor(kRed);
  sim1_sfON_W->SetLineColor(kRed);
  sim1_sfON_W->SetMarkerSize(0.75);
  sim1_sfON_W->SetMarkerStyle(20);
  data1_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  data1_sfON_W->SetMaximum(data1_sfON_W->GetMaximum()*1.2);
  data1_sfON_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim1_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  data1_sfON_W->Draw("HISTE0");
  sim1_sfON_W->Draw("SAMEE0");
  

  c2_sfON_W->cd(3);

  
  data2_sfON_W->SetMarkerColor(kBlue);
  data2_sfON_W->SetMarkerStyle(22);
  data2_sfON_W->SetMarkerSize(0.75);
  data2_sfON_W->SetLineColor(kBlue);
  //data2_sfON_W->SetFillStyle(3002);
  //data2_sfON_W->SetFillColor(kBlue);
  data2_sfON_W->SetLineWidth(3);
  data2_sfON_W->SetMinimum(0.);
  sim2_sfON_W->SetMarkerColor(kRed);
  sim2_sfON_W->SetLineColor(kRed);
  sim2_sfON_W->SetMarkerSize(0.75);
  sim2_sfON_W->SetMarkerStyle(20);
  data2_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  data2_sfON_W->SetMaximum(data2_sfON_W->GetMaximum()*1.2);
  data2_sfON_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim2_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  data2_sfON_W->Draw("HISTE0");
  sim2_sfON_W->Draw("SAMEE0");


  c2_sfON_W->cd(4);


  data3_sfON_W->SetMarkerColor(kBlue);
  data3_sfON_W->SetMarkerStyle(22);
  data3_sfON_W->SetMarkerSize(0.75);
  data3_sfON_W->SetLineColor(kBlue);
  //data3_sfON_W->SetFillStyle(3002);
  //data3_sfON_W->SetFillColor(kBlue);
  data3_sfON_W->SetLineWidth(3);
  data3_sfON_W->SetMinimum(0.);
  sim3_sfON_W->SetMarkerColor(kRed);
  sim3_sfON_W->SetLineColor(kRed);
  sim3_sfON_W->SetMarkerSize(0.75);
  sim3_sfON_W->SetMarkerStyle(20);
  data3_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  data3_sfON_W->SetMaximum(data3_sfON_W->GetMaximum()*1.2);
  data3_sfON_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim3_sfON_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  data3_sfON_W->Draw("HISTE0");
  sim3_sfON_W->Draw("SAMEE0");


  pdfFile = TString::Format("spectraComp_%i-%i_Type%s_%0.0f-%0.0f_sfON_W.pdf",octetStart,octetEnd,normType.Data(), normLow, normHigh);
  //TString pdfFileStart = pdfFile + TString("(");
  //TString pdfFileEnd = pdfFile + TString
  c1_sfON_W->Print(TString::Format("%s(",pdfFile.Data()));
  c2_sfON_W->Print(TString::Format("%s)",pdfFile.Data()));
}
  
