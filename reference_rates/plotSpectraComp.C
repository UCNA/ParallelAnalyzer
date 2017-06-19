

{
  gStyle->SetOptStat(0);

  bool withEvis = true;

  int octetStart=60;
  int octetEnd=121;

  TString normType = "0";

  Double_t normLow = 0.;
  Double_t normHigh = 200.;
  
  Double_t xAxisMax = 1200.;

  //Storing event fractions for data, E/W and sfON/OFF
  Double_t t0E_sfON, t0E_sfOFF, t1E_sfON, t1E_sfOFF, t23E_sfON, t23E_sfOFF, t0W_sfON, t0W_sfOFF, t1W_sfON, t1W_sfOFF, t23W_sfON, t23W_sfOFF;
  t0E_sfON = t0E_sfOFF = t1E_sfON = t1E_sfOFF = t23E_sfON = t23E_sfOFF = t0W_sfON = t0W_sfOFF = t1W_sfON = t1W_sfOFF = t23W_sfON = t23W_sfOFF = 0.;
  
  TFile *fg_file = new TFile(TString::Format("ForegroundSpectra_%i-%i%s.root",octetStart,octetEnd,(withEvis?"_Evis":"")),"READ");
  TFile *bg_file = new TFile(TString::Format("BackgroundSpectra_%i-%i%s.root",octetStart,octetEnd,(withEvis?"_Evis":"")),"READ");
  TFile *sim_file = new TFile(TString::Format("SIM_ForegroundSpectra_%i-%i%s.root",octetStart,octetEnd,(withEvis?"_Evis":"")),"READ");


  //Load the histograms

  TH1D *fg0E_sfOFF = (TH1D*)fg_file->Get("Type0E_sfOFF");
  TH1D *fg1E_sfOFF = (TH1D*)fg_file->Get("Type1E_sfOFF");
  TH1D *fg23E_sfOFF = (TH1D*)fg_file->Get("Type23E_sfOFF");

  TH1D *fg0W_sfOFF = (TH1D*)fg_file->Get("Type0W_sfOFF");
  TH1D *fg1W_sfOFF = (TH1D*)fg_file->Get("Type1W_sfOFF");
  TH1D *fg23W_sfOFF = (TH1D*)fg_file->Get("Type23W_sfOFF");

  TH1D *bg0E_sfOFF = (TH1D*)bg_file->Get("Type0E_sfOFF");
  TH1D *bg1E_sfOFF = (TH1D*)bg_file->Get("Type1E_sfOFF");
  TH1D *bg23E_sfOFF = (TH1D*)bg_file->Get("Type23E_sfOFF");

  TH1D *bg0W_sfOFF = (TH1D*)bg_file->Get("Type0W_sfOFF");
  TH1D *bg1W_sfOFF = (TH1D*)bg_file->Get("Type1W_sfOFF");
  TH1D *bg23W_sfOFF = (TH1D*)bg_file->Get("Type23W_sfOFF");

  TH1D *sim0E_sfOFF = (TH1D*)sim_file->Get("Type0E_sfOFF");
  TH1D *sim1E_sfOFF = (TH1D*)sim_file->Get("Type1E_sfOFF");
  TH1D *sim23E_sfOFF = (TH1D*)sim_file->Get("Type23E_sfOFF");

  TH1D *sim0W_sfOFF = (TH1D*)sim_file->Get("Type0W_sfOFF");
  TH1D *sim1W_sfOFF = (TH1D*)sim_file->Get("Type1W_sfOFF");
  TH1D *sim23W_sfOFF = (TH1D*)sim_file->Get("Type23W_sfOFF");
  
  TH1D *fg0E_sfON = (TH1D*)fg_file->Get("Type0E_sfON");
  TH1D *fg1E_sfON = (TH1D*)fg_file->Get("Type1E_sfON");
  TH1D *fg23E_sfON = (TH1D*)fg_file->Get("Type23E_sfON");

  TH1D *fg0W_sfON = (TH1D*)fg_file->Get("Type0W_sfON");
  TH1D *fg1W_sfON = (TH1D*)fg_file->Get("Type1W_sfON");
  TH1D *fg23W_sfON = (TH1D*)fg_file->Get("Type23W_sfON");

  TH1D *bg0E_sfON = (TH1D*)bg_file->Get("Type0E_sfON");
  TH1D *bg1E_sfON = (TH1D*)bg_file->Get("Type1E_sfON");
  TH1D *bg23E_sfON = (TH1D*)bg_file->Get("Type23E_sfON");

  TH1D *bg0W_sfON = (TH1D*)bg_file->Get("Type0W_sfON");
  TH1D *bg1W_sfON = (TH1D*)bg_file->Get("Type1W_sfON");
  TH1D *bg23W_sfON = (TH1D*)bg_file->Get("Type23W_sfON");

  TH1D *sim0E_sfON = (TH1D*)sim_file->Get("Type0E_sfON");
  TH1D *sim1E_sfON = (TH1D*)sim_file->Get("Type1E_sfON");
  TH1D *sim23E_sfON = (TH1D*)sim_file->Get("Type23E_sfON");

  TH1D *sim0W_sfON = (TH1D*)sim_file->Get("Type0W_sfON");
  TH1D *sim1W_sfON = (TH1D*)sim_file->Get("Type1W_sfON");
  TH1D *sim23W_sfON = (TH1D*)sim_file->Get("Type23W_sfON");


  //subtract bg from fg
  fg0E_sfON->Add(bg0E_sfON,-1.);
  fg0W_sfON->Add(bg0W_sfON,-1.);
  fg1E_sfON->Add(bg1E_sfON,-1.);
  fg1W_sfON->Add(bg1W_sfON,-1.);
  fg23E_sfON->Add(bg23E_sfON,-1.);
  fg23W_sfON->Add(bg23W_sfON,-1.);
  fg0E_sfOFF->Add(bg0E_sfOFF,-1.);
  fg0W_sfOFF->Add(bg0W_sfOFF,-1.);
  fg1E_sfOFF->Add(bg1E_sfOFF,-1.);
  fg1W_sfOFF->Add(bg1W_sfOFF,-1.);
  fg23E_sfOFF->Add(bg23E_sfOFF,-1.);
  fg23W_sfOFF->Add(bg23W_sfOFF,-1.);

  t0E_sfOFF = fg0E_sfOFF->Integral(18,79);
  t0E_sfON = fg0E_sfON->Integral(18,79);
  t1E_sfOFF = fg1E_sfOFF->Integral(18,79);
  t1E_sfON = fg1E_sfON->Integral(18,79);
  t23E_sfOFF = fg23E_sfOFF->Integral(18,79);
  t23E_sfON = fg23E_sfON->Integral(18,79);
  t0W_sfOFF = fg0W_sfOFF->Integral(18,79);
  t0W_sfON = fg0W_sfON->Integral(18,79);
  t1W_sfOFF = fg1W_sfOFF->Integral(18,79);
  t1W_sfON = fg1W_sfON->Integral(18,79);
  t23W_sfOFF = fg23W_sfOFF->Integral(18,79);
  t23W_sfON = fg23W_sfON->Integral(18,79);

  double totalE_sfOFF = t0E_sfOFF + t1E_sfOFF + t23E_sfOFF;
  double totalE_sfON = t0E_sfON + t1E_sfON + t23E_sfON;
  double totalW_sfOFF = t0W_sfOFF + t1W_sfOFF + t23W_sfOFF;
  double totalW_sfON = t0W_sfON + t1W_sfON + t23W_sfON;

  
  std::cout << "********************************************\n"
	    << "          Event ratios (180-780 keV)\n\n"
	    << "\t\tsfON\t\t\tsfOFF\n"
	    << "Type 0E:\t" << t0E_sfON/totalE_sfON << "\t\t" << t0E_sfOFF/totalE_sfOFF << "\n"
	    << "Type 0W:\t" << t0W_sfON/totalW_sfON << "\t\t" << t0W_sfOFF/totalW_sfOFF << "\n\n"
	    << "Type 1E:\t" << t1E_sfON/totalE_sfON << "\t\t" << t1E_sfOFF/totalE_sfOFF << "\n"
	    << "Type 1W:\t" << t1W_sfON/totalW_sfON << "\t\t" << t1W_sfOFF/totalW_sfOFF << "\n\n"
	    << "Type 23E:\t" << t23E_sfON/totalE_sfON << "\t\t" << t23E_sfOFF/totalE_sfOFF << "\n"
	    << "Type 23W:\t" << t23W_sfON/totalW_sfON << "\t\t" << t23W_sfOFF/totalW_sfOFF << "\n\n"
	    << "********************************************\n";
    
  //make copies of sfON and then add sfOFF
  
  TH1D *uk0_E = (TH1D*)fg0E_sfON->Clone("uk0_E");
  TH1D *sim0_E = (TH1D*)sim0E_sfON->Clone("sim0_E");
  TH1D *uk1_E = (TH1D*)fg1E_sfON->Clone("uk1_E");
  TH1D *sim1_E = (TH1D*)sim1E_sfON->Clone("sim1_E");
  TH1D *uk23_E = (TH1D*)fg23E_sfON->Clone("uk23_E");
  TH1D *sim23_E = (TH1D*)sim23E_sfON->Clone("sim23_E");
  TH1D *uk0_W = (TH1D*)fg0W_sfON->Clone("uk0_W");
  TH1D *sim0_W = (TH1D*)sim0W_sfON->Clone("sim0_W");
  TH1D *uk1_W = (TH1D*)fg1W_sfON->Clone("uk1_W");
  TH1D *sim1_W = (TH1D*)sim1W_sfON->Clone("sim1_W");
  TH1D *uk23_W = (TH1D*)fg23W_sfON->Clone("uk23_W");
  TH1D *sim23_W = (TH1D*)sim23W_sfON->Clone("sim23_W");

  uk0_E->Add(fg0E_sfOFF,1.);
  uk0_W->Add(fg0W_sfOFF,1.);
  uk1_E->Add(fg1E_sfOFF,1.);
  uk1_W->Add(fg1W_sfOFF,1.);
  uk23_E->Add(fg23E_sfOFF,1.);
  uk23_W->Add(fg23W_sfOFF,1.);
  
  sim0_E->Add(sim0E_sfOFF,1.);
  sim0_W->Add(sim0W_sfOFF,1.);
  sim1_E->Add(sim1E_sfOFF,1.);
  sim1_W->Add(sim1W_sfOFF,1.);
  sim23_E->Add(sim23E_sfOFF,1.);
  sim23_W->Add(sim23W_sfOFF,1.);

  TH1D *ukALL_E = (TH1D*)uk0_E->Clone("ukALL_E");
  TH1D *simALL_E = (TH1D*)sim0_E->Clone("simALL_E");
  TH1D *ukALL_W = (TH1D*)uk0_W->Clone("ukALL_W");
  TH1D *simALL_W = (TH1D*)sim0_W->Clone("simALL_W");
  
  ukALL_E->Add(uk1_E,1.); ukALL_E->Add(uk23_E,1.);
  ukALL_W->Add(uk1_W,1.); ukALL_W->Add(uk23_W,1.);
  simALL_E->Add(sim1_E,1.); simALL_E->Add(sim23_E,1.);
  simALL_W->Add(sim1_W,1.); simALL_W->Add(sim23_W,1.);
    
  //Normalize
  Double_t normFactorE = 1.;
  Double_t normFactorW = 1.;
  
  if (normType==TString("ALL")) {
    normFactorE = ukALL_E->Integral(ukALL_E->GetXaxis()->FindFixBin(normLow),ukALL_E->GetXaxis()->FindFixBin(normHigh))/simALL_E->Integral(simALL_E->GetXaxis()->FindFixBin(normLow),simALL_E->GetXaxis()->FindFixBin(normHigh));
    normFactorW = ukALL_W->Integral(ukALL_W->GetXaxis()->FindFixBin(normLow),ukALL_W->GetXaxis()->FindFixBin(normHigh))/simALL_W->Integral(simALL_W->GetXaxis()->FindFixBin(normLow),simALL_W->GetXaxis()->FindFixBin(normHigh));
  }

  if (normType==TString("0")) {
    normFactorE = uk0_E->Integral(uk0_E->GetXaxis()->FindFixBin(normLow),uk0_E->GetXaxis()->FindFixBin(normHigh))/sim0_E->Integral(sim0_E->GetXaxis()->FindFixBin(normLow),sim0_E->GetXaxis()->FindFixBin(normHigh));
    normFactorW = uk0_W->Integral(uk0_W->GetXaxis()->FindFixBin(normLow),uk0_W->GetXaxis()->FindFixBin(normHigh))/sim0_W->Integral(sim0_W->GetXaxis()->FindFixBin(normLow),sim0_W->GetXaxis()->FindFixBin(normHigh));
  }
  
  if (normType==TString("1")) {
    normFactorE = uk1_E->Integral(uk1_E->GetXaxis()->FindFixBin(normLow),uk1_E->GetXaxis()->FindFixBin(normHigh))/sim1_E->Integral(sim1_E->GetXaxis()->FindFixBin(normLow),sim1_E->GetXaxis()->FindFixBin(normHigh));
    normFactorW = uk1_W->Integral(uk1_W->GetXaxis()->FindFixBin(normLow),uk1_W->GetXaxis()->FindFixBin(normHigh))/sim1_W->Integral(sim1_W->GetXaxis()->FindFixBin(normLow),sim1_W->GetXaxis()->FindFixBin(normHigh));
  }
  

  if (normType==TString("23")) {
    normFactorE = uk23_E->Integral(uk23_E->GetXaxis()->FindFixBin(normLow),uk23_E->GetXaxis()->FindFixBin(normHigh))/sim23_E->Integral(sim23_E->GetXaxis()->FindFixBin(normLow),sim23_E->GetXaxis()->FindFixBin(normHigh));
    normFactorW = uk23_W->Integral(uk23_W->GetXaxis()->FindFixBin(normLow),uk23_W->GetXaxis()->FindFixBin(normHigh))/sim23_W->Integral(sim23_W->GetXaxis()->FindFixBin(normLow),sim23_W->GetXaxis()->FindFixBin(normHigh));
  }
  
  
  simALL_E->Scale(normFactorE);
  simALL_W->Scale(normFactorW);

  //Residuals...
  TCanvas *c1E = new TCanvas("c1E","c1E", 1600., 600.);
  c1E->Divide(2,1);
  c1E->cd(2);
  Int_t nBins = ukALL_E->GetNbinsX();
  Double_t Min = ukALL_E->GetXaxis()->GetBinLowEdge(ukALL_E->GetXaxis()->GetFirst());
  Double_t Max = ukALL_E->GetXaxis()->GetBinUpEdge(ukALL_E->GetXaxis()->GetLast());

  //residual
  TH1D *residE = new TH1D("residE","East Residuals: MC-Data", nBins, Min, Max);
  residE->Add(simALL_E,ukALL_E,1,-1);
  residE->GetXaxis()->SetRangeUser(0., xAxisMax);
  residE->GetYaxis()->SetTitle("event rate (mHz/keV)");
  //residE->SetMaximum(2.);
  //residE->SetMinimum(-2.);
  residE->SetLineWidth(2);
  residE->SetMarkerColor(kBlue);
  residE->SetMarkerStyle(34);
  residE->SetMarkerSize(1.);

  TH1D *perc_residE = new TH1D("perc_residE","East Fractional Residuals: (MC-Data)/MC", nBins, Min, Max);
  perc_residE->Add(simALL_E,ukALL_E,1,-1);
  perc_residE->Divide(perc_residE,simALL_E,1.,1.);
  perc_residE->GetXaxis()->SetRangeUser(0., xAxisMax);
  //perc_resid->GetYaxis()->SetTitle("event rate (mHz/keV)");
  perc_residE->SetMaximum(.1);
  perc_residE->SetMinimum(-.1);
  perc_residE->SetLineWidth(2);
  


  //perc_residE->Draw();
  residE->Draw("P0");
  c1E->Update();
  TLine *l = new TLine(Min, 0., xAxisMax, 0.);
  l->SetLineStyle(8);
  l->Draw();

  c1E->cd(1);
  
  ukALL_E->SetTitle("East: All Event Types");
  ukALL_E->SetMarkerColor(kBlue);
  ukALL_E->SetMarkerStyle(22);
  ukALL_E->SetMarkerSize(0.75);
  ukALL_E->SetLineColor(kBlue);
  //ukALL_E->SetFillStyle(3002);
  //ukALL_E->SetFillColor(kBlue);
  ukALL_E->SetLineWidth(3);
  simALL_E->SetMarkerColor(kRed);
  simALL_E->SetLineColor(kRed);
  simALL_E->SetMarkerSize(0.75);
  simALL_E->SetMarkerStyle(20);
  ukALL_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  ukALL_E->SetMaximum(ukALL_E->GetMaximum()*1.2);
  ukALL_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  simALL_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  ukALL_E->Draw("HIST");
  simALL_E->Draw("SAMEP0");


  TCanvas *c2E = new TCanvas("c2E", "c2E", 1600., 600.);
  c2E->Divide(3,1);
  c2E->cd(1);
  
  

  uk0_E->SetMarkerColor(kBlue);
  uk0_E->SetMarkerStyle(22);
  uk0_E->SetMarkerSize(0.75);
  uk0_E->SetLineColor(kBlue);
  //uk0_E->SetFillStyle(3002);
  //uk0_E->SetFillColor(kBlue);
  uk0_E->SetLineWidth(3);
  sim0_E->SetMarkerColor(kRed);
  sim0_E->SetLineColor(kRed);
  sim0_E->SetMarkerSize(0.75);
  sim0_E->SetMarkerStyle(20);
  uk0_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  uk0_E->SetMaximum(uk0_E->GetMaximum()*1.2);
  uk0_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim0_E->Scale(normFactorE);
  sim0_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  uk0_E->Draw("HIST");
  sim0_E->Draw("SAMEP0");
  
  c2E->cd(2);
  

  uk1_E->SetMarkerColor(kBlue);
  uk1_E->SetMarkerStyle(22);
  uk1_E->SetMarkerSize(0.75);
  uk1_E->SetLineColor(kBlue);
  //uk1_E->SetFillStyle(3002);
  //uk1_E->SetFillColor(kBlue);
  uk1_E->SetLineWidth(3);
  uk1_E->SetMinimum(0.);
  sim1_E->SetMarkerColor(kRed);
  sim1_E->SetLineColor(kRed);
  sim1_E->SetMarkerSize(0.75);
  sim1_E->SetMarkerStyle(20);
  uk1_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  uk1_E->SetMaximum(uk1_E->GetMaximum()*1.2);
  uk1_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim1_E->Scale(normFactorE);
  sim1_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  uk1_E->Draw("HIST");
  sim1_E->Draw("SAMEP0");
  

  c2E->cd(3);

  
  uk23_E->SetMarkerColor(kBlue);
  uk23_E->SetMarkerStyle(22);
  uk23_E->SetMarkerSize(0.75);
  uk23_E->SetLineColor(kBlue);
  //uk23_E->SetFillStyle(3002);
  //uk23_E->SetFillColor(kBlue);
  uk23_E->SetLineWidth(3);
  uk23_E->SetMinimum(0.);
  sim23_E->SetMarkerColor(kRed);
  sim23_E->SetLineColor(kRed);
  sim23_E->SetMarkerSize(0.75);
  sim23_E->SetMarkerStyle(20);
  uk23_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  uk23_E->SetMaximum(uk23_E->GetMaximum()*1.2);
  uk23_E->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim23_E->Scale(normFactorE);
  sim23_E->GetXaxis()->SetRangeUser(0., xAxisMax);
  uk23_E->Draw("HIST");
  sim23_E->Draw("SAMEP0");


  //West
  TCanvas *c1W = new TCanvas("c1W","c1W", 1600., 600.);
  c1W->Divide(2,1);
  c1W->cd(2);
  nBins = ukALL_W->GetNbinsX();
  Min = ukALL_W->GetXaxis()->GetBinLowEdge(ukALL_W->GetXaxis()->GetFirst());
  Max = ukALL_W->GetXaxis()->GetBinUpEdge(ukALL_W->GetXaxis()->GetLast());

  //residual
  TH1D *residW = new TH1D("residW","West Residuals: MC-Data", nBins, Min, Max);
  residW->Add(simALL_W,ukALL_W,1.,-1.);
  residW->GetXaxis()->SetRangeUser(0., xAxisMax);
  residW->GetYaxis()->SetTitle("event rate (mHz/keV)");
  //residW->SetMaximum(2.);
  //residW->SetMinimum(-2.);
  residW->SetLineWidth(2);
  residW->SetMarkerColor(kBlue);
  residW->SetMarkerStyle(34);
  residW->SetMarkerSize(1.);

  TH1D *perc_residW = new TH1D("perc_residW","West Fractional Residuals: (MC-Data)/MC", nBins, Min, Max);
  perc_residW->Add(simALL_W,ukALL_W,1.,-1.);
  perc_residW->Divide(perc_residW,simALL_W,1.,1.);
  perc_residW->GetXaxis()->SetRangeUser(0., xAxisMax);
  //perc_resid->GetYaxis()->SetTitle("event rate (mHz/keV)");
  perc_residW->SetMaximum(.1);
  perc_residW->SetMinimum(-.1);
  perc_residW->SetLineWidth(2);


  //perc_residW->Draw();
  residW->Draw("P0");
  c1W->Update();
  TLine *l2 = new TLine(Min, 0., xAxisMax, 0.);
  l2->SetLineStyle(8);
  l2->Draw();

  c1W->cd(1);
  
  ukALL_W->SetTitle("West: All Event Types");
  ukALL_W->SetMarkerColor(kBlue);
  ukALL_W->SetMarkerStyle(22);
  ukALL_W->SetMarkerSize(0.75);
  ukALL_W->SetLineColor(kBlue);
  //ukALL_W->SetFillStyle(3002);
  //ukALL_W->SetFillColor(kBlue);
  ukALL_W->SetLineWidth(3);
  simALL_W->SetMarkerColor(kRed);
  simALL_W->SetLineColor(kRed);
  simALL_W->SetMarkerSize(0.75);
  simALL_W->SetMarkerStyle(20);
  ukALL_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  ukALL_W->SetMaximum(ukALL_W->GetMaximum()*1.2);
  ukALL_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  simALL_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  ukALL_W->Draw("HIST");
  simALL_W->Draw("SAMEP0");


  TCanvas *c2W = new TCanvas("c2W", "c2W", 1600., 600.);
  c2W->Divide(3,1);
  c2W->cd(1);
  
  

  uk0_W->SetMarkerColor(kBlue);
  uk0_W->SetMarkerStyle(22);
  uk0_W->SetMarkerSize(0.75);
  uk0_W->SetLineColor(kBlue);
  //uk0_W->SetFillStyle(3002);
  //uk0_W->SetFillColor(kBlue);
  uk0_W->SetLineWidth(3);
  sim0_W->SetMarkerColor(kRed);
  sim0_W->SetLineColor(kRed);
  sim0_W->SetMarkerSize(0.75);
  sim0_W->SetMarkerStyle(20);
  uk0_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  uk0_W->SetMaximum(uk0_W->GetMaximum()*1.2);
  uk0_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim0_W->Scale(normFactorW);
  sim0_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  uk0_W->Draw("HIST");
  sim0_W->Draw("SAMEP0");
  
  c2W->cd(2);
  

  uk1_W->SetMarkerColor(kBlue);
  uk1_W->SetMarkerStyle(22);
  uk1_W->SetMarkerSize(0.75);
  uk1_W->SetLineColor(kBlue);
  //uk1_W->SetFillStyle(3002);
  //uk1_W->SetFillColor(kBlue);
  uk1_W->SetLineWidth(3);
  uk1_W->SetMinimum(0.);
  sim1_W->SetMarkerColor(kRed);
  sim1_W->SetLineColor(kRed);
  sim1_W->SetMarkerSize(0.75);
  sim1_W->SetMarkerStyle(20);
  uk1_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  uk1_W->SetMaximum(uk1_W->GetMaximum()*1.2);
  uk1_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim1_W->Scale(normFactorW);
  sim1_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  uk1_W->Draw("HIST");
  sim1_W->Draw("SAMEP0");
  

  c2W->cd(3);

  
  uk23_W->SetMarkerColor(kBlue);
  uk23_W->SetMarkerStyle(22);
  uk23_W->SetMarkerSize(0.75);
  uk23_W->SetLineColor(kBlue);
  //uk23_W->SetFillStyle(3002);
  //uk23_W->SetFillColor(kBlue);
  uk23_W->SetLineWidth(3);
  uk23_W->SetMinimum(0.);
  sim23_W->SetMarkerColor(kRed);
  sim23_W->SetLineColor(kRed);
  sim23_W->SetMarkerSize(0.75);
  sim23_W->SetMarkerStyle(20);
  uk23_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  uk23_W->SetMaximum(uk23_W->GetMaximum()*1.2);
  uk23_W->GetYaxis()->SetTitle("event rate (mHz/keV)");
  sim23_W->Scale(normFactorW);
  sim23_W->GetXaxis()->SetRangeUser(0., xAxisMax);
  uk23_W->Draw("HIST");
  sim23_W->Draw("SAMEP0");



  TString pdfFile = TString::Format("spectraComp_%i-%i_Type%s_%0.0f-%0.0f%s.pdf",
				    octetStart,octetEnd,normType.Data(), normLow, normHigh,(withEvis?"_Evis":""));
  //TString pdfFileStart = pdfFile + TString("(");
  //TString pdfFileEnd = pdfFile + TString
  c1E->Print(TString::Format("%s(",pdfFile.Data()));
  c2E->Print(pdfFile);
  c1W->Print(pdfFile);
  c2W->Print(TString::Format("%s)",pdfFile.Data()));
  
}
