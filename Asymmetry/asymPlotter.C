{

  bool removeOct10 = true;
   
  gStyle->SetTitleSize(0.08,"t");
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleYSize(0.06);
  gStyle->SetTitleYOffset(0.65);
  gStyle->SetTitleXSize(0.06);
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetLabelSize(0.05,"xyz");
  
  gStyle->SetOptFit(1111);
  // gStyle->SetTitleX(0.25);
  gStyle->SetStatX(0.85);
  gStyle->SetStatY(0.85);
  
  std::vector <Double_t> oct_NCSU;
  std::vector <Double_t> asym_NCSU;
  std::vector <Double_t> err_NCSU;

  std::vector <Double_t> oct_UK;
  std::vector <Double_t> asym_UK;
  std::vector <Double_t> err_UK;
  
  ifstream infile("octval.dat");

  string hold[8];
  double octet, asym, err, pull;
  
  infile >> hold[0] >> hold[1] >> hold[2] >> hold[3] >> hold[4] >> hold[5] >> hold[6] >> hold[7];

  while (infile >> octet >> asym >> err >> pull) {
    if (removeOct10 && octet==10) continue;
    oct_NCSU.push_back(octet);
    asym_NCSU.push_back(asym);
    err_NCSU.push_back(err);
  }

  infile.close();

  infile.open("octvalUK.dat");

  infile >> hold[0] >> hold[1] >> hold[2] >> hold[3];
    
  while (infile >> octet >> asym >> err >> pull) {
    if (removeOct10 && octet==10) continue;
    oct_UK.push_back(octet);
    asym_UK.push_back(asym);
    err_UK.push_back(err);
  }

  std::vector <Double_t> diff;

  TCanvas *c1 = new TCanvas("c1","c1",750,600);

  c1->Divide(1,2);

  c1->cd(1);

  TGraphErrors *UK = new TGraphErrors(oct_UK.size(), &oct_UK[0], &asym_UK[0], 0, &err_UK[0]);
  UK->SetTitle("UK Asymmetry");
  UK->GetXaxis()->SetTitle("Octet Number");
  UK->GetYaxis()->SetTitle("Raw Asymmetry");
  UK->GetXaxis()->CenterTitle();
  UK->GetYaxis()->CenterTitle();
  UK->SetMarkerColor(kBlue);
  UK->SetMarkerStyle(20);
  UK->SetLineWidth(2);
  UK->SetLineColor(kBlue);
  UK->SetMinimum(0.04);
  UK->SetMaximum(0.06);
  UK->GetXaxis()->SetLimits(oct_UK[0]-2., oct_UK[oct_UK.size()-1]+2.);
  UK->Draw("AP");

  TF1 *fit = new TF1("fit","[0]",oct_UK[0]-1., oct_UK[oct_UK.size()-1]+1.);
  fit->SetLineColor(1);
  fit->SetLineWidth(3);
  fit->SetParameter(0,0.05);
  
  UK->Fit("fit","R");

  cout << "UK\t" << fit->GetParameter(0) << "\t+/-\t" << fit->GetParError(0) << std::endl;


  c1->cd(2);

  TGraphErrors *NCSU = new TGraphErrors(oct_NCSU.size(), &oct_NCSU[0], &asym_NCSU[0], 0, &err_NCSU[0]);
  NCSU->SetTitle("NCSU Asymmetry");
  NCSU->GetXaxis()->SetTitle("Octet Number");
  NCSU->GetYaxis()->SetTitle("Raw Asymmetry");
  NCSU->GetXaxis()->CenterTitle();
  NCSU->GetYaxis()->CenterTitle();
  NCSU->SetMarkerColor(kRed);
  NCSU->SetMarkerStyle(20);
  NCSU->SetLineWidth(2);
  NCSU->SetLineColor(kRed);
  NCSU->SetMinimum(0.04);
  NCSU->SetMaximum(0.06);
  NCSU->GetXaxis()->SetLimits(oct_NCSU[0]-2., oct_UK[oct_NCSU.size()-1]+2.);
  NCSU->Draw("AP");

  TF1 *fit2 = new TF1("fit2","[0]",oct_NCSU[0]-1., oct_NCSU[oct_NCSU.size()-1]+1.);
  fit2->SetLineColor(01);
  fit2->SetLineWidth(3);
  fit2->SetParameter(0,0.05);
  
  NCSU->Fit("fit2","R");

  cout << "NCSU\t" << fit2->GetParameter(0) << "\t+/-\t" << fit2->GetParError(0) << std::endl;

  c1->Update();


  ///////////////////////// PLOT Residual /////////////////////////

  TCanvas *c2 = new TCanvas("c2","c2",750,600);
 
  c2->Divide(1,2);

  c2->cd(1);

  std::vector <Double_t> res(oct_UK.size(),0.);
  std::vector <Double_t> resErr(oct_UK.size(),0.);
  
  for ( UInt_t i=0; i<oct_UK.size(); i++ ) {

    res[i] = asym_UK[i] - asym_NCSU[i];
    resErr[i] = TMath::Sqrt(err_UK[i]*err_UK[i] + err_NCSU[i]*err_NCSU[i]);

  }
  
  TGraphErrors *resid = new TGraphErrors(oct_UK.size(), &oct_UK[0], &res[0], 0, &resErr[0]);
  resid->SetTitle("Residual ( UK - NCSU )");
  resid->GetXaxis()->SetTitle("Octet Number");
  resid->GetYaxis()->SetTitle("Residual");
  resid->GetXaxis()->CenterTitle();
  resid->GetYaxis()->CenterTitle();
  resid->SetMarkerColor(6);
  resid->SetMarkerStyle(20);
  resid->SetLineWidth(2);
  resid->SetLineColor(6);
  resid->SetMinimum(-0.01);
  resid->SetMaximum(0.01);
  resid->GetXaxis()->SetLimits(oct_UK[0]-2., oct_UK[oct_UK.size()-1]+2.);
  resid->Draw("AP");

  TF1 *fit3 = new TF1("fit3","[0]",oct_UK[0]-1., oct_UK[oct_UK.size()-1]+1.);
  fit3->SetLineColor(1);
  fit3->SetLineWidth(3);
  fit3->SetParameter(0,0.05);
  
  resid->Fit("fit3","R");

  cout << "resid\t" << fit3->GetParameter(0) << "\t+/-\t" << fit3->GetParError(0) << std::endl;


  c2->cd(2);
  
  std::vector <Double_t> err_res(oct_UK.size(),0.);
  std::vector <Double_t> err_resErr(oct_UK.size(),0.);
  
  for ( UInt_t i=0; i<oct_UK.size(); i++ ) {

    //err_res[i] = ( err_UK[i]/asym_UK[i] ) / ( err_NCSU[i]/asym_NCSU[i] );
    err_res[i] = ( err_UK[i] ) / ( err_NCSU[i] );
    err_resErr[i] = 0.;//TMath::Sqrt(err_UK[i]*err_UK[i] + err_NCSU[i]*err_NCSU[i]);

  }
  
  TGraphErrors *err_resid = new TGraphErrors(oct_UK.size(), &oct_UK[0], &err_res[0], 0, &err_resErr[0]);
  err_resid->SetTitle("Residuals on Errors ( errUK / errNCSU )");
  err_resid->GetXaxis()->SetTitle("Octet Number");
  err_resid->GetYaxis()->SetTitle("Residual");
  err_resid->GetXaxis()->CenterTitle();
  err_resid->GetYaxis()->CenterTitle();
  err_resid->SetMarkerColor(3);
  err_resid->SetMarkerStyle(22);
  err_resid->SetLineWidth(2);
  err_resid->SetLineColor(6);
  //err_resid->SetMinimum(-0.00005);
  //err_resid->SetMaximum(0.00005);
  err_resid->GetXaxis()->SetLimits(oct_UK[0]-2., oct_UK[oct_UK.size()-1]+2.);
  
  err_resid->Draw("AP");

  /*TF1 *fit4 = new TF1("fit4","[0]",oct_UK[0]-1., oct_UK[oct_UK.size()-1]+1.);
  fit4->SetLineColor(1);
  fit4->SetLineWidth(3);
  fit4->SetParameter(0,0.05);
  
  
  err_resid->Fit("fit4","R");

  cout << "err_resid\t" << fit4->GetParameter(0) << "\t+/-\t" << fit4->GetParError(0) << std::endl;*/

  TLine *l = new TLine(oct_UK[0]-2.,1.,oct_UK[oct_UK.size()-1]+2.,1.);
  l->SetLineStyle(2);

  l->Draw();
  
  
  

  
  
  TString fnamePNG_both = "UK_NCSU_asymmetries.png";
  TString fnamePNG_comp = "UK_NCSU_residuals.png";
  TString fnamePDF = "UK_NCSU_asymmetryComparison.pdf";

  c1->Print(fnamePNG_both.Data());
  c1->Print(TString::Format("%s(",fnamePDF.Data()));

  c2->Print(fnamePNG_comp.Data());
  c2->Print(TString::Format("%s)",fnamePDF.Data()));

}
  
		  
