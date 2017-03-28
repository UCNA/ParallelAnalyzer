{
  #include "../include/sourcePeaks.h"

  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  // Style options
  gROOT->SetStyle("Plain");
  //gStyle->SetOptStat(11);
  gStyle->SetOptStat(0);
  gStyle->SetStatFontSize(0.020);
  gStyle->SetOptFit(0); // 1111
  gStyle->SetOptTitle(1);
  gStyle->SetTitleFontSize(0.05);
  //gStyle->SetTitleX(0.17);
  //gStyle->SetTitleAlign(13);
  gStyle->SetTitleOffset(1.20, "x");
  gStyle->SetTitleOffset(1.10, "y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetNdivisions(510,"X");
  gStyle->SetNdivisions(510,"Y");
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.12);

  //gStyle->SetCanvasColor(-1); 
  //gStyle->SetPadColor(-1); 
  //gStyle->SetFrameFillColor(-1); 
  //gStyle->SetHistFillColor(-1); 
  //gStyle->SetTitleFillColor(-1); 
  //gStyle->SetFillColor(-1); 
  gStyle->SetFillStyle(0000); 
  gStyle->SetStatStyle(0); 
  gStyle->SetTitleStyle(0); 
  gStyle->SetCanvasBorderSize(0); 
  gStyle->SetFrameBorderSize(0); 
  gStyle->SetLegendBorderSize(0); 
  gStyle->SetStatBorderSize(0); 
  gStyle->SetTitleBorderSize(0);
  
  // East residuals
  double En[5];
  En[0] = peakCe;//98.2;
  En[1] = peakIn;
  En[2] = peakSn;//331.2;
  En[3] = peakBiLow; //443.0;
  En[4] = peakBiHigh;//928.0;
 
  double dEn[5] = {0.};

  //2011
  double res2011[5];
  res2011[0] = 0.25;
  res2011[1] = -.4;
  res2011[2] =  -3.0;
  res2011[3] = 0.75;
  res2011[4] = -1.69;

  double sig2011[5];
  sig2011[0] = 1.95;
  sig2011[1] = 1.9;
  sig2011[2] = 3.9;
  sig2011[3] = 5.7;
  sig2011[4] = 8.;
  

  //2012
  double res2012[5];
  res2012[0] = -0.34;
  res2012[1] = 2.8;
  res2012[2] = 0.64;
  res2012[3] = -1.79;
  res2012[4] = -3.9;

  double sig2012[5];
  sig2012[0] = 1.72;
  sig2012[1] = 2.1;
  sig2012[2] = 2.93;
  sig2012[3] = 3.8;
  sig2012[4] = 5.52;



  // Plot
  c1 = new TCanvas("c1", "canvas");
  c1->SetLogy(0);

  TMultiGraph *mg = new TMultiGraph();

  TGraphErrors *gr2011 = new TGraphErrors(5,En,res2011,dEn,sig2011);
  gr2011->SetTitle("2011-2012");
  gr2011->SetMarkerColor(3);
  gr2011->SetLineColor(3);
  gr2011->SetLineWidth(2);
  gr2011->SetMarkerStyle(21);
  gr2011->SetMarkerSize(1.5);
  gr2011->SetFillStyle(0);

  TGraphErrors *gr2012 = new TGraphErrors(5,En,res2012,dEn,sig2012);
  gr2012->SetTitle("2012-2013");
  gr2012->SetMarkerColor(4);
  gr2012->SetLineColor(4);
  gr2012->SetLineWidth(2);
  gr2012->SetMarkerStyle(22);
  gr2012->SetMarkerSize(1.5);
  gr2012->SetFillStyle(0);
  

  
  mg->Add(gr2011,"P");
  //mg->Draw("A");
  mg->Add(gr2012,"P");
 
  mg->Draw("A");
  mg->GetXaxis()->SetTitle("E_{recon} [keV]");
  mg->GetXaxis()->SetTitleOffset(1.2);
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->SetTitle("Calibration Residual [keV]");
  mg->GetYaxis()->SetTitleOffset(1.2);
  mg->GetYaxis()->CenterTitle();
  c1->BuildLegend();

  mg->GetXaxis()->SetLimits(0.0,1200.0);
  mg->SetMinimum(-30.0);
  mg->SetMaximum( 30.0);

  

 Int_t n = 2;
  Double_t x[n] = {0, 1200};
  Double_t y[n] = {0.0, 0.0};

  TGraph *gr0 = new TGraph(n,x,y);
  gr0->Draw("Same");
  gr0->SetLineWidth(2);
  gr0->SetLineColor(1);
  gr0->SetLineStyle(2);

  Int_t nn = 5;
  Double_t perc1=0.017;
  Double_t perc2=0.014;
  Double_t perc3=0.007;
  Double_t perc4=0.007;
  Double_t x2[nn] = {0., peakCe, peakSn, peakBiHigh, 1200.};
  Double_t percent[nn] = {1., perc1, perc2, perc3, peakBiHigh/1200.*perc4};
  Double_t y_upper[nn], y_lower[nn];
  for (int i=1; i<nn; i++) {
    Double_t val = x2[i]*percent[i];
    y_upper[i]=val;
    y_lower[i]=-val;
  }
  y_upper[0]=y_upper[1];
  y_lower[0]=y_lower[1];
    

  TGraph *env_upper = new TGraph(nn,x2,y_upper);
  env_upper->Draw("Same");
  env_upper->SetLineWidth(3);
  env_upper->SetLineColor(2);
  env_upper->SetLineStyle(8);
  
  TGraph *env_lower = new TGraph(nn,x2,y_lower);
  env_lower->Draw("Same");
  env_lower->SetLineWidth(3);
  env_lower->SetLineColor(2);
  env_lower->SetLineStyle(8); 
}
