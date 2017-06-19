
//Should make an array of all the calibration periods and the
// runs they apply to, and then put vertical lines at each of these
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>

#include <TString.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TStyle.h>
#include <TF1.h>
#include <TList.h>
#include <TPaveStats.h>
#include <TLegend.h>

const static double peakCe = 130.3;// 131.5956;//80.5;
const static double peakIn = 174.35;
const static double peakSn = 368.4938;//317.8;
const static double peakBiLow = 498.;//501.420;//448.8;
const static double peakBiHigh = 993.789;//926.;

std::vector < Double_t > getOctetAsym(Int_t octet) {

  std::vector <Double_t> ret;

  TString path = TString::Format("%s/Octet_%i/OctetAsymmetry/UnCorr_withPOL_FittedAsymmetry_Octet%i_AnaChA_220-680.dat", getenv("ANALYSIS_RESULTS"),octet,octet);
  std::cout << path << std::endl;

  std::ifstream infile(path.Data());  
  std::string txt = "";
  Double_t Asym = 0., AsymError = 0.;

  if (infile.is_open()) {
    infile >> txt >> Asym >> AsymError;
    ret.push_back(Asym);
    ret.push_back(AsymError);
    infile.close();
  }
  else { 
    std::cout << "Couldn't open Asymmetry File!!\n";
    ret.push_back(0.); 
    ret.push_back(0.);
  }
  
  return ret;
}

//This returns a vector or RMS values, where the order is Ce,In,Sn,Bi2,Bi1
std::vector < std::vector < Double_t > > getXePeriodEnvelope(Int_t XePeriod) {

  
  std::vector < std::vector <Double_t> > ret(2,std::vector<Double_t>(0));
  Int_t lowCal, highCal;
  lowCal = highCal = 0;

  if ( XePeriod == 2 ) lowCal = 1, highCal = 4;
  else if (XePeriod == 3 ) lowCal = 5, highCal = 5;
  else if (XePeriod == 4 ) lowCal = 6, highCal = 6;
  else if (XePeriod == 5 ) lowCal = 7, highCal = 8;
  else if (XePeriod == 7 ) lowCal = 9, highCal = 12;
  else if (XePeriod == 8 ) lowCal = 16, highCal = 17;
  else if (XePeriod == 9 ) lowCal = 18, highCal = 20;
  else if (XePeriod == 10 ) lowCal = 21, highCal = 24;
  else std::cout << "Bad Xe Period\n", exit(0);

  TString path = TString::Format("%s/error_envelope/error_envelope_calPeriods_%i-%i.dat", getenv("ANALYSIS_CODE"),lowCal,highCal);
  std::cout << path << std::endl;

  std::ifstream infile(path.Data());  
  std::string txt = "";

  if (infile.is_open()) {
    for ( UInt_t i = 1; i<16; ++i ) {
      infile >> txt >> txt >> txt;
      if ( (i+2)%3 == 0 && (i+2)!=6 ) ret[1].push_back(atof(txt.c_str()));
      if ( i%3 == 0 && i!=6 ) ret[0].push_back(atof(txt.c_str()));
    }
  }
  else {
    std::cout << "Couldn't open error envelope file!!\n";
    exit(0);
  }
  infile.close();
  return ret;
}

//This returns mean, sigma, and rms, where the order is Ce,In,Sn,Bi2,Bi1
std::vector < std::vector < Double_t > > getGeomEnvelope(TString geom) {

  
  std::vector < std::vector <Double_t> > ret(3,std::vector<Double_t>(0));
  Int_t lowCal, highCal;
  
  if ( geom==TString("2011-2012") ) lowCal=1, highCal=12;
  else if ( geom==TString("2012-2013") ) lowCal=16, highCal=24;
  
  else std::cout << "Bad geometry passed to getGeomEnvelope\n", exit(0);

  TString path = TString::Format("%s/error_envelope/error_envelope_calPeriods_%i-%i.dat", getenv("ANALYSIS_CODE"),lowCal,highCal);
  std::cout << path << std::endl;

  std::ifstream infile(path.Data());  
  std::string txt = "";

  if (infile.is_open()) {
    for ( UInt_t i = 1; i<16; ++i ) {
      infile >> txt >> txt >> txt;
      if ( i%3 == 0 && i!=6 ) ret[2].push_back(atof(txt.c_str())); // RMS
      if ( (i+1)%3 == 0 && i!=5 ) ret[1].push_back(atof(txt.c_str())); // Err
      if ( (i+2)%3 == 0 && i!=4 ) ret[0].push_back(atof(txt.c_str())); // mean
     
    }
  }
  else {
    std::cout << "Couldn't open error envelope file!!\n";
    exit(0);
  }
  infile.close();
  return ret;
}


void asymm_vs_cal(TString year) {

  gStyle->SetTitleSize(0.08,"t");
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleYSize(0.05);
  gStyle->SetTitleYOffset(0.8);
  gStyle->SetTitleXSize(0.06);
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetOptFit(1111);
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

  std::vector <Int_t> badOct;
  Int_t octs[] = {7,9,59,60,61,62,63,64,65,66,70,92};
  badOct.assign(octs, octs+12);

  std::vector <Int_t> XePeriodOctetEnd; 
  std::vector <Int_t> XePeriod;
  std::vector <Int_t> SrcPeriodOctetEnd; 
  std::vector <Int_t> SrcPeriod;

  Int_t octMin, octMax;

  if (year==TString("2011-2012") ) {
    octMin = 0;
    octMax = 59;
    
    XePeriodOctetEnd.push_back(14); XePeriod.push_back(2);
    XePeriodOctetEnd.push_back(23); XePeriod.push_back(3);
    XePeriodOctetEnd.push_back(31); XePeriod.push_back(4);
    XePeriodOctetEnd.push_back(46); XePeriod.push_back(5);
    XePeriodOctetEnd.push_back(59); XePeriod.push_back(7);

    SrcPeriodOctetEnd.push_back(4); SrcPeriod.push_back(1);
    SrcPeriodOctetEnd.push_back(6); SrcPeriod.push_back(2);
    SrcPeriodOctetEnd.push_back(9); SrcPeriod.push_back(3);
    SrcPeriodOctetEnd.push_back(14); SrcPeriod.push_back(4);
    SrcPeriodOctetEnd.push_back(23); SrcPeriod.push_back(5);
    SrcPeriodOctetEnd.push_back(31); SrcPeriod.push_back(6);
    SrcPeriodOctetEnd.push_back(39); SrcPeriod.push_back(7);
    SrcPeriodOctetEnd.push_back(46); SrcPeriod.push_back(8);
    SrcPeriodOctetEnd.push_back(50); SrcPeriod.push_back(9);
    SrcPeriodOctetEnd.push_back(59); SrcPeriod.push_back(11);
  }
  else {
    octMin = 60;
    octMax = 121;

    XePeriodOctetEnd.push_back(79); XePeriod.push_back(8);
    XePeriodOctetEnd.push_back(95); XePeriod.push_back(9); 
    XePeriodOctetEnd.push_back(121); XePeriod.push_back(10);

    SrcPeriodOctetEnd.push_back(71); SrcPeriod.push_back(16);
    SrcPeriodOctetEnd.push_back(79); SrcPeriod.push_back(17);
    SrcPeriodOctetEnd.push_back(85); SrcPeriod.push_back(18);
    SrcPeriodOctetEnd.push_back(91); SrcPeriod.push_back(19);
    SrcPeriodOctetEnd.push_back(95); SrcPeriod.push_back(20);
    SrcPeriodOctetEnd.push_back(105); SrcPeriod.push_back(22);
    SrcPeriodOctetEnd.push_back(121); SrcPeriod.push_back(23);
    
  }

  


  // making all of the TLines to reperesent the Xe calibration periods

  std::vector <TLine*> linesXe(XePeriodOctetEnd.size(),0);
  std::vector <TLine*> linesSrc(SrcPeriodOctetEnd.size(),0);

  std::vector <Double_t> dOctetList;
  std::vector < std::vector <Double_t> > AsymAndError(2,std::vector<Double_t>(0));

  //Read in the Asymmetries

  for ( Int_t oct = octMin; oct<=octMax; ++oct ) {

    if ( std::find(badOct.begin(), badOct.end(),oct) != badOct.end() ) continue;  //Checking if octet should be ignored for data quality reasons
    
    std::vector <Double_t> vec = getOctetAsym(oct); 
    
    if ( vec[1]>0. ) {
      dOctetList.push_back(oct);
      AsymAndError[0].push_back(vec[0]);
      AsymAndError[1].push_back(vec[1]);
    }
  }
    
  TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
  c1->Divide(1,2);
  c1->cd(1);

  //Plotting Asymmetries
  TGraphErrors *g = new TGraphErrors(dOctetList.size(), &dOctetList[0], &AsymAndError[0][0], 0, &AsymAndError[1][0]);
  TString title = TString::Format("%s Raw Measured Asymmetry 220-680 keV Window",year.Data());
  g->SetTitle(title);
  g->SetMarkerStyle(20);
  g->SetLineWidth(2);
  g->GetXaxis()->SetLimits(dOctetList[0]-2., dOctetList[dOctetList.size()-1]+2.);
  g->GetXaxis()->SetTitle("Number");
  g->GetYaxis()->SetTitle("Raw Asymmetry");
  g->GetXaxis()->CenterTitle();
  g->GetYaxis()->CenterTitle();
  g->Draw("APZ");
  
  c1->Update();
  TF1 *fit = new TF1("fit","[0]",dOctetList[0], dOctetList[dOctetList.size()-1]);
  fit->SetLineColor(kRed);
  fit->SetLineWidth(3);
  fit->SetParameter(0,0.05);
  
  g->Fit("fit","R");

  /* TPaveStats *ps = (TPaveStats *)g->GetListOfFunctions()->FindObject("stats");
  ps->SetX1NDC(0.75);
  ps->SetX2NDC(0.95);

  c1->Update();
  c1->Modified();*/

  Double_t min = 0., max = 0.;
  min = gPad->GetUymin();
  max = gPad->GetUymax();

  
  for (UInt_t i=0; i<SrcPeriodOctetEnd.size(); i++) {
    
    linesSrc[i] = new TLine((double)SrcPeriodOctetEnd[i]+0.5,min,(double)SrcPeriodOctetEnd[i]+0.5,max);
    linesSrc[i]->SetLineStyle(2);
    //linesSrc[i]->SetLineColor(kRed);
    linesSrc[i]->SetLineWidth(1);
    if ( SrcPeriodOctetEnd[i] > octMin && SrcPeriodOctetEnd[i] < octMax ) {
      linesSrc[i]->Draw();
    }
  } 

  for (UInt_t i=0; i<XePeriodOctetEnd.size(); i++) {
    
    linesXe[i] = new TLine((double)XePeriodOctetEnd[i]+0.5,min,(double)XePeriodOctetEnd[i]+0.5,max);
    linesXe[i]->SetLineStyle(2);
    linesXe[i]->SetLineColor(kBlue);
    linesXe[i]->SetLineWidth(3);
    if ( XePeriodOctetEnd[i] > octMin && XePeriodOctetEnd[i] < octMax ) {
      linesXe[i]->Draw();
    }
  } 

  c1->Update();

  c1->cd(2);
  
  gPad->Divide(XePeriod.size(),1);

  // Now making the error envelope plots
  
  Double_t xval[] = {peakCe, peakSn, peakBiLow, peakBiHigh};
  Double_t yval[] = {0,0,0,0};
  
  std::vector <TGraphErrors*> gEnv(XePeriod.size());
  std::vector <TGraph*> gEnv_mean(XePeriod.size());
  
  Double_t xmin = 0.0;
  Double_t xmax = 1200.;
  TLine *zeroLine = new TLine(xmin,0.,xmax,0.);

  for ( UInt_t ii=0; ii<XePeriod.size(); ++ii ) {

    c1->cd(2); gPad->cd(ii+1);
    
    std::vector < std::vector <Double_t> > env = getXePeriodEnvelope(XePeriod[ii]);

    gEnv[ii] = new TGraphErrors(4, xval, yval, 0, &env[0][0]);
    title = TString::Format("Error Envelope Xe Period %i",XePeriod[ii]);
    gEnv[ii]->SetTitle(title);
    gEnv[ii]->SetMarkerStyle(21);
    gEnv[ii]->SetMarkerSize(0);
    gEnv[ii]->SetMarkerColor(kBlue);
    gEnv[ii]->SetLineWidth(3);
    gEnv[ii]->SetLineColor(kBlue);
    gEnv[ii]->SetMinimum(-12.);
    gEnv[ii]->SetMaximum(12.);
    gEnv[ii]->GetXaxis()->SetLimits(xmin,xmax);
    gEnv[ii]->GetXaxis()->SetTitle("Energy (keV)");
    gEnv[ii]->GetYaxis()->SetTitle("Residual (keV)");
    gEnv[ii]->GetXaxis()->CenterTitle();
    gEnv[ii]->GetYaxis()->CenterTitle();
    //gEnv[ii]->SetFillStyle(3002);
    gEnv[ii]->Draw("APZ");
    
    c1->Update();


    gEnv_mean[ii] = new TGraph(4, xval, &env[1][0]);
    title = TString::Format("Error Envelope Xe Period %i",XePeriod[ii]);
    //gEnv_mean[ii]->SetTitle(title);
    gEnv_mean[ii]->SetMarkerStyle(21);
    gEnv_mean[ii]->SetMarkerSize(1);
    gEnv_mean[ii]->SetMarkerColor(kRed);
    //gEnv_mean[ii]->SetLineWidth(3);
    //gEnv_mean[ii]->SetLineColor(kBlue);
    //gEnv_mean[ii]->SetMinimum(-12.);
    //gEnv_mean[ii]->SetMaximum(12.);
    //gEnv_mean[ii]->GetXaxis()->SetLimits(xmin,xmax);
    //gEnv_mean[ii]->GetXaxis()->SetTitle("Energy (keV)");
    //gEnv_mean[ii]->GetYaxis()->SetTitle("Residual (keV)");
    //gEnv_mean[ii]->GetXaxis()->CenterTitle();
    //gEnv_mean[ii]->GetYaxis()->CenterTitle();
    //gEnv_mean[ii]->SetFillStyle(3002);
    gEnv_mean[ii]->Draw("PSAME");
    
    c1->Update();

    

    zeroLine->Draw();

  }

  c1->Update();

  TString filename = TString::Format("errEnv_vs_XePeriod_%s",year.Data());
  c1->Print(TString::Format("%s.pdf(",filename.Data()));
 

  ////////////////////////////////////////////////////////////////////////////////
  //Now do the total error envelope
  ////////////////////////////////////////////////////////////////////////////////


  TCanvas *c2 = new TCanvas("c2","c2",1600,1200);
  
  double En[4];
  En[0] = peakCe;//98.2;
  //  En[1] = peakIn;
  En[1] = peakSn;//331.2;
  En[2] = peakBiLow; //443.0;
  En[3] = peakBiHigh;//928.0;
 
  double dEn[4] = {0.};
  
  double rmsY[4] = {0.};

  std::vector < std::vector<Double_t> > env = getGeomEnvelope(year);

  TMultiGraph *mg = new TMultiGraph();

  TGraphErrors *RMS = new TGraphErrors(4,En,rmsY,dEn,&env[2][0]);
  RMS->SetTitle(TString::Format("RMS"));
  RMS->SetMarkerColor(8);
  RMS->SetLineColor(8);
  RMS->SetLineWidth(15);
  RMS->SetMarkerStyle(21);
  RMS->SetMarkerSize(0);
  RMS->SetFillStyle(0);

  TGraphErrors *gr = new TGraphErrors(4,En,&env[0][0],dEn,&env[1][0]);
  gr->SetTitle(TString::Format("Mean & Sigma"));
  gr->SetMarkerColor(kBlue);
  gr->SetLineColor(kBlue);
  gr->SetLineWidth(2);
  gr->SetMarkerStyle(21);
  gr->SetMarkerSize(1.25);
  gr->SetFillStyle(0);
  

  
  mg->Add(RMS,"PZ");
  //mg->Draw("A");
  mg->Add(gr,"P");
 
  mg->Draw("A");
  mg->SetTitle(TString::Format("Error Envelope %s",year.Data()));
  mg->GetXaxis()->SetTitle("E_{recon} [keV]");
  mg->GetXaxis()->SetTitleOffset(1.2);
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->SetTitle("Calibration Residual [keV]");
  // mg->GetYaxis()->SetTitleOffset(1.2);
  mg->GetYaxis()->CenterTitle();
  

  mg->GetXaxis()->SetLimits(0.0,1200.0);
  mg->SetMinimum(-30.0);
  mg->SetMaximum( 30.0);

  TLegend *leg = new TLegend(0.65,0.70,0.875,0.875);
  leg->AddEntry(RMS,0,"l");
  leg->AddEntry(gr,0,"lp");
  leg->Draw();

  

  const Int_t n = 2;
  Double_t x[n] = {0, 1200};
  Double_t y[n] = {0.0, 0.0};

  TGraph *gr0 = new TGraph(n,x,y);
  gr0->Draw("Same");
  gr0->SetLineWidth(2);
  gr0->SetLineColor(1);
  gr0->SetLineStyle(2);

  const Int_t nn = 5;
  Double_t perc1=0.017;
  Double_t perc2=0.010;
  Double_t perc3=0.0065;
  Double_t perc4=0.0065;
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

  c2->Update();

  c2->Print(TString::Format("%s.pdf)",filename.Data()));
 
}
