
#include <string>
#include <vector>

// year = 2011 or 2012
void TwoDim_ErrorPlotter(Int_t year) {

  gStyle->SetTitleSize(0.08,"t");
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleYSize(0.04);//0.5
  gStyle->SetTitleYOffset(1.9);//1.3
  gStyle->SetTitleXSize(0.04);
  gStyle->SetTitleXOffset(1.9);
  gStyle->SetTitleSize(0.04,"Z");
  gStyle->SetTitleOffset(1.6,"Z");
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetOptFit(1111);
  gStyle->SetStatY(0.85);
  gStyle->SetStatX(0.975);
  gStyle->SetStatW(.09);

  
  

  gStyle->SetFillStyle(0000); 
  gStyle->SetStatStyle(0);
  gStyle->SetOptStat(0);
  //gStyle->SetTitleStyle(0); 
  //gStyle->SetCanvasBorderSize(0); 
  //gStyle->SetFrameBorderSize(0); 
  gStyle->SetLegendBorderSize(0); 
  //gStyle->SetStatBorderSize(0); 
  //gStyle->SetTitleBorderSize(0);

  std::vector < Double_t > asymm;
  std::vector < Double_t > energy;
  std::vector < Double_t > statistics;
  std::vector < Double_t > sim;
  std::vector < Double_t > total;

  std::vector <Int_t> xBinVal;
  std::vector <Int_t> yBinVal;
 
  std::vector <Double_t> en150;
  std::vector <Double_t> en200;
  std::vector <Double_t> en250;
  std::vector <Double_t> en300;
 
  std::vector <Double_t> a150;
  std::vector <Double_t> a200;
  std::vector <Double_t> a250;
  std::vector <Double_t> a300;

  std::vector <Double_t> a150_err;
  std::vector <Double_t> a200_err;
  std::vector <Double_t> a250_err;
  std::vector <Double_t> a300_err;
  
  Double_t as_minval=1000., as_maxval=-1000.;
  Double_t en_minval=1000., en_maxval=-1000.;
  Double_t stat_minval=1000., stat_maxval=-1000.;
  Double_t mc_minval=1000., mc_maxval=-1000.;
  Double_t tot_minval=1000., tot_maxval=-1000.;

  ifstream infile(TString::Format("%i_Uncertainties.txt",year).Data());

  std::string stars = "";
  std::string windowStr = "";
  std::string txtHold = "";
  Double_t asymm_hold,en_hold,stat_hold,mc_hold,tot_hold;
  
  while ( infile >> stars 
	  >> txtHold >> asymm_hold
	  >> windowStr
	  >> txtHold >> tot_hold  
	  >> txtHold >> en_hold
	  >> txtHold >> stat_hold
	  >> txtHold >> mc_hold )  {

    Double_t upperWindow = atof(windowStr.substr(4,3).c_str());
    Double_t lowerWindow = atof(windowStr.substr(0,3).c_str());

    
    if ( lowerWindow==150. ) {
      std::cout << lowerWindow << endl;
      
      en150.push_back(upperWindow);
      a150.push_back(asymm_hold);
      a150_err.push_back(TMath::Abs(asymm_hold*stat_hold));
    }
    if ( lowerWindow==200. ) {
      en200.push_back(upperWindow);
      a200.push_back(asymm_hold);
      a200_err.push_back(TMath::Abs(asymm_hold*stat_hold));
    }
    if ( lowerWindow==250. ) {
      en250.push_back(upperWindow);
      a250.push_back(asymm_hold);
      a250_err.push_back(TMath::Abs(asymm_hold*stat_hold));
    }
    if ( lowerWindow==300. ) {
      en300.push_back(upperWindow);
      a300.push_back(asymm_hold);
      a300_err.push_back(TMath::Abs(asymm_hold*stat_hold));
    }
    
    tot_hold*=100.;
    en_hold*=100.;
    mc_hold*=100.;
    stat_hold*=100.;

    if (as_minval>asymm_hold) as_minval=asymm_hold;
    if (as_maxval<asymm_hold) as_maxval=asymm_hold;
    if (en_minval>en_hold) en_minval=en_hold;
    if (en_maxval<en_hold) en_maxval=en_hold;
    if (stat_minval>stat_hold) stat_minval=stat_hold;
    if (stat_maxval<stat_hold) stat_maxval=stat_hold;
    if (mc_minval>mc_hold) mc_minval=mc_hold;
    if (mc_maxval<mc_hold) mc_maxval=mc_hold;
    if (tot_minval>tot_hold) tot_minval=tot_hold;
    if (tot_maxval<tot_hold) tot_maxval=tot_hold;
    
    xBinVal.push_back( upperWindow );
    yBinVal.push_back( lowerWindow );

    asymm.push_back(asymm_hold);
    energy.push_back(en_hold);
    statistics.push_back(stat_hold);
    sim.push_back(mc_hold);
    total.push_back(tot_hold);
    
  }
  Double_t ybinLow = 95.;//495.;
  Double_t ybinHigh = 495.;//805.;
  Int_t nbinsY = (ybinHigh-ybinLow)/10.;

  Double_t xbinLow = 495.;//95.;
  Double_t xbinHigh = 805.;//495;
  Int_t nbinsX = (xbinHigh-xbinLow)/10.;
  
  //Two dimensional plots for each error in minimization
  TH2D *as = new TH2D("as","Asymmetry vs. Analysis Window",
		      nbinsX,xbinLow,xbinHigh,nbinsY,ybinLow,ybinHigh);
  
  TH2D *en = new TH2D("en","Energy Uncertainty vs. Analysis Window",
		      nbinsX,xbinLow,xbinHigh,nbinsY,ybinLow,ybinHigh);

  TH2D *stat = new TH2D("stat","Statistical Uncertainty vs. Analysis Window",
		      nbinsX,xbinLow,xbinHigh,nbinsY,ybinLow,ybinHigh);

  TH2D *mc = new TH2D("mc","Monte Carlo Uncertainty vs. Analysis Window",
		      nbinsX,xbinLow,xbinHigh,nbinsY,ybinLow,ybinHigh);

  TH2D *tot = new TH2D("tot","Total Uncertainty vs. Analysis Window",
		      nbinsX,xbinLow,xbinHigh,nbinsY,ybinLow,ybinHigh);

  for (UInt_t bin=0; bin<xBinVal.size(); ++bin) {

    Int_t xbin = tot->GetXaxis()->FindBin(xBinVal[bin]);
    Int_t ybin = tot->GetYaxis()->FindBin(yBinVal[bin]);

    as->SetBinContent(xbin,ybin,-asymm[bin]);
    as->SetBinError(xbin,ybin,-asymm[bin]*statistics[bin]/100.);
    en->SetBinContent(xbin,ybin,energy[bin]);
    stat->SetBinContent(xbin,ybin,statistics[bin]);
    mc->SetBinContent(xbin,ybin,sim[bin]);
    tot->SetBinContent(xbin,ybin,total[bin]);

  }


  /*const Int_t Number=4;
  Double_t Red[Number]    = { 1.00, 0.00, 1.00, 1.00};
  Double_t Green[Number]   = { 1.00, 0.00, 1.00, 0.00};
  Double_t Blue[Number]  = { 1.00, 1.00, 1.00, 0.00};
  Int_t nb=500;
  
  Double_t Length[Number] = {0.0, en_minval, 0.50*(en_maxval+en_minval), en_maxval}; 
  TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);*/

  gStyle->SetPalette(104);
  
  TCanvas *cEn = new TCanvas("cEn","cEn");
  en->Draw("LEGO20");
  en->GetXaxis()->SetTitle("Maximum Energy (keV)");
  en->GetYaxis()->SetTitle("Minimum Energy (keV)");
  en->GetZaxis()->SetRangeUser(en_minval*0.99,en_maxval);
  en->GetZaxis()->SetTitle("% Error");
  en->GetXaxis()->CenterTitle();
  en->GetYaxis()->CenterTitle();
  en->GetZaxis()->CenterTitle();  
    
  TCanvas *cStat = new TCanvas("cStat","cStat");
  stat->Draw("lego20");
  stat->GetXaxis()->SetTitle("Maximum Energy (keV)");
  stat->GetYaxis()->SetTitle("Minimum Energy (keV)");
  stat->GetZaxis()->SetRangeUser(stat_minval*0.99,stat_maxval);
  stat->GetZaxis()->SetTitle("% Error");
  stat->GetXaxis()->CenterTitle();
  stat->GetYaxis()->CenterTitle();
  stat->GetZaxis()->CenterTitle();
  
  TCanvas *cMc = new TCanvas("cMc","cMc");
  mc->Draw("lego20");
  mc->GetXaxis()->SetTitle("Maximum Energy (keV)");
  mc->GetYaxis()->SetTitle("Minimum Energy (keV)");
  mc->GetZaxis()->SetRangeUser(mc_minval*0.99,mc_maxval);
  mc->GetZaxis()->SetTitle("% Error");
  mc->GetXaxis()->CenterTitle();
  mc->GetYaxis()->CenterTitle();
  mc->GetZaxis()->CenterTitle();


  /*Length[1] = tot_minval;
  Length[2] = 0.50*(tot_maxval+tot_minval);
  Length[3] = tot_maxval; 
  TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
  */
  
 
  
  TCanvas *cTot = new TCanvas("cTot","cTot");
  tot->Draw("lego20");
  tot->GetXaxis()->SetTitle("Maximum Energy (keV)");
  tot->GetYaxis()->SetTitle("Minimum Energy (keV)");
  tot->GetZaxis()->SetRangeUser(tot_minval*0.99,tot_maxval);
  tot->GetZaxis()->SetTitle("% Error");
  tot->GetXaxis()->CenterTitle();
  tot->GetYaxis()->CenterTitle();
  tot->GetZaxis()->CenterTitle();

  gStyle->SetErrorX(0);

  TCanvas *cAs = new TCanvas("cAs","cAs");
  as->Draw("E");
  as->GetXaxis()->SetTitle("Maximum Energy (keV)");
  as->GetYaxis()->SetTitle("Minimum Energy (keV)");
  as->SetMarkerStyle(25);
  as->SetMarkerColor(kBlue);
  as->GetZaxis()->SetRangeUser(-as_maxval*0.99,-as_minval*1.01);
  //as->SetMinimum(0.120);
  //as->SetMaximum(0.125);
  as->GetZaxis()->SetTitle("|A|");
  as->GetXaxis()->CenterTitle();
  as->GetYaxis()->CenterTitle();
  as->GetZaxis()->CenterTitle();
  cAs->Update();

  TCanvas *c1 = new TCanvas("c1","c1",1200,1200);
  c1->Divide(2,2);

  c1->cd(1);

  for (UInt_t i=0; i<en150.size(); ++i) {
    std::cout << en150[i] << " " << a150[i] << " " << a150_err[i] << endl;
   
  }
  
  TGraphErrors *g150 = new TGraphErrors(en150.size(),&en150[0],&a150[0],
					0,&a150_err[0]);
  g150->SetMarkerStyle(kOpenSquare);
  g150->SetLineWidth(2);
  g150->SetTitle("Lower Window Edge at 150 keV");
  g150->GetXaxis()->SetTitle("Upper Window Edge (keV)");
  g150->GetYaxis()->SetTitle("Asymmetry");
  
  g150->Draw("AP");

  c1->cd(2);

  TGraphErrors *g200 = new TGraphErrors(en200.size(),&en200[0],&a200[0],
					0,&a200_err[0]);
  g200->SetMarkerStyle(kOpenSquare);
  g200->SetLineWidth(2);
  g200->GetXaxis()->SetTitle("Upper Window Edge (keV)");
  g200->GetYaxis()->SetTitle("Asymmetry");
  g200->SetTitle("Lower Window Edge at 200 keV");

  
  g200->Draw("AP");

  c1->cd(3);

  TGraphErrors *g250 = new TGraphErrors(en250.size(),&en250[0],&a250[0],
					0,&a250_err[0]);
  g250->SetMarkerStyle(kOpenSquare);
  g250->SetLineWidth(2);
  g250->GetXaxis()->SetTitle("Upper Window Edge (keV)");
  g250->GetYaxis()->SetTitle("Asymmetry");
  g250->SetTitle("Lower Window Edge at 250 keV");

  g250->Draw("AP");

  c1->cd(4);

  TGraphErrors *g300 = new TGraphErrors(en300.size(),&en300[0],&a300[0],
					0,&a300_err[0]);
  g300->SetMarkerStyle(kOpenSquare);
  g300->SetLineWidth(2);
  g300->GetXaxis()->SetTitle("Upper Window Edge (keV)");
  g300->GetYaxis()->SetTitle("Asymmetry");
  g300->SetTitle("Lower Window Edge at 300 keV");

  g300->Draw("AP");
  
}
