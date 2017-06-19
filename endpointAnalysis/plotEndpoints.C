#include <TGraphErrors.h>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <TCanvas.h>

std::vector < std::vector <Double_t> > readEndpoints(Int_t octet, bool data) {

  std::vector < std::vector <Double_t> > ep(2,std::vector<Double_t>(2,0.));
  
  std::ifstream infile(TString::Format("%s/FinalEndpoints/%s_finalEndpoints_Octet_%i.txt",getenv("ENDPOINT_ANALYSIS"),data?"UK":"SIM",octet).Data());

  std::string hold;
  Double_t val, valErr;
  Int_t inc=0;
  
  if ( infile.is_open() ) {
    while ( infile >> hold >> val >> hold >> valErr ) {
      ep[inc][0] = val;
      ep[inc][1] = valErr;
      inc++;
    }
  }
  else {
    std::cout << "Couldn't read endpoints for octet " << octet << "\n";
    ep[0][0] = 780.; ep[0][1] = 1000.;
    ep[1][0] = 780.; ep[1][1] = 1000.;
  }

  return ep;

};

void plotEndpoints(TString geometry) {

  gStyle->SetTitleSize(0.08,"t");
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetTitleYSize(0.05);
  gStyle->SetTitleYOffset(1.);
  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetOptStat(0);
  //gStyle->SetOptFit(1111);
  //gStyle->SetStatY(0.85);
  //gStyle->SetStatX(0.975);
  //gStyle->SetStatW(.09);

  //gStyle->SetFillStyle(0000); 
  //gStyle->SetStatStyle(0); 
  //gStyle->SetTitleStyle(0); 
  //gStyle->SetCanvasBorderSize(0); 
  //gStyle->SetFrameBorderSize(0); 
  //gStyle->SetLegendBorderSize(0); 
  //gStyle->SetStatBorderSize(0); 
  //gStyle->SetTitleBorderSize(0);

  //gStyle->SetGridStyle(4);
  //gStyle->SetGridColor(kBlue);

  Int_t badocts[] = {7,60,61,62,63,64,65,66,67,91,93,101,107,121};
  

  Int_t octMin, octMax;

  if (geometry==TString("2011-2012") ) {
    octMin = 0;
    octMax = 59;
  }
  else {
    octMin = 60; 
    octMax = 121;
  }

  //Vectors for putting endpoint information
  std::vector <Double_t> DATA_east_ep;
  std::vector <Double_t> DATA_east_epErr;
  std::vector <Double_t> DATA_west_ep;
  std::vector <Double_t> DATA_west_epErr;

  std::vector <Double_t> SIM_east_ep;
  std::vector <Double_t> SIM_east_epErr;
  std::vector <Double_t> SIM_west_ep;
  std::vector <Double_t> SIM_west_epErr;

  std::vector <Double_t> octets;

  //// Read in all necessary endpoints
  for ( int octet=octMin ; octet<=octMax ; ++octet ) {

    if ( std::find(badocts, badocts+14,octet) != (badocts+14) ) continue;  //Checking if octet should be ignored for data quality reasons

    octets.push_back(octet);

    std::vector < std::vector <Double_t> > DATA_ep = readEndpoints(octet,true);
    DATA_east_ep.push_back( DATA_ep[0][0] );
    DATA_east_epErr.push_back( DATA_ep[0][1] );
    DATA_west_ep.push_back( DATA_ep[1][0] );
    DATA_west_epErr.push_back( DATA_ep[1][1] );

    std::vector < std::vector <Double_t> > SIM_ep = readEndpoints(octet,false);
    SIM_east_ep.push_back( SIM_ep[0][0] );
    SIM_east_epErr.push_back( SIM_ep[0][1] );
    SIM_west_ep.push_back( SIM_ep[1][0] );
    SIM_west_epErr.push_back( SIM_ep[1][1] );

  }

  TGraphErrors *DATA_eastGraph = new TGraphErrors(octets.size(),&octets[0],
				       &DATA_east_ep[0],0,&DATA_east_epErr[0]);
  DATA_eastGraph->SetMarkerColor(kBlue);
  DATA_eastGraph->SetLineColor(kBlue);
  DATA_eastGraph->SetMarkerStyle(21);
  DATA_eastGraph->SetTitle(TString::Format("%s East Beta Decay Endpoints",geometry.Data()));

  TGraphErrors *DATA_westGraph = new TGraphErrors(octets.size(),&octets[0],
				       &DATA_west_ep[0],0,&DATA_west_epErr[0]);
  DATA_westGraph->SetMarkerColor(kBlue);
  DATA_westGraph->SetLineColor(kBlue);
  DATA_westGraph->SetMarkerStyle(21);
  DATA_westGraph->SetTitle(TString::Format("%s West Beta Decay Endpoints",geometry.Data()));
  

  TGraphErrors *SIM_eastGraph = new TGraphErrors(octets.size(),&octets[0],
				       &SIM_east_ep[0],0,&SIM_east_epErr[0]);
  SIM_eastGraph->SetMarkerColor(kRed);
  SIM_eastGraph->SetLineColor(kRed);
  SIM_eastGraph->SetMarkerStyle(22);
  SIM_eastGraph->SetTitle("East Beta Decay Endpoints");

  TGraphErrors *SIM_westGraph = new TGraphErrors(octets.size(),&octets[0],
				       &SIM_west_ep[0],0,&SIM_west_epErr[0]);
  SIM_westGraph->SetMarkerColor(kRed);
  SIM_westGraph->SetLineColor(kRed);
  SIM_westGraph->SetMarkerStyle(22);
  SIM_westGraph->SetTitle("West Beta Decay Endpoints");

  TCanvas *c = new TCanvas("c","c",1200,1200);
  c->SetGrid();
  c->Divide(1,2);

  c->cd(1);
  gPad->SetGridy(1);
  DATA_eastGraph->Draw("AP");
  DATA_eastGraph->SetMinimum(750.);
  DATA_eastGraph->SetMaximum(815.);
  DATA_eastGraph->GetXaxis()->SetRangeUser(octets[0]-10.,octets[octets.size()-1]+10.); 
  DATA_eastGraph->GetXaxis()->SetTitle("Octet number"); 
  DATA_eastGraph->GetYaxis()->SetTitle("Kinetic Energy Endpoint (keV)"); 
  SIM_eastGraph->Draw("PSAME");

  c->cd(2);
  gPad->SetGridy(1);
  DATA_westGraph->Draw("AP");
  DATA_westGraph->SetMinimum(750.);
  DATA_westGraph->SetMaximum(815.);
  DATA_westGraph->GetXaxis()->SetRangeUser(octets[0]-10.,octets[octets.size()-1]+10.);  
  DATA_westGraph->GetXaxis()->SetTitle("Octet number"); 
  DATA_westGraph->GetYaxis()->SetTitle("Kinetic Energy Endpoint (keV)"); 
  SIM_westGraph->Draw("PSAME");

  c->Update();
  c->Modified();

  c->Print(TString::Format("%s_endpoints.pdf(",geometry.Data()));

  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);
  
  
  TH1D *epDataE = new TH1D("epDataE",TString::Format("%s East Beta Decay Endpoint Distribution",geometry.Data()), 150, 700., 820.);
  TH1D *epDataW = new TH1D("epDataW",TString::Format("%s West Beta Decay Endpoint Distribution",geometry.Data()), 150, 700., 820.);
  TH1D *epSimE = new TH1D("epSimE",TString::Format("%s East Beta Decay Endpoint Distribution",geometry.Data()), 150, 700., 820.);
  TH1D *epSimW = new TH1D("epSimW",TString::Format("%s West Beta Decay Endpoint Distribution",geometry.Data()), 150, 700., 820.);

  epSimE->SetLineColor(kRed);
  epSimW->SetLineColor(kRed);

  for ( UInt_t i=0; i<octets.size(); ++i ) {
    
    epDataE->Fill(DATA_east_ep[i]);
    epDataW->Fill(DATA_west_ep[i]);
    epSimE->Fill(SIM_east_ep[i]);
    epSimW->Fill(SIM_west_ep[i]);

  }

  std::cout << "endpoint Mean and RMS East: " << epSimE->GetMean() - epDataE->GetMean() << " +/- " << epDataE->GetRMS() << std::endl;
  std::cout << "endpoint Mean and RMS West: " << epSimW->GetMean() - epDataW->GetMean() << " +/- " << epDataW->GetRMS() << std::endl;

  c2->cd(1);
  epSimE->Draw();
  epDataE->Draw("SAME");
  

  c2->cd(2); 
  epSimW->Draw();
  epDataW->Draw("SAME");

  c2->Print(TString::Format("%s_endpoints.pdf)",geometry.Data()));

}

  
