// Usage: root[0] .x fit_source_positions.C("runNumber")

#include <iostream>
#include <fstream>

#include "TSpectrum2.h"
#include "TH2.h"
#include "TF2.h"
#include "TMath.h"

void fit_source_positions(TString runNumber)
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  // Style options
  gROOT->SetStyle("Plain");
  //gStyle->SetOptSta=t(11);
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

  // Input ntuple
  TString filenameIn;
  filenameIn  = TString(getenv("REPLAY_PASS2"))+TString("/replay_pass2_");
  filenameIn += runNumber;
  filenameIn += ".root";
  cout << "Processing ... " << filenameIn << endl;
  TFile *filein = new TFile(filenameIn);

  // Read source list
  int nSources;
  string sourceName[3];
  TString filenameList;
  filenameList  = TString(getenv("SOURCE_LIST"))+TString("/source_list_");
  filenameList += runNumber;
  filenameList += ".dat";
  cout << "Reading sources from ... " << filenameList << endl;

  ifstream fileList(filenameList);
  fileList >> nSources;
  cout << "... Number of sources: " << nSources << endl;
  for (int j=0; j<nSources; j++) {
    fileList >> sourceName[j];
    cout << sourceName[j] << endl;
  }

  // Histograms
  TH2F *hisexy = new TH2F("hisexy", "", 100,-50.0,50.0, 100,-50.0,50.0);
  TH2F *hiswxy = new TH2F("hiswxy", "", 100,-50.0,50.0, 100,-50.0,50.0);

  // Cuts
  TCut *type = new TCut("(Type == 0)");
  TCut *east = new TCut("(Side == 0)");
  TCut *west = new TCut("(Side == 1)");

  // Project ntuple data into histograms
  pass2->Draw("yE.center:xE.center >> hisexy", *type && *east);
  pass2->Draw("yW.center:xW.center >> hiswxy", *type && *west);

  // Find East peaks
  c1 = new TCanvas("c1","c1");
  hisexy->Draw();
  TSpectrum2 *specEast = new TSpectrum2(nSources, 1.0); // (max peaks, resolution)
  Int_t nFoundEast = specEast->Search(hisexy, 2.0, "nomarkov");
  Float_t *xpeaksEast = specEast->GetPositionX();
  Float_t *ypeaksEast = specEast->GetPositionY();
  for (Int_t p=0; p<nFoundEast; p++) {
    cout << p << " " << xpeaksEast[p] << " " << ypeaksEast[p] << endl;
  }

  // East peaks mean and sigma
  double xEast[3]={100.,100.,100.}, yEast[3]={100.,100.,100.}, sigmaEast[3];
  for (Int_t p=0; p<nFoundEast; p++) {
    TF2 *fit2DGauss = new TF2("fit2DGauss", "[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))*exp(-(y-[3])*(y-[3])/(2.*[4]*[4]))", xpeaksEast[p]-5., xpeaksEast[p]+5., ypeaksEast[p]-5., ypeaksEast[p]+5.);
    fit2DGauss->SetParameter(0, 100.0);
    fit2DGauss->SetParameter(1, xpeaksEast[p]);
    fit2DGauss->SetParameter(2, 2.0);
    fit2DGauss->SetParameter(3, ypeaksEast[p]);
    fit2DGauss->SetParameter(4, 2.0);

    hisexy->Fit("fit2DGauss", "R0");
    xEast[p] = fit2DGauss->GetParameter(1);
    yEast[p] = fit2DGauss->GetParameter(3);
    sigmaEast[p] = (fabs(fit2DGauss->GetParameter(2)) + fabs(fit2DGauss->GetParameter(4)))/2.;

    cout << xEast[p] << " " << yEast[p] << " " << sigmaEast[p] << endl;
  }

  // Find West peaks
  c2 = new TCanvas("c2","c2");
  hiswxy->Draw();
  TSpectrum2 *specWest = new TSpectrum2(nSources, 1.0); // (max peaks, resolution)
  Int_t nFoundWest = specWest->Search(hiswxy,2.0, "nomarkov");
  Float_t *xpeaksWest = specWest->GetPositionX();
  Float_t *ypeaksWest = specWest->GetPositionY();
  for (Int_t p=0; p<nFoundWest; p++) {
    cout << p << " " << xpeaksWest[p] << " " << ypeaksWest[p] << endl;
  }

  // West peaks mean and sigma
  double xWest[3], yWest[3], sigmaWest[3];
  for (Int_t p=0; p<nFoundWest; p++) {
    TF2 *fit2DGauss = new TF2("fit2DGauss", "[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))*exp(-(y-[3])*(y-[3])/(2.*[4]*[4]))", xpeaksWest[p]-5., xpeaksWest[p]+5., ypeaksWest[p]-5., ypeaksWest[p]+5.);
    fit2DGauss->SetParameter(0, 100.0);
    fit2DGauss->SetParameter(1, xpeaksWest[p]);
    fit2DGauss->SetParameter(2, 2.0);
    fit2DGauss->SetParameter(3, ypeaksWest[p]);
    fit2DGauss->SetParameter(4, 2.0);

    hiswxy->Fit("fit2DGauss", "R0");
    xWest[p] = fit2DGauss->GetParameter(1);
    yWest[p] = fit2DGauss->GetParameter(3);
    sigmaWest[p] = (fabs(fit2DGauss->GetParameter(2)) + fabs(fit2DGauss->GetParameter(4)))/2.;

    cout << xWest[p] << " " << yWest[p] << " " << sigmaWest[p] << endl;
  }

  // Write fit results to file
  TString filenameOut;
  filenameOut  = TString(getenv("SOURCE_POSITIONS"))+TString("/source_positions_");
  filenameOut += runNumber;
  filenameOut += ".dat";
  ofstream outFit(filenameOut);
  cout << "Writing fitted (x,y) positions to ... " << filenameOut << endl;

  if (nSources == 1) { outFit << xEast[0] << " " << yEast[0] << " " << sigmaEast[0] << " "
		       << xWest[0] << " " << yWest[0] << " " << sigmaWest[0] << endl;
  }

  if (nSources == 2) {
    int pEastMin, pWestMin;
    int pEastMax, pWestMax;
    if (xEast[0] < xEast[1]) {pEastMin = 0; pEastMax = 1;}
    else {pEastMin = 1; pEastMax = 0;}
    if (xWest[0] < xWest[1]) {pWestMin = 0; pWestMax = 1;}
    else {pWestMin = 1; pWestMax = 0;}

    outFit << xEast[pEastMin] << " " << yEast[pEastMin] << " " << sigmaEast[pEastMin] << " "
           << xWest[pWestMin] << " " << yWest[pWestMin] << " " << sigmaWest[pWestMin] << endl;
    outFit << xEast[pEastMax] << " " << yEast[pEastMax] << " " << sigmaEast[pEastMax] << " "
           << xWest[pWestMax] << " " << yWest[pWestMax] << " " << sigmaWest[pWestMax] << endl;

      }

  if (nSources == 3) {
    double xEastMin =  1.0e99;
    double xEastMax = -1.0e99;
    double xWestMin =  1.0e99;
    double xWestMax = -1.0e99;
    int pEastMin, pEastMax, pEastMid;
    int pWestMin, pWestMax, pWestMid;

    for (int p=0; p<nSources; p++) {
      if (xEast[p] < xEastMin) {xEastMin = xEast[p]; pEastMin = p;}
      if (xEast[p] > xEastMax) {xEastMax = xEast[p]; pEastMax = p;}
    }
    for (int p=0; p<nSources; p++) {
      if ((p != pEastMin) && (p != pEastMax)) pEastMid = p;
    }

    for (int p=0; p<nSources; p++) {
      if (xWest[p] < xWestMin) {xWestMin = xWest[p]; pWestMin = p;}
      if (xWest[p] > xWestMax) {xWestMax = xWest[p]; pWestMax = p;}
    }
    for (int p=0; p<nSources; p++) {
      if ((p != pWestMin) && (p != pWestMax)) pWestMid = p;
    }

    outFit << xEast[pEastMin] << " " << yEast[pEastMin] << " " << sigmaEast[pEastMin] << " "
           << xWest[pWestMin] << " " << yWest[pWestMin] << " " << sigmaWest[pWestMin] << endl;
    outFit << xEast[pEastMid] << " " << yEast[pEastMid] << " " << sigmaEast[pEastMid] << " "
           << xWest[pWestMid] << " " << yWest[pWestMid] << " " << sigmaWest[pWestMid] << endl;
    outFit << xEast[pEastMax] << " " << yEast[pEastMax] << " " << sigmaEast[pEastMax] << " "
    << xWest[pWestMax] << " " << yWest[pWestMax] << " " << sigmaWest[pWestMax] << endl;

  }

}
