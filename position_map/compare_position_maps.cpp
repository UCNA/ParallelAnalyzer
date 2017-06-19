

#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ios>
#include "positionMapHandler.hh"

#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TPaveText.h>
#include <TMath.h>

using namespace std;

void plot_position_map(TString fitType1, TString fitType2, int iRunPeriod, double binWidth);


int main(int argc, char* argv[]) {
  TString fitType1 = TString(argv[1]);
  TString fitType2 = TString(argv[2]);
  int runPeriod = atoi(argv[3]);
  double binWidth = atof(argv[4]);

  if ( argc!=5 ) {
    std::cout << "USAGE: ./compare_position_maps.exe [fitType1] [fitType2] [XePeriod] [binWidth]\n"
	      << "fitType = {ave, peak, endpoint}\n";
    exit(0);
  }

  plot_position_map(fitType1, fitType2, runPeriod, binWidth);
  

}

void plot_position_map( TString fitType1, TString fitType2, int XePeriod, double binWidth=5.)
{

  if ( fitType1!=TString("ave") && fitType1!=TString("peak") && fitType1!=TString("endpoint") ) { 
    std::cout << "bad fit1 choice\n"; 
    exit(0); 
  }
  if ( fitType2!=TString("ave") && fitType2!=TString("peak") && fitType2!=TString("endpoint") ) { 
    std::cout << "bad fit2 choice\n"; 
    exit(0); 
  }
			 

  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  // Style options
  //gROOT->SetStyle("Plain");
  Int_t palette[10];
  palette[0] = 0;
  palette[1] = 19;
  palette[2] = 18;
  palette[3] = 17;
  palette[4] = 16;
  palette[5] = 15;
  palette[6] = 14;
  palette[7] = 13;
  palette[8] = 12;
  palette[9] = 1;
  gStyle->SetPalette(10,palette);  // z-axis color scale for 2D histograms
  //gStyle->SetPalette(2);  // z-axis color scale for 2D histograms
  //gStyle->SetOptStat(11);
  gStyle->SetOptStat(0);
  gStyle->SetStatFontSize(0.030);
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFontSize(0.05);
  //gStyle->SetTitleX(0.17);
  //gStyle->SetTitleAlign(13);
  gStyle->SetTitleOffset(0.80, "x");
  gStyle->SetTitleOffset(1.30, "y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetNdivisions(510,"X");
  gStyle->SetNdivisions(510,"Y");
  gStyle->SetNdivisions(9,"Z");
  gStyle->SetPadLeftMargin(0.13); // 0.13
  gStyle->SetPadRightMargin(0.15); // 0.04
  //gStyle->SetPadBottomMargin(0.37); // 0.30
  gStyle->SetLabelSize(0.045, "X");
  gStyle->SetLabelSize(0.045, "Y");
  gStyle->SetLabelSize(0.045, "Z");
  gStyle->SetLabelOffset(0.00, "X");
  gStyle->SetLabelOffset(0.01, "Y");
  gStyle->SetTitleSize(0.050, "X");
  gStyle->SetTitleSize(0.050, "Y");
  //gROOT->ForceStyle();

  
  PositionMap posmap1(binWidth, 50.);
  posmap1.setRCflag(false); // plotting non RC maps
  posmap1.readPositionMap(XePeriod,fitType1);

  PositionMap posmap2(binWidth, 50.);
  posmap2.setRCflag(false); // plotting non RC maps
  posmap2.readPositionMap(XePeriod,fitType2);

  
  PositionMap plotMap(0.5,50.);
  plotMap.setRCflag(false); // plotting non RC maps

  int nPosBinsXY = plotMap.getNbinsXY();
  double xyMin = plotMap.getBinLower(0);
  double xyMax = plotMap.getBinUpper(nPosBinsXY-1);  

  cout << nPosBinsXY << " " << xyMin << " " << xyMax << endl;

  // Histograms
  TH2D *hisE0 = new TH2D("E0", "", nPosBinsXY,xyMin,xyMax, nPosBinsXY,xyMin,xyMax);
  TH2D *hisE1 = new TH2D("E1", "", nPosBinsXY,xyMin,xyMax, nPosBinsXY,xyMin,xyMax);
  TH2D *hisE2 = new TH2D("E2", "", nPosBinsXY,xyMin,xyMax, nPosBinsXY,xyMin,xyMax);
  TH2D *hisE3 = new TH2D("E3", "", nPosBinsXY,xyMin,xyMax, nPosBinsXY,xyMin,xyMax);
  TH2D *hisW0 = new TH2D("W0", "", nPosBinsXY,xyMin,xyMax, nPosBinsXY,xyMin,xyMax);
  TH2D *hisW1 = new TH2D("W1", "", nPosBinsXY,xyMin,xyMax, nPosBinsXY,xyMin,xyMax);
  TH2D *hisW2 = new TH2D("W2", "", nPosBinsXY,xyMin,xyMax, nPosBinsXY,xyMin,xyMax);
  TH2D *hisW3 = new TH2D("W3", "", nPosBinsXY,xyMin,xyMax, nPosBinsXY,xyMin,xyMax);

  //Fill all the bins
  for (int xb=0; xb<nPosBinsXY; xb++) {

    for (int yb=0; yb<nPosBinsXY; yb++) {

      if (TMath::Sqrt(plotMap.getBinCenter(xb)*plotMap.getBinCenter(xb)+plotMap.getBinCenter(yb)*plotMap.getBinCenter(yb)) <= 50.) {

	std::vector <Double_t> eta1 = posmap1.getInterpolatedEta(plotMap.getBinCenter(xb),plotMap.getBinCenter(yb),plotMap.getBinCenter(xb),plotMap.getBinCenter(yb));
	std::vector <Double_t> eta2 = posmap2.getInterpolatedEta(plotMap.getBinCenter(xb),plotMap.getBinCenter(yb),plotMap.getBinCenter(xb),plotMap.getBinCenter(yb));
	//std::vector <Double_t> eta = posmap.getInterpolatedEta(0., 0., 0., 0.);

	//cout << eta[0] << endl;
        hisE0->Fill(plotMap.getBinCenter(xb),plotMap.getBinCenter(yb), eta2[0] / eta1[0]);
        hisE1->Fill(plotMap.getBinCenter(xb),plotMap.getBinCenter(yb), eta2[1] / eta1[1]);
        hisE2->Fill(plotMap.getBinCenter(xb),plotMap.getBinCenter(yb), eta2[2] / eta1[2]);
        hisE3->Fill(plotMap.getBinCenter(xb),plotMap.getBinCenter(yb), eta2[3] / eta1[3]);
        hisW0->Fill(plotMap.getBinCenter(xb),plotMap.getBinCenter(yb), eta2[4] / eta1[4]);
        hisW1->Fill(plotMap.getBinCenter(xb),plotMap.getBinCenter(yb), eta2[5] / eta1[5]);
        hisW2->Fill(plotMap.getBinCenter(xb),plotMap.getBinCenter(yb), eta2[6] / eta1[6]);
        hisW3->Fill(plotMap.getBinCenter(xb),plotMap.getBinCenter(yb), eta2[7] / eta1[7]);
      }
    }
  }
  

  // Output PDF file
  TString filenameOut;
  Char_t runPeriodString[10];
  sprintf(runPeriodString,"%i",XePeriod);
  filenameOut = "posmapComp_"+TString(runPeriodString)+TString::Format("_%0.1fmm_%s_vs_%s.pdf",binWidth,fitType1.Data(),fitType2.Data());
  
  /*if (iRunPeriod == 1) filenameOut  = "position_map_1_RC_123.pdf";
  if (iRunPeriod == 2) filenameOut  = "position_map_2_RC_123.pdf";
  if (iRunPeriod == 3) filenameOut  = "position_map_3_RC_123.pdf";
  if (iRunPeriod == 4) filenameOut  = "position_map_4_RC_123.pdf";
  if (iRunPeriod == 5) filenameOut  = "position_map_5_RC_123.pdf";
  if (iRunPeriod == 6) filenameOut  = "position_map_6_RC_123.pdf";
  if (iRunPeriod == 7) filenameOut  = "position_map_7_RC_123.pdf";*/

  TString filenameOutFirst;
  filenameOutFirst  = filenameOut;
  filenameOutFirst += "(";

  TString filenameOutLast;
  filenameOutLast  = filenameOut;
  filenameOutLast += ")";

  TCanvas * c0 = new TCanvas("c0", "canvas", 1000, 900);
  c0->Divide(2,2);

  // East 3
  c0->cd(1);
  hisE3->Draw("colz");
  hisE3->SetXTitle("x [mm]");
  hisE3->SetYTitle("y [mm]");
  hisE3->GetXaxis()->CenterTitle();
  hisE3->GetYaxis()->CenterTitle();
  hisE3->GetYaxis()->CenterTitle();
  hisE3->SetAxisRange(0.75,1.25,"Z");

  /*
  TEllipse *ell = new TEllipse(0,0,45,45);
  ell->SetFillStyle(0);
  ell->SetLineColor(1);
  ell->SetLineStyle(2);
  ell->SetLineWidth(3);
  ell->Draw("same");
  */

  Double_t x1_text =  -47;
  Double_t y1_text =  45;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.052);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("E3");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  // East 2
  c0->cd(2);
  hisE2->Draw("colz");
  hisE2->SetXTitle("x [mm]");
  hisE2->SetYTitle("y [mm]");
  hisE2->GetXaxis()->CenterTitle();
  hisE2->GetYaxis()->CenterTitle();
  hisE2->GetYaxis()->CenterTitle();
  hisE2->SetAxisRange(0.75,1.25,"Z");

  /*
  TEllipse *ell = new TEllipse(0,0,45,45);
  ell->SetFillStyle(0);
  ell->SetLineColor(1);
  ell->SetLineStyle(2);
  ell->SetLineWidth(3);
  ell->Draw("same");
  */

  x1_text = 40;
  y1_text = 45;

  pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.052);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("E2");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  //c0->Print(filenameOutFirst);

  //TCanvas * c1 = new TCanvas("c1", "canvas");
  //c1->Divide(2,1);
  //c1->SetLogy(0);

  // East 0
  c0->cd(3);
  hisE0->Draw("colz");
  hisE0->SetXTitle("x [mm]");
  hisE0->SetYTitle("y [mm]");
  hisE0->GetXaxis()->CenterTitle();
  hisE0->GetYaxis()->CenterTitle();
  hisE0->GetYaxis()->CenterTitle();
  hisE0->SetAxisRange(0.75,1.25,"Z");

  /*
  TEllipse *ell = new TEllipse(0,0,45,45);
  ell->SetFillStyle(0);
  ell->SetLineColor(1);
  ell->SetLineStyle(2);
  ell->SetLineWidth(3);
  ell->Draw("same");
  */

  x1_text = -47;
  y1_text = -45;

  pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.052);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("E0");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  // East 1
  c0->cd(4);
  hisE1->Draw("colz");
  hisE1->SetXTitle("x [mm]");
  hisE1->SetYTitle("y [mm]");
  hisE1->GetXaxis()->CenterTitle();
  hisE1->GetYaxis()->CenterTitle();
  hisE1->GetYaxis()->CenterTitle();
  hisE1->SetAxisRange(0.75,1.25,"Z");

  /*
  TEllipse *ell = new TEllipse(0,0,45,45);
  ell->SetFillStyle(0);
  ell->SetLineColor(1);
  ell->SetLineStyle(2);
  ell->SetLineWidth(3);
  ell->Draw("same");
  */

  x1_text =  40;
  y1_text = -45;

  pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.052);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("E1");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  c0->Print(filenameOutFirst);

  TCanvas * c2 = new TCanvas("c2", "canvas", 1000, 900);
  c2->Divide(2,2);

  // West 1
  c2->cd(1);
  hisW1->Draw("colz");
  hisW1->SetXTitle("x [mm]");
  hisW1->SetYTitle("y [mm]");
  hisW1->GetXaxis()->CenterTitle();
  hisW1->GetYaxis()->CenterTitle();
  hisW1->GetYaxis()->CenterTitle();
  hisW1->SetAxisRange(0.75,1.25,"Z");

  /*
  TEllipse *ell = new TEllipse(0,0,45,45);
  ell->SetFillStyle(0);
  ell->SetLineColor(1);
  ell->SetLineStyle(2);
  ell->SetLineWidth(3);
  ell->Draw("same");
  */

  x1_text = -47;
  y1_text = 45;

  pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.052);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("W1");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  // West 0
  c2->cd(2);
  hisW0->Draw("colz");
  hisW0->SetXTitle("x [mm]");
  hisW0->SetYTitle("y [mm]");
  hisW0->GetXaxis()->CenterTitle();
  hisW0->GetYaxis()->CenterTitle();
  hisW0->GetYaxis()->CenterTitle();
  hisW0->SetAxisRange(0.75,1.25,"Z");

  /*
  TEllipse *ell = new TEllipse(0,0,45,45);
  ell->SetFillStyle(0);
  ell->SetLineColor(1);
  ell->SetLineStyle(2);
  ell->SetLineWidth(3);
  ell->Draw("same");
  */

  x1_text =  40;
  y1_text = 45;

  pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.052);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("W0");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  //c2->Print(filenameOut);

  //TCanvas * c3 = new TCanvas("c3", "canvas");
  //c3->Divide(2,1);
  //c3->SetLogy(0);

  // West 2
  c2->cd(3);
  hisW2->Draw("colz");
  hisW2->SetXTitle("x [mm]");
  hisW2->SetYTitle("y [mm]");
  hisW2->GetXaxis()->CenterTitle();
  hisW2->GetYaxis()->CenterTitle();
  hisW2->GetYaxis()->CenterTitle();
  hisW2->SetAxisRange(0.75,1.25,"Z");

  /*
  TEllipse *ell = new TEllipse(0,0,45,45);
  ell->SetFillStyle(0);
  ell->SetLineColor(1);
  ell->SetLineStyle(2);
  ell->SetLineWidth(3);
  ell->Draw("same");
  */

  x1_text =  -47;
  y1_text =  -45;

  pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.052);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("W2");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  // West 3
  c2->cd(4);
  hisW3->Draw("colz");
  hisW3->SetXTitle("x [mm]");
  hisW3->SetYTitle("y [mm]");
  hisW3->GetXaxis()->CenterTitle();
  hisW3->GetYaxis()->CenterTitle();
  hisW3->GetYaxis()->CenterTitle();
  hisW3->SetAxisRange(0.75,1.25,"Z");

  /*
  TEllipse *ell = new TEllipse(0,0,45,45);
  ell->SetFillStyle(0);
  ell->SetLineColor(1);
  ell->SetLineStyle(2);
  ell->SetLineWidth(3);
  ell->Draw("same");
  */

  x1_text = 40;
  y1_text = -45;

  pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.052);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("W3");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  // Redraw axis covered up by gray band
  //gPad->RedrawAxis();

  c2->Print(filenameOutLast);

  delete pt1; delete c0;  delete c2; //delete c3;
  delete hisE0; delete hisE1; delete hisE2; delete hisE3;
  delete hisW0; delete hisW1; delete hisW2; delete hisW3;
  
}
