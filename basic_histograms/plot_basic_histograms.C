// Usage: root[0] .x plot_basic_histograms.C("runNumber")
#include <cstdlib>

void plot_basic_histograms(TString runNumber)
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  // Style options
  gROOT->SetStyle("Plain");
  //gStyle->SetOptStat(11);
  gStyle->SetOptStat(0);
  gStyle->SetStatFontSize(0.030);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFontSize(0.05);
  //gStyle->SetTitleX(0.17);
  //gStyle->SetTitleAlign(13);
  gStyle->SetTitleOffset(0.90, "x");
  gStyle->SetTitleOffset(1.25, "y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetNdivisions(510,"X");
  gStyle->SetNdivisions(510,"Y");
  gStyle->SetPadLeftMargin(0.15); // 0.13
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.03);
  gStyle->SetPadBottomMargin(0.15); // 0.12
  gStyle->SetLabelSize(0.052,"X");
  gStyle->SetLabelSize(0.052,"Y");
  gStyle->SetTitleSize(0.052,"X");
  gStyle->SetTitleSize(0.052,"Y");
  gROOT->ForceStyle();
  gStyle->SetErrorX(0);

  // Open input file
  TString filenameIn;
  filenameIn  = TString(getenv("BASIC_HISTOGRAMS"))+"/basic_histograms_";
  filenameIn += runNumber;
  filenameIn += ".root";
  cout << "Processing ... " << filenameIn << endl;
  TFile *filein = new TFile(filenameIn);

  // Output file
  TString filenameOut;
  filenameOut  = TString(getenv("BASIC_HISTOGRAMS"))+"/basic_histograms_";
  filenameOut += runNumber;
  filenameOut += ".pdf";

  TString filenameOutFirst;
  filenameOutFirst  = filenameOut;
  filenameOutFirst += "(";

  TString filenameOutLast;
  filenameOutLast  = filenameOut;
  filenameOutLast += ")";

  c1 = new TCanvas("c1","c1");
  c1->Divide(2,2);

  c1_1->cd();
  c1_1->SetLogy(1);
  his1->SetXTitle("East MWPC Anode PADC");
  his1->GetXaxis()->CenterTitle();
  his1->GetXaxis()->SetRangeUser(100.0,500.0);
  his1->Draw();

  c1_2->cd();
  c1_2->SetLogy(1);
  his2->SetXTitle("West MWPC Anode PADC");
  his2->GetXaxis()->CenterTitle();
  his2->GetXaxis()->SetRangeUser(50.0,450.0);
  his2->Draw();

  c1_3->cd();
  c1_3->SetLogy(1);
  his3->SetXTitle("East Two-Fold Timing [ns]");
  his3->GetXaxis()->CenterTitle();
  his3->GetXaxis()->SetRangeUser(100.0,160.0);
  his3->Draw();

  c1_4->cd();
  c1_4->SetLogy(1);
  his4->SetXTitle("West Two-Fold Timing [ns]");
  his4->GetXaxis()->CenterTitle();
  his4->GetXaxis()->SetRangeUser(100.0,160.0);
  his4->Draw();

  c1->Print(filenameOutFirst);

  c2 = new TCanvas("c2","c2");
  c2->Divide(2,2);

  c2_1->cd();
  c2_1->SetLogy(1);
  his11->SetXTitle("East PMT #1 QADC Scintillator Triggers");
  his11->GetXaxis()->CenterTitle();
  his11->GetXaxis()->SetRangeUser(0.0,1000.0);
  his11->Draw();

  c2_2->cd();
  c2_2->SetLogy(1);
  his12->SetXTitle("East PMT #2 QADC Scintillator Triggers");
  his12->GetXaxis()->CenterTitle();
  his12->GetXaxis()->SetRangeUser(0.0,1000.0);
  his12->Draw();

  c2_3->cd();
  c2_3->SetLogy(1);
  his13->SetXTitle("East PMT #3 QADC Scintillator Triggers");
  his13->GetXaxis()->CenterTitle();
  his13->GetXaxis()->SetRangeUser(0.0,1000.0);
  his13->Draw();

  c2_4->cd();
  c2_4->SetLogy(1);
  his14->SetXTitle("East PMT #4 QADC Scintillator Triggers");
  his14->GetXaxis()->CenterTitle();
  his14->GetXaxis()->SetRangeUser(0.0,1000.0);
  his14->Draw();

  c2->Print(filenameOut);

  c3 = new TCanvas("c3","c3");
  c3->Divide(2,2);

  c3_1->cd();
  c3_1->SetLogy(1);
  his15->SetXTitle("West PMT #1 QADC Scintillator Triggers");
  his15->GetXaxis()->CenterTitle();
  his15->GetXaxis()->SetRangeUser(0.0,1000.0);
  his15->Draw();

  c3_2->cd();
  c3_2->SetLogy(1);
  his16->SetXTitle("West PMT #2 QADC Scintillator Triggers");
  his16->GetXaxis()->CenterTitle();
  his16->GetXaxis()->SetRangeUser(0.0,1000.0);
  his16->Draw();

  c3_3->cd();
  c3_3->SetLogy(1);
  his17->SetXTitle("West PMT #3 QADC Scintillator Triggers");
  his17->GetXaxis()->CenterTitle();
  his17->GetXaxis()->SetRangeUser(0.0,1000.0);
  his17->Draw();

  c3_4->cd();
  c3_4->SetLogy(1);
  his18->SetXTitle("West PMT #4 QADC Scintillator Triggers");
  his18->GetXaxis()->CenterTitle();
  his18->GetXaxis()->SetRangeUser(0.0,1000.0);
  his18->Draw();

  c3->Print(filenameOut);

  c4 = new TCanvas("c4","c4");
  c4->Divide(2,2);

  c4_1->cd();
  c4_1->SetLogy(1);
  his21->SetXTitle("East PMT #1 QADC Bi Pulser Triggers");
  his21->GetXaxis()->CenterTitle();
  his21->GetXaxis()->SetRangeUser(500.0,4000.0);
  his21->Draw();

  c4_2->cd();
  c4_2->SetLogy(1);
  his22->SetXTitle("East PMT #2 QADC Bi Pulser Triggers");
  his22->GetXaxis()->CenterTitle();
  his22->GetXaxis()->SetRangeUser(500.0,4000.0);
  his22->Draw();

  c4_3->cd();
  c4_3->SetLogy(1);
  his23->SetXTitle("East PMT #3 QADC Bi Pulser Triggers");
  his23->GetXaxis()->CenterTitle();
  his23->GetXaxis()->SetRangeUser(500.0,4000.0);
  his23->Draw();

  c4_4->cd();
  c4_4->SetLogy(1);
  his24->SetXTitle("East PMT #4 QADC Bi Pulser Triggers");
  his24->GetXaxis()->CenterTitle();
  his24->GetXaxis()->SetRangeUser(500.0,4000.0);
  his24->Draw();

  c4->Print(filenameOut);

  c5 = new TCanvas("c5","c5");
  c5->Divide(2,2);

  c5_1->cd();
  c5_1->SetLogy(1);
  his25->SetXTitle("West PMT #1 QADC Bi Pulser Triggers");
  his25->GetXaxis()->CenterTitle();
  his25->GetXaxis()->SetRangeUser(500.0,4000.0);
  his25->Draw();

  c5_2->cd();
  c5_2->SetLogy(1);
  his26->SetXTitle("West PMT #2 QADC Bi Pulser Triggers");
  his26->GetXaxis()->CenterTitle();
  his26->GetXaxis()->SetRangeUser(500.0,4000.0);
  his26->Draw();

  c5_3->cd();
  c5_3->SetLogy(1);
  his27->SetXTitle("West PMT #3 QADC Bi Pulser Triggers");
  his27->GetXaxis()->CenterTitle();
  his27->GetXaxis()->SetRangeUser(500.0,4000.0);
  his27->Draw();

  c5_4->cd();
  c5_4->SetLogy(1);
  his28->SetXTitle("West PMT #4 QADC Bi Pulser Triggers");
  his28->GetXaxis()->CenterTitle();
  his28->GetXaxis()->SetRangeUser(500.0,4000.0);
  his28->Draw();

  c5->Print(filenameOut);

  c6 = new TCanvas("c6","c6");
  c6->Divide(2,2);

  c6_1->cd();
  c6_1->SetLogy(1);
  his31->SetXTitle("East PMT #1 QADC LED Triggers");
  his31->GetXaxis()->CenterTitle();
  his31->GetXaxis()->SetRangeUser(500.0,4000.0);
  his31->Draw();

  c6_2->cd();
  c6_2->SetLogy(1);
  his32->SetXTitle("East PMT #2 QADC LED Triggers");
  his32->GetXaxis()->CenterTitle();
  his32->GetXaxis()->SetRangeUser(500.0,4000.0);
  his32->Draw();

  c6_3->cd();
  c6_3->SetLogy(1);
  his33->SetXTitle("East PMT #3 QADC LED Triggers");
  his33->GetXaxis()->CenterTitle();
  his33->GetXaxis()->SetRangeUser(500.0,4000.0);
  his33->Draw();

  c6_4->cd();
  c6_4->SetLogy(1);
  his34->SetXTitle("East PMT #4 QADC LED Triggers");
  his34->GetXaxis()->CenterTitle();
  his34->GetXaxis()->SetRangeUser(500.0,4000.0);
  his34->Draw();

  c6->Print(filenameOut);

  c7 = new TCanvas("c7","c7");
  c7->Divide(2,2);

  c7_1->cd();
  c7_1->SetLogy(1);
  his35->SetXTitle("West PMT #1 QADC LED Triggers");
  his35->GetXaxis()->CenterTitle();
  his35->GetXaxis()->SetRangeUser(500.0,4000.0);
  his35->Draw();

  c7_2->cd();
  c7_2->SetLogy(1);
  his36->SetXTitle("West PMT #2 QADC LED Triggers");
  his36->GetXaxis()->CenterTitle();
  his36->GetXaxis()->SetRangeUser(500.0,4000.0);
  his36->Draw();

  c7_3->cd();
  c7_3->SetLogy(1);
  his37->SetXTitle("West PMT #3 QADC LED Triggers");
  his37->GetXaxis()->CenterTitle();
  his37->GetXaxis()->SetRangeUser(500.0,4000.0);
  his37->Draw();

  c7_4->cd();
  c7_4->SetLogy(1);
  his38->SetXTitle("West PMT #4 QADC LED Triggers");
  his38->GetXaxis()->CenterTitle();
  his38->GetXaxis()->SetRangeUser(500.0,4000.0);
  his38->Draw();

  c7->Print(filenameOut);

  c8 = new TCanvas("c8","c8");
  c8->Divide(4,4);

  c8_1->cd();
  c8_1->SetLogy(1);
  his101->SetXTitle("East MWPC y = +76.20 mm");
  his101->GetXaxis()->CenterTitle();
  his101->Draw();

  c8_2->cd();
  c8_2->SetLogy(1);
  his102->SetXTitle("East MWPC y = +66.04 mm");
  his102->GetXaxis()->CenterTitle();
  his102->Draw();

  c8_3->cd();
  c8_3->SetLogy(1);
  his103->SetXTitle("East MWPC y = +55.88 mm");
  his103->GetXaxis()->CenterTitle();
  his103->Draw();

  c8_4->cd();
  c8_4->SetLogy(1);
  his104->SetXTitle("East MWPC y = +45.72 mm");
  his104->GetXaxis()->CenterTitle();
  his104->Draw();

  c8_5->cd();
  c8_5->SetLogy(1);
  his105->SetXTitle("East MWPC y = +35.56 mm");
  his105->GetXaxis()->CenterTitle();
  his105->Draw();

  c8_6->cd();
  c8_6->SetLogy(1);
  his106->SetXTitle("East MWPC y = +25.40 mm");
  his106->GetXaxis()->CenterTitle();
  his106->Draw();

  c8_7->cd();
  c8_7->SetLogy(1);
  his107->SetXTitle("East MWPC y = +15.24 mm");
  his107->GetXaxis()->CenterTitle();
  his107->Draw();

  c8_8->cd();
  c8_8->SetLogy(1);
  his108->SetXTitle("East MWPC y = +5.08 mm");
  his108->GetXaxis()->CenterTitle();
  his108->Draw();

  c8_9->cd();
  c8_9->SetLogy(1);
  his109->SetXTitle("East MWPC y = -5.08 mm");
  his109->GetXaxis()->CenterTitle();
  his109->Draw();

  c8_10->cd();
  c8_10->SetLogy(1);
  his110->SetXTitle("East MWPC y = -15.24 mm");
  his110->GetXaxis()->CenterTitle();
  his110->Draw();

  c8_11->cd();
  c8_11->SetLogy(1);
  his111->SetXTitle("East MWPC y = -25.40 mm");
  his111->GetXaxis()->CenterTitle();
  his111->Draw();

  c8_12->cd();
  c8_12->SetLogy(1);
  his112->SetXTitle("East MWPC y = -35.56 mm");
  his112->GetXaxis()->CenterTitle();
  his112->Draw();

  c8_13->cd();
  c8_13->SetLogy(1);
  his113->SetXTitle("East MWPC y = -45.72 mm");
  his113->GetXaxis()->CenterTitle();
  his113->Draw();

  c8_14->cd();
  c8_14->SetLogy(1);
  his114->SetXTitle("East MWPC y = -55.88 mm");
  his114->GetXaxis()->CenterTitle();
  his114->Draw();

  c8_15->cd();
  c8_15->SetLogy(1);
  his115->SetXTitle("East MWPC y = -66.04 mm");
  his115->GetXaxis()->CenterTitle();
  his115->Draw();

  c8_16->cd();
  c8_16->SetLogy(1);
  his116->SetXTitle("East MWPC y = -76.20 mm");
  his116->GetXaxis()->CenterTitle();
  his116->Draw();

  c8->Print(filenameOut);

  c9 = new TCanvas("c9","c9");
  c9->Divide(4,4);

  c9_1->cd();
  c9_1->SetLogy(1);
  his117->SetXTitle("East MWPC x = +76.20 mm");
  his117->GetXaxis()->CenterTitle();
  his117->Draw();

  c9_2->cd();
  c9_2->SetLogy(1);
  his118->SetXTitle("East MWPC x = +66.04 mm");
  his118->GetXaxis()->CenterTitle();
  his118->Draw();

  c9_3->cd();
  c9_3->SetLogy(1);
  his119->SetXTitle("East MWPC x = +55.88 mm");
  his119->GetXaxis()->CenterTitle();
  his119->Draw();

  c9_4->cd();
  c9_4->SetLogy(1);
  his120->SetXTitle("East MWPC x = +45.72 mm");
  his120->GetXaxis()->CenterTitle();
  his120->Draw();

  c9_5->cd();
  c9_5->SetLogy(1);
  his121->SetXTitle("East MWPC x = +35.56 mm");
  his121->GetXaxis()->CenterTitle();
  his121->Draw();

  c9_6->cd();
  c9_6->SetLogy(1);
  his122->SetXTitle("East MWPC x = +25.40 mm");
  his122->GetXaxis()->CenterTitle();
  his122->Draw();

  c9_7->cd();
  c9_7->SetLogy(1);
  his123->SetXTitle("East MWPC x = +15.24 mm");
  his123->GetXaxis()->CenterTitle();
  his123->Draw();

  c9_8->cd();
  c9_8->SetLogy(1);
  his124->SetXTitle("East MWPC x = +5.08 mm");
  his124->GetXaxis()->CenterTitle();
  his124->Draw();

  c9_9->cd();
  c9_9->SetLogy(1);
  his125->SetXTitle("East MWPC x = -5.08 mm");
  his125->GetXaxis()->CenterTitle();
  his125->Draw();

  c9_10->cd();
  c9_10->SetLogy(1);
  his126->SetXTitle("East MWPC x = -15.24 mm");
  his126->GetXaxis()->CenterTitle();
  his126->Draw();

  c9_11->cd();
  c9_11->SetLogy(1);
  his127->SetXTitle("East MWPC x = -25.40 mm");
  his127->GetXaxis()->CenterTitle();
  his127->Draw();

  c9_12->cd();
  c9_12->SetLogy(1);
  his128->SetXTitle("East MWPC x = -35.56 mm");
  his128->GetXaxis()->CenterTitle();
  his128->Draw();

  c9_13->cd();
  c9_13->SetLogy(1);
  his129->SetXTitle("East MWPC x = -45.72 mm");
  his129->GetXaxis()->CenterTitle();
  his129->Draw();

  c9_14->cd();
  c9_14->SetLogy(1);
  his130->SetXTitle("East MWPC x = -55.88 mm");
  his130->GetXaxis()->CenterTitle();
  his130->Draw();

  c9_15->cd();
  c9_15->SetLogy(1);
  his131->SetXTitle("East MWPC x = -66.04 mm");
  his131->GetXaxis()->CenterTitle();
  his131->Draw();

  c9_16->cd();
  c9_16->SetLogy(1);
  his132->SetXTitle("East MWPC x = -76.20 mm");
  his132->GetXaxis()->CenterTitle();
  his132->Draw();

  c9->Print(filenameOut);

  c10 = new TCanvas("c10","c10");
  c10->Divide(4,4);

  c10_1->cd();
  c10_1->SetLogy(1);
  his201->SetXTitle("West MWPC y = +76.20 mm");
  his201->GetXaxis()->CenterTitle();
  his201->Draw();

  c10_2->cd();
  c10_2->SetLogy(1);
  his202->SetXTitle("West MWPC y = +66.04 mm");
  his202->GetXaxis()->CenterTitle();
  his202->Draw();

  c10_3->cd();
  c10_3->SetLogy(1);
  his203->SetXTitle("West MWPC y = +55.88 mm");
  his203->GetXaxis()->CenterTitle();
  his203->Draw();

  c10_4->cd();
  c10_4->SetLogy(1);
  his204->SetXTitle("West MWPC y = +45.72 mm");
  his204->GetXaxis()->CenterTitle();
  his204->Draw();

  c10_5->cd();
  c10_5->SetLogy(1);
  his205->SetXTitle("West MWPC y = +35.56 mm");
  his205->GetXaxis()->CenterTitle();
  his205->Draw();

  c10_6->cd();
  c10_6->SetLogy(1);
  his206->SetXTitle("West MWPC y = +25.40 mm");
  his206->GetXaxis()->CenterTitle();
  his206->Draw();

  c10_7->cd();
  c10_7->SetLogy(1);
  his207->SetXTitle("West MWPC y = +15.24 mm");
  his207->GetXaxis()->CenterTitle();
  his207->Draw();

  c10_8->cd();
  c10_8->SetLogy(1);
  his208->SetXTitle("West MWPC y = +5.08 mm");
  his208->GetXaxis()->CenterTitle();
  his208->Draw();

  c10_9->cd();
  c10_9->SetLogy(1);
  his209->SetXTitle("West MWPC y = -5.08 mm");
  his209->GetXaxis()->CenterTitle();
  his209->Draw();

  c10_10->cd();
  c10_10->SetLogy(1);
  his210->SetXTitle("West MWPC y = -15.24 mm");
  his210->GetXaxis()->CenterTitle();
  his210->Draw();

  c10_11->cd();
  c10_11->SetLogy(1);
  his211->SetXTitle("West MWPC y = -25.40 mm");
  his211->GetXaxis()->CenterTitle();
  his211->Draw();

  c10_12->cd();
  c10_12->SetLogy(1);
  his212->SetXTitle("West MWPC y = -35.56 mm");
  his212->GetXaxis()->CenterTitle();
  his212->Draw();

  c10_13->cd();
  c10_13->SetLogy(1);
  his213->SetXTitle("West MWPC y = -45.72 mm");
  his213->GetXaxis()->CenterTitle();
  his213->Draw();

  c10_14->cd();
  c10_14->SetLogy(1);
  his214->SetXTitle("West MWPC y = -55.88 mm");
  his214->GetXaxis()->CenterTitle();
  his214->Draw();

  c10_15->cd();
  c10_15->SetLogy(1);
  his215->SetXTitle("West MWPC y = -66.04 mm");
  his215->GetXaxis()->CenterTitle();
  his215->Draw();

  c10_16->cd();
  c10_16->SetLogy(1);
  his216->SetXTitle("West MWPC y = -76.20 mm");
  his216->GetXaxis()->CenterTitle();
  his216->Draw();

  c10->Print(filenameOut);

  c11 = new TCanvas("c11","c11");
  c11->Divide(4,4);

  c11_1->cd();
  c11_1->SetLogy(1);
  his217->SetXTitle("West MWPC x = -76.20 mm");
  his217->GetXaxis()->CenterTitle();
  his217->Draw();

  c11_2->cd();
  c11_2->SetLogy(1);
  his218->SetXTitle("West MWPC x = -66.04 mm");
  his218->GetXaxis()->CenterTitle();
  his218->Draw();

  c11_3->cd();
  c11_3->SetLogy(1);
  his219->SetXTitle("West MWPC x = -55.88 mm");
  his219->GetXaxis()->CenterTitle();
  his219->Draw();

  c11_4->cd();
  c11_4->SetLogy(1);
  his220->SetXTitle("West MWPC x = -45.72 mm");
  his220->GetXaxis()->CenterTitle();
  his220->Draw();

  c11_5->cd();
  c11_5->SetLogy(1);
  his221->SetXTitle("West MWPC x = -35.56 mm");
  his221->GetXaxis()->CenterTitle();
  his221->Draw();

  c11_6->cd();
  c11_6->SetLogy(1);
  his222->SetXTitle("West MWPC x = -25.40 mm");
  his222->GetXaxis()->CenterTitle();
  his222->Draw();

  c11_7->cd();
  c11_7->SetLogy(1);
  his223->SetXTitle("West MWPC x = -15.24 mm");
  his223->GetXaxis()->CenterTitle();
  his223->Draw();

  c11_8->cd();
  c11_8->SetLogy(1);
  his224->SetXTitle("West MWPC x = -5.08 mm");
  his224->GetXaxis()->CenterTitle();
  his224->Draw();

  c11_9->cd();
  c11_9->SetLogy(1);
  his225->SetXTitle("West MWPC x = +5.08 mm");
  his225->GetXaxis()->CenterTitle();
  his225->Draw();

  c11_10->cd();
  c11_10->SetLogy(1);
  his226->SetXTitle("West MWPC x = +15.24 mm");
  his226->GetXaxis()->CenterTitle();
  his226->Draw();

  c11_11->cd();
  c11_11->SetLogy(1);
  his227->SetXTitle("West MWPC x = +25.40 mm");
  his227->GetXaxis()->CenterTitle();
  his227->Draw();

  c11_12->cd();
  c11_12->SetLogy(1);
  his228->SetXTitle("West MWPC x = +35.56 mm");
  his228->GetXaxis()->CenterTitle();
  his228->Draw();

  c11_13->cd();
  c11_13->SetLogy(1);
  his229->SetXTitle("West MWPC x = +45.72 mm");
  his229->GetXaxis()->CenterTitle();
  his229->Draw();

  c11_14->cd();
  c11_14->SetLogy(1);
  his230->SetXTitle("West MWPC x = +55.88 mm");
  his230->GetXaxis()->CenterTitle();
  his230->Draw();

  c11_15->cd();
  c11_15->SetLogy(1);
  his231->SetXTitle("West MWPC x = +66.04 mm");
  his231->GetXaxis()->CenterTitle();
  his231->Draw();

  c11_16->cd();
  c11_16->SetLogy(1);
  his232->SetXTitle("West MWPC x = +76.20 mm");
  his232->GetXaxis()->CenterTitle();
  his232->Draw();

  c11->Print(filenameOut);

  c12 = new TCanvas("c12","c12");
  c12->Divide(2,2);

  c12_1->cd();
  c12_1->SetLogy(0);
  his301->SetXTitle("Scintillator Triggers Timing [s]");
  his301->GetXaxis()->CenterTitle();
  //his301->GetXaxis()->SetRangeUser(100.0,300.0);
  his301->Draw();

  c12_2->cd();
  c12_2->SetLogy(0);
  his302->SetXTitle("Gate Valve UCN Monitor Triggers Timing [s]");
  his302->GetXaxis()->CenterTitle();
  //his302->GetXaxis()->SetRangeUser(50.0,250.0);
  his302->Draw();

  c12_3->cd();
  c12_3->SetLogy(0);
  his304->SetXTitle("AFP Fe Foil UCN Monitor Triggers Timing [s]");
  his304->GetXaxis()->CenterTitle();
  //his304->GetXaxis()->SetRangeUser(100.0,160.0);
  his304->Draw();

  c12_4->cd();
  c12_4->SetLogy(0);
  his305->SetXTitle("SCS UCN Monitor Triggers Timing [s]");
  his305->GetXaxis()->CenterTitle();
  //his305->GetXaxis()->SetRangeUser(100.0,160.0);
  his305->Draw();

  c12->Print(filenameOut);

  c13 = new TCanvas("c13","c13");
  c13->Divide(2,2);

  c13_1->cd();
  c13_1->SetLogy(1);
  his401->SetXTitle("East Top Veto QADC");
  his401->GetXaxis()->CenterTitle();
  //his401->GetXaxis()->SetRangeUser(100.0,300.0);
  his401->Draw();

  c13_2->cd();
  c13_2->SetLogy(1);
  his402->SetXTitle("East Top Veto TDC");
  his402->GetXaxis()->CenterTitle();
  //his402->GetXaxis()->SetRangeUser(50.0,250.0);
  his402->Draw();

  c13_3->cd();
  c13_3->SetLogy(1);
  his403->SetXTitle("East Drift Tube Veto TAC");
  his403->GetXaxis()->CenterTitle();
  //his403->GetXaxis()->SetRangeUser(100.0,160.0);
  his403->Draw();

  c13_4->cd();
  c13_4->SetLogy(1);
  his404->SetXTitle("West Drift Tube Veto TAC");
  his404->GetXaxis()->CenterTitle();
  //his404->GetXaxis()->SetRangeUser(100.0,160.0);
  his404->Draw();

  c13->Print(filenameOut);

  c14 = new TCanvas("c14","c14");
  c14->Divide(2,2);

  c14_1->cd();
  c14_1->SetLogy(1);
  his405->SetXTitle("East Backing Veto QADC");
  his405->GetXaxis()->CenterTitle();
  his405->GetXaxis()->SetRangeUser(0.0,1500.0);
  his405->Draw();

  c14_2->cd();
  c14_2->SetLogy(1);
  his406->SetXTitle("East Backing Veto TDC");
  his406->GetXaxis()->CenterTitle();
  //his406->GetXaxis()->SetRangeUser(50.0,250.0);
  his406->Draw();

  c14_3->cd();
  c14_3->SetLogy(1);
  his407->SetXTitle("West Backing Veto QADC");
  his407->GetXaxis()->CenterTitle();
  his407->GetXaxis()->SetRangeUser(0.0,1500.0);
  his407->Draw();

  c14_4->cd();
  c14_4->SetLogy(1);
  his408->SetXTitle("West Backing Veto TDC");
  his408->GetXaxis()->CenterTitle();
  //his408->GetXaxis()->SetRangeUser(100.0,160.0);
  his408->Draw();

  c14->Print(filenameOutLast);

}
