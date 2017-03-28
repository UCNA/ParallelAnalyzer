// Usage: root[0] .x plot_source_peaks.C("runNumber")

void plot_source_peaks(TString runNumber)
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

  // Open input ntuple
  TString filenameIn;
  filenameIn  = TString(getenv("SOURCE_PEAKS"))+TString("/source_peaks_");
  filenameIn += runNumber;
  filenameIn += ".root";
  cout << "Processing ... " << filenameIn << endl;
  TFile *filein = new TFile(filenameIn);

  // Open source list file
  int nSources;
  string sourceName[3];
  bool useSource[3] = {false,false,false};
  TString filenameList;
  filenameList  = TString(getenv("SOURCE_LIST"))+TString("/source_list_");
  filenameList += runNumber;
  filenameList += ".dat";
  ifstream fileList(filenameList);
  fileList >> nSources;
  cout << "... nSources: " << nSources << endl;
  for (int n=0; n<nSources;n++) {
    fileList >> sourceName[n];
    if (sourceName[n]=="Ce" || sourceName[n]=="Sn" || sourceName[n]=="Bi") useSource[n]=true;
  }

  // Output file
  TString filenameOut;
  filenameOut  = TString(getenv("SOURCE_PEAKS"))+TString("/source_peaks_");
  filenameOut += runNumber;
  filenameOut += ".pdf";

  TString filenameOutFirst;
  filenameOutFirst  = filenameOut;
  filenameOutFirst += "(";

  TString filenameOutLast;
  filenameOutLast  = filenameOut;
  filenameOutLast += ")";

  // Source #1: East
  if (useSource[0]) {
    double upperRange = 2000.;
    if (sourceName[0]=="Bi") upperRange = 3000.;
    c1 = new TCanvas("c1","c1");
    c1->Divide(2,2);

    c1_1->cd();
    c1_1->SetLogy(0);
    his1_E0->SetXTitle("East PMT 1");
    his1_E0->GetXaxis()->CenterTitle();
    his1_E0->GetXaxis()->SetRangeUser(-50.0,upperRange);
    his1_E0->Draw();
    
    c1_2->cd();
    c1_2->SetLogy(0);
    his1_E1->SetXTitle("East PMT 2");
    his1_E1->GetXaxis()->CenterTitle();
    his1_E1->GetXaxis()->SetRangeUser(-50.0,upperRange);
    his1_E1->Draw();
    
    c1_3->cd();
    c1_3->SetLogy(0);
    his1_E2->SetXTitle("East PMT 3");
    his1_E2->GetXaxis()->CenterTitle();
    his1_E2->GetXaxis()->SetRangeUser(-50.0,upperRange);
    his1_E2->Draw();
    
    c1_4->cd();
    c1_4->SetLogy(0);
    his1_E3->SetXTitle("East PMT 4");
    his1_E3->GetXaxis()->CenterTitle();
    his1_E3->GetXaxis()->SetRangeUser(-50.0,upperRange);
    his1_E3->Draw();
    
    c1->Print(filenameOutFirst);
    
    // Source #1: West
    c2 = new TCanvas("c2","c2");
    c2->Divide(2,2);
    
    c2_1->cd();
    c2_1->SetLogy(0);
    his1_W0->SetXTitle("West PMT 1");
    his1_W0->GetXaxis()->CenterTitle();
    his1_W0->GetXaxis()->SetRangeUser(-50.0,upperRange);
    his1_W0->Draw();
    
    c2_2->cd();
    c2_2->SetLogy(0);
    his1_W1->SetXTitle("West PMT 2");
    his1_W1->GetXaxis()->CenterTitle();
    his1_W1->GetXaxis()->SetRangeUser(-50.0,upperRange);
    his1_W1->Draw();
    
    c2_3->cd();
    c2_3->SetLogy(0);
    his1_W2->SetXTitle("West PMT 3");
    his1_W2->GetXaxis()->CenterTitle();
    his1_W2->GetXaxis()->SetRangeUser(-50.0,upperRange);
    his1_W2->Draw();

    c2_4->cd();
    c2_4->SetLogy(0);
    his1_W3->SetXTitle("West PMT 4");
    his1_W3->GetXaxis()->CenterTitle();
    his1_W3->GetXaxis()->SetRangeUser(-50.0,upperRange+980.);
    his1_W3->Draw();
    
    if (nSources > 1 && (useSource[1] || useSource[2])) c2->Print(filenameOut);
    else c2->Print(filenameOutLast);
  }
  if (nSources > 1 && useSource[1]) {
    double upperRange = 2000.;
    if (sourceName[1]=="Bi") upperRange = 3000.;
      
    // Source #2: East
    c3 = new TCanvas("c3","c3");
    c3->Divide(2,2);
    
    c3_1->cd();
    c3_1->SetLogy(0);
    his2_E0->SetXTitle("East PMT 1");
    his2_E0->GetXaxis()->CenterTitle();
    his2_E0->GetXaxis()->SetRangeUser(-50.0,upperRange);
    his2_E0->Draw();
    
    c3_2->cd();
    c3_2->SetLogy(0);
    his2_E1->SetXTitle("East PMT 2");
    his2_E1->GetXaxis()->CenterTitle();
    his2_E1->GetXaxis()->SetRangeUser(-50.0,upperRange);
    his2_E1->Draw();
    
    c3_3->cd();
    c3_3->SetLogy(0);
    his2_E2->SetXTitle("East PMT 3");
    his2_E2->GetXaxis()->CenterTitle();
    his2_E2->GetXaxis()->SetRangeUser(-50.0,upperRange);
    his2_E2->Draw();
    
    c3_4->cd();
    c3_4->SetLogy(0);
    his2_E3->SetXTitle("East PMT 4");
    his2_E3->GetXaxis()->CenterTitle();
    his2_E3->GetXaxis()->SetRangeUser(-50.0,upperRange);
    his2_E3->Draw();
    
    if (useSource[0]) c3->Print(filenameOut);
    else c3->Print(filenameOutFirst);
    
    // Source #2: West
    c4 = new TCanvas("c4","c4");
    c4->Divide(2,2);
    
    c4_1->cd();
    c4_1->SetLogy(0);
    his2_W0->SetXTitle("West PMT 1");
    his2_W0->GetXaxis()->CenterTitle();
    his2_W0->GetXaxis()->SetRangeUser(-50.0,upperRange);
    his2_W0->Draw();
    
    c4_2->cd();
    c4_2->SetLogy(0);
    his2_W1->SetXTitle("West PMT 2");
    his2_W1->GetXaxis()->CenterTitle();
    his2_W1->GetXaxis()->SetRangeUser(-50.0,upperRange);
    his2_W1->Draw();
    
    c4_3->cd();
    c4_3->SetLogy(0);
    his2_W2->SetXTitle("West PMT 3");
    his2_W2->GetXaxis()->CenterTitle();
    his2_W2->GetXaxis()->SetRangeUser(-50.0,upperRange);
    his2_W2->Draw();
    
    c4_4->cd();
    c4_4->SetLogy(0);
    his2_W3->SetXTitle("West PMT 4");
    his2_W3->GetXaxis()->CenterTitle();
    his2_W3->GetXaxis()->SetRangeUser(-50.0,upperRange+980.);
    his2_W3->Draw();
    
    if (nSources > 2 && useSource[2]) c4->Print(filenameOut);
    else c4->Print(filenameOutLast);

  }

  if (nSources > 2 && useSource[2]) {
    double upperRange = 2000.;
    if (sourceName[2]=="Bi") upperRange = 3000.;

    // Source #3: East
    c5 = new TCanvas("c5","c5");
    c5->Divide(2,2);
    
    c5_1->cd();
    c5_1->SetLogy(0);
    his3_E0->SetXTitle("East PMT 1");
    his3_E0->GetXaxis()->CenterTitle();
    his3_E0->GetXaxis()->SetRangeUser(-50.0,upperRange);
    his3_E0->Draw();
    
    c5_2->cd();
    c5_2->SetLogy(0);
    his3_E1->SetXTitle("East PMT 2");
    his3_E1->GetXaxis()->CenterTitle();
    his3_E1->GetXaxis()->SetRangeUser(-50.0,upperRange);
    his3_E1->Draw();
    
    c5_3->cd();
    c5_3->SetLogy(0);
    his3_E2->SetXTitle("East PMT 3");
    his3_E2->GetXaxis()->CenterTitle();
    his3_E2->GetXaxis()->SetRangeUser(-50.0,upperRange);
    his3_E2->Draw();
    
    c5_4->cd();
    c5_4->SetLogy(0);
    his3_E3->SetXTitle("East PMT 4");
    his3_E3->GetXaxis()->CenterTitle();
    his3_E3->GetXaxis()->SetRangeUser(-50.0,upperRange);
    his3_E3->Draw();
    
    if (useSource[0] || useSource[1]) c5->Print(filenameOut);
    else c5->Print(filenameOutFirst);
    
    // Source #3: West
    c6 = new TCanvas("c6","c6");
    c6->Divide(2,2);
    
    c6_1->cd();
    c6_1->SetLogy(0);
    his3_W0->SetXTitle("West PMT 1");
    his3_W0->GetXaxis()->CenterTitle();
    his3_W0->GetXaxis()->SetRangeUser(-50.0,upperRange);
    his3_W0->Draw();
    
    c6_2->cd();
    c6_2->SetLogy(0);
    his3_W1->SetXTitle("West PMT 2");
    his3_W1->GetXaxis()->CenterTitle();
    his3_W1->GetXaxis()->SetRangeUser(-50.0,upperRange);
    his3_W1->Draw();
    
    c6_3->cd();
    c6_3->SetLogy(0);
    his3_W2->SetXTitle("West PMT 3");
    his3_W2->GetXaxis()->CenterTitle();
    his3_W2->GetXaxis()->SetRangeUser(-50.0,upperRange);
    his3_W2->Draw();
    
    c6_4->cd();
    c6_4->SetLogy(0);
    his3_W3->SetXTitle("West PMT 4");
    his3_W3->GetXaxis()->CenterTitle();
    his3_W3->GetXaxis()->SetRangeUser(-50.0,upperRange+980.);
    his3_W3->Draw();
    
    c6->Print(filenameOutLast);
    
  }

}
