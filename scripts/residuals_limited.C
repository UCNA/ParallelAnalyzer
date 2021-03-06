{
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

  Int_t runLow = 18745;
  Int_t runHigh = 19365;

  // Read East data file
  char tempEast[500];
  sprintf(tempEast, "residuals_East_%i-%i.dat", runLow, runHigh);
  ifstream fileEast(tempEast);

  size_t N = 500;
  TString sourceEast[N];
  Int_t runEast[N];
  Double_t resEast[N];
  Double_t resCeEast[N], resSnEast[N], resBiEast[N];

  Int_t i = 0;
  Int_t nCeEast = 0;
  Int_t nSnEast = 0;
  Int_t nBiEast = 0;
  while (!fileEast.eof()) {
    fileEast >> sourceEast[i] >> runEast[i] >> resEast[i];
    if (sourceEast[i] == "Ce_East") {
      resCeEast[nCeEast] = resEast[i];
      nCeEast++;
    }
    if (sourceEast[i] == "Sn_East") {
      resSnEast[nSnEast] = resEast[i];
      nSnEast++;
    }
    if (sourceEast[i] == "Bi_East") {
      resBiEast[nBiEast] = resEast[i];
      nBiEast++;
    }
    if (fileEast.fail()) break;
    i++;
  }
  cout << nCeEast << " " << nSnEast << " " << nBiEast << endl;

  // Read West data file
  char tempWest[500];
  sprintf(tempWest, "residuals_West_%i-%i.dat", runLow, runHigh);
  ifstream fileWest(tempWest);

  TString sourceWest[N];
  Int_t runWest[N];
  Double_t resWest[N];
  Double_t resCeWest[N], resSnWest[N], resBiWest[N];

  Int_t i = 0;
  Int_t nCeWest = 0;
  Int_t nSnWest = 0;
  Int_t nBiWest = 0;
  while (!fileWest.eof()) {
    fileWest >> sourceWest[i] >> runWest[i] >> resWest[i];
    if (sourceWest[i] == "Ce_West") {
      resCeWest[nCeWest] = resWest[i];
      nCeWest++;
    }
    if (sourceWest[i] == "Sn_West") {
      resSnWest[nSnWest] = resWest[i];
      nSnWest++;
    }
    if (sourceWest[i] == "Bi_West") {
      resBiWest[nBiWest] = resWest[i];
      nBiWest++;
    }
    if (fileWest.fail()) break;
    i++;
  }
  cout << nCeWest << " " << nSnWest << " " << nBiWest << endl;

  // Histograms
  TH1F *hisCeEast = new TH1F("hisCeEast", "", 60, -30.0, 30.0);
  TH1F *hisSnEast = new TH1F("hisSnEast", "", 30, -30.0, 30.0);
  TH1F *hisBiEast = new TH1F("hisBiEast", "", 30, -30.0, 30.0);
  TH1F *hisCeWest = new TH1F("hisCeWest", "", 60, -30.0, 30.0);
  TH1F *hisSnWest = new TH1F("hisSnWest", "", 30, -30.0, 30.0);
  TH1F *hisBiWest = new TH1F("hisBiWest", "", 30, -30.0, 30.0);

  // Fill histograms
  for (int j=0; j<nCeEast; j++) {
    hisCeEast->Fill(resCeEast[j]);
  }
  for (int j=0; j<nSnEast; j++) {
    hisSnEast->Fill(resSnEast[j]);
  }
  for (int j=0; j<nBiEast; j++) {
    hisBiEast->Fill(resBiEast[j]);
  }
  for (int j=0; j<nCeWest; j++) {
    hisCeWest->Fill(resCeWest[j]);
  }
  for (int j=0; j<nSnWest; j++) {
    hisSnWest->Fill(resSnWest[j]);
  }
  for (int j=0; j<nBiWest; j++) {
    hisBiWest->Fill(resBiWest[j]);
  }
 
  // Calculate mean and standard deviation
  double meanCeEast = 0.;
  for (int j=0; j<nCeEast; j++) {
    meanCeEast += resCeEast[j];
  }
  meanCeEast = meanCeEast / (double) nCeEast;

  double sigmaCeEast = 0.;
  for (int j=0; j<nCeEast; j++) {
    sigmaCeEast += (resCeEast[j] - meanCeEast)*(resCeEast[j] - meanCeEast);
  }
  sigmaCeEast = sigmaCeEast / (double) nCeEast;
  sigmaCeEast = sqrt( sigmaCeEast );

  cout << "meanCeEast = " << meanCeEast << endl;
  cout << "     sigma = " << sigmaCeEast << endl;

  double meanSnEast = 0.;
  for (int j=0; j<nSnEast; j++) {
    meanSnEast += resSnEast[j];
  }
  meanSnEast = meanSnEast / (double) nSnEast;

  double sigmaSnEast = 0.;
  for (int j=0; j<nSnEast; j++) {
    sigmaSnEast += (resSnEast[j] - meanSnEast)*(resSnEast[j] - meanSnEast);
  }
  sigmaSnEast = sigmaSnEast / (double) nSnEast;
  sigmaSnEast = sqrt( sigmaSnEast );

  cout << "meanSnEast = " << meanSnEast << endl;
  cout << "     sigma = " << sigmaSnEast << endl;

  double meanBiEast = 0.;
  for (int j=0; j<nBiEast; j++) {
    meanBiEast += resBiEast[j];
  }
  meanBiEast = meanBiEast / (double) nBiEast;

  double sigmaBiEast = 0.;
  for (int j=0; j<nBiEast; j++) {
    sigmaBiEast += (resBiEast[j] - meanBiEast)*(resBiEast[j] - meanBiEast);
  }
  sigmaBiEast = sigmaBiEast / (double) nBiEast;
  sigmaBiEast = sqrt( sigmaBiEast );

  cout << "meanBiEast = " << meanBiEast << endl;
  cout << "     sigma = " << sigmaBiEast << endl;

  // Ce East
  c1 = new TCanvas("c1", "c1");
  c1->SetLogy(0);

  hisCeEast->SetXTitle("East Ce E_{Q} Error [keV]");
  hisCeEast->GetXaxis()->CenterTitle();
  hisCeEast->SetLineColor(1);
  hisCeEast->Draw();
  //hisCeEast->Fit("gaus", "", "", -5.0, 5.0);

  // Sn East
  c2 = new TCanvas("c2", "c2");
  c2->SetLogy(0);

  hisSnEast->SetXTitle("East Sn E_{Q} Error [keV]");
  hisSnEast->GetXaxis()->CenterTitle();
  hisSnEast->SetLineColor(1);
  hisSnEast->Draw();
  //hisSnEast->Fit("gaus", "", "", -10.0, 10.0);

  // Bi East
  c3 = new TCanvas("c3", "c3");
  c3->SetLogy(0);

  hisBiEast->SetXTitle("East Bi E_{Q} Error [keV]");
  hisBiEast->GetXaxis()->CenterTitle();
  hisBiEast->SetLineColor(1);
  hisBiEast->Draw();
  //hisBiEast->Fit("gaus", "");

  // Ce West
  c4 = new TCanvas("c4", "c4");
  c4->SetLogy(0);

  hisCeWest->SetXTitle("West Ce E_{Q} Error [keV]");
  hisCeWest->GetXaxis()->CenterTitle();
  hisCeWest->SetLineColor(1);
  hisCeWest->Draw();
  //hisCeWest->Fit("gaus", "");

  // Sn West
  c5 = new TCanvas("c5", "c5");
  c5->SetLogy(0);

  hisSnWest->SetXTitle("West Sn E_{Q} Error [keV]");
  hisSnWest->GetXaxis()->CenterTitle();
  hisSnWest->SetLineColor(1);
  hisSnWest->Draw();
  //hisSnWest->Fit("gaus", "");

  // Bi West
  c6 = new TCanvas("c6", "c6");
  c6->SetLogy(0);

  hisBiWest->SetXTitle("West Bi E_{Q} Error [keV]");
  hisBiWest->GetXaxis()->CenterTitle();
  hisBiWest->SetLineColor(1);
  hisBiWest->Draw();
  //hisBiWest->Fit("gaus", "");

}
