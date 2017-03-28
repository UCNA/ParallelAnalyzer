/*

Study the edge effects from purely geometric effects in the 
simulated geometry.

Plotted will be all peaks for each ring at a radius r
*/
#include <vector>
#include "../revCalSim.h"

void srcPositionStudy(Int_t binWidth, TString source, TString geometry) {

  int numFiles = 200; 

  TString simLocation;
  TChain *chain = new TChain("anaTree");

  if (geometry==TString("2011-2012")) simLocation = TString(getenv("SIM_2011_2012"));
  else if (geometry==TString("2012-2013")) simLocation = TString(getenv("SIM_2012_2013"));
  else if (geometry==TString("2012-2013_ISOBUTANE")) simLocation = TString(getenv("SIM_2012_2013_ISOBUTANE"));
  else { std::cout << "BAD GEOMETRY\n"; exit(0); }

  for (int i=0; i<numFiles; i++) {
    chain->AddFile(TString::Format("%s/%s/analyzed_%i.root",simLocation.Data(),source.Data(),i));
  }
  

  Double_t maxEn = 1500.;
  Double_t minEn = 0.;
  Double_t fidMax = 55.;
  Int_t nHists = (int)(fidMax/binWidth);
  
  
  
  std::vector <TH1D*> histsE(nHists, 0);
  std::vector <TH1D*> histsW(nHists, 0);
  std::vector <Double_t> rmin(nHists,0);
  std::vector <Double_t> rmax(nHists,0);
  std::vector <Double_t> rmid(nHists,0);

  //final means and errors
  std::vector < Double_t > EastMeans(nHists,0.);
  std::vector < Double_t > WestMeans(nHists,0.);
  
  std::vector < Double_t > EastMeanErrors(nHists,0.);
  std::vector < Double_t > WestMeanErrors(nHists,0.);

  for (Int_t i=0; i<nHists; i++) {

    rmin[i] = i*binWidth;
    rmax[i] = (i+1)*binWidth;
    rmid[i] = (double)rmin[i] + (double)binWidth/2.;

    histsE[i] = new TH1D(TString::Format("his%iE",i), 
			 TString::Format("%s %i mm radius bins East", source.Data(),binWidth),
			550, 0., 1100.);
    histsW[i] = new TH1D(TString::Format("his%iW",i), 
			 TString::Format("%s %i mm radius bins West", source.Data(),binWidth),
			550, 0., 1100.);

    histsE[i]->SetLineColor(i+1);
    histsW[i]->SetLineColor(i+1);
  }

  Double_t primPos[4];

  // Set the addresses of the information read in from the simulation file
  chain->SetBranchAddress("MWPCEnergy",&mwpcE);
  chain->SetBranchAddress("time",&Time);
  chain->SetBranchAddress("Edep",&edep);
  chain->SetBranchAddress("EdepQ",&edepQ);
  chain->SetBranchAddress("MWPCPos",&mwpc_pos);
  chain->SetBranchAddress("ScintPos",&scint_pos);
  chain->SetBranchAddress("primKE",&primKE);
  chain->SetBranchAddress("primTheta",&primTheta);
  chain->SetBranchAddress("primPos",&primPos);


  //Get total number of events in TChain
  UInt_t nevents = chain->GetEntries();
  cout << "events = " << nevents << endl;

  for (Int_t i=0; i<nevents; i++) {
    
    chain->GetEvent(i);

    Int_t nBin = primPos[3]*1000./binWidth;

    if (edepQ.EdepQE>0. && primKE<maxEn && primKE>minEn && mwpcE.MWPCEnergyE>0.1 && primTheta>TMath::Pi()/2.) histsE[nBin]->Fill(edepQ.EdepQE);
    if (edepQ.EdepQW>0. && primKE<maxEn && primKE>minEn && mwpcE.MWPCEnergyW>0.1 && primTheta<TMath::Pi()/2.) histsW[nBin]->Fill(edepQ.EdepQW);

    if (i%100000==0) std::cout << "*";
  }
  
  std::cout << std::endl;

  TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
  c1->Divide(2,2);
  //TCanvas *c2 = new TCanvas("c2");
  
  //histsE[1]->Draw("SAME");
  
  Double_t refMeanE = 0., refMeanW = 0.; //This will hold the mean of the center pixel

  std::cout << nHists << endl;  
  for (Int_t i=0; i<nHists; i++) {
    c1->cd(1);
    histsE[i]->Draw("SAME");
    
    EastMeans[i] = histsE[i]->GetMean();
    EastMeanErrors[i] = EastMeans[i]>0. ? EastMeans[i]/sqrt(histsE[i]->GetEntries()) : 0.;
    //cout << EastMeans[i][0] << " " << EastMeans[i][1] << endl;

    c1->cd(2);
    histsW[i]->Draw("SAME");
    
    WestMeans[i] = histsW[i]->GetMean();
    WestMeanErrors[i] = WestMeans[i]>0. ? WestMeans[i]/sqrt(histsW[i]->GetEntries()) : 0.;

    if (i==0) { refMeanE = EastMeans[i], refMeanW = WestMeans[i]; }
  }

  

  //TCanvas *c3 = new TCanvas("c3");
  c1->cd(3);

  std::vector <Double_t> xerr(nHists,0.);
  
  TGraphErrors *gEast = new TGraphErrors(11, &rmax[0],&EastMeans[0],&xerr[0],&EastMeanErrors[0]);
  gEast->SetMarkerStyle(21);
  gEast->SetTitle(TString::Format("East %s Mean EQ vs. position",source.Data()));
  gEast->GetXaxis()->SetTitle("Position Bin Edge (mm)");
  gEast->GetYaxis()->SetTitle("Energy (keV)");
  gEast->Draw("AP");

  TLine *eastLine = new TLine(gEast->GetXaxis()->GetXmin(),refMeanE,gEast->GetXaxis()->GetXmax(),refMeanE);
  eastLine->SetLineStyle(2);
  eastLine->Draw();
 

  //TCanvas *c4 = new TCanvas("c4");
  c1->cd(4);

  TGraphErrors *gWest = new TGraphErrors(11, &rmax[0],&WestMeans[0],&xerr[0],&WestMeanErrors[0]);
  gWest->SetMarkerStyle(21);
  gWest->SetTitle(TString::Format("West %s Mean EQ vs. position",source.Data()));
  gWest->GetXaxis()->SetTitle("Position Bin Edge (mm)");
  gWest->GetYaxis()->SetTitle("Energy (keV)");
  gWest->Draw("AP");

  TLine *westLine = new TLine(gWest->GetXaxis()->GetXmin(),refMeanW,gWest->GetXaxis()->GetXmax(),refMeanW);
  westLine->SetLineStyle(2);
  westLine->Draw();

  
    


}
