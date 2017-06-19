#include <vector>
#include <algorithm>

std::vector < Float_t> readOctetFileForFGruns(int octet) {
  
  std::vector <Float_t> vec;
  
  char fileName[200];
  sprintf(fileName,"%s/All_Octets/octet_list_%i.dat",getenv("OCTET_LIST"),octet);
  std::ifstream infile(fileName);
  std::string runTypeHold;
  Float_t runNumberHold;
  
  // Populate the map with the runType and runNumber from this octet
  while (infile >> runTypeHold >> runNumberHold) {

    if ( runTypeHold=="A2" || runTypeHold=="A10" || 
	 runTypeHold=="B5" || runTypeHold=="B7" || 
	 runTypeHold=="B2" || runTypeHold=="B10" 
	 || runTypeHold=="A5" || runTypeHold=="A7" )  {
      
      vec.push_back(runNumberHold);
    
    }

  }

  infile.close();
 
  return vec;
};

void plotRates(TString year, TString anaCh, int emin, int emax) {

  gStyle->SetTitleSize(0.06,"t");
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleYSize(0.05);
  gStyle->SetTitleYOffset(0.8);
  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetLabelSize(0.04,"xyz");
  //gStyle->SetOptFit(1111);
  gStyle->SetStatY(0.85);
  gStyle->SetStatX(0.975);
  gStyle->SetStatW(.09);

  gStyle->SetFillStyle(0000); 
  gStyle->SetStatStyle(0); 
  //gStyle->SetTitleStyle(0); 
  //gStyle->SetCanvasBorderSize(0); 
  //gStyle->SetFrameBorderSize(0); 
  //gStyle->SetLegendBorderSize(0); 
  //gStyle->SetStatBorderSize(0); 
  //gStyle->SetTitleBorderSize(0);

  

  //This will plot the rates as a function of Beta Runs

  Int_t binMin = emin/10;
  Int_t binMax = emax/10; // Not including this bin, so the sum goes from binMin->binMax-1

  
  std::vector <Int_t> badOct;// {7,60,61,62,63,64,65,66};
  Int_t octs[] = {7,59,60,61,62,63,64,65,66,67,91,93,101,107,121};

  for (int j=0; j<8; ++j) badOct.push_back(octs[j]);

  Int_t octMin, octMax;

  if (year==TString("2011-2012") ) {
    octMin = 0;
    octMax = 59;
  }

  else {
    octMin = 60;
    octMax = 121;
  }

  

  std::vector < Float_t > runNumber;
  std::vector < std::vector<Float_t> > fgrate(2,std::vector<Float_t>(0));
  std::vector < std::vector<Float_t> > bgrate(2,std::vector<Float_t>(0));
  std::vector < std::vector<Float_t> > fgrateErr(2,std::vector<Float_t>(0));
  std::vector < std::vector<Float_t> > bgrateErr(2,std::vector<Float_t>(0));

  std::vector < Float_t > octetNumber;
  std::vector < std::vector<Float_t> > octet_fgrate(2,std::vector<Float_t>(0));
  std::vector < std::vector<Float_t> > octet_bgrate(2,std::vector<Float_t>(0));
  std::vector < std::vector<Float_t> > octet_fgrateErr(2,std::vector<Float_t>(0));
  std::vector < std::vector<Float_t> > octet_bgrateErr(2,std::vector<Float_t>(0));


  double eastMax=0., westMax=0., eastMin=100., westMin=100.;
  double octet_eastMax=0., octet_westMax=0., octet_eastMin=100., octet_westMin=100.;

  
  for ( Int_t oct=octMin; oct<octMax+1; ++oct ) {
    
    if (std::find(badOct.begin(),badOct.end(),oct)!=badOct.end()) continue;
    std::vector<Float_t> runs = readOctetFileForFGruns(oct);

    octetNumber.push_back(oct);
    
    Float_t octet_eastfgsum = 0., octet_eastbgsum = 0.;
    Float_t octet_easttimeFG = 0., octet_easttimeBG = 0.;
    Float_t octet_westfgsum = 0., octet_westbgsum = 0.;
    Float_t octet_westtimeFG = 0., octet_westtimeBG = 0.;

    for (unsigned i=0; i<runs.size(); ++i) {

      cout << "Processing Run : " << runs[i] << endl;
    
      // These are all double runs of one runtype. They are already accounted for in the previous 
      // run's rate.
      if (runs[i]==17167 || runs[i]==18963 || runs[i]==22200 || runs[i]==22823 || runs[i]==23073) continue;
      runNumber.push_back(runs[i]);

      Float_t fgsum = 0;
      Float_t bgsum = 0;
      
      TString fname = TString::Format("%s/Asymmetry/BinByBinComparison/UK_Erun%0.0f_anaCh%s.dat",
				      getenv("ANALYSIS_CODE"),runs[i],anaCh.Data());
      
      std::ifstream infile(fname.Data());
      
      std::string hold;
      Float_t timeFG, timeBG;
      
      infile >> hold >> hold >> hold >> timeFG
	     >> hold >> hold >> hold >> timeBG;
      
      octet_easttimeFG+=timeFG, octet_easttimeBG+=timeBG;
      
      int bin=0;
      
      Float_t en, fg, bg;
      
      while ( infile >> en >> fg >> bg ) {	
	if ( bin>=binMin && bin<binMax ) {
	  //cout << " " << bin << " " << en << " " << fg << " " << bg << endl;
	  fgsum+=fg;
	  bgsum+=bg;
	  
	}
	bin++;
      }
      
      octet_eastfgsum+=fgsum*timeFG, octet_eastbgsum+=bgsum*timeBG;
      
      eastMax = fgsum>eastMax?fgsum:eastMax;
      eastMin = bgsum<eastMin?bgsum:eastMin;
      
      fgrate[0].push_back(fgsum);
      bgrate[0].push_back(bgsum);
      
      fgrateErr[0].push_back(TMath::Sqrt(fgsum/timeFG));
      bgrateErr[0].push_back(TMath::Sqrt(bgsum/timeBG));
      
      infile.close();
      
      
      //West
      fgsum = 0;
      bgsum = 0;
      
      fname = TString::Format("%s/Asymmetry/BinByBinComparison/UK_Wrun%0.0f_anaCh%s.dat",
			      getenv("ANALYSIS_CODE"),runs[i],anaCh.Data());
      
      infile.open(fname.Data());
      
      infile >> hold >> hold >> hold >> timeFG
	     >> hold >> hold >> hold >> timeBG;
      
      octet_westtimeFG+=timeFG, octet_westtimeBG+=timeBG;
      
      
      bin=0;
      
      while ( infile >> en >> fg >> bg ) {	
	if ( bin>=binMin && bin<binMax ) {
	  //cout << " " << bin << " " << en << " " << fg << " " << bg << endl;
	  fgsum+=fg;
	  bgsum+=bg;
	  
	}
	bin++;
      }
      
      octet_westfgsum+=fgsum*timeFG, octet_westbgsum+=bgsum*timeBG;
      
      westMax = fgsum>westMax?fgsum:westMax;
      westMin = bgsum<westMin?bgsum:westMin;
      
      fgrate[1].push_back(fgsum);
      bgrate[1].push_back(bgsum);
      
      fgrateErr[1].push_back(TMath::Sqrt(fgsum/timeFG));
      bgrateErr[1].push_back(TMath::Sqrt(bgsum/timeBG));
    }
    
    Float_t octet_eastFGRate = octet_eastfgsum/octet_easttimeFG;
    Float_t octet_eastBGRate = octet_eastbgsum/octet_easttimeBG;
    Float_t octet_westFGRate = octet_westfgsum/octet_westtimeFG;
    Float_t octet_westBGRate = octet_westbgsum/octet_westtimeBG;

    octet_eastMax = octet_eastFGRate>octet_eastMax?octet_eastFGRate:octet_eastMax;
    octet_eastMin = octet_eastBGRate<octet_eastMin?octet_eastBGRate:octet_eastMin;
    
    octet_fgrate[0].push_back(octet_eastFGRate);
    octet_bgrate[0].push_back(octet_eastBGRate);
    
    octet_fgrateErr[0].push_back(TMath::Sqrt(octet_eastFGRate/octet_easttimeFG));
    octet_bgrateErr[0].push_back(TMath::Sqrt(octet_eastBGRate/octet_easttimeBG));

    octet_westMax = octet_westFGRate>octet_westMax?octet_westFGRate:octet_westMax;
    octet_westMin = octet_westBGRate<octet_westMin?octet_westBGRate:octet_westMin;
    
    octet_fgrate[1].push_back(octet_westFGRate);
    octet_bgrate[1].push_back(octet_westBGRate);
    
    octet_fgrateErr[1].push_back(TMath::Sqrt(octet_westFGRate/octet_westtimeFG));
    octet_bgrateErr[1].push_back(TMath::Sqrt(octet_westBGRate/octet_westtimeBG));

    
  }

  //Plotting

  //East
  TCanvas *cEast = new TCanvas("cEast","East",1200, 600);
  //cEast->Divide(1,2);
  //cEast->cd(1);
  
  TGraphErrors *eastFG = new TGraphErrors(runNumber.size(),&runNumber[0], &fgrate[0][0],0,&fgrateErr[0][0]);
  eastFG->SetTitle(TString::Format("%s East %0.0f-%0.0f keV Integrated Rates",year.Data(),binMin*10., binMax*10.));
  eastFG->GetYaxis()->SetTitle("Event Rate (Hz)");
  eastFG->GetXaxis()->SetTitle("Run Number");
  eastFG->SetMarkerStyle(20);
  eastFG->SetMarkerColor(kBlack);
  eastFG->SetMinimum(0.75*eastMin);
  eastFG->SetMaximum(1.5*eastMax);
  eastFG->Draw("AP");

  TF1 *eastFGfit = new TF1("eastFGfit","[0]",runNumber[0],runNumber[runNumber.size()-1]);
  eastFGfit->SetParameter(0,10.);
  eastFGfit->SetLineColor(kBlack);
  eastFGfit->SetLineStyle(2);
  eastFGfit->SetLineWidth(3);
  eastFG->Fit(eastFGfit,"R");

  //East->cd(2);

  TGraphErrors *eastBG = new TGraphErrors(runNumber.size(),&runNumber[0], &bgrate[0][0], 0, &bgrateErr[0][0]);
  eastBG->SetTitle(TString::Format("%s East %0.0f-%0.0f keV Integrated Rates",year.Data(),binMin*10., binMax*10.));
  eastBG->GetYaxis()->SetTitle("Event Rate (Hz)");
  eastBG->GetXaxis()->SetTitle("Run Number");
  eastBG->SetMarkerStyle(26);
  eastBG->SetMarkerColor(kBlue);
  eastBG->SetLineColor(kBlue);
  eastBG->Draw("PSAME");

  TF1 *eastBGfit = new TF1("eastBGfit","[0]",runNumber[0],runNumber[runNumber.size()-1]);
  eastBGfit->SetParameter(0,0.1);
  eastBGfit->SetLineColor(kBlue);
  eastBGfit->SetLineStyle(6);
  eastBGfit->SetLineWidth(3);
  eastBG->Fit(eastBGfit,"R");

  TLegend *eastLeg = new TLegend(0.55,0.5,0.8,0.66);
  eastLeg->AddEntry(eastFG,"Foreground","ep");
  eastLeg->AddEntry(eastFGfit,TString::Format("Ave Rate = %0.2f Hz",eastFGfit->GetParameter(0)),"l");
  eastLeg->AddEntry(eastBG,"Background","ep");
  eastLeg->AddEntry(eastBGfit,TString::Format("Ave Rate = %0.2f Hz",eastBGfit->GetParameter(0)),"l");
  eastLeg->Draw("SAME");

  gPad->SetLogy();
  cEast->Update();


  //West
  TCanvas *cWest = new TCanvas("cWest","West",1200, 600);
  //cWest->Divide(1,2);
  //cWest->cd(1);
  
  TGraphErrors *westFG = new TGraphErrors(runNumber.size(),&runNumber[0], &fgrate[1][0],0,&fgrateErr[1][0]);
  westFG->SetTitle(TString::Format("%s West %0.0f-%0.0f keV Integrated Rates",year.Data(),binMin*10., binMax*10.));  westFG->GetYaxis()->SetTitle("Event Rate (Hz)");
  westFG->GetXaxis()->SetTitle("Run Number");
  westFG->SetMarkerStyle(20);
  westFG->SetMarkerColor(kBlack);
  westFG->SetMinimum(0.75*westMin);
  westFG->SetMaximum(1.5*westMax);
  westFG->Draw("AP");

  TF1 *westFGfit = new TF1("westFGfit","[0]",runNumber[0],runNumber[runNumber.size()-1]);
  westFGfit->SetParameter(0,10.);
  westFGfit->SetLineColor(kBlack);
  westFGfit->SetLineStyle(2);
  westFGfit->SetLineWidth(3);
  westFG->Fit(westFGfit,"R");

  //West->cd(2);

  TGraphErrors *westBG = new TGraphErrors(runNumber.size(),&runNumber[0], &bgrate[1][0], 0, &bgrateErr[1][0]);
  westBG->SetTitle(TString::Format("%s West %0.0f-%0.0f keV Integrated Rates",year.Data(),binMin*10., binMax*10.));
  westBG->GetYaxis()->SetTitle("Event Rate (Hz)");
  westBG->GetXaxis()->SetTitle("Run Number");
  westBG->SetMarkerStyle(26);
  westBG->SetMarkerColor(kBlue);
  westBG->SetLineColor(kBlue);
  westBG->Draw("PSAME");

  TF1 *westBGfit = new TF1("westBGfit","[0]",runNumber[0],runNumber[runNumber.size()-1]);
  westBGfit->SetParameter(0,0.1);
  westBGfit->SetLineColor(kBlue);
  westBGfit->SetLineStyle(6);
  westBGfit->SetLineWidth(3);
  westBG->Fit(westBGfit,"R");

  TLegend *westLeg = new TLegend(0.55,0.5,0.8,0.66);
  westLeg->AddEntry(westFG,"Foreground","ep");
  westLeg->AddEntry(westFGfit,TString::Format("Ave Rate = %0.2f Hz",westFGfit->GetParameter(0)),"l");
  westLeg->AddEntry(westBG,"Background","ep");
  westLeg->AddEntry(westBGfit,TString::Format("Ave Rate = %0.2f Hz",westBGfit->GetParameter(0)),"l");
  westLeg->Draw("SAME");

  gPad->SetLogy();
  cWest->Update();

  
  /////////////////////////////////// Octets //////////////////////////////////////////////////////
  //East Octets
  TCanvas *cEast_Octets = new TCanvas("cEast_Octets","East_Octets",1200, 600);
  //cEast_Octets->Divide(1,2);
  //cEast_Octets->cd(1);
  
  TGraphErrors *east_octetsFG = new TGraphErrors(octetNumber.size(),&octetNumber[0], &octet_fgrate[0][0],0,&octet_fgrateErr[0][0]);
  east_octetsFG->SetTitle(TString::Format("%s East Octets %0.0f-%0.0f keV Integrated Rates",year.Data(),binMin*10., binMax*10.));
  east_octetsFG->GetYaxis()->SetTitle("Event Rate (Hz)");
  east_octetsFG->GetXaxis()->SetTitle("Octet Number");
  east_octetsFG->SetMarkerStyle(20);
  east_octetsFG->SetMarkerColor(kBlack);
  east_octetsFG->SetMinimum(0.75*octet_eastMin);
  east_octetsFG->SetMaximum(1.5*octet_eastMax);
  east_octetsFG->Draw("AP");

  TF1 *east_octetsFGfit = new TF1("east_octetsFGfit","[0]",octetNumber[0],octetNumber[octetNumber.size()-1]);
  east_octetsFGfit->SetParameter(0,10.);
  east_octetsFGfit->SetLineColor(kBlack);
  east_octetsFGfit->SetLineStyle(2);
  east_octetsFGfit->SetLineWidth(3);
  east_octetsFG->Fit(east_octetsFGfit,"R");

  //East_Octets->cd(2);

  TGraphErrors *east_octetsBG = new TGraphErrors(octetNumber.size(),&octetNumber[0], &octet_bgrate[0][0], 0, &octet_bgrateErr[0][0]);
  east_octetsBG->SetTitle(TString::Format("%s East Octets %0.0f-%0.0f keV Integrated Rates",year.Data(),binMin*10., binMax*10.));
  east_octetsBG->GetYaxis()->SetTitle("Event Rate (Hz)");
  east_octetsBG->GetXaxis()->SetTitle("Octet Number");
  east_octetsBG->SetMarkerStyle(26);
  east_octetsBG->SetMarkerColor(kBlue);
  east_octetsBG->SetLineColor(kBlue);
  east_octetsBG->Draw("PSAME");

  TF1 *east_octetsBGfit = new TF1("east_octetsBGfit","[0]",octetNumber[0],octetNumber[octetNumber.size()-1]);
  east_octetsBGfit->SetParameter(0,0.1);
  east_octetsBGfit->SetLineColor(kBlue);
  east_octetsBGfit->SetLineStyle(6);
  east_octetsBGfit->SetLineWidth(3);
  east_octetsBG->Fit(east_octetsBGfit,"R");

  TLegend *east_octetsLeg = new TLegend(0.55,0.5,0.8,0.66);
  east_octetsLeg->AddEntry(east_octetsFG,"Foreground","ep");
  east_octetsLeg->AddEntry(east_octetsFGfit,TString::Format("Ave Rate = %0.2f Hz",east_octetsFGfit->GetParameter(0)),"l");
  east_octetsLeg->AddEntry(east_octetsBG,"Background","ep");
  east_octetsLeg->AddEntry(east_octetsBGfit,TString::Format("Ave Rate = %0.2f Hz",east_octetsBGfit->GetParameter(0)),"l");
  east_octetsLeg->Draw("SAME");

  gPad->SetLogy();
  cEast_Octets->Update();


  //West Octets
  TCanvas *cWest_Octets = new TCanvas("cWest_Octets","West_Octets",1200, 600);
  //cWest_Octets->Divide(1,2);
  //cWest_Octets->cd(1);
  
  TGraphErrors *west_octetsFG = new TGraphErrors(octetNumber.size(),&octetNumber[0], &octet_fgrate[1][0],0,&octet_fgrateErr[1][0]);
  west_octetsFG->SetTitle(TString::Format("%s West Octets %0.0f-%0.0f keV Integrated Rates",year.Data(),binMin*10., binMax*10.));
  west_octetsFG->GetYaxis()->SetTitle("Event Rate (Hz)");
  west_octetsFG->GetXaxis()->SetTitle("Octet Number");
  west_octetsFG->SetMarkerStyle(20);
  west_octetsFG->SetMarkerColor(kBlack);
  west_octetsFG->SetMinimum(0.75*octet_westMin);
  west_octetsFG->SetMaximum(1.5*octet_westMax);
  west_octetsFG->Draw("AP");

  TF1 *west_octetsFGfit = new TF1("west_octetsFGfit","[0]",octetNumber[0],octetNumber[octetNumber.size()-1]);
  west_octetsFGfit->SetParameter(0,10.);
  west_octetsFGfit->SetLineColor(kBlack);
  west_octetsFGfit->SetLineStyle(2);
  west_octetsFGfit->SetLineWidth(3);
  west_octetsFG->Fit(west_octetsFGfit,"R");

  //West_Octets->cd(2);

  TGraphErrors *west_octetsBG = new TGraphErrors(octetNumber.size(),&octetNumber[0], &octet_bgrate[1][0], 0, &octet_bgrateErr[1][0]);
  west_octetsBG->SetTitle(TString::Format("%s West Octets %0.0f-%0.0f keV Integrated Rates",year.Data(),binMin*10., binMax*10.));
  west_octetsBG->GetYaxis()->SetTitle("Event Rate (Hz)");
  west_octetsBG->GetXaxis()->SetTitle("Octet Number");
  west_octetsBG->SetMarkerStyle(26);
  west_octetsBG->SetMarkerColor(kBlue);
  west_octetsBG->SetLineColor(kBlue);
  west_octetsBG->Draw("PSAME");

  TF1 *west_octetsBGfit = new TF1("west_octetsBGfit","[0]",octetNumber[0],octetNumber[octetNumber.size()-1]);
  west_octetsBGfit->SetParameter(0,0.1);
  west_octetsBGfit->SetLineColor(kBlue);
  west_octetsBGfit->SetLineStyle(6);
  west_octetsBGfit->SetLineWidth(3);
  west_octetsBG->Fit(west_octetsBGfit,"R");

  TLegend *west_octetsLeg = new TLegend(0.55,0.5,0.8,0.66);
  west_octetsLeg->AddEntry(west_octetsFG,"Foreground","ep");
  west_octetsLeg->AddEntry(west_octetsFGfit,TString::Format("Ave Rate = %0.2f Hz",west_octetsFGfit->GetParameter(0)),"l");
  west_octetsLeg->AddEntry(west_octetsBG,"Background","ep");
  west_octetsLeg->AddEntry(west_octetsBGfit,TString::Format("Ave Rate = %0.2f Hz",west_octetsBGfit->GetParameter(0)),"l");
  west_octetsLeg->Draw("SAME");

  gPad->SetLogy();
  cWest_Octets->Update();


  
  //////////////////////////////////////////////////

  cEast->Print(TString::Format("integratedRatesByRun_%s_anaCh%s_%0.0f-%0.0f.pdf(",year.Data(),anaCh.Data(),binMin*10.,binMax*10.));
  cWest->Print(TString::Format("integratedRatesByRun_%s_anaCh%s_%0.0f-%0.0f.pdf",year.Data(),anaCh.Data(),binMin*10.,binMax*10.));
  cEast_Octets->Print(TString::Format("integratedRatesByRun_%s_anaCh%s_%0.0f-%0.0f.pdf",year.Data(),anaCh.Data(),binMin*10.,binMax*10.));
  cWest_Octets->Print(TString::Format("integratedRatesByRun_%s_anaCh%s_%0.0f-%0.0f.pdf)",year.Data(),anaCh.Data(),binMin*10.,binMax*10.));

}
