//Plots the pedestals as a function of run number (time)

//Should make an array of all the calibration periods and the
// runs they apply to, and then put vertical lines at each of these
#include <vector>

std::vector <Int_t> getPMTQuality(Int_t runNumber) {
  //Read in PMT quality file
  //cout << "Reading in PMT Quality file ...\n";
  vector <Int_t>  pmtQuality (8,0);
  Char_t temp[200];
  sprintf(temp,"%s/residuals/PMT_runQuality_master.dat",getenv("ANALYSIS_CODE")); 
  ifstream pmt;
  //std::cout << temp << std::endl;
  pmt.open(temp);
  Int_t run_hold;
  while (pmt >> run_hold >> pmtQuality[0] >> pmtQuality[1] >> pmtQuality[2]
	 >> pmtQuality[3] >> pmtQuality[4] >> pmtQuality[5]
	 >> pmtQuality[6] >> pmtQuality[7]) {
    if (run_hold==runNumber) break;
    if (pmt.fail()) break;
  }
  pmt.close();
  if (run_hold!=runNumber) {
    cout << "Run not found in PMT quality file!" << endl;
    exit(0);
  }
  return pmtQuality;
}

void plotPedestalTimeDependence(TString year) {

  gStyle->SetTitleSize(0.08,"t");
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleYSize(0.06);
  gStyle->SetTitleYOffset(0.6);
  gStyle->SetTitleXSize(0.06);
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetLabelSize(0.06,"xyz");
  
  std::vector <Int_t> calPeriodRunEnd; // TO BE IMPLEMENTED LATER
  std::vector <Int_t> calPeriodRefRun;



  Int_t runMin, runMax;

  if (year==TString("2011-2012") ) {
    runMin = 17000;
    runMax = 20000;
    
    calPeriodRefRun.push_back(17238); calPeriodRunEnd.push_back(17298);
    calPeriodRefRun.push_back(17370); calPeriodRunEnd.push_back(17440);
    calPeriodRefRun.push_back(17521); calPeriodRunEnd.push_back(17735);
    calPeriodRefRun.push_back(17892); calPeriodRunEnd.push_back(17956);
    calPeriodRefRun.push_back(18361); calPeriodRunEnd.push_back(18387);
    calPeriodRefRun.push_back(18621); calPeriodRunEnd.push_back(18684);
    calPeriodRefRun.push_back(18749); calPeriodRunEnd.push_back(18995);
    calPeriodRefRun.push_back(19232); calPeriodRunEnd.push_back(19240);
    calPeriodRefRun.push_back(19359); calPeriodRunEnd.push_back(19545);
    calPeriodRefRun.push_back(19857); calPeriodRunEnd.push_back(20000);
  }
  else {
    runMin = 21000;
    runMax = 23300;

    calPeriodRefRun.push_back(20519); calPeriodRunEnd.push_back(20742);
    calPeriodRefRun.push_back(20820); calPeriodRunEnd.push_back(20838);
    calPeriodRefRun.push_back(21091); calPeriodRunEnd.push_back(21238);
    calPeriodRefRun.push_back(21315); calPeriodRunEnd.push_back(21606);
    calPeriodRefRun.push_back(21683); calPeriodRunEnd.push_back(21864);
    calPeriodRefRun.push_back(21918); calPeriodRunEnd.push_back(22119);
    calPeriodRefRun.push_back(22219); calPeriodRunEnd.push_back(22239);
    calPeriodRefRun.push_back(22441); calPeriodRunEnd.push_back(22631);
    calPeriodRefRun.push_back(22771); calPeriodRunEnd.push_back(23174);
  }

  std::vector < std::vector <Double_t> > runNumber( 8, std::vector<Double_t>(4000) );
  std::vector < std::vector <Double_t> > eastPedestal( 4, std::vector<Double_t>(4000) );
  std::vector < std::vector <Double_t> > westPedestal( 4, std::vector<Double_t>(4000) );
  std::vector < std::vector <Double_t> > eastPedestalRMS( 4, std::vector<Double_t>(4000) );
  std::vector < std::vector <Double_t> > westPedestalRMS( 4, std::vector<Double_t>(4000) );
  std::vector < Int_t > numRuns(8,0);

  //These are hold the reference runs specifically to be drawn in a different color
  std::vector < std::vector <Double_t> > runNumber_ref( 8, std::vector<Double_t>(20) );
  std::vector < std::vector <Double_t> > eastPedestal_ref( 4, std::vector<Double_t>(20) );
  std::vector < std::vector <Double_t> > westPedestal_ref( 4, std::vector<Double_t>(20) );
  std::vector < std::vector <Double_t> > eastPedestalRMS_ref( 4, std::vector<Double_t>(20) );
  std::vector < std::vector <Double_t> > westPedestalRMS_ref( 4, std::vector<Double_t>(20) );
  std::vector < Int_t > numRuns_ref( 8,0 );
  

  double e1=0., e2=0., e3=0., e4=0., w1=0., w2=0., w3=0., w4=0.;
  double e1_RMS=0., e2_RMS=0., e3_RMS=0., e4_RMS=0., w1_RMS=0., w2_RMS=0., w3_RMS=0., w4_RMS=0.;
  int run = 0;
  ifstream pedestalFile;

  for (Int_t rn=runMin; rn<runMax; rn++) {
    
    if (rn%500==0) std::cout << "NOW AT RUN " << rn << std::endl;

    TString filename = TString::Format("%s/PMT_pedestals_%i.dat",getenv("PEDESTALS"),rn);
    pedestalFile.open(filename.Data());
    
    if (pedestalFile.is_open()) {
      if ( rn==17588 || (rn>19583 && rn<19616) || (rn>=16000 && rn<=17078) ) { pedestalFile.close(); continue; } //Groups of runs to ignore due to bad crap
      // 17588 is empty, the ranges are unused Xe runs which weren't processed, 
      // (rn>20515 && rn<21087) is all the garbage which isn't really usable from the beginning of 2012-2013

      std::vector <Int_t> pmtQuality = getPMTQuality(rn);
      
      //std::cout << "Opened " << filename << std::endl;
      pedestalFile >> run >> e1 >> e1_RMS
	       >> run >> e2 >> e2_RMS
	       >> run >> e3 >> e3_RMS
	       >> run >> e4 >> e4_RMS
	       >> run >> w1 >> w1_RMS
	       >> run >> w2 >> w2_RMS
	       >> run >> w3 >> w3_RMS
	       >> run >> w4 >> w4_RMS;
      
      if (pmtQuality[0]) { 
	runNumber[0][numRuns[0]]=rn; eastPedestal[0][numRuns[0]]=e1; eastPedestalRMS[0][numRuns[0]]=e1_RMS; numRuns[0]++; 
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    runNumber_ref[0][numRuns_ref[0]]=rn; eastPedestal_ref[0][numRuns_ref[0]]=e1; eastPedestalRMS_ref[0][numRuns_ref[0]]=e1_RMS; numRuns_ref[0]++; 
	    continue;
	  }
	}
      }
      if (pmtQuality[1]) { 
	runNumber[1][numRuns[1]]=rn; eastPedestal[1][numRuns[1]]=e2; eastPedestalRMS[1][numRuns[1]]=e2_RMS; numRuns[1]++; 
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    runNumber_ref[1][numRuns_ref[1]]=rn; eastPedestal_ref[1][numRuns_ref[1]]=e2; eastPedestalRMS_ref[1][numRuns_ref[1]]=e2_RMS; numRuns_ref[1]++; 
	    continue;
	  }
	}
      }
      if (pmtQuality[2]) {
	runNumber[2][numRuns[2]]=rn; eastPedestal[2][numRuns[2]]=e3; eastPedestalRMS[2][numRuns[2]]=e3_RMS; numRuns[2]++; 
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    runNumber_ref[2][numRuns_ref[2]]=rn; eastPedestal_ref[2][numRuns_ref[2]]=e3; eastPedestalRMS_ref[2][numRuns_ref[2]]=e3_RMS; numRuns_ref[2]++; 
	    continue;
	  }
	}
      }
      if (pmtQuality[3]) { 
	runNumber[3][numRuns[3]]=rn; eastPedestal[3][numRuns[3]]=e4; eastPedestalRMS[3][numRuns[3]]=e4_RMS; numRuns[3]++; 
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    runNumber_ref[3][numRuns_ref[3]]=rn; eastPedestal_ref[3][numRuns_ref[3]]=e4; eastPedestalRMS_ref[3][numRuns_ref[3]]=e4_RMS; numRuns_ref[3]++;
	    continue;; 
	  }
	}
      }
      if (pmtQuality[4]) { 
	runNumber[4][numRuns[4]]=rn; westPedestal[0][numRuns[4]]=w1; westPedestalRMS[0][numRuns[4]]=w1_RMS; numRuns[4]++; 
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    runNumber_ref[4][numRuns_ref[4]]=rn; westPedestal_ref[0][numRuns_ref[4]]=w1; westPedestalRMS_ref[0][numRuns_ref[4]]=w1_RMS; numRuns_ref[4]++; 
	    continue;
	  }
	}
      }
      else w1=1000.; //so it doesn't trigger the check below
      if (pmtQuality[5]) { 
	runNumber[5][numRuns[5]]=rn; westPedestal[1][numRuns[5]]=w2; westPedestalRMS[1][numRuns[5]]=w2_RMS; numRuns[5]++; 
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    runNumber_ref[5][numRuns_ref[5]]=rn; westPedestal_ref[1][numRuns_ref[5]]=w2; westPedestalRMS_ref[1][numRuns_ref[5]]=w2_RMS; numRuns_ref[5]++; 
	    continue;
	  }
	}
      }
      if (pmtQuality[6]) { 
	runNumber[6][numRuns[6]]=rn; westPedestal[2][numRuns[6]]=w3; westPedestalRMS[2][numRuns[6]]=w3_RMS; numRuns[6]++; 
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    runNumber_ref[6][numRuns_ref[6]]=rn; westPedestal_ref[2][numRuns_ref[6]]=w3; westPedestalRMS_ref[2][numRuns_ref[6]]=w3_RMS; numRuns_ref[6]++;
	    continue;
	  }
	}
      }
      if (pmtQuality[7]) { 
	runNumber[7][numRuns[7]]=rn; westPedestal[3][numRuns[7]]=w4; westPedestalRMS[3][numRuns[7]]=w4_RMS; numRuns[7]++; 
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    runNumber_ref[7][numRuns_ref[7]]=rn; westPedestal_ref[3][numRuns_ref[7]]=w4; westPedestalRMS_ref[3][numRuns_ref[7]]=w4_RMS; numRuns_ref[7]++; 
	    continue;
	  }
	}
      }
      
      
      pedestalFile.close();
    }
  }


  // making all of the TLines to reperesent the calibration periods

  std::vector <TLine*> linesE1(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesE2(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesE3(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesE4(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesW1(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesW2(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesW3(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesW4(calPeriodRunEnd.size(),0);

  std::vector <TLine*> linesE1_fac(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesE2_fac(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesE3_fac(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesE4_fac(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesW1_fac(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesW2_fac(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesW3_fac(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesW4_fac(calPeriodRunEnd.size(),0);

  Double_t min = 0., max = 0.;
  

  //EAST SIDE
  TCanvas *cEast = new TCanvas("cEast","East PMTs",1000, 1200);
  cEast->Divide(1,4);

  cEast->cd(1);

  TGraphErrors *east1 = new TGraphErrors(numRuns[0],&runNumber[0][0],&eastPedestal[0][0],0,&eastPedestalRMS[0][0]);
  east1->SetTitle("PMT East 1");
  east1->GetXaxis()->SetTitle("Run Number");
  east1->GetYaxis()->SetTitle("ADC Channels");
  east1->GetXaxis()->SetLimits(runMin, runMax);
  east1->SetMarkerStyle(8);
  // east1->SetMinimum(0.);
  //east1->SetMaximum(4096.);
  east1->Draw("AP");
  
  cEast->Update();
  
  min = gPad->GetUymin();
  max = gPad->GetUymax();

  for (UInt_t i=0; i<calPeriodRunEnd.size(); i++) {
    
    linesE1[i] = new TLine((double)calPeriodRunEnd[i],min,(double)calPeriodRunEnd[i],max);
    linesE1[i]->SetLineStyle(2);
    linesE1[i]->SetLineColor(kRed);
    linesE1[i]->SetLineWidth(2);
    if ( (double)calPeriodRunEnd[i] > runNumber[0][0] && (double)calPeriodRunEnd[i] < runNumber[0][numRuns[0]-1] ) {
      linesE1[i]->Draw();
    }
  } 
  
  TGraphErrors *east1_ref = new TGraphErrors(numRuns_ref[0],&runNumber_ref[0][0],&eastPedestal_ref[0][0],0,&eastPedestalRMS_ref[0][0]);
  east1_ref->SetTitle("PMT East 1");
  east1_ref->GetXaxis()->SetTitle("Run Number");
  east1_ref->GetYaxis()->SetTitle("ADC Channels");
  east1_ref->SetMarkerStyle(8);
  east1_ref->SetMarkerColor(kRed);
  // east1_ref->SetMinimum(0.);
  //east1_ref->SetMaximum(4096.);
  east1_ref->Draw("P SAME");

  

  

  cEast->cd(2);

  TGraphErrors *east2 = new TGraphErrors(numRuns[1],&runNumber[1][0],&eastPedestal[1][0],0,&eastPedestalRMS[1][0]);
  east2->SetTitle("PMT East 2");
  east2->GetXaxis()->SetTitle("Run Number");
  east2->GetYaxis()->SetTitle("ADC Channels");
  east2->GetXaxis()->SetLimits(runMin, runMax);
  east2->SetMarkerStyle(8);
  //east2->SetMinimum(0.);
  //east2->SetMaximum(4096.);
  east2->Draw("AP");

  cEast->Update();
  
  min = gPad->GetUymin();
  max = gPad->GetUymax();

  for (UInt_t i=0; i<calPeriodRunEnd.size(); i++) {
    
    linesE2[i] = new TLine((double)calPeriodRunEnd[i],min,(double)calPeriodRunEnd[i],max);
    linesE2[i]->SetLineStyle(2);
    linesE2[i]->SetLineColor(kRed);
    linesE2[i]->SetLineWidth(2);
    if ( (double)calPeriodRunEnd[i] > runNumber[1][0] && (double)calPeriodRunEnd[i] < runNumber[1][numRuns[1]-1] ) {
      linesE2[i]->Draw();
    }
  }

  TGraphErrors *east2_ref = new TGraphErrors(numRuns_ref[1],&runNumber_ref[1][0],&eastPedestal_ref[1][0],0,&eastPedestalRMS_ref[1][0]);
  east2_ref->SetTitle("PMT East 2");
  east2_ref->GetXaxis()->SetTitle("Run Number");
  east2_ref->GetYaxis()->SetTitle("ADC Channels");
  east2_ref->SetMarkerStyle(8);
  east2_ref->SetMarkerColor(kRed);
  // east2_ref->SetMinimum(0.);
  //east2_ref->SetMaximum(4096.);
  east2_ref->Draw("P SAME");


  cEast->cd(3);

  TGraphErrors *east3 = new TGraphErrors(numRuns[2],&runNumber[2][0],&eastPedestal[2][0],0,&eastPedestalRMS[2][0]);
  east3->SetTitle("PMT East 3");
  east3->GetXaxis()->SetTitle("Run Number");
  east3->GetYaxis()->SetTitle("ADC Channels");
  east3->GetXaxis()->SetLimits(runMin, runMax);
  east3->SetMarkerStyle(8);
  //east3->SetMinimum(0.);
  //east3->SetMaximum(4096.);
  east3->Draw("AP");

  cEast->Update();
  
  min = gPad->GetUymin();
  max = gPad->GetUymax();

  for (UInt_t i=0; i<calPeriodRunEnd.size(); i++) {
    
    linesE3[i] = new TLine((double)calPeriodRunEnd[i],min,(double)calPeriodRunEnd[i],max);
    linesE3[i]->SetLineStyle(2);
    linesE3[i]->SetLineColor(kRed);
    linesE3[i]->SetLineWidth(2);
    if ( (double)calPeriodRunEnd[i] > runNumber[2][0] && (double)calPeriodRunEnd[i] < runNumber[2][numRuns[2]-1] ) {
      linesE3[i]->Draw();
    }
  }

  TGraphErrors *east3_ref = new TGraphErrors(numRuns_ref[2],&runNumber_ref[2][0],&eastPedestal_ref[2][0],0,&eastPedestalRMS_ref[2][0]);
  east3_ref->SetTitle("PMT East 3");
  east3_ref->GetXaxis()->SetTitle("Run Number");
  east3_ref->GetYaxis()->SetTitle("ADC Channels");
  east3_ref->SetMarkerStyle(8);
  east3_ref->SetMarkerColor(kRed);
  // east3_ref->SetMinimum(0.);
  //east3_ref->SetMaximum(4096.);
  east3_ref->Draw("P SAME");
  

  cEast->cd(4);

  TGraphErrors *east4 = new TGraphErrors(numRuns[3],&runNumber[3][0],&eastPedestal[3][0],0,&eastPedestalRMS[3][0]);
  east4->SetTitle("PMT East 4");
  east4->GetXaxis()->SetTitle("Run Number");
  east4->GetYaxis()->SetTitle("ADC Channels");
  east4->GetXaxis()->SetLimits(runMin, runMax);
  east4->SetMarkerStyle(8);
  //east4->SetMinimum(0.);
  //east4->SetMaximum(4096.);
  east4->Draw("AP");

  cEast->Update();
  
  min = gPad->GetUymin();
  max = gPad->GetUymax();

  for (UInt_t i=0; i<calPeriodRunEnd.size(); i++) {
    
    linesE4[i] = new TLine((double)calPeriodRunEnd[i],min,(double)calPeriodRunEnd[i],max);
    linesE4[i]->SetLineStyle(2);
    linesE4[i]->SetLineColor(kRed);
    linesE4[i]->SetLineWidth(2);
    if ( (double)calPeriodRunEnd[i] > runNumber[3][0] && (double)calPeriodRunEnd[i] < runNumber[3][numRuns[3]-1] ) {
      linesE4[i]->Draw();
    }
  }

  TGraphErrors *east4_ref = new TGraphErrors(numRuns_ref[3],&runNumber_ref[3][0],&eastPedestal_ref[3][0],0,&eastPedestalRMS_ref[3][0]);
  east4_ref->SetTitle("PMT East 4");
  east4_ref->GetXaxis()->SetTitle("Run Number");
  east4_ref->GetYaxis()->SetTitle("ADC Channels");
  east4_ref->SetMarkerStyle(8);
  east4_ref->SetMarkerColor(kRed);
  // east4_ref->SetMinimum(0.);
  //east4_ref->SetMaximum(4096.);
  east4_ref->Draw("P SAME");



  //WEST SIDE
  TCanvas *cWest = new TCanvas("cWest","West PMTs",1000, 1200);
  cWest->Divide(1,4);

  cWest->cd(1);

  TGraphErrors *west1 = new TGraphErrors(numRuns[4],&runNumber[4][0],&westPedestal[0][0],0,&westPedestalRMS[0][0]);
  west1->SetTitle("PMT West 1");
  west1->GetXaxis()->SetTitle("Run Number");
  west1->GetYaxis()->SetTitle("ADC Channels");
  west1->GetXaxis()->SetLimits(runMin, runMax);
  west1->SetMarkerStyle(8);
  //west1->SetMinimum(0.);
  //west1->SetMaximum(4096.);
  west1->Draw("AP");

  cWest->Update();
  
  min = gPad->GetUymin();
  max = gPad->GetUymax();

  for (UInt_t i=0; i<calPeriodRunEnd.size(); i++) {
    
    linesW1[i] = new TLine((double)calPeriodRunEnd[i],min,(double)calPeriodRunEnd[i],max);
    linesW1[i]->SetLineStyle(2);
    linesW1[i]->SetLineColor(kRed);
    linesW1[i]->SetLineWidth(2);
    if ( (double)calPeriodRunEnd[i] > runNumber[4][0] && (double)calPeriodRunEnd[i] < runNumber[4][numRuns[4]-1] ) {
      linesW1[i]->Draw();
    }
  }

  TGraphErrors *west1_ref = new TGraphErrors(numRuns_ref[4],&runNumber_ref[4][0],&westPedestal_ref[0][0],0,&westPedestalRMS_ref[0][0]);
  west1_ref->SetTitle("PMT West 1");
  west1_ref->GetXaxis()->SetTitle("Run Number");
  west1_ref->GetYaxis()->SetTitle("ADC Channels");
  west1_ref->SetMarkerStyle(8);
  west1_ref->SetMarkerColor(kRed);
  // west1_ref->SetMinimum(0.);
  //west1_ref->SetMaximum(4096.);
  west1_ref->Draw("P SAME");



  cWest->cd(2);

  TGraphErrors *west2 = new TGraphErrors(numRuns[5],&runNumber[5][0],&westPedestal[1][0],0,&westPedestalRMS[1][0]);
  west2->SetTitle("PMT West 2");
  west2->GetXaxis()->SetTitle("Run Number");
  west2->GetYaxis()->SetTitle("ADC Channels");
  west2->GetXaxis()->SetLimits(runMin, runMax);
  west2->SetMarkerStyle(8);
  //west2->SetMinimum(0.);
  //west2->SetMaximum(4096.);
  west2->Draw("AP");

  cWest->Update();
  
  min = gPad->GetUymin();
  max = gPad->GetUymax();

  for (UInt_t i=0; i<calPeriodRunEnd.size(); i++) {
    
    linesW2[i] = new TLine((double)calPeriodRunEnd[i],min,(double)calPeriodRunEnd[i],max);
    linesW2[i]->SetLineStyle(2);
    linesW2[i]->SetLineColor(kRed);
    linesW2[i]->SetLineWidth(2);
    if ( (double)calPeriodRunEnd[i] > runNumber[5][0] && (double)calPeriodRunEnd[i] < runNumber[5][numRuns[5]-1] ) {
      linesW2[i]->Draw();
    }
  }

  TGraphErrors *west2_ref = new TGraphErrors(numRuns_ref[5],&runNumber_ref[5][0],&westPedestal_ref[1][0],0,&westPedestalRMS_ref[1][0]);
  west2_ref->SetTitle("PMT West 2");
  west2_ref->GetXaxis()->SetTitle("Run Number");
  west2_ref->GetYaxis()->SetTitle("ADC Channels");
  west2_ref->SetMarkerStyle(8);
  west2_ref->SetMarkerColor(kRed);
  // west2_ref->SetMinimum(0.);
  //west2_ref->SetMaximum(4096.);
  west2_ref->Draw("P SAME");



  cWest->cd(3);

  TGraphErrors *west3 = new TGraphErrors(numRuns[6],&runNumber[6][0],&westPedestal[2][0],0,&westPedestalRMS[2][0]);
  west3->SetTitle("PMT West 3");
  west3->GetXaxis()->SetTitle("Run Number");
  west3->GetYaxis()->SetTitle("ADC Channels");
  west3->GetXaxis()->SetLimits(runMin, runMax);
  west3->SetMarkerStyle(8);
  //west3->SetMinimum(0.);
  //west3->SetMaximum(4096.);
  west3->Draw("AP");

  cWest->Update();
  
  min = gPad->GetUymin();
  max = gPad->GetUymax();

  for (UInt_t i=0; i<calPeriodRunEnd.size(); i++) {
    
    linesW3[i] = new TLine((double)calPeriodRunEnd[i],min,(double)calPeriodRunEnd[i],max);
    linesW3[i]->SetLineStyle(2);
    linesW3[i]->SetLineColor(kRed);
    linesW3[i]->SetLineWidth(2);
    if ( (double)calPeriodRunEnd[i] > runNumber[6][0] && (double)calPeriodRunEnd[i] < runNumber[6][numRuns[6]-1] ) {
      linesW3[i]->Draw();
    }
  }

  TGraphErrors *west3_ref = new TGraphErrors(numRuns_ref[6],&runNumber_ref[6][0],&westPedestal_ref[2][0],0,&westPedestalRMS_ref[2][0]);
  west3_ref->SetTitle("PMT West 3");
  west3_ref->GetXaxis()->SetTitle("Run Number");
  west3_ref->GetYaxis()->SetTitle("ADC Channels");
  west3_ref->SetMarkerStyle(8);
  west3_ref->SetMarkerColor(kRed);
  // west3_ref->SetMinimum(0.);
  //west3_ref->SetMaximum(4096.);
  west3_ref->Draw("P SAME");



  if (numRuns[7]>0) {
    cWest->cd(4);
    
    
    TGraphErrors *west4 = new TGraphErrors(numRuns[7],&runNumber[7][0],&westPedestal[3][0],0,&westPedestalRMS[3][0]);
    west4->SetTitle("PMT West 4");
    west4->GetXaxis()->SetTitle("Run Number");
    west4->GetYaxis()->SetTitle("ADC Channels");
    west4->GetXaxis()->SetLimits(runMin, runMax);
    west4->SetMarkerStyle(8);
    //west4->SetMinimum(0.);
    //west4->SetMaximum(4096.);
    west4->Draw("AP");

    cWest->Update();
  
    min = gPad->GetUymin();
    max = gPad->GetUymax();
    
    for (UInt_t i=0; i<calPeriodRunEnd.size(); i++) {
      
      linesW4[i] = new TLine((double)calPeriodRunEnd[i],min,(double)calPeriodRunEnd[i],max);
      linesW4[i]->SetLineStyle(2);
      linesW4[i]->SetLineColor(kRed);
      linesW4[i]->SetLineWidth(2);
      if ( (double)calPeriodRunEnd[i] > runNumber[7][0] && (double)calPeriodRunEnd[i] < runNumber[7][numRuns[7]-1] ) {
	linesW4[i]->Draw();
      }
    }

    TGraphErrors *west4_ref = new TGraphErrors(numRuns_ref[7],&runNumber_ref[7][0],&westPedestal_ref[3][0],0,&westPedestalRMS_ref[3][0]);
    west4_ref->SetTitle("PMT West 4");
    west4_ref->GetXaxis()->SetTitle("Run Number");
    west4_ref->GetYaxis()->SetTitle("ADC Channels");
    west4_ref->SetMarkerStyle(8);
    west4_ref->SetMarkerColor(kRed);
    // west4_ref->SetMinimum(0.);
    //west4_ref->SetMaximum(4096.);
    west4_ref->Draw("P SAME");

    
  }
    

  TString fnamebase = TString::Format("%s_pedestals.pdf",year.Data());
  cEast->Print(TString::Format("%s(",fnamebase.Data()));
  cWest->Print(TString::Format("%s)",fnamebase.Data()));
  
  

}
