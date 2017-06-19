
#include <TMinuit.h>
#include <vector>

std::vector <Int_t> getPMTQuality(Int_t runNumber) {
  //Read in PMT quality file
  cout << "Reading in PMT Quality file ...\n";
  vector <Int_t>  pmtQuality (8,0);
  Char_t temp[200];
  sprintf(temp,"%s/residuals/PMT_runQuality_master.dat",getenv("ANALYSIS_CODE")); 
  ifstream pmt;
  std::cout << temp << std::endl;
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
};

vector <vector <double> > returnSourcePosition (Int_t runNumber, string src) {
  Char_t temp[500];
//cout << "made it here in returnSourcePosition\n";   
  sprintf(temp,"%s/source_list_%i.dat",getenv("SOURCE_LIST"),runNumber);
  ifstream file(temp);
  cout << src << endl;
  int num = 0;
  file >> num;
  cout << num << endl;
  int srcNum = 0;
  string src_check;
  for (int i=0; i<num;srcNum++,i++) {
    file >> src_check;
    cout << src_check << endl;
    if (src_check==src) break;   
  }
  cout << "The source Number is: " << srcNum << endl;
  if (srcNum==num) {
    cout << "Didn't find source in that run\n"; exit(0);
  }
  file.close();
  
  sprintf(temp,"%s/source_positions_%i.dat",getenv("SOURCE_POSITIONS"),runNumber);
  file.open(temp);
  
  vector < vector < double > > srcPos;
  srcPos.resize(2,vector <double> (3,0.));
  
  for (int i=0; i<srcNum+1; i++) {
    for (int j=0; j<2; j++) {
      for (int jj=0; jj<3; jj++) {
	file >> srcPos[j][jj];
      }
    }
  }
  return srcPos;
}

bool isSourceInFidCut(Int_t run, string src, Double_t fidCut, Int_t side) {
//cout << "made it here in isSourceInFidCut\n";   
  std::vector < std::vector <Double_t> > pos = returnSourcePosition(run,src);

  return ( fidCut*fidCut > ( pos[side][0]*pos[side][0] + pos[side][1]*pos[side][1] ) ) ? true : false;

}


void width_fitter(Int_t calPeriod)
{

  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetTitleYOffset(1.5);

  //Read in sim and data widths

  std::vector < std::vector < Int_t> > simWidthRun(8,std::vector <Int_t>(0));
  std::vector < std::vector < Int_t> > dataWidthRun(8,std::vector <Int_t> (0) );
  std::vector < std::vector < Int_t> > simWidthErrorRun(8,std::vector <Int_t> (0) );
  std::vector < std::vector < Int_t> > dataWidthErrorRun(8,std::vector <Int_t> (0) );
  std::vector < std::vector < std::string> > simWidthSrc(8,std::vector <std::string> (0) );
  std::vector < std::vector < std::string> > dataWidthSrc(8,std::vector <std::string> (0) );
  std::vector < std::vector < std::string> > simWidthErrorSrc(8,std::vector <std::string> (0) );
  std::vector < std::vector < std::string> > dataWidthErrorSrc(8,std::vector <std::string> (0) );
  std::vector < std::vector < Double_t > > simWidths(8,std::vector <Double_t> (0));
  std::vector < std::vector < Double_t > > dataWidths(8,std::vector <Double_t> (0));
  std::vector < std::vector < Double_t > > simWidthErrors(8,std::vector <Double_t> (0));
  std::vector < std::vector < Double_t > > dataWidthErrors(8,std::vector <Double_t> (0));
 
  std::vector < Double_t > old_k(8,0.);
  std::vector < Double_t > new_k(8,0.);
  std::vector < Double_t > slope(8,0.);

  Char_t tempfile[200];
  sprintf(tempfile, "%s/simulation_comparison/nPE_per_keV/nPE_per_keV_%i.dat",getenv("ANALYSIS_CODE"),calPeriod);
  ifstream kfilein(tempfile);
 
  Int_t k=0;
  std::cout << "OLD k VALUES:\n";
  while (kfilein >> old_k[k]) {
    std::cout << old_k[k] << std::endl;
    k++;
  } 
  kfilein.close();


  // Read in all residuals files
  Double_t p[8] ={0.}; 
  Int_t run=0;
  std::string nm = "";

  sprintf(tempfile,"%s/residuals/source_runs_EvisWidth_RunPeriod_%i.dat",getenv("ANALYSIS_CODE"),calPeriod);  
  ifstream dataWidthsFile(tempfile); 
  if (dataWidthsFile.is_open()) cout << tempfile << endl;

  while ( dataWidthsFile >> run >> nm >> p[0] >> p[1] >> p[2] >> p[3] >> p[4] >> p[5] >> p[6] >> p[7] ) {
    
    if ( nm==std::string("In") ) continue;
    if ( nm==std::string("Cd") ) continue;
    if ( nm==std::string("Cs") ) continue;
    if ( nm==std::string("Bi2") ) continue;

      for ( UInt_t i=0; i<8; ++i ) {  dataWidthRun[i].push_back(run); dataWidthSrc[i].push_back(nm); dataWidths[i].push_back(p[i]); }
    
  }
	   

  sprintf(tempfile,"%s/residuals/source_runs_EvisWidthError_RunPeriod_%i.dat",getenv("ANALYSIS_CODE"),calPeriod);  
  ifstream dataWidthErrorsFile(tempfile);
  if (dataWidthErrorsFile.is_open()) cout << tempfile << endl;
  while ( dataWidthErrorsFile >> run >> nm >> p[0] >> p[1] >> p[2] >> p[3] >> p[4] >> p[5] >> p[6] >> p[7] ) {
    if ( nm==std::string("In") ) continue;
    if ( nm==std::string("Cd") ) continue;
    if ( nm==std::string("Cs") ) continue;
    if ( nm==std::string("Bi2") ) continue;

    for ( UInt_t i=0; i<8; ++i ) {  dataWidthErrorRun[i].push_back(run); dataWidthErrorSrc[i].push_back(nm); dataWidthErrors[i].push_back(p[i]); }
      
  }

  sprintf(tempfile,"%s/residuals/SIM_source_runs_EvisWidth_RunPeriod_%i.dat",getenv("ANALYSIS_CODE"),calPeriod);  
  ifstream simWidthsFile(tempfile);
  if (simWidthsFile.is_open()) cout << tempfile << endl;
  while ( simWidthsFile >> run >> nm >> p[0] >> p[1] >> p[2] >> p[3] >> p[4] >> p[5] >> p[6] >> p[7] ) {
    if ( nm==std::string("In") ) continue;
    if ( nm==std::string("Cd") ) continue;
    if ( nm==std::string("Cs") ) continue;
    if ( nm==std::string("Bi2") ) continue;

for ( UInt_t i=0; i<8; ++i ) {  simWidthRun[i].push_back(run); simWidthSrc[i].push_back(nm); simWidths[i].push_back(p[i]); }
     
  }

  sprintf(tempfile,"%s/residuals/SIM_source_runs_EvisWidthError_RunPeriod_%i.dat",getenv("ANALYSIS_CODE"),calPeriod);  
  ifstream simWidthErrorsFile(tempfile); 
  if (simWidthErrorsFile.is_open()) cout << tempfile << endl;
  while ( simWidthErrorsFile >> run >> nm >> p[0] >> p[1] >> p[2] >> p[3] >> p[4] >> p[5] >> p[6] >> p[7] ) {
    if ( nm==std::string("In") ) continue;
    if ( nm==std::string("Cd") ) continue;
    if ( nm==std::string("Cs") ) continue;
    if ( nm==std::string("Bi2") ) continue;

    for ( UInt_t i=0; i<8; ++i ) {  simWidthErrorRun[i].push_back(run); simWidthErrorSrc[i].push_back(nm); simWidthErrors[i].push_back(p[i]); }
    
  }
  
  dataWidthsFile.close();
  simWidthsFile.close();
  dataWidthErrorsFile.close();
  simWidthErrorsFile.close();
  
  std::vector < std::vector < Double_t > > finalSimWidths(8,std::vector <Double_t> (0));
  std::vector < std::vector < Double_t > > finalDataWidths(8,std::vector <Double_t> (0));
  std::vector < std::vector < Double_t > > finalSimWidthErrors(8,std::vector <Double_t> (0));
  std::vector < std::vector < Double_t > > finalDataWidthErrors(8,std::vector <Double_t> (0));

  
  for ( UInt_t i=0; i<dataWidthRun[0].size(); ++i ) {
    
    Double_t fiducialCut = 35.;
   
    bool passFidCutEast = isSourceInFidCut(dataWidthRun[0][i], dataWidthSrc[0][i].substr(0,2), fiducialCut, 0);
    bool passFidCutWest = isSourceInFidCut(dataWidthRun[0][i], dataWidthSrc[0][i].substr(0,2), fiducialCut, 1);
    

    for (UInt_t j=0; j<8; j++) {
      
      if ( j<4 && !passFidCutEast ) { 
	std::cout << "didn't pass east fiducial cut\n";
	continue;
      }

      else if ( j>3 && ( !passFidCutWest || ( j==5 && dataWidthRun[j][i]>16983 && dataWidthRun[j][i]<=17249 ) ) ) { 
	std::cout << "didn't pass west fiducial cut or pmt is bad\n";
	continue;  
      }
      
      if ( (TMath::Abs(simWidths[j][i]-dataWidths[j][i])/dataWidths[j][i] < 1.5) ) { //Checking for outliers
	std::cout << dataWidthRun[j][i] << " PMT " << j << std::endl; 
	finalDataWidths[j].push_back(dataWidths[j][i]);
	finalSimWidths[j].push_back(simWidths[j][i]);
	finalDataWidthErrors[j].push_back(dataWidthErrors[j][i]);
	finalSimWidthErrors[j].push_back(simWidthErrors[j][i]);
      }
    }
  }
  

	    

  

  cout << i << endl;

  TCanvas *c1 = new TCanvas("c1", "c1", 1400, 1000);

  TPad *p0 = new TPad("p0","East 1", 0.0, 0.5, 0.25, 1.0);
  TPad *p1 = new TPad("p1","East 2", 0.25, 0.5, 0.5, 1.0);
  TPad *p2 = new TPad("p2","East 3", 0.5, 0.5, 0.75, 1.0);
  TPad *p3 = new TPad("p3","East 4", 0.75, 0.5, 1.0, 1.0);
  TPad *p4 = new TPad("p4","West 1", 0.0, 0., 0.25, .5);
  TPad *p5 = new TPad("p5","West 2", 0.25, 0., 0.5, .5);
  TPad *p6 = new TPad("p6","West 3", 0.5, 0., 0.75, .5);
  TPad *p7 = new TPad("p7","West 4", 0.75, 0., 1.0, .5);
  p0->Draw();
  p1->Draw();
  p2->Draw();
  p3->Draw();
  p4->Draw(); 
  p5->Draw();
  p6->Draw(); 
  p7->Draw(); 


  TF1 *f1 = new TF1("f1","[0]*x",0., 300.);
  f1->SetParameter(0,1.);
  f1->SetParLimits(0,0.5,1.5);
  f1->SetLineStyle(2);
  f1->SetLineWidth(1);
  gStyle->SetOptFit();
  
  TString status = " ";
  Int_t nAttempts = 5;
  Int_t attempt = 0;

  Double_t maxRange = 200.;

  
  p0->cd();

  if (finalDataWidths[0].size()>0) {

    TGraphErrors *pmt0 = new TGraphErrors(finalDataWidths[0].size(), &finalSimWidths[0][0], &finalDataWidths[0][0],
					  &finalSimWidthErrors[0][0], &finalDataWidthErrors[0][0]);
    pmt0->SetTitle("PMT East 1");
    pmt0->SetMarkerColor(1);
    pmt0->SetLineColor(1);
    pmt0->SetMarkerStyle(20);
    pmt0->SetMarkerSize(0.25);
    pmt0->GetXaxis()->SetTitle("Simulated Width (keV)");
    pmt0->GetYaxis()->SetTitle("Actual Width (keV)");
    pmt0->SetMinimum(0.0);
    pmt0->SetMaximum(TMath::MaxElement(finalDataWidths[0].size(),&finalDataWidths[0][0])+50.);
    pmt0->GetXaxis()->SetLimits(0.0,TMath::MaxElement(finalDataWidths[0].size(),&finalSimWidths[0][0])+50.);   
    pmt0->Draw("AP");
    
    while (status!=TString("CONVERGED ") && attempt!=nAttempts) {    
      pmt0->Fit("f1");
      status = gMinuit->fCstatu;
      attempt++;
    }
    if (status!=TString("CONVERGED ")) slope[0] = 1.;
    else slope[0] = f1->GetParameter(0);
    
  }
  
  else slope[0]=1.;

  status = " ";
  attempt = 0;
  f1->SetParameter(0,1.);

  p1->cd();

  if (finalDataWidths[1].size()>0) {
    TGraphErrors *pmt1 = new TGraphErrors(finalDataWidths[1].size(), &finalSimWidths[1][0], &finalDataWidths[1][0],
					  &finalSimWidthErrors[1][0], &finalDataWidthErrors[1][0]);
    pmt1->SetTitle("PMT East 2");
    pmt1->SetMarkerColor(1);
    pmt1->SetLineColor(1);
    pmt1->SetMarkerStyle(20);
    pmt1->SetMarkerSize(0.25);
    pmt1->GetXaxis()->SetTitle("Simulated Width (keV)");
    pmt1->GetYaxis()->SetTitle("Actual Width (keV)");
    pmt1->SetMinimum(0.0);
    pmt1->SetMaximum(TMath::MaxElement(finalDataWidths[1].size(),&finalDataWidths[1][0])+50.);
    pmt1->GetXaxis()->SetLimits(0.0,TMath::MaxElement(finalDataWidths[1].size(),&finalSimWidths[1][0])+50.);   
    pmt1->Draw("AP");
    
    while (status!=TString("CONVERGED ") && attempt!=nAttempts) {    
      pmt1->Fit("f1");
      status = gMinuit->fCstatu;
      attempt++;
    }
    if (status!=TString("CONVERGED ")) slope[1] = 1.;
    else slope[1] = f1->GetParameter(0);
  }

  else slope[1]=1.;

  status = " ";
  attempt = 0;
  f1->SetParameter(0,1.);

  p2->cd();


  if (finalDataWidths[2].size()>0) {
    TGraphErrors *pmt2 = new TGraphErrors(finalDataWidths[2].size(), &finalSimWidths[2][0], &finalDataWidths[2][0],
					  &finalSimWidthErrors[2][0], &finalDataWidthErrors[2][0]);   
    pmt2->SetTitle("PMT East 3");
    pmt2->SetMarkerColor(1);
    pmt2->SetLineColor(1);
    pmt2->SetMarkerStyle(20);
    pmt2->SetMarkerSize(0.25);
    pmt2->GetXaxis()->SetTitle("Simulated Width (keV)");
    pmt2->GetYaxis()->SetTitle("Actual Width (keV)");  
    pmt2->SetMinimum(0.0);
    pmt2->SetMaximum(TMath::MaxElement(finalDataWidths[2].size(),&finalDataWidths[2][0])+50.);
    pmt2->GetXaxis()->SetLimits(0.0,TMath::MaxElement(finalDataWidths[2].size(),&finalSimWidths[2][0])+50.);   
    pmt2->Draw("AP");
    
    while (status!=TString("CONVERGED ") && attempt!=nAttempts) {    
      pmt2->Fit("f1");
      status = gMinuit->fCstatu;
      attempt++;
    }
    if (status!=TString("CONVERGED ")) slope[2] = 1.;
    else slope[2] = f1->GetParameter(0);
  }

  else slope[2]=1.;

  status = " ";
  attempt = 0;
  f1->SetParameter(0,1.);


  p3->cd();

  if (finalDataWidths[3].size()>0) {

    TGraphErrors *pmt3 = new TGraphErrors(finalDataWidths[3].size(), &finalSimWidths[3][0], &finalDataWidths[3][0],
					  &finalSimWidthErrors[3][0], &finalDataWidthErrors[3][0]);
    pmt3->SetTitle("PMT East 4");
    pmt3->SetMarkerColor(1);
    pmt3->SetLineColor(1);
    pmt3->SetMarkerStyle(20);
    pmt3->SetMarkerSize(0.25);
    pmt3->GetXaxis()->SetTitle("Simulated Width (keV)");
    pmt3->GetYaxis()->SetTitle("Actual Width (keV)");
    pmt3->SetMinimum(0.0);
    pmt3->SetMaximum(TMath::MaxElement(finalDataWidths[3].size(),&finalDataWidths[3][0])+50.);
    pmt3->GetXaxis()->SetLimits(0.0,TMath::MaxElement(finalDataWidths[3].size(),&finalSimWidths[3][0])+50.);
    pmt3->Draw("AP");
    
    while (status!=TString("CONVERGED ") && attempt!=nAttempts) {    
      pmt3->Fit("f1");
      status = gMinuit->fCstatu;
      attempt++;
    }
    if (status!=TString("CONVERGED ")) slope[3] = 1.;
    else slope[3] = f1->GetParameter(0);

  }
  
  else slope[3]=1.;

  status = " ";
  attempt = 0;
  f1->SetParameter(0,1.);


  p4->cd();

  if (finalDataWidths[4].size()>0) {
    
    TGraphErrors *pmt4 = new TGraphErrors(finalDataWidths[4].size(), &finalSimWidths[4][0], &finalDataWidths[4][0],
					  &finalSimWidthErrors[4][0], &finalDataWidthErrors[4][0]);
    pmt4->SetTitle("PMT West 1");
    pmt4->SetMarkerColor(1);
    pmt4->SetLineColor(1);
    pmt4->SetMarkerStyle(20);
    pmt4->SetMarkerSize(0.25);
    pmt4->GetXaxis()->SetTitle("Simulated Width (keV)");
    pmt4->GetYaxis()->SetTitle("Actual Width (keV)");
    pmt4->SetMinimum(0.0);
    pmt4->SetMaximum(TMath::MaxElement(finalDataWidths[4].size(),&finalDataWidths[4][0])+50.);
    pmt4->GetXaxis()->SetLimits(0.0,TMath::MaxElement(finalDataWidths[4].size(),&finalSimWidths[4][0])+50.);   
    pmt4->Draw("AP");
    
    while (status!=TString("CONVERGED ") && attempt!=nAttempts) {    
      pmt4->Fit("f1");
      status = gMinuit->fCstatu;
      attempt++;
    }
    if (status!=TString("CONVERGED ")) slope[4] = 1.;
    else slope[4] = f1->GetParameter(0);
  }

  else slope[4]=1.;

  status = " ";
  attempt = 0;
  f1->SetParameter(0,1.);


  p5->cd();

  if (finalDataWidths[5].size()>0) {

    TGraphErrors *pmt5 = new TGraphErrors(finalDataWidths[5].size(), &finalSimWidths[5][0], &finalDataWidths[5][0],
					  &finalSimWidthErrors[5][0], &finalDataWidthErrors[5][0]);
    pmt5->SetTitle("PMT West 2");
    pmt5->SetMarkerColor(1);
    pmt5->SetLineColor(1);
    pmt5->SetMarkerStyle(20);
    pmt5->SetMarkerSize(0.25);
    pmt5->GetXaxis()->SetTitle("Simulated Width (keV)");
    pmt5->GetYaxis()->SetTitle("Actual Width (keV)");
    pmt5->SetMinimum(0.0);
    pmt5->SetMaximum(TMath::MaxElement(finalDataWidths[5].size(),&finalDataWidths[5][0])+50.);
    pmt5->GetXaxis()->SetLimits(0.0,TMath::MaxElement(finalDataWidths[5].size(),&finalSimWidths[5][0])+50.);   
    pmt5->Draw("AP");
    
    while (status!=TString("CONVERGED ") && attempt!=nAttempts) {    
      pmt5->Fit("f1");
      status = gMinuit->fCstatu;
      attempt++;
    }
    if (status!=TString("CONVERGED ")) slope[5] = 1.;
    else slope[5] = f1->GetParameter(0);
  }

  else slope[5]=1.;

  status = " ";
  attempt = 0;
  f1->SetParameter(0,1.);


  p6->cd();

  if (finalDataWidths[6].size()>0) {

    TGraphErrors *pmt6 = new TGraphErrors(finalDataWidths[6].size(), &finalSimWidths[6][0], &finalDataWidths[6][0],
					  &finalSimWidthErrors[6][0], &finalDataWidthErrors[6][0]);
    pmt6->SetTitle("PMT West 3");
    pmt6->SetMarkerColor(1);
    pmt6->SetLineColor(1);
    pmt6->SetMarkerStyle(20);
    pmt6->SetMarkerSize(0.25);
    pmt6->GetXaxis()->SetTitle("Simulated Width (keV)");
    pmt6->GetYaxis()->SetTitle("Actual Width (keV)");  
    pmt6->SetMinimum(0.0);
    pmt6->SetMaximum(TMath::MaxElement(finalDataWidths[6].size(),&finalDataWidths[6][0])+50.);
    pmt6->GetXaxis()->SetLimits(0.0,TMath::MaxElement(finalDataWidths[6].size(),&finalSimWidths[6][0])+50.);   
    pmt6->Draw("AP");
    
    while (status!=TString("CONVERGED ") && attempt!=nAttempts) {    
    pmt6->Fit("f1");
    status = gMinuit->fCstatu;
    attempt++;
    }
    if (status!=TString("CONVERGED ")) slope[6] = 1.;
    else slope[6] = f1->GetParameter(0);
  }
  
  else slope[6]=1.;
  
  status = " ";
  attempt = 0;
  f1->SetParameter(0,1.);


  p7->cd();

  if (finalDataWidths[7].size()>0) {
    
    TGraphErrors *pmt7 = new TGraphErrors(finalDataWidths[7].size(), &finalSimWidths[7][0], &finalDataWidths[7][0],
					  &finalSimWidthErrors[7][0], &finalDataWidthErrors[7][0]);
    pmt7->SetTitle("PMT West 4");
    pmt7->SetMarkerColor(1);
    pmt7->SetLineColor(1);
    pmt7->SetMarkerStyle(20);
    pmt7->SetMarkerSize(0.25);
    pmt7->GetXaxis()->SetTitle("Simulated Width (keV)");
    pmt7->GetYaxis()->SetTitle("Actual Width (keV)");
    pmt7->SetMinimum(0.0);
    pmt7->SetMaximum(TMath::MaxElement(finalDataWidths[7].size(),&finalDataWidths[7][0])+50.);
    pmt7->GetXaxis()->SetLimits(0.0,TMath::MaxElement(finalDataWidths[7].size(),&finalSimWidths[7][0])+50.);   
    pmt7->Draw("AP");
    
    while (status!=TString("CONVERGED ") && attempt!=nAttempts) {    
      pmt7->Fit("f1");
      status = gMinuit->fCstatu;
      attempt++;
    }
    if (status!=TString("CONVERGED ")) slope[7] = 1.;
    else slope[7] = f1->GetParameter(0);
  }
  
  else slope[7]=1.;


  TString pdffile = TString::Format("%s/simulation_comparison/nPE_per_keV/width_comp_%i.pdf",getenv("ANALYSIS_CODE"),calPeriod);
  c1->Print(pdffile);


  //Calculate new k 

  sprintf(tempfile,"%s/simulation_comparison/nPE_per_keV/nPE_per_keV_%i.dat",getenv("ANALYSIS_CODE"),calPeriod);
  ofstream new_k_fileout(tempfile);

  sprintf(tempfile,"%s/simulation_comparison/nPE_per_keV/prev_nPE_per_keV_%i.dat",getenv("ANALYSIS_CODE"),calPeriod);
  ofstream old_k_fileout(tempfile);

  for (UInt_t jj=0; jj<8; jj++) {
    if (jj!=0) {
      new_k_fileout << std::endl;
      old_k_fileout << std::endl;
    }
    std::cout << old_k[jj] << " " << slope[jj] << std::endl;

    new_k[jj] = old_k[jj]/(slope[jj]*slope[jj]);
    new_k_fileout << new_k[jj];
    old_k_fileout << old_k[jj];
  }
  old_k_fileout.close();
  new_k_fileout.close();
  
}
