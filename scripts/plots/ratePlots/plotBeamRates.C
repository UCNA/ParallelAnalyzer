// Script to plot the beam rates both before and after making the beam cuts. 
// Also will print out the total time cut out, and how many runs this affected

#include <vector>

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

void plotBeamRates(TString year) {

  //read in all beta decay runs and their UCN rates
  
  std::vector <Int_t> badOct {7,9,59,60,61,62,63,64,65,66,70,92};;
  //Int_t octs[] = {7,9,59,60,61,62,63,64,65,66,70,92};
  //badOct.assign(octs, octs+12);
  
  Int_t octMin, octMax;

  if (year==TString("2011-2012") ) {
    octMin = 0;
    octMax = 59;
  }

  else {
    octMin = 67;
    octMax = 121;
  }

  std::vector <Float_t> runs;
  
  for ( Int_t oct=octMin; oct<octMax+1; ++oct ) {
    
    std::vector<Float_t> octruns = readOctetFileForFGruns(oct);
    runs.insert(runs.end(),octruns.begin(),octruns.end());

  }
  
  // Read in the beam cut files for all runs
  std::vector <Float_t> ratesBefore; 
  std::vector <Float_t> ratesAfter; 
  std::vector <Int_t> runsWithCuts; //all runs with a beam drop

  for ( UInt_t i = 0; i<runs.size(); i++ ) {

    ifstream ifile(TString::Format("%s/beamCuts_%0.0f.dat",getenv("CUTS"),runs[i]).Data());
    
    std::string holdTxt;
    std::string rateBefore, rateAfter;
    Int_t numCuts;

    ifile >> holdTxt >> rateBefore >> holdTxt
	  >> holdTxt >> rateAfter  >> holdTxt
	  >> holdTxt >> numCuts;
    
    ifile.close();
    
    ratesBefore.push_back(rateBefore!="-nan"?atof(rateBefore.c_str()):0.);
    ratesAfter.push_back(rateAfter!="-nan"?atof(rateAfter.c_str()):0.);
    
    if ( numCuts ) runsWithCuts.push_back(runs[i]);

  }

  TCanvas *c1 = new TCanvas("c1","c1", 1600, 1000);
  c1->Divide(1,2);

  c1->cd(1);

  TGraph *g1 = new TGraph(ratesBefore.size(), &runs[0], &ratesBefore[0]);
  g1->SetTitle("UCN Monitor 1 Before Cutting Beam Drops");
  g1->GetXaxis()->SetTitle("run number");
  g1->GetYaxis()->SetTitle("Average Rate (Hz)");
  g1->SetMinimum(0.);
  g1->SetMaximum(28.);
  g1->SetMarkerStyle(20);
  g1->SetMarkerColor(kBlue);
  g1->SetMarkerSize(0.75);
  g1->Draw("AP");
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0001);

  TF1 *f1 = new TF1("f1","[0]",runs[0]-5., runs[runs.size()-1]+5.);
  f1->SetParameter(0,16.);

  g1->Fit(f1,"R");

  c1->cd(2);

  TGraph *g2 = new TGraph(ratesBefore.size(), &runs[0], &ratesAfter[0]);
  g2->SetTitle("UCN Monitor 1 After Cutting Beam Drops");
  g2->GetXaxis()->SetTitle("run number");
  g2->GetYaxis()->SetTitle("Average Rate (Hz)");
  g2->SetMinimum(0.);
  g2->SetMaximum(28.);
  g2->SetMarkerStyle(20);
  g2->SetMarkerColor(kBlue);
  g2->SetMarkerSize(0.75);
  g2->Draw("AP");
  
  TF1 *f2 = new TF1("f2","[0]",runs[0]-5., runs[runs.size()-1]+5.);
  f2->SetParameter(0,16.);

  g2->Fit(f2,"R");
  
  std::cout << "Runs affected by cut: " << runsWithCuts.size() 
	    << " out of " << runs.size() << std::endl;

};
