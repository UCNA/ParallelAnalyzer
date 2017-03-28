#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "DataTree.hh"

// ROOT libraries
#include "TRandom3.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TH1D.h>



int main(int argc, char *argv[]) {

  Int_t runNumber = atoi(argv[1]);
  Float_t rateCut = 2.; // Hz

  std::cout << "Calculating beam drops for run " << runNumber << std::endl;

  
  DataTree *t = new DataTree();
  
  TString infile = TString::Format("%s/replay_pass1_%i.root",
				   getenv("REPLAY_PASS1"),runNumber);
  t->setupInputTree(infile.Data(),"pass1");

  // Look through histogram and find ranges that need to be cut.

  std::vector < std::vector <Float_t> > ranges;

  Int_t nBins = t->UCN_Mon_1_Rate->GetNbinsX();
  Int_t binCounter = 0;

  // Variables to calculate average before and after cut
  Float_t aveBefore = 0.;
  Float_t aveAfter = 0.;
  Int_t binsAfter = 0;

  for ( Int_t b=1; b<=nBins-2; b++ ) { // The minus 2 removes the last 2 empty bins

    Float_t rate = t->UCN_Mon_1_Rate->GetBinContent(b);

    aveBefore+=rate;
    
    if ( rate < rateCut ) {
      binCounter++;
    }

    if ( rate > rateCut || b == (nBins-3) ) {

      if (binCounter<3) {
	aveAfter += rate;
        binsAfter += 1;
      }
      
      if ( binCounter > 1 ) { // We'll record this chunk

	// We go 3 bins back from first bin under rateCut because the rate seems
	// to fall slowly. It then rises quickly, so we use the upper edge of the first
	// good rate bin to come back on
	Float_t lowCut = ( (b-binCounter)-3 ) > 0 ? t->UCN_Mon_1_Rate->GetXaxis()->GetBinLowEdge((b-binCounter)-3) : 0.;
	Float_t UpCut = t->UCN_Mon_1_Rate->GetXaxis()->GetBinUpEdge(b);
	std::vector <Float_t> minMax {lowCut, UpCut};
	ranges.push_back(minMax);
      }
      binCounter = 0;
    }
  }

  aveBefore = aveBefore/(nBins-2);
  aveAfter = aveAfter/binsAfter;
  
  delete t;

  std::cout << "Average Before cuts: " << aveBefore << " Hz\n";
  std::cout << "Average After cuts: " << aveAfter << " Hz\n\n";
  std::cout << "Beam Cuts:\n";

  for ( UInt_t i=0; i<ranges.size(); ++i ) {
    std::cout << ranges[i][0] << " -> " << ranges[i][1] << std::endl;
  }

  std::cout << "//////////////////////////////////\n\n";

  TString outPath = TString::Format("%s/beamCuts_%i.dat",getenv("CUTS"),runNumber);
  std::ofstream ofile(outPath.Data());

  ofile << "Average_Before_cuts: " << aveBefore << " Hz\n";
  ofile << "Average_After_cuts: " << aveAfter << " Hz\n";
  ofile << "Num_Beam_Cuts: " << ranges.size() << "\n";

  for ( UInt_t i=0; i<ranges.size(); ++i ) {
    ofile << ranges[i][0] << " " << ranges[i][1] << "\n";
  }

  ofile.close();
  return 0;
};
    
	
	

      

    
    
