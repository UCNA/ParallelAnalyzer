#include "MBUtils.hh"
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

#include <TMath.h>
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TPDF.h>
#include <TLegend.h>
#include <TPaveText.h>

std::vector <Int_t> badOct {7,60,61,62,63,64,65,66}; 

int main(int argc, char *argv[])
{
  if (argc!=6) {
    std::cout << "Usage: ./PositionComparisons.exe [octet start] [octet end] [bool simulation] [bool quadrants=true (false will do radial)] [analysisChoice=C] \n";
    std::cout << "The code will produce comparisons for every octet in the range given,\non an octet-by-octet basis, and as a whole, using the Super-Sum\n";
    exit(0);
  }

  gStyle->SetOptStat(0);

  int octetNumStart = atoi(argv[1]);
  int octetNumEnd = atoi(argv[2]);
  bool sim = ( std::string(argv[3])=="true" ) ? true : false;
  std::string prefix = sim ? "SIM" : "UK";
  bool quad = ( std::string(argv[4])=="true" ) ? true : false;
  std::string asymmType = quad ? "Quadrants" : "Radial";

  std::string anaCh(argc==6?argv[5]:"C");
  
  //double enWinLow = (double) atof(argv[3]);
  //double enWinHigh = (double) atof(argv[4]);
  
  bool asymOn = false;
  
  double enBinWidth = 10.;
  int numBins = 1200/enBinWidth;
  

  //set up energy binning vector
  std::vector <double> enBins(numBins,0.);
  
  for (int bin=0; bin<numBins; bin++) {
    enBins[bin] = ((double)bin+.5)*(double)enBinWidth;
  }

  Double_t totalTimeON_E, totalTimeON_W, totalTimeOFF_E, totalTimeOFF_W;
  totalTimeON_E = totalTimeON_W = totalTimeOFF_E = totalTimeOFF_W = 0.;

  Double_t totalBGTimeON_E, totalBGTimeON_W, totalBGTimeOFF_E, totalBGTimeOFF_W;
  totalBGTimeON_E = totalBGTimeON_W = totalBGTimeOFF_E = totalBGTimeOFF_W = 0.;
 
  //Count vectors for holding all of the total (BG subtracted) counts for each spin state

  
  UInt_t numAsymms = quad?4:6;

  std::vector <std::vector <double> > W_TotalCountsON(numAsymms, std::vector <double>(numBins,0.)); 
  std::vector <std::vector <double> > W_TotalErrorON(numAsymms, std::vector <double>(numBins,0.));
  std::vector <std::vector <double> > W_TotalCountsOFF(numAsymms, std::vector <double>(numBins,0.)); 
  std::vector <std::vector <double> > W_TotalErrorOFF(numAsymms, std::vector <double>(numBins,0.));

  std::vector <std::vector <double> > E_TotalCountsON(numAsymms, std::vector <double>(numBins,0.)); 
  std::vector <std::vector <double> > E_TotalErrorON(numAsymms, std::vector <double>(numBins,0.));
  std::vector <std::vector <double> > E_TotalCountsOFF(numAsymms, std::vector <double>(numBins,0.)); 
  std::vector <std::vector <double> > E_TotalErrorOFF(numAsymms, std::vector <double>(numBins,0.));

  std::vector <std::vector <double> > W_BGTotalCountsON(numAsymms, std::vector <double>(numBins,0.)); 
  std::vector <std::vector <double> > W_BGTotalErrorON(numAsymms, std::vector <double>(numBins,0.));
  std::vector <std::vector <double> > W_BGTotalCountsOFF(numAsymms, std::vector <double>(numBins,0.)); 
  std::vector <std::vector <double> > W_BGTotalErrorOFF(numAsymms, std::vector <double>(numBins,0.));

  std::vector <std::vector <double> > E_BGTotalCountsON(numAsymms, std::vector <double>(numBins,0.)); 
  std::vector <std::vector <double> > E_BGTotalErrorON(numAsymms, std::vector <double>(numBins,0.));
  std::vector <std::vector <double> > E_BGTotalCountsOFF(numAsymms, std::vector <double>(numBins,0.)); 
  std::vector <std::vector <double> > E_BGTotalErrorOFF(numAsymms, std::vector <double>(numBins,0.));

  
  

  for (int octetNum=octetNumStart; octetNum<octetNumEnd+1; octetNum++) {
    
    if (std::find(badOct.begin(), badOct.end(),octetNum) != badOct.end()) { continue; } //Checking if octet should be ignored for data quality reasons

    // read in runs in octet and then open the rate files and fill the vectors according
    // to the type of run
    
    std::ifstream octFile(TString::Format("%s/All_Octets/octet_list_%i.dat",
					  getenv("OCTET_LIST"),octetNum));
    Double_t en, fg, bg, time, bg_time;
    std::string txt_hold;
    Int_t bin = 0;
    
    Int_t runNumber = 0;
    std::string runType = "";
    
    while ( octFile >> runType >> runNumber ) {
      
      std::ifstream rateFile;
      
      /////////////// Spin flipper off///////////////////
      if ( runType == "A2" || runType == "A10" || runType == "B5" || runType == "B7" ) {	
	
	//////////////// EAST RATES //////////////////
	
	// All event types
	rateFile.open(TString::Format("../BinByBinComparison/%s_Erun%i_anaCh%s_%s.dat", prefix.c_str(),runNumber, anaCh.c_str(),asymmType.c_str()));
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
	bin = 0; 
	
	std::vector <Double_t> fgrate(numAsymms,0);
	std::vector <Double_t> bgrate(numAsymms,0);
		
	while ( rateFile >> en )  {
	  std::cout << en << "\t";
	  for (UInt_t i=0; i<numAsymms; ++i) {
	    rateFile >> fg >> bg;
	    std::cout << fg << "\t" << bg << "\t";
	    E_TotalCountsOFF[i][bin] += (fg - bg)*time;
	    E_BGTotalCountsOFF[i][bin] += (bg)*bg_time;
	  }
	  std::cout << std::endl;
	  ++bin;
	}

	totalTimeOFF_E += time;
	totalBGTimeOFF_E += bg_time;
	rateFile.close();

	//////////////// WEST RATES //////////////////
	
	// All event types
	rateFile.open(TString::Format("../BinByBinComparison/%s_Wrun%i_anaCh%s_%s.dat", prefix.c_str(),runNumber, anaCh.c_str(),asymmType.c_str()));
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
	bin = 0; 
	
		
	while ( rateFile >> en )  {
	  std::cout << en << "\t";
	  for (UInt_t i=0; i<numAsymms; ++i) {
	    rateFile >> fg >> bg;
	    std::cout << fg << "\t" << bg << "\t";
	    W_TotalCountsOFF[i][bin] += (fg - bg)*time;
	    W_BGTotalCountsOFF[i][bin] += (bg)*bg_time;
	  }
	  std::cout << std::endl;
	  ++bin;
	}

	totalTimeOFF_W += time;
	totalBGTimeOFF_W += bg_time;
	rateFile.close();
	
      }
    
    
      /////////////// Spin flipper on///////////////////
      if ( runType == "B2" || runType == "B10" || runType == "A5" || runType == "A7" ) {	
	
	//////////////// EAST RATES //////////////////
	
	// All event types
	rateFile.open(TString::Format("../BinByBinComparison/%s_Erun%i_anaCh%s_%s.dat", prefix.c_str(),runNumber, anaCh.c_str(),asymmType.c_str()));
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
	bin = 0; 
	
	std::vector <Double_t> fgrate(numAsymms,0);
	std::vector <Double_t> bgrate(numAsymms,0);
	
	while ( rateFile >> en )  {
	  std::cout << en << "\t";
	  for (UInt_t i=0; i<numAsymms; ++i) {
	    rateFile >> fg >> bg;
	    std::cout << fg << "\t" << bg << "\t";
	    E_TotalCountsON[i][bin] += (fg - bg)*time;
	    E_BGTotalCountsON[i][bin] += (bg)*bg_time;
	  }
	  std::cout << std::endl;
	  ++bin;
	}
	
	totalTimeON_E += time;
	totalBGTimeON_E += bg_time;
	rateFile.close();
	
	//////////////// WEST RATES //////////////////
	
	// All event types
	rateFile.open(TString::Format("../BinByBinComparison/%s_Wrun%i_anaCh%s_%s.dat", prefix.c_str(),runNumber, anaCh.c_str(),asymmType.c_str()));
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
	bin = 0; 
	
	
	while ( rateFile >> en )  {
	  std::cout << en << "\t";
	  for (UInt_t i=0; i<numAsymms; ++i) {
	    rateFile >> fg >> bg;
	    std::cout << fg << "\t" << bg << "\t";
	    W_TotalCountsON[i][bin] += (fg - bg)*time;
	    W_BGTotalCountsON[i][bin] += (bg)*bg_time;
	  }
	  std::cout << std::endl;
	  ++bin;
	}
	
	totalTimeON_W += time;
	totalBGTimeON_W += bg_time;
	rateFile.close();
      }
    }
  }
  

  // Calculate the simple errors
  for ( unsigned int i = 0; i<numAsymms; ++i ) {
    for ( unsigned int b = 0; b<120; ++b ) {
      
      W_TotalErrorON[i][b] = sqrt( fabs(W_TotalCountsON[i][b]) );
      W_TotalErrorOFF[i][b] = sqrt( fabs(W_TotalCountsOFF[i][b]) );
      E_TotalErrorON[i][b] = sqrt( fabs(E_TotalCountsON[i][b]) );
      E_TotalErrorOFF[i][b] = sqrt( fabs(E_TotalCountsOFF[i][b]) );
      
      W_BGTotalErrorON[i][b] = sqrt( fabs(W_BGTotalCountsON[i][b]) );
      W_BGTotalErrorOFF[i][b] = sqrt( fabs(W_BGTotalCountsOFF[i][b]) );
      E_BGTotalErrorON[i][b] = sqrt( fabs(E_BGTotalCountsON[i][b]) );
      E_BGTotalErrorOFF[i][b] = sqrt( fabs(E_BGTotalCountsOFF[i][b]) );
      
    }
  }
  
  // Calculate the super-sum in every bin
  std::vector <std::vector <double> > SuperSum(numAsymms, std::vector <double>(numBins,0.)); 
  std::vector <std::vector <double> > SuperSumError(numAsymms, std::vector <double>(numBins,0.)); 

  std::vector <std::vector <double> > bg_SuperSum(numAsymms, std::vector <double>(numBins,0.)); 
  std::vector <std::vector <double> > bg_SuperSumError(numAsymms, std::vector <double>(numBins,0.)); 
  

  // First for data
  for ( unsigned int i=0; i<numAsymms; ++i ) {
    for ( unsigned int b = 0; b<120; ++b ) {
      
      double sfON[2]={0.};
      double sfOFF[2]={0.};
      double sfON_err[2]={0.};
      double sfOFF_err[2]={0.};
      
      
      // AFP Off
      sfOFF[0] = E_TotalCountsOFF[i][b] / totalTimeOFF_E;
      sfOFF_err[0] = E_TotalErrorOFF[i][b] / totalTimeOFF_E;
      
      sfOFF[1] = W_TotalCountsOFF[i][b] / totalTimeOFF_W;
      sfOFF_err[1] = W_TotalErrorOFF[i][b] / totalTimeOFF_W;
      
      // AFP ON
      sfON[0] = E_TotalCountsON[i][b] / totalTimeON_E;
      sfON_err[0] = E_TotalErrorON[i][b] / totalTimeON_E;
      
      sfON[1] = W_TotalCountsON[i][b] / totalTimeON_W;
      sfON_err[1] = W_TotalErrorON[i][b] / totalTimeON_W;
      
      
      
      //The geometric mean as defined by http://www.arpapress.com/Volumes/Vol11Issue3/IJRRAS_11_3_08.pdf
      double R1 = (sfON[1]>0. && sfOFF[0]>0.) ? sqrt(sfOFF[0]*sfON[1]) : (sfON[1]<0. && sfOFF[0]<0.) ? -sqrt(sfOFF[0]*sfON[1]) : 0.5*(sfOFF[0]+sfON[1]); 
      double R2 = (sfON[0]>0. && sfOFF[1]>0.) ? sqrt(sfOFF[1]*sfON[0]) : (sfON[0]<0. && sfOFF[1]<0.) ? -sqrt(sfOFF[1]*sfON[0]) : 0.5*(sfOFF[1]+sfON[0]);
      double deltaR1 = ( ( (sfON[1]>0. && sfOFF[0]>0.) || (sfON[1]<0. && sfOFF[0]<0.) ) ? 0.5*sqrt((power(sfOFF[0]*sfON_err[1],2)+power(sfON[1]*sfOFF_err[0],2))/(sfON[1]*sfOFF[0])) : 
			 0.5*sqrt( power(sfON_err[1],2) + power(sfOFF_err[0],2) ) );
      double deltaR2 = ((sfON[0]>0. && sfOFF[1]>0.) || (sfON[0]<0. && sfOFF[1]<0.)) ? 0.5*sqrt((power(sfOFF[1]*sfON_err[0],2)+power(sfON[0]*sfOFF_err[1],2))/(sfON[0]*sfOFF[1])) : 
	0.5*sqrt(power(sfON_err[0],2) + power(sfOFF_err[1],2));
      
      SuperSum[i][b] = R1 + R2;
      //std::cout << sfOFF[0] << " " << sfOFF_err[0] << std::endl;
      SuperSumError[i][b] = sqrt(power(deltaR1,2) + power(deltaR2,2));
      
      std::cout << SuperSum[i][b] << "\t" << SuperSumError[i][b] << "\n";
      
    }
  }

  // Now for bg_

  if ( !sim ) {
    for ( unsigned int i=0; i<numAsymms; ++i ) {
      for ( unsigned int b = 0; b<120; ++b ) {
	
	double sfON[2]={0.};
	double sfOFF[2]={0.};
	double sfON_err[2]={0.};
	double sfOFF_err[2]={0.};
	
	
	// AFP Off
	sfOFF[0] = E_BGTotalCountsOFF[i][b] / totalBGTimeOFF_E;
	sfOFF_err[0] = E_BGTotalErrorOFF[i][b] / totalBGTimeOFF_E;
	
	sfOFF[1] = W_BGTotalCountsOFF[i][b] / totalBGTimeOFF_W;
	sfOFF_err[1] = W_BGTotalErrorOFF[i][b] / totalBGTimeOFF_W;
	
	// AFP ON
	sfON[0] = E_BGTotalCountsON[i][b] / totalBGTimeON_E;
	sfON_err[0] = E_BGTotalErrorON[i][b] / totalBGTimeON_E;
	
	sfON[1] = W_BGTotalCountsON[i][b] / totalBGTimeON_W;
	sfON_err[1] = W_BGTotalErrorON[i][b] / totalBGTimeON_W;
	
	
	
	//The geometric mean as defined by http://www.arpapress.com/Volumes/Vol11Issue3/IJRRAS_11_3_08.pdf
	double R1 = (sfON[1]>0. && sfOFF[0]>0.) ? sqrt(sfOFF[0]*sfON[1]) : (sfON[1]<0. && sfOFF[0]<0.) ? -sqrt(sfOFF[0]*sfON[1]) : 0.5*(sfOFF[0]+sfON[1]); 
	double R2 = (sfON[0]>0. && sfOFF[1]>0.) ? sqrt(sfOFF[1]*sfON[0]) : (sfON[0]<0. && sfOFF[1]<0.) ? -sqrt(sfOFF[1]*sfON[0]) : 0.5*(sfOFF[1]+sfON[0]);
	double deltaR1 = ( ( (sfON[1]>0. && sfOFF[0]>0.) || (sfON[1]<0. && sfOFF[0]<0.) ) ? 0.5*sqrt((power(sfOFF[0]*sfON_err[1],2)+power(sfON[1]*sfOFF_err[0],2))/(sfON[1]*sfOFF[0])) : 
			   0.5*sqrt( power(sfON_err[1],2) + power(sfOFF_err[0],2) ) );
	double deltaR2 = ((sfON[0]>0. && sfOFF[1]>0.) || (sfON[0]<0. && sfOFF[1]<0.)) ? 0.5*sqrt((power(sfOFF[1]*sfON_err[0],2)+power(sfON[0]*sfOFF_err[1],2))/(sfON[0]*sfOFF[1])) : 
	  0.5*sqrt(power(sfON_err[0],2) + power(sfOFF_err[1],2));
	
	bg_SuperSum[i][b] = R1 + R2;
	//std::cout << sfOFF[0] << " " << sfOFF_err[0] << std::endl;
	bg_SuperSumError[i][b] = sqrt(power(deltaR1,2) + power(deltaR2,2));
	
	//std::cout << SuperSum[i][b] << "\t" << SuperSumError[i][b] << "\n";
	
      }
    }
  }
  
  TFile *f2 = new TFile(TString::Format("Octets_%i-%i_%s_anaCh%s_%s.root", octetNumStart, octetNumEnd, sim?"SIM":"DATA", anaCh.c_str(), asymmType.c_str()),"RECREATE");

  std::vector <TH1D*> histsFG;
  std::vector <TH1D*> histsBG;
 

  for (UInt_t i=0; i<numAsymms; ++i) {
    //TH1D * hold = new TH1D(TString::Format("%s_%i",asymmType.c_str(),i),TString::Format("%s_%i",asymmType.c_str(),i),numBins,0.,1200.);
    histsFG.push_back(new TH1D(TString::Format("FG_%s_%i",asymmType.c_str(),i),TString::Format("FG_%s_%i",asymmType.c_str(),i),numBins,0.,1200.));
    histsFG[i]->GetXaxis()->SetTitle("Energy (keV)");
    
    histsBG.push_back(new TH1D(TString::Format("BG_%s_%i",asymmType.c_str(),i),TString::Format("BG_%s_%i",asymmType.c_str(),i),numBins,0.,1200.));
    histsBG[i]->GetXaxis()->SetTitle("Energy (keV)");
  }
  
  for (int i=0; i<numAsymms; ++i) {
    for (int bin=0; bin<numBins; bin++) {
      histsFG[i]->SetBinContent(bin+1,100.*SuperSum[i][bin]);      
      if (!sim) histsBG[i]->SetBinContent(bin+1,100.*bg_SuperSum[i][bin]);
      
    }
    
  }
    
  f2->Write();
  f2->Close();
  
  
  return 0;
  
}
      
