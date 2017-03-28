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

std::vector <Int_t> badOct = {7,9,59,60,61,62,63,64,65,66}; 

int main(int argc, char *argv[])
{
  if (argc!=4) {
    std::cout << "Usage: ./EastWestComparisons.exe [octet start] [octet end] [bool simulation]\n";
    std::cout << "The code will produce comparisons for every octet in the range given,\non an octet-by-octet basis, and as a whole, using the Super-Sum\n";
    exit(0);
  }

  gStyle->SetOptStat(0);

  int octetNumStart = atoi(argv[1]);
  int octetNumEnd = atoi(argv[2]);
  bool sim = ( std::string(argv[3])=="true" ) ? true : false;
  std::string prefix = sim ? "SIM" : "UK";
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
  
  //Count vectors for holding all of the total (BG subtracted) counts for each spin state

  //type0=0, type1=1, type23=2, ALL=3
  std::vector <std::vector <double> > W_TotalCountsON(5, std::vector <double>(numBins,0.)); 
  std::vector <std::vector <double> > W_TotalErrorON(5, std::vector <double>(numBins,0.));
  std::vector <std::vector <double> > W_TotalCountsOFF(5, std::vector <double>(numBins,0.)); 
  std::vector <std::vector <double> > W_TotalErrorOFF(5, std::vector <double>(numBins,0.));

  std::vector <std::vector <double> > E_TotalCountsON(5, std::vector <double>(numBins,0.)); 
  std::vector <std::vector <double> > E_TotalErrorON(5, std::vector <double>(numBins,0.));
  std::vector <std::vector <double> > E_TotalCountsOFF(5, std::vector <double>(numBins,0.)); 
  std::vector <std::vector <double> > E_TotalErrorOFF(5, std::vector <double>(numBins,0.));
  

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
	rateFile.open(TString::Format("../BinByBinComparison/%s_Erun%i_anaChC.dat", prefix.c_str(),runNumber));
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
	bin = 0; 
	
	while ( rateFile >> en >> fg >> bg )  { 
	  E_TotalCountsOFF[4][bin] += (fg - bg)*time;
	  ++bin;
	}

	totalTimeOFF_E += time;
	rateFile.close();
	
	
	// Type 0 only
	rateFile.open(TString::Format("../BinByBinComparison/%s_Erun%i_anaChD.dat", prefix.c_str(), runNumber));
	
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
        bin = 0;
	
	while ( rateFile >> en >> fg >> bg )  { 
	  E_TotalCountsOFF[0][bin] += (fg - bg)*time;
	  ++bin;
	}
	
	rateFile.close();
	
	// Type 1 only
	rateFile.open(TString::Format("../BinByBinComparison/%s_Erun%i_anaChF.dat", prefix.c_str(), runNumber));
	
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
        bin = 0;
	
	while ( rateFile >> en >> fg >> bg )  { 
	  E_TotalCountsOFF[1][bin] += (fg - bg)*time;
	  ++bin;
	}
	
	rateFile.close();
	
	// Type 2 only
	rateFile.open(TString::Format("../BinByBinComparison/%s_Erun%i_anaChJ.dat", prefix.c_str(), runNumber));
	
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
        bin = 0;
	
	while ( rateFile >> en >> fg >> bg )  { 
	  E_TotalCountsOFF[2][bin] += (fg - bg)*time;
	  ++bin;
	}
	
	rateFile.close();

	// Type 3 only
	rateFile.open(TString::Format("../BinByBinComparison/%s_Erun%i_anaChK.dat", prefix.c_str(), runNumber));
	
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
        bin = 0;
	
	while ( rateFile >> en >> fg >> bg )  { 
	  E_TotalCountsOFF[3][bin] += (fg - bg)*time;
	  ++bin;
	}
	
	rateFile.close();
      
      
      
	///////////////// WEST RATES ///////////////////
	
	// All event types
	rateFile.open(TString::Format("../BinByBinComparison/%s_Wrun%i_anaChC.dat", prefix.c_str(), runNumber));
	
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
	bin = 0; 
	
	while ( rateFile >> en >> fg >> bg )  { 
	  W_TotalCountsOFF[4][bin] += (fg - bg)*time;
	  ++bin;
	}
	
	totalTimeOFF_W += time;
	rateFile.close();
	
	
	// Type 0 only
	rateFile.open(TString::Format("../BinByBinComparison/%s_Wrun%i_anaChD.dat", prefix.c_str(), runNumber));
	
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
        bin = 0;
	
	while ( rateFile >> en >> fg >> bg )  { 
	  W_TotalCountsOFF[0][bin] += (fg - bg)*time;
	  ++bin;
	}
	
	rateFile.close();
	
	// Type 1 only
	rateFile.open(TString::Format("../BinByBinComparison/%s_Wrun%i_anaChF.dat", prefix.c_str(), runNumber));
	
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
        bin = 0;
	
	while ( rateFile >> en >> fg >> bg )  { 
	  W_TotalCountsOFF[1][bin] += (fg - bg)*time;
	  ++bin;
	}
	
	rateFile.close();
	
	// Type 2 only
	rateFile.open(TString::Format("../BinByBinComparison/%s_Wrun%i_anaChJ.dat", prefix.c_str(), runNumber));
	
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
        bin = 0;
	
	while ( rateFile >> en >> fg >> bg )  { 
	  W_TotalCountsOFF[2][bin] += (fg - bg)*time;
	  ++bin;
	}
	
	rateFile.close();
      

      // Type 3 only
	rateFile.open(TString::Format("../BinByBinComparison/%s_Wrun%i_anaChK.dat", prefix.c_str(), runNumber));
	
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
        bin = 0;
	
	while ( rateFile >> en >> fg >> bg )  { 
	  W_TotalCountsOFF[3][bin] += (fg - bg)*time;
	  ++bin;
	}
	
	rateFile.close();
      }
      
      
      /////////////// Spin flipper ON ///////////////////
      if ( runType == "B2" || runType == "B10" || runType == "A5" || runType == "A7" ) {	
	
	//////////////// EAST RATES //////////////////
	
	// All event types
	rateFile.open(TString::Format("../BinByBinComparison/%s_Erun%i_anaChC.dat", prefix.c_str(), runNumber));
	
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
	bin = 0; 

	while ( rateFile >> en >> fg >> bg )  { 
	  E_TotalCountsON[4][bin] += (fg - bg)*time;
	  ++bin;
	}

	totalTimeON_E += time;
	rateFile.close();


	// Type 0 only
	rateFile.open(TString::Format("../BinByBinComparison/%s_Erun%i_anaChD.dat", prefix.c_str(), runNumber));
	
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;

        bin = 0;

	while ( rateFile >> en >> fg >> bg )  { 
	  E_TotalCountsON[0][bin] += (fg - bg)*time;
	  ++bin;
	}

	rateFile.close();

	// Type 1 only
	rateFile.open(TString::Format("../BinByBinComparison/%s_Erun%i_anaChF.dat", prefix.c_str(), runNumber));
	
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;

        bin = 0;

	while ( rateFile >> en >> fg >> bg )  { 
	  E_TotalCountsON[1][bin] += (fg - bg)*time;
	  ++bin;
	}

	rateFile.close();

	// Type 2 only
	rateFile.open(TString::Format("../BinByBinComparison/%s_Erun%i_anaChJ.dat", prefix.c_str(), runNumber));
	
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;

        bin = 0;

	while ( rateFile >> en >> fg >> bg )  { 
	  E_TotalCountsON[2][bin] += (fg - bg)*time;
	  ++bin;
	}

	rateFile.close();

	// Type 3 only
	rateFile.open(TString::Format("../BinByBinComparison/%s_Erun%i_anaChK.dat", prefix.c_str(), runNumber));
	
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;

        bin = 0;

	while ( rateFile >> en >> fg >> bg )  { 
	  E_TotalCountsON[3][bin] += (fg - bg)*time;
	  ++bin;
	}

	rateFile.close();
   

      
	///////////////// WEST RATES ///////////////////
	
	// All event types
	rateFile.open(TString::Format("../BinByBinComparison/%s_Wrun%i_anaChC.dat", prefix.c_str(), runNumber));
	
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
	bin = 0; 
	
	while ( rateFile >> en >> fg >> bg )  { 
	  W_TotalCountsON[4][bin] += (fg - bg)*time;
	  ++bin;
	}
	
	totalTimeON_W += time;
	rateFile.close();
	
	
	// Type 0 only
	rateFile.open(TString::Format("../BinByBinComparison/%s_Wrun%i_anaChD.dat", prefix.c_str(), runNumber));
	
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
	bin = 0;
	
	while ( rateFile >> en >> fg >> bg )  { 
	  W_TotalCountsON[0][bin] += (fg - bg)*time;
	  ++bin;
	}
	
	rateFile.close();
	
	// Type 1 only
	rateFile.open(TString::Format("../BinByBinComparison/%s_Wrun%i_anaChF.dat", prefix.c_str(), runNumber));
	
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
	bin = 0;
	
	while ( rateFile >> en >> fg >> bg )  { 
	  W_TotalCountsON[1][bin] += (fg - bg)*time;
	  ++bin;
	}
	
	rateFile.close();
	
	// Type 2 only
	rateFile.open(TString::Format("../BinByBinComparison/%s_Wrun%i_anaChJ.dat", prefix.c_str(), runNumber));
	
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
	bin = 0;
	
	while ( rateFile >> en >> fg >> bg )  { 
	  W_TotalCountsON[2][bin] += (fg - bg)*time;
	  ++bin;
	}
	
	rateFile.close();

	// Type 3 only
	rateFile.open(TString::Format("../BinByBinComparison/%s_Wrun%i_anaChK.dat", prefix.c_str(), runNumber));
	
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
	bin = 0;
	
	while ( rateFile >> en >> fg >> bg )  { 
	  W_TotalCountsON[3][bin] += (fg - bg)*time;
	  ++bin;
	}
	
	rateFile.close();
      }

    }
  }

  // Calculate the simple errors
  for ( unsigned int i = 0; i<5; ++i ) {
    for ( unsigned int b = 0; b<120; ++b ) {
      
      W_TotalErrorON[i][b] = sqrt( fabs(W_TotalCountsON[i][b]) );
      W_TotalErrorOFF[i][b] = sqrt( fabs(W_TotalCountsOFF[i][b]) );
      E_TotalErrorON[i][b] = sqrt( fabs(E_TotalCountsON[i][b]) );
      E_TotalErrorOFF[i][b] = sqrt( fabs(E_TotalCountsOFF[i][b]) );
      
    }
  }
  
  // Calculate the super-sum in every bin
  std::vector <std::vector <double> > data_SuperSum(5, std::vector <double>(numBins,0.)); 
  std::vector <std::vector <double> > data_SuperSumError(5, std::vector <double>(numBins,0.)); 
  
  for ( unsigned int i=0; i<5; ++i ) {
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
      
      data_SuperSum[i][b] = R1 + R2;
      //std::cout << sfOFF[0] << " " << sfOFF_err[0] << std::endl;
      data_SuperSumError[i][b] = sqrt(power(deltaR1,2) + power(deltaR2,2));
      
      std::cout << data_SuperSum[i][b] << "\t" << data_SuperSumError[i][b] << "\n";
      
    }
  }
  
  TFile *f2 = new TFile(TString::Format("Octets_%i-%i_%s.root", octetNumStart, octetNumEnd, sim?"SIM":"DATA"),"RECREATE");
  
  TH1D Erecon0("Erecon0","Type 0 Super-Sum",numBins,0.,1200.);
  Erecon0.GetXaxis()->SetTitle("Energy (keV)");
  TH1D Erecon1("Erecon1","Type 1 Super-Sum",numBins,0.,1200.);
  Erecon1.GetXaxis()->SetTitle("Energy (keV)");
  TH1D Erecon2("Erecon2","Type 2 Super-Sum",numBins,0.,1200.);
  Erecon2.GetXaxis()->SetTitle("Energy (keV)");
  TH1D Erecon3("Erecon3","Type 3 Super-Sum",numBins,0.,1200.);
  Erecon3.GetXaxis()->SetTitle("Energy (keV)");
  TH1D EreconALL("EreconALL","Type ALL Super-Sum",numBins,0.,1200.);
  EreconALL.GetXaxis()->SetTitle("Energy (keV)");
  
  for (int bin=0; bin<numBins; bin++) {
    Erecon0.SetBinContent(bin+1,100.*data_SuperSum[0][bin]);
    Erecon1.SetBinContent(bin+1,100.*data_SuperSum[1][bin]);
    Erecon2.SetBinContent(bin+1,100.*data_SuperSum[2][bin]);
    Erecon3.SetBinContent(bin+1,100.*data_SuperSum[3][bin]);
    EreconALL.SetBinContent(bin+1,100.*data_SuperSum[4][bin]);
    
    Erecon0.SetBinError(bin+1,100.*data_SuperSumError[0][bin]);
    Erecon1.SetBinError(bin+1,100.*data_SuperSumError[1][bin]);
    Erecon2.SetBinError(bin+1,100.*data_SuperSumError[2][bin]);
    Erecon3.SetBinError(bin+1,100.*data_SuperSumError[3][bin]);
    EreconALL.SetBinError(bin+1,100.*data_SuperSumError[4][bin]);
    
  }

  TH1D Erecon0_sfOFF_E("Erecon0_sfOFF_E","East Type 0 Spin Flipper OFF",numBins,0.,1200.);
  Erecon0_sfOFF_E.GetXaxis()->SetTitle("Energy (keV)");
  TH1D Erecon1_sfOFF_E("Erecon1_sfOFF_E","East Type 1 Spin Flipper OFF",numBins,0.,1200.);
  Erecon1_sfOFF_E.GetXaxis()->SetTitle("Energy (keV)");
  TH1D Erecon2_sfOFF_E("Erecon2_sfOFF_E","East Type 2 Spin Flipper OFF",numBins,0.,1200.);
  Erecon2_sfOFF_E.GetXaxis()->SetTitle("Energy (keV)");
  TH1D Erecon3_sfOFF_E("Erecon3_sfOFF_E","East Type 3 Spin Flipper OFF",numBins,0.,1200.);
  Erecon3_sfOFF_E.GetXaxis()->SetTitle("Energy (keV)");
  TH1D EreconALL_sfOFF_E("EreconALL_sfOFF_E","East Type ALL Spin Flipper OFF",numBins,0.,1200.);
  EreconALL_sfOFF_E.GetXaxis()->SetTitle("Energy (keV)");
  
  for (int bin=0; bin<numBins; bin++) {
    Erecon0_sfOFF_E.SetBinContent(bin+1,100.*E_TotalCountsOFF[0][bin]/totalTimeOFF_E);
    Erecon1_sfOFF_E.SetBinContent(bin+1,100.*E_TotalCountsOFF[1][bin]/totalTimeOFF_E);
    Erecon2_sfOFF_E.SetBinContent(bin+1,100.*E_TotalCountsOFF[2][bin]/totalTimeOFF_E);
    Erecon3_sfOFF_E.SetBinContent(bin+1,100.*E_TotalCountsOFF[3][bin]/totalTimeOFF_E);
    EreconALL_sfOFF_E.SetBinContent(bin+1,100.*E_TotalCountsOFF[4][bin]/totalTimeOFF_E);
    
    Erecon0_sfOFF_E.SetBinError(bin+1,100.*E_TotalErrorOFF[0][bin]/totalTimeOFF_E);
    Erecon1_sfOFF_E.SetBinError(bin+1,100.*E_TotalErrorOFF[1][bin]/totalTimeOFF_E);
    Erecon2_sfOFF_E.SetBinError(bin+1,100.*E_TotalErrorOFF[2][bin]/totalTimeOFF_E);
    Erecon3_sfOFF_E.SetBinError(bin+1,100.*E_TotalErrorOFF[3][bin]/totalTimeOFF_E);
    EreconALL_sfOFF_E.SetBinError(bin+1,100.*E_TotalErrorOFF[4][bin]/totalTimeOFF_E);
    
  }

  TH1D Erecon0_sfOFF_W("Erecon0_sfOFF_W","West Type 0 Spin Flipper OFF",numBins,0.,1200.);
  Erecon0_sfOFF_W.GetXaxis()->SetTitle("Energy (keV)");
  TH1D Erecon1_sfOFF_W("Erecon1_sfOFF_W","West Type 1 Spin Flipper OFF",numBins,0.,1200.);
  Erecon1_sfOFF_W.GetXaxis()->SetTitle("Energy (keV)");
  TH1D Erecon2_sfOFF_W("Erecon2_sfOFF_W","West Type 2 Spin Flipper OFF",numBins,0.,1200.);
  Erecon2_sfOFF_W.GetXaxis()->SetTitle("Energy (keV)");
  TH1D Erecon3_sfOFF_W("Erecon3_sfOFF_W","West Type 3 Spin Flipper OFF",numBins,0.,1200.);
  Erecon3_sfOFF_W.GetXaxis()->SetTitle("Energy (keV)");
  TH1D EreconALL_sfOFF_W("EreconALL_sfOFF_W","West Type ALL Spin Flipper OFF",numBins,0.,1200.);
  EreconALL_sfOFF_W.GetXaxis()->SetTitle("Energy (keV)");
  
  for (int bin=0; bin<numBins; bin++) {
    Erecon0_sfOFF_W.SetBinContent(bin+1,100.*W_TotalCountsOFF[0][bin]/totalTimeOFF_W);
    Erecon1_sfOFF_W.SetBinContent(bin+1,100.*W_TotalCountsOFF[1][bin]/totalTimeOFF_W);
    Erecon2_sfOFF_W.SetBinContent(bin+1,100.*W_TotalCountsOFF[2][bin]/totalTimeOFF_W);
    Erecon3_sfOFF_W.SetBinContent(bin+1,100.*W_TotalCountsOFF[3][bin]/totalTimeOFF_W);
    EreconALL_sfOFF_W.SetBinContent(bin+1,100.*W_TotalCountsOFF[4][bin]/totalTimeOFF_W);
    
    Erecon0_sfOFF_W.SetBinError(bin+1,100.*W_TotalErrorOFF[0][bin]/totalTimeOFF_W);
    Erecon1_sfOFF_W.SetBinError(bin+1,100.*W_TotalErrorOFF[1][bin]/totalTimeOFF_W);
    Erecon2_sfOFF_W.SetBinError(bin+1,100.*W_TotalErrorOFF[2][bin]/totalTimeOFF_W);
    Erecon3_sfOFF_W.SetBinError(bin+1,100.*W_TotalErrorOFF[3][bin]/totalTimeOFF_W);
    EreconALL_sfOFF_W.SetBinError(bin+1,100.*W_TotalErrorOFF[4][bin]/totalTimeOFF_W);
    
  }

  TH1D Erecon0_sfON_E("Erecon0_sfON_E","East Type 0 Spin Flipper ON",numBins,0.,1200.);
  Erecon0_sfON_E.GetXaxis()->SetTitle("Energy (keV)");
  TH1D Erecon1_sfON_E("Erecon1_sfON_E","East Type 1 Spin Flipper ON",numBins,0.,1200.);
  Erecon1_sfON_E.GetXaxis()->SetTitle("Energy (keV)");
  TH1D Erecon2_sfON_E("Erecon2_sfON_E","East Type 2 Spin Flipper ON",numBins,0.,1200.);
  Erecon2_sfON_E.GetXaxis()->SetTitle("Energy (keV)");
  TH1D Erecon3_sfON_E("Erecon3_sfON_E","East Type 3 Spin Flipper ON",numBins,0.,1200.);
  Erecon3_sfON_E.GetXaxis()->SetTitle("Energy (keV)");
  TH1D EreconALL_sfON_E("EreconALL_sfON_E","East Type ALL Spin Flipper ON",numBins,0.,1200.);
  EreconALL_sfON_E.GetXaxis()->SetTitle("Energy (keV)");
  
  for (int bin=0; bin<numBins; bin++) {
    Erecon0_sfON_E.SetBinContent(bin+1,100.*E_TotalCountsON[0][bin]/totalTimeON_E);
    Erecon1_sfON_E.SetBinContent(bin+1,100.*E_TotalCountsON[1][bin]/totalTimeON_E);
    Erecon2_sfON_E.SetBinContent(bin+1,100.*E_TotalCountsON[2][bin]/totalTimeON_E);
    Erecon3_sfON_E.SetBinContent(bin+1,100.*E_TotalCountsON[3][bin]/totalTimeON_E);
    EreconALL_sfON_E.SetBinContent(bin+1,100.*E_TotalCountsON[4][bin]/totalTimeON_E);
    
    Erecon0_sfON_E.SetBinError(bin+1,100.*E_TotalErrorON[0][bin]/totalTimeON_E);
    Erecon1_sfON_E.SetBinError(bin+1,100.*E_TotalErrorON[1][bin]/totalTimeON_E);
    Erecon2_sfON_E.SetBinError(bin+1,100.*E_TotalErrorON[2][bin]/totalTimeON_E);
    Erecon3_sfON_E.SetBinError(bin+1,100.*E_TotalErrorON[3][bin]/totalTimeON_E);
    EreconALL_sfON_E.SetBinError(bin+1,100.*E_TotalErrorON[4][bin]/totalTimeON_E);
    
  }

  TH1D Erecon0_sfON_W("Erecon0_sfON_W","West Type 0 Spin Flipper ON",numBins,0.,1200.);
  Erecon0_sfON_W.GetXaxis()->SetTitle("Energy (keV)");
  TH1D Erecon1_sfON_W("Erecon1_sfON_W","West Type 1 Spin Flipper ON",numBins,0.,1200.);
  Erecon1_sfON_W.GetXaxis()->SetTitle("Energy (keV)");
  TH1D Erecon2_sfON_W("Erecon2_sfON_W","West Type 2 Spin Flipper ON",numBins,0.,1200.);
  Erecon2_sfON_W.GetXaxis()->SetTitle("Energy (keV)");
  TH1D Erecon3_sfON_W("Erecon3_sfON_W","West Type 3 Spin Flipper ON",numBins,0.,1200.);
  Erecon3_sfON_W.GetXaxis()->SetTitle("Energy (keV)");
  TH1D EreconALL_sfON_W("EreconALL_sfON_W","West Type ALL Spin Flipper ON",numBins,0.,1200.);
  EreconALL_sfON_W.GetXaxis()->SetTitle("Energy (keV)");
  
  for (int bin=0; bin<numBins; bin++) {
    Erecon0_sfON_W.SetBinContent(bin+1,100.*W_TotalCountsON[0][bin]/totalTimeON_W);
    Erecon1_sfON_W.SetBinContent(bin+1,100.*W_TotalCountsON[1][bin]/totalTimeON_W);
    Erecon2_sfON_W.SetBinContent(bin+1,100.*W_TotalCountsON[2][bin]/totalTimeON_W);
    Erecon3_sfON_W.SetBinContent(bin+1,100.*W_TotalCountsON[3][bin]/totalTimeON_W);
    EreconALL_sfON_W.SetBinContent(bin+1,100.*W_TotalCountsON[4][bin]/totalTimeON_W);
    
    Erecon0_sfON_W.SetBinError(bin+1,100.*W_TotalErrorON[0][bin]/totalTimeON_W);
    Erecon1_sfON_W.SetBinError(bin+1,100.*W_TotalErrorON[1][bin]/totalTimeON_W);
    Erecon2_sfON_W.SetBinError(bin+1,100.*W_TotalErrorON[2][bin]/totalTimeON_W);
    Erecon3_sfON_W.SetBinError(bin+1,100.*W_TotalErrorON[3][bin]/totalTimeON_W);
    EreconALL_sfON_W.SetBinError(bin+1,100.*W_TotalErrorON[4][bin]/totalTimeON_W);
    
  }
    
  f2->Write();
  f2->Close();
  
  
  return 0;
  
}
      
