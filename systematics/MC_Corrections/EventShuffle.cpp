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
#include <TF1.h>

#include "BetaSpectrum.hh"

TString anaChoices[10] = {"A","B","C","D","E","F","G","H","J","K"};

//Returns the values of the systematic corrections and statistical error on this for now
std::vector < std::vector <Double_t> >  LoadOctetSystematics(Int_t octet, std::string anaChoice, std::vector <Double_t> enBinMidpoint);

//Returns a vector containing all the theory corrections to A0 for a particular bin
std::vector <Double_t> LoadTheoryCorrections(std::vector <Double_t> enBinMidpoint);

// Function to return beta when given the kinetic energy of an electron            
Double_t returnBeta(Double_t En) {
  Double_t me = 510.998928; //rest mass energy of electron in keV     
  return sqrt(En*En+2.*En*me)/(En+me);
};


void calcBinByBinSuperRatio(std::vector< std::vector<double> > &A2, 
			    std::vector< std::vector<double> > &A5,
			    std::vector< std::vector<double> > &A7, 
			    std::vector< std::vector<double> > &A10,
			    std::vector< std::vector<double> > &B2, 
			    std::vector< std::vector<double> > &B5,
			    std::vector< std::vector<double> > &B7,
			    std::vector< std::vector<double> > &B10,
			    std::vector< std::vector<double> > &A2err, 
			    std::vector< std::vector<double> > &A5err,
			    std::vector< std::vector<double> > &A7err, 
			    std::vector< std::vector<double> > &A10err,
			    std::vector< std::vector<double> > &B2err, 
			    std::vector< std::vector<double> > &B5err,
			    std::vector< std::vector<double> > &B7err, 
			    std::vector< std::vector<double> > &B10err,
			    std::vector<double> &asymmetry, std::vector<double> &asymmetryError) {

//arrays to hold the weighted averages of the rates for each spin state and side
  double sfON[2]={0.};
  double sfOFF[2]={0.};
  double sfON_err[2]={0.};
  double sfOFF_err[2]={0.};

  for (unsigned int bin=0; bin<asymmetry.size(); bin++) {

    double R = 0.;
    double deltaR = 0.;

    for (unsigned int side=0; side<2; side++) {

      double weightsum=0.; // Holds the sum of the weights for the denominator of the weighted av
      
      //Start with flipper off... Only use rate if it has an error to avoid dividing by 0
      sfOFF[side] = ( ( A2err[side][bin] > 0. ? A2[side][bin]/power(A2err[side][bin],2)   : 0. ) +
		      ( A10err[side][bin]> 0. ? A10[side][bin]/power(A10err[side][bin],2) : 0. ) +
		      (	B5err[side][bin] > 0. ? B5[side][bin]/power(B5err[side][bin],2)   : 0. ) +
		      ( B7err[side][bin] > 0. ? B7[side][bin]/power(B7err[side][bin],2)   : 0. ) );

      weightsum = ( ( A2err[side][bin] > 0. ? 1./power(A2err[side][bin],2)  : 0. ) + 
		    ( A10err[side][bin]> 0. ? 1./power(A10err[side][bin],2) : 0. ) +
		    ( B5err[side][bin] > 0. ? 1./power(B5err[side][bin],2)  : 0. ) +
		    ( B7err[side][bin] > 0. ? 1./power(B7err[side][bin],2)  : 0. ) );

      sfOFF[side] = weightsum>0. ? sfOFF[side]/weightsum : 0.;
      sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;


      weightsum=0.;

      // Now for flipper on
      sfON[side] = ( ( B2err[side][bin] > 0. ? B2[side][bin]/power(B2err[side][bin],2)   : 0. ) +
		      ( B10err[side][bin]> 0. ? B10[side][bin]/power(B10err[side][bin],2) : 0. ) +
		      (	A5err[side][bin] > 0. ? A5[side][bin]/power(A5err[side][bin],2)   : 0. ) +
		      ( A7err[side][bin] > 0. ? A7[side][bin]/power(A7err[side][bin],2)   : 0. ) );

      weightsum = ( ( B2err[side][bin] > 0. ? 1./power(B2err[side][bin],2)  : 0. ) + 
		    ( B10err[side][bin]> 0. ? 1./power(B10err[side][bin],2) : 0. ) +
		    ( A5err[side][bin] > 0. ? 1./power(A5err[side][bin],2)  : 0. ) +
		    ( A7err[side][bin] > 0. ? 1./power(A7err[side][bin],2)  : 0. ) );
      
      sfON[side] = weightsum>0. ? sfON[side]/weightsum : 0.;
      sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
      

      
    }
    
    //Calculate super-ratio, R. Note that any of the rates that go into this could be negative 
    R = (sfOFF[1]*sfON[0]) / (sfON[1]*sfOFF[0]) ;
    
    //Only use rates which are real and not zero to avoid nan (couldn't get TMath::IsNaN() to work properly :( ...)
    if ( R!=0. && sfON[1]!=0. && sfOFF[0]!=0.) { 
     
      deltaR = sqrt( power( sfOFF_err[1]*sfON[0]/(sfOFF[0]*sfON[1]), 2 ) + 
		     power( sfON_err[0]*sfOFF[1]/(sfOFF[0]*sfON[1]), 2 ) +
		     power( sfOFF_err[0]*sfON[0]*sfOFF[1]/(sfOFF[0]*sfOFF[0]*sfON[1]), 2 ) +
		     power( sfON_err[1]*sfON[0]*sfOFF[1]/(sfOFF[0]*sfON[1]*sfON[1]), 2 ) );
      //deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
      asymmetry[bin] = (1.-sqrt(std::abs(R)))/(1+sqrt(std::abs(R)));
      asymmetryError[bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
    }
    else {
      asymmetry[bin] = 0.;
      asymmetryError[bin] = 1.;
    }

  } 
};


std::vector <Int_t> badOct {7,9,59,60,61,62,63,64,65,66,67,91,93,101,107,121};  

int main(int argc, char *argv[])
{
  if (argc!=9) {
    std::cout << "Usage: ./EventShuffle.exe [octet start] [octet end] [Emax] [Emin] [T1 percent shuffle] [T2 percent shuffle] [T3 percent Shuffle] [Lost event percent shuffle]\n";
    exit(0);
  }

  gStyle->SetOptStat(0);

  int octetNumStart = atoi(argv[1]);
  int octetNumEnd = atoi(argv[2]);
  double emin = atof(argv[3]);
  double emax = atof(argv[4]);  
  double percentT1 = atof(argv[5]);
  double percentT2 = atof(argv[6]);
  double percentT3 = atof(argv[7]);  
  double percentT4 = atof(argv[8]);  

  double enBinWidth = 10.;
  int numBins = 1200/enBinWidth;
  

  //set up energy binning vector
  std::vector <double> enBins(numBins,0.);
  
  for (int bin=0; bin<numBins; bin++) {
    enBins[bin] = ((double)bin+.5)*(double)enBinWidth;
  }

  

  // vectors to hold the final weighted average asymmetries in every bin for both
  // the Type0 asymm, and the asymm with an event shuffle
  std::vector <std::vector<Double_t> > Asymm(2,std::vector<Double_t>(numBins,0));
  std::vector <std::vector<Double_t> > AsymmWithShuffle(2,std::vector<Double_t>(numBins,0));
  
  

  for (int octetNum=octetNumStart; octetNum<octetNumEnd+1; octetNum++) {
    
    if (std::find(badOct.begin(), badOct.end(),octetNum) != badOct.end()) { continue; } //Checking if octet should be ignored for data quality reasons

    // read in runs in octet and then open the rate files and fill the vectors according
    // to the type of run
    
    std::ifstream octFile(TString::Format("%s/All_Octets/octet_list_%i.dat",
					  getenv("OCTET_LIST"),octetNum));
    Double_t en, fg, bg, time, bg_time;
    std::string txt_hold;
   

    std::vector < Int_t > runNumber;
    std::vector <std::string> runType;

    Int_t rn=0;
    std::string rt="";
    
    while ( octFile >> rt >> rn )   runNumber.push_back(rn), runType.push_back(rt);
   
    if ( runNumber.size()<8 ) continue;

    std::vector < std::vector<Double_t> > A2(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > B2(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > A5(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > B5(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > A7(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > B7(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > A10(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > B10(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > A2err(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > B2err(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > A5err(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > B5err(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > A7err(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > B7err(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > A10err(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > B10err(2,std::vector <double> (numBins,0.));

    std::vector < std::vector<Double_t> > A2_shift(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > B2_shift(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > A5_shift(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > B5_shift(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > A7_shift(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > B7_shift(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > A10_shift(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > B10_shift(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > A2err_shift(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > B2err_shift(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > A5err_shift(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > B5err_shift(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > A7err_shift(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > B7err_shift(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > A10err_shift(2,std::vector <double> (numBins,0.));
    std::vector < std::vector<Double_t> > B10err_shift(2,std::vector <double> (numBins,0.));
    
    for (UInt_t i=0;i<runNumber.size(); ++i) {
      
      std::ifstream rateFile;


      // East Type 0
      rateFile.open(TString::Format("%s/Asymmetry/BinByBinComparison/SIM_Erun%i_anaChD.dat", getenv("ANALYSIS_CODE"),runNumber[i]));
      rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
      rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
      
      int bin = 0; 
	
      while ( rateFile >> en >> fg >> bg )  { 
      
	if ( runType[i] == "A2" ) { 
	  A2[0][bin]=fg, A2err[0][bin]=TMath::Sqrt( fg/time );
	  A2_shift[0][bin]=fg, A2err_shift[0][bin]+= fg/time;
	}
	else if ( runType[i] == "A5" ) { 
	  A5[0][bin]=fg, A5err[0][bin]=TMath::Sqrt( fg/time );
	  A5_shift[0][bin]=fg, A5err_shift[0][bin]+= fg/time;
	}
	else if ( runType[i] == "A7" ) { 
	  A7[0][bin]=fg, A7err[0][bin]=TMath::Sqrt( fg/time );
	  A7_shift[0][bin]=fg, A7err_shift[0][bin]+= fg/time;
	}
	else if ( runType[i] == "A10" ) { 
	  A10[0][bin]=fg, A10err[0][bin]=TMath::Sqrt( fg/time );
	  A10_shift[0][bin]=fg, A10err_shift[0][bin]+= fg/time;
	}
	else if ( runType[i] == "B2" ) { 
	  B2[0][bin]=fg, B2err[0][bin]=TMath::Sqrt( fg/time );
	  B2_shift[0][bin]=fg, B2err_shift[0][bin]+= fg/time;
	}
	else if ( runType[i] == "B5" ) { 
	  B5[0][bin]=fg, B5err[0][bin]=TMath::Sqrt( fg/time );
	  B5_shift[0][bin]=fg, B5err_shift[0][bin]+= fg/time;
	}
	else if ( runType[i] == "B7" ) { 
	  B7[0][bin]=fg, B7err[0][bin]=TMath::Sqrt( fg/time );
	  B7_shift[0][bin]=fg, B7err_shift[0][bin]+= fg/time;
	}
	else if ( runType[i] == "B10" ) { 
	  B10[0][bin]=fg, B10err[0][bin]=TMath::Sqrt( fg/time );
	  B10_shift[0][bin]=fg, B10err_shift[0][bin]+= fg/time;
	}
	bin++;
      }
      rateFile.close();
      
      // West Type 0
      rateFile.open(TString::Format("%s/Asymmetry/BinByBinComparison/SIM_Wrun%i_anaChD.dat", getenv("ANALYSIS_CODE"),runNumber[i]));
      rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
      rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
      
      bin = 0; 
	
      while ( rateFile >> en >> fg >> bg )  { 
      
	if ( runType[i] == "A2" ) { 
	  A2[1][bin]=fg, A2err[1][bin]=TMath::Sqrt( fg/time );
	  A2_shift[1][bin]=fg, A2err_shift[1][bin]+= fg/time;
	}
	else if ( runType[i] == "A5" ) { 
	  A5[1][bin]=fg, A5err[1][bin]=TMath::Sqrt( fg/time );
	  A5_shift[1][bin]=fg, A5err_shift[1][bin]+= fg/time;
	}
	else if ( runType[i] == "A7" ) { 
	  A7[1][bin]=fg, A7err[1][bin]=TMath::Sqrt( fg/time );
	  A7_shift[1][bin]=fg, A7err_shift[1][bin]+= fg/time;
	}
	else if ( runType[i] == "A10" ) { 
	  A10[1][bin]=fg, A10err[1][bin]=TMath::Sqrt( fg/time );
	  A10_shift[1][bin]=fg, A10err_shift[1][bin]+= fg/time;
	}
	else if ( runType[i] == "B2" ) { 
	  B2[1][bin]=fg, B2err[1][bin]=TMath::Sqrt( fg/time );
	  B2_shift[1][bin]=fg, B2err_shift[1][bin]+= fg/time;
	}
	else if ( runType[i] == "B5" ) { 
	  B5[1][bin]=fg, B5err[1][bin]=TMath::Sqrt( fg/time );
	  B5_shift[1][bin]=fg, B5err_shift[1][bin]+= fg/time;
	}
	else if ( runType[i] == "B7" ) { 
	  B7[1][bin]=fg, B7err[1][bin]=TMath::Sqrt( fg/time );
	  B7_shift[1][bin]=fg, B7err_shift[1][bin]+= fg/time;
	}
	else if ( runType[i] == "B10" ) { 
	  B10[1][bin]=fg, B10err[1][bin]=TMath::Sqrt( fg/time );
	  B10_shift[1][bin]=fg, B10err_shift[1][bin]+= fg/time;
	}
	bin++;
      }

      rateFile.close();


      // Now check to see if we need to do shuffling of Type 1

      if ( percentT1!=0. ) {
	
	rateFile.open(TString::Format("%s/Asymmetry/BinByBinComparison/SIM_Erun%i_anaChF.dat", getenv("ANALYSIS_CODE"),runNumber[i]));
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
	bin = 0; 
	
	while ( rateFile >> en >> fg >> bg )  { 
	  
	  if ( runType[i] == "A2" ) { 
	    A2_shift[0][bin]+=percentT1*fg, A2err_shift[0][bin]+=percentT1*fg/time;	}
	  else if ( runType[i] == "A5" ) { 
	    A5_shift[0][bin]+=percentT1*fg, A5err_shift[0][bin]+=percentT1*fg/time;	}
	  else if ( runType[i] == "A7" ) { 
	    A7_shift[0][bin]+=percentT1*fg, A7err_shift[0][bin]+=percentT1*fg/time;	}
	  else if ( runType[i] == "A10" ) { 
	    A10_shift[0][bin]+=percentT1*fg, A10err_shift[0][bin]+=percentT1*fg/time;
	  }
	  else if ( runType[i] == "B2" ) { 
	    B2_shift[0][bin]+=percentT1*fg, B2err_shift[0][bin]+=percentT1*percentT1*fg/time;	}
	  else if ( runType[i] == "B5" ) { 
	    B5_shift[0][bin]+=percentT1*fg, B5err_shift[0][bin]+=percentT1*fg/time;	}
	  else if ( runType[i] == "B7" ) { 
	    B7_shift[0][bin]+=percentT1*fg, B7err_shift[0][bin]+=percentT1*fg/time;	}
	  else if ( runType[i] == "B10" ) { 
	    B10_shift[0][bin]+=percentT1*fg, B10err_shift[0][bin]+=percentT1*fg/time;
	  }
	  bin++;
	}
	rateFile.close();


	rateFile.open(TString::Format("%s/Asymmetry/BinByBinComparison/SIM_Wrun%i_anaChF.dat", getenv("ANALYSIS_CODE"),runNumber[i]));
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
	bin = 0; 
	
	while ( rateFile >> en >> fg >> bg )  { 
	  
	  if ( runType[i] == "A2" ) { 
	    A2_shift[1][bin]+=percentT1*fg, A2err_shift[1][bin]+=percentT1*fg/time;	}
	  else if ( runType[i] == "A5" ) { 
	    A5_shift[1][bin]+=percentT1*fg, A5err_shift[1][bin]+=percentT1*fg/time;	}
	  else if ( runType[i] == "A7" ) { 
	    A7_shift[1][bin]+=percentT1*fg, A7err_shift[1][bin]+=percentT1*fg/time;	}
	  else if ( runType[i] == "A10" ) { 
	    A10_shift[1][bin]+=percentT1*fg, A10err_shift[1][bin]+=percentT1*fg/time;
	  }
	  else if ( runType[i] == "B2" ) { 
	    B2_shift[1][bin]+=percentT1*fg, B2err_shift[1][bin]+=percentT1*fg/time;	}
	  else if ( runType[i] == "B5" ) { 
	    B5_shift[1][bin]+=percentT1*fg, B5err_shift[1][bin]+=percentT1*fg/time;	}
	  else if ( runType[i] == "B7" ) { 
	    B7_shift[1][bin]+=percentT1*fg, B7err_shift[1][bin]+=percentT1*fg/time;	}
	  else if ( runType[i] == "B10" ) { 
	    B10_shift[1][bin]+=percentT1*fg, B10err_shift[1][bin]+=percentT1*fg/time;
	  }
	  bin++;
	}
	
	rateFile.close();

      }


      // Now check to see if we need to do shuffling of Type 1

      if ( percentT2!=0. ) {
	// Note that for T2, we shuffle the west events to the east, as this would mimic a wirechamber not triggering on
	// the side that had a scintillator not trigger. Otherwise, the event would not be classified as an electron like event
	rateFile.open(TString::Format("%s/Asymmetry/BinByBinComparison/SIM_Wrun%i_anaChJ.dat", getenv("ANALYSIS_CODE"),runNumber[i]));
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
	bin = 0; 
	
	while ( rateFile >> en >> fg >> bg )  { 
	  
	  if ( runType[i] == "A2" ) { 
	   
	    A2_shift[0][bin]+=percentT2*fg, A2err_shift[0][bin]+=percentT2*fg/time;	}
	  else if ( runType[i] == "A5" ) { 
	    A5_shift[0][bin]+=percentT2*fg, A5err_shift[0][bin]+=percentT2*fg/time;	}
	  else if ( runType[i] == "A7" ) { 
	    A7_shift[0][bin]+=percentT2*fg, A7err_shift[0][bin]+=percentT2*fg/time;	}
	  else if ( runType[i] == "A10" ) { 
	    A10_shift[0][bin]+=percentT2*fg, A10err_shift[0][bin]+=percentT2*fg/time;
	  }
	  else if ( runType[i] == "B2" ) { 
	    B2_shift[0][bin]+=percentT2*fg, B2err_shift[0][bin]+=percentT2*fg/time;	}
	  else if ( runType[i] == "B5" ) { 
	    B5_shift[0][bin]+=percentT2*fg, B5err_shift[0][bin]+=percentT2*fg/time;	}
	  else if ( runType[i] == "B7" ) { 
	    B7_shift[0][bin]+=percentT2*fg, B7err_shift[0][bin]+=percentT2*fg/time;	}
	  else if ( runType[i] == "B10" ) { 
	    B10_shift[0][bin]+=percentT2*fg, B10err_shift[0][bin]+=percentT2*fg/time;
	  }
	  bin++;
	}
	rateFile.close();

	rateFile.open(TString::Format("%s/Asymmetry/BinByBinComparison/SIM_Erun%i_anaChJ.dat", getenv("ANALYSIS_CODE"),runNumber[i]));
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
	bin = 0; 
	
	while ( rateFile >> en >> fg >> bg )  { 
	  
	  if ( runType[i] == "A2" ) { 
	    A2_shift[1][bin]+=percentT2*fg, A2err_shift[1][bin]+=percentT2*fg/time;	}
	  else if ( runType[i] == "A5" ) { 
	    A5_shift[1][bin]+=percentT2*fg, A5err_shift[1][bin]+=percentT2*fg/time;	}
	  else if ( runType[i] == "A7" ) { 
	    A7_shift[1][bin]+=percentT2*fg, A7err_shift[1][bin]+=percentT2*fg/time;	}
	  else if ( runType[i] == "A10" ) { 
	    A10_shift[1][bin]+=percentT2*fg, A10err_shift[1][bin]+=percentT2*fg/time;
	  }
	  else if ( runType[i] == "B2" ) { 
	    B2_shift[1][bin]+=percentT2*fg, B2err_shift[1][bin]+=percentT2*fg/time;	}
	  else if ( runType[i] == "B5" ) { 
	    B5_shift[1][bin]+=percentT2*fg, B5err_shift[1][bin]+=percentT2*fg/time;	}
	  else if ( runType[i] == "B7" ) { 
	    B7_shift[1][bin]+=percentT2*fg, B7err_shift[1][bin]+=percentT2*fg/time;	}
	  else if ( runType[i] == "B10" ) { 
	    B10_shift[1][bin]+=percentT2*fg, B10err_shift[1][bin]+=percentT2*fg/time;
	  }
	  bin++;
	}
	rateFile.close();

      }


      // Now check to see if we need to do shuffling of Type 1

      if ( percentT3!=0. ) {

	rateFile.open(TString::Format("%s/Asymmetry/BinByBinComparison/SIM_Erun%i_anaChK.dat", getenv("ANALYSIS_CODE"),runNumber[i]));
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
	bin = 0; 
	
	while ( rateFile >> en >> fg >> bg )  { 
	  
	  if ( runType[i] == "A2" ) { 
	    A2_shift[0][bin]+=percentT3*fg, A2err_shift[0][bin]+=percentT3*fg/time;	}
	  else if ( runType[i] == "A5" ) { 
	    A5_shift[0][bin]+=percentT3*fg, A5err_shift[0][bin]+=percentT3*fg/time;	}
	  else if ( runType[i] == "A7" ) { 
	    A7_shift[0][bin]+=percentT3*fg, A7err_shift[0][bin]+=percentT3*fg/time;	}
	  else if ( runType[i] == "A10" ) { 
	    A10_shift[0][bin]+=percentT3*fg, A10err_shift[0][bin]+=percentT3*fg/time;
	  }
	  else if ( runType[i] == "B2" ) { 
	    B2_shift[0][bin]+=percentT3*fg, B2err_shift[0][bin]+=percentT3*fg/time;	}
	  else if ( runType[i] == "B5" ) { 
	    B5_shift[0][bin]+=percentT3*fg, B5err_shift[0][bin]+=percentT3*fg/time;	}
	  else if ( runType[i] == "B7" ) { 
	    B7_shift[0][bin]+=percentT3*fg, B7err_shift[0][bin]+=percentT3*fg/time;	}
	  else if ( runType[i] == "B10" ) { 
	    B10_shift[0][bin]+=percentT3*fg, B10err_shift[0][bin]+=percentT3*fg/time;
	  }
	  bin++;
	}
	rateFile.close();

	rateFile.open(TString::Format("%s/Asymmetry/BinByBinComparison/SIM_Wrun%i_anaChK.dat", getenv("ANALYSIS_CODE"),runNumber[i]));
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
	bin = 0; 
	
	while ( rateFile >> en >> fg >> bg )  { 
	  
	  if ( runType[i] == "A2" ) { 
	    A2_shift[1][bin]+=percentT3*fg, A2err_shift[1][bin]+=percentT3*fg/time;	}
	  else if ( runType[i] == "A5" ) { 
	    A5_shift[1][bin]+=percentT3*fg, A5err_shift[1][bin]+=percentT3*fg/time;	}
	  else if ( runType[i] == "A7" ) { 
	    A7_shift[1][bin]+=percentT3*fg, A7err_shift[1][bin]+=percentT3*fg/time;	}
	  else if ( runType[i] == "A10" ) { 
	    A10_shift[1][bin]+=percentT3*fg, A10err_shift[1][bin]+=percentT3*fg/time;
	  }
	  else if ( runType[i] == "B2" ) { 
	    B2_shift[1][bin]+=percentT3*fg, B2err_shift[1][bin]+=percentT3*fg/time;	}
	  else if ( runType[i] == "B5" ) { 
	    B5_shift[1][bin]+=percentT3*fg, B5err_shift[1][bin]+=percentT3*fg/time;	}
	  else if ( runType[i] == "B7" ) { 
	    B7_shift[1][bin]+=percentT3*fg, B7err_shift[1][bin]+=percentT3*fg/time;	}
	  else if ( runType[i] == "B10" ) { 
	    B10_shift[1][bin]+=percentT3*fg, B10err_shift[1][bin]+=percentT3*fg/time;
	  }
	  bin++;
	}
	rateFile.close();

      }
    
      ////////////////////////////////////////////////////////////////////////////////////////////////
      // Shuffling lost events.....
      if ( percentT4!=0. ) {

	TFile *rfile = new TFile(TString::Format("%s/systematics/MC_Corrections/lost_events_%s.root",
						 getenv("ANALYSIS_CODE"),rn<20000?"2011-2012":"2012-2013"),"READ");

	TH1D *h = (TH1D*)rfile->Get("percLost");

	std::vector <Double_t> perc(120,0.);

	for (int j=0; j<h->GetNbinsX(); ++j) {
	  perc[j] = h->GetBinContent(j+1);
	  std::cout << perc[j] << std::endl;
	}

	for (UInt_t bin=0; bin<perc.size(); ++bin) {

	  if ( runType[i] == "A2" ) { 
	    A2_shift[0][bin]+=percentT4*perc[bin]*A2[0][bin], A2err_shift[0][bin]+=percentT4*perc[bin]*A2[0][bin]/time;	}
	  else if ( runType[i] == "A5" ) { 
	    A5_shift[0][bin]+=percentT4*perc[bin]*A5[0][bin], A5err_shift[0][bin]+=percentT4*perc[bin]*A5[0][bin]/time;	}
	  else if ( runType[i] == "A7" ) { 
	    A7_shift[0][bin]+=percentT4*perc[bin]*A7[0][bin], A7err_shift[0][bin]+=percentT4*perc[bin]*A7[0][bin]/time;	}
	  else if ( runType[i] == "A10" ) { 
	    A10_shift[0][bin]+=percentT4*perc[bin]*A10[0][bin], A10err_shift[0][bin]+=percentT4*perc[bin]*A10[0][bin]/time; }
	  else if ( runType[i] == "B2" ) { 
	    B2_shift[0][bin]+=percentT4*perc[bin]*B2[0][bin], B2err_shift[0][bin]+=percentT4*perc[bin]*B2[0][bin]/time;	}
	  else if ( runType[i] == "B5" ) { 
	    B5_shift[0][bin]+=percentT4*perc[bin]*B5[0][bin], B5err_shift[0][bin]+=percentT4*perc[bin]*B5[0][bin]/time;	}
	  else if ( runType[i] == "B7" ) { 
	    B7_shift[0][bin]+=percentT4*perc[bin]*B7[0][bin], B7err_shift[0][bin]+=percentT4*perc[bin]*B7[0][bin]/time;	}
	  else if ( runType[i] == "B10" ) { 
	    B10_shift[0][bin]+=percentT4*perc[bin]*B10[0][bin], B10err_shift[0][bin]+=percentT4*perc[bin]*B10[0][bin]/time; }

	  if ( runType[i] == "A2" ) { 
	    A2_shift[1][bin]+=percentT4*perc[bin]*A2[1][bin], A2err_shift[1][bin]+=percentT4*perc[bin]*A2[1][bin]/time;	}
	  else if ( runType[i] == "A5" ) { 
	    A5_shift[1][bin]+=percentT4*perc[bin]*A5[1][bin], A5err_shift[1][bin]+=percentT4*perc[bin]*A5[1][bin]/time;	}
	  else if ( runType[i] == "A7" ) { 
	    A7_shift[1][bin]+=percentT4*perc[bin]*A7[1][bin], A7err_shift[1][bin]+=percentT4*perc[bin]*A7[1][bin]/time;	}
	  else if ( runType[i] == "A10" ) { 
	    A10_shift[1][bin]+=percentT4*perc[bin]*A10[1][bin], A10err_shift[1][bin]+=percentT4*perc[bin]*A10[1][bin]/time; }
	  else if ( runType[i] == "B2" ) { 
	    B2_shift[1][bin]+=percentT4*perc[bin]*B2[1][bin], B2err_shift[1][bin]+=percentT4*perc[bin]*B2[1][bin]/time;	}
	  else if ( runType[i] == "B5" ) { 
	    B5_shift[1][bin]+=percentT4*perc[bin]*B5[1][bin], B5err_shift[1][bin]+=percentT4*perc[bin]*B5[1][bin]/time;	}
	  else if ( runType[i] == "B7" ) { 
	    B7_shift[1][bin]+=percentT4*perc[bin]*B7[1][bin], B7err_shift[1][bin]+=percentT4*perc[bin]*B7[1][bin]/time;	}
	  else if ( runType[i] == "B10" ) { 
	    B10_shift[1][bin]+=percentT4*perc[bin]*B10[1][bin], B10err_shift[1][bin]+=percentT4*perc[bin]*B10[1][bin]/time; }

	
	}
	delete rfile;
      }	
  
  //////////////////////////////////////////////////////////////////////////

    }
    
    // take the sq root of the shifted errors since they would be added in quadrature
    for (int j=0; j<2; ++j) {
      for (int k=0; k<numBins; ++k) {
	A2err_shift[j][k] = TMath::Sqrt( A2err_shift[j][k] );
	A5err_shift[j][k] = TMath::Sqrt( A5err_shift[j][k] );
	A7err_shift[j][k] = TMath::Sqrt( A7err_shift[j][k] );
	A10err_shift[j][k] = TMath::Sqrt( A10err_shift[j][k] );
	
	B2err_shift[j][k] = TMath::Sqrt( B2err_shift[j][k] );
	B5err_shift[j][k] = TMath::Sqrt( B5err_shift[j][k] );
	B7err_shift[j][k] = TMath::Sqrt( B7err_shift[j][k] );
	B10err_shift[j][k] = TMath::Sqrt( B10err_shift[j][k] );
      }
    }
    // Now calculate the asymmetries
    std::vector <double> asymmetry(numBins,0);
    std::vector <double> asymmetryError(numBins,0);
    
    std::vector <double> asymmetry_shift(numBins,0);
    std::vector <double> asymmetryError_shift(numBins,0);
    
    calcBinByBinSuperRatio(A2,A5,A7,A10,B2,B5,B7,B10,
			   A2err,A5err,A7err,A10err,B2err,B5err,B7err,B10err,
			   asymmetry,asymmetryError);
    
    calcBinByBinSuperRatio(A2_shift,A5_shift,A7_shift,A10_shift,B2_shift,B5_shift,B7_shift,B10_shift,
			   A2err_shift,A5err_shift,A7err_shift,A10err_shift,B2err_shift,B5err_shift,B7err_shift,B10err_shift,
			   asymmetry_shift,asymmetryError_shift);
    
   
    // TODO: Add in the application of the corrections
    std::vector < std::vector <Double_t> > deltaSys = LoadOctetSystematics(octetNum,"D",enBins);
    std::vector <Double_t> theoryCorr(numBins,1.);
    theoryCorr = LoadTheoryCorrections(enBins);


    for ( int j=0; j<numBins; ++j ) {
      asymmetry[j] = asymmetry[j]*deltaSys[j][0]/theoryCorr[j];
      asymmetryError[j] = asymmetryError[j]*deltaSys[j][0]/theoryCorr[j];        
    // Since Asym is multiplied by delta, AsymError would be as well...
      asymmetry_shift[j] = asymmetry_shift[j]*deltaSys[j][0]/theoryCorr[j];
      asymmetryError_shift[j] = asymmetryError_shift[j]*deltaSys[j][0]/theoryCorr[j];        

    }

    // Add the bin by bin asymmetries to the total asymmetry in each bin, weighted by its own asymmetry
    for ( int j=0; j<numBins; ++j ) {
      Asymm[0][j]+= ( asymmetryError[j]!=0. ? ( asymmetry[j] / (asymmetryError[j]*asymmetryError[j]) ) : 0. );
      Asymm[1][j]+= ( asymmetryError[j]!=0. ? 1./(asymmetryError[j]*asymmetryError[j]) : 0. );
      AsymmWithShuffle[0][j]+= ( asymmetryError_shift[j]!=0. ? 
				 ( asymmetry_shift[j] / (asymmetryError_shift[j]*asymmetryError_shift[j]) ) : 0. );
      AsymmWithShuffle[1][j]+= ( asymmetryError_shift[j]!=0. ? 1./(asymmetryError_shift[j]*asymmetryError_shift[j]): 0. );
    }
    
    
  }

  // Finish calculation of asymm
  for ( int j=0; j<numBins; ++j ) {
    Asymm[0][j] /= ( Asymm[1][j]!=0. ? Asymm[1][j] : 1. );
    Asymm[1][j] = Asymm[1][j]!=0. ? 1. / TMath::Sqrt( Asymm[1][j] ) : 0.;
    AsymmWithShuffle[0][j] /= ( AsymmWithShuffle[1][j]!=0. ? AsymmWithShuffle[1][j] : 1. );
    AsymmWithShuffle[1][j] = AsymmWithShuffle[1][j]!=0. ? 1./TMath::Sqrt( AsymmWithShuffle[1][j] ) : 0.;
      
  }

  //Divide by Beta/2.
  for ( int j=0; j<numBins; ++j ) {
    Double_t betaOver2 = returnBeta(enBins[j])/2.;
    Asymm[0][j] /= betaOver2;
    Asymm[1][j] /= betaOver2;
    AsymmWithShuffle[0][j] /= betaOver2;
    AsymmWithShuffle[1][j] /= betaOver2;
      
  }

  //Now we need to fit for the asymmetries
  TF1 *fit = new TF1("fit","[0]",emin,emax);
  fit->SetParameter(0,-0.12);
  
  TGraphErrors *a = new TGraphErrors(numBins,&enBins[0],&Asymm[0][0],0,&Asymm[1][0]);
  a->Fit(fit,"R");
  
  double asym = fit->GetParameter(0);
  double asymErr = fit->GetParError(0);
  
  TGraphErrors *aWithShuffle = new TGraphErrors(numBins,&enBins[0],&AsymmWithShuffle[0][0],0,&AsymmWithShuffle[1][0]);
  aWithShuffle->Fit(fit,"R");
  
  double asymWithShuffle = fit->GetParameter(0);
  double asymErrWithShuffle = fit->GetParError(0);

  std::cout << "percentT1 = " << percentT1 << std::endl
	    << "percentT2 = " << percentT2 << std::endl
	    << "percentT3 = " << percentT3 << std::endl
	    << "asym = " << asym << " +/- " << asymErr << std::endl
	    << "asymWithShuffle = " << asymWithShuffle << " +/- " << asymErrWithShuffle << std::endl
	    << "Percent Difference deltaA/A = " << (asym-asymWithShuffle)/asym << std::endl;

  std::ofstream resultFile("resultHolderEventShuffle.txt");
  resultFile << "asym=\t" << asym << "\n" << "asymWithShuffle=\t" << asymWithShuffle << "\n"
	   << "percentDifference=\t" << (asym-asymWithShuffle)/asym;
  resultFile.close();

  return 0;
  
}
      

std::vector < std::vector <Double_t> > LoadOctetSystematics(Int_t octet, std::string anaChoice, std::vector <Double_t> enBinMidpoint) {

  Int_t iAnaChoice;

  for (UInt_t i = 0; i<10; i++) {
    
    if (anaChoices[i]==anaChoice) iAnaChoice = i+1;

  }
  
  //TString filename = TString::Format("%s/Octet_%i/OctetAsymmetry/Systematics/ThOverProc_Octet-%i_Analysis-%i.txt",getenv("ANALYSIS_RESULTS"),octet,octet,iAnaChoice);
  TString filename = TString::Format("%s/systematics/MC_Corrections/DeltaExp_OctetByOctetCorrections/ThOverProc_Octet-%i_Analysis-%i.txt",getenv("ANALYSIS_CODE"),octet,iAnaChoice);
  //std::cout << filename.Data() << std::endl;
  std::vector < std::vector <Double_t> > syst(enBinMidpoint.size(), std::vector<Double_t>(2,1.));
  
  std::ifstream infile(filename.Data());

  if (!infile.is_open()) throw "Couldn't open file in LoadOctetSystematics!";

  //Read in the header crap
  //  std::string firstline[8];
  //for (int i=0; i<8; i++) {
  //  infile >> firstline[i];
  //  std::cout << firstline[i] << " " ;
  // }
  std::cout << "\n";

  //Read in the systematics
  Double_t mid; //midpoint of bin
  Double_t ratio;  //systematics correction (Apure/Aproc)
  Double_t midErr; //BinWidth/2
  Double_t ratioErr; //Systematic error on ratio

  Int_t it = 0;

  while (infile >> mid >> ratio >> ratioErr) {
    if (mid==enBinMidpoint[it]) {
      syst[it][0] = ratio!=0. ? ratio : 1.;
      syst[it][1] = ratioErr;
      std::cout << mid << " " << ratio << " " << ratioErr << "\n";
      it++;
    }
  }

  std::cout << "Loaded Systematic Errors for Octet " << octet << " ...\n";

  return syst;
};

std::vector <Double_t> LoadTheoryCorrections(std::vector <Double_t> enBinMidpoint) {

  std::vector <Double_t> syst(enBinMidpoint.size(), 1.);

  for (UInt_t i=0; i<syst.size(); i++) {
   
    syst[i] = asymmetryCorrectionFactor(enBinMidpoint[i]); //As defined in BetaSpectrum.hh by MPM

  }

  return syst;

};
