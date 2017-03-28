/*
 Compares the BG subtracted data spectra to simulation data which has 
 detector effects already taken into account
*/

#include "Asymmetries.hh"
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
  if (argc!=5) {
    std::cout << "Usage: ./SuperSumComparison.exe [octet start] [octet end] [energy bin width] [BOOL Sim Asymmetry Weight On]\n";
    std::cout << "The code will produce comparisons for every octet in the range given,\non an octet-by-octet basis, and as a whole, using the Super-Sum\n";
    exit(0);
  }

  gStyle->SetOptStat(0);

  int octetNumStart = atoi(argv[1]);
  int octetNumEnd = atoi(argv[2]);
  //double enWinLow = (double) atof(argv[3]);
  //double enWinHigh = (double) atof(argv[4]);
  int enBinWidth = atoi(argv[3]);
  
  bool asymOn = false;
  
  if (std::string(argv[4])=="true" || std::string(argv[4])=="1") asymOn=true;
  
  int numBins = 1200/enBinWidth;

  //set up energy binning vector
  std::vector <double> enBins(numBins,0.);
  
  for (int bin=0; bin<numBins; bin++) {
    enBins[bin] = ((double)bin+.5)*(double)enBinWidth;
  }
  
  std::vector <std::vector <double> > superSumTotal_uk(4, std::vector <double>(numBins,0.)); //type0=0, type1=1, type23=2, ALL=3
  std::vector <std::vector <double> > superSumTotalError_uk(4, std::vector <double>(numBins,0.));
  std::vector <std::vector <double> > superSumTotal_sim(4, std::vector <double>(numBins,0.));
  std::vector <std::vector <double> > superSumTotalError_sim(4, std::vector <double>(numBins,0.));


  for (int octetNum=octetNumStart; octetNum<octetNumEnd+1; octetNum++) {
    
    if (std::find(badOct.begin(), badOct.end(),octetNum) != badOct.end()) { continue; } //Checking if octet should be ignored for data quality reasons

    std::vector < std::vector <double> > superSum_uk(4, std::vector <double>(numBins,0.));//,std::vector <double>(numBins,0.));
    std::vector < std::vector <double> > superSumError_uk(4, std::vector <double>(numBins,0.));//,std::vector <double>(numBins,0.));
    std::vector < std::vector <double> > superSum_sim(4, std::vector <double>(numBins,0.));//,std::vector <double>(numBins,0.));
    std::vector < std::vector <double> > superSumError_sim(4, std::vector <double>(numBins,0.));//,std::vector <double>(numBins,0.));
   
 
    try {
      
      //double normFactor = 0., ukIntegral=0., simIntegral=0.;

      ifstream infile;
      double binLow, rate, rateErr;
      int inc;

      { //Setting scope for first OctetAsymmetry so that it will be deleted to clear memory
	
	infile.open(TString::Format("%s/Octet_%i/OctetAsymmetry/superSum_Octet%i_AnaChA.dat",getenv("ANALYSIS_RESULTS"),octetNum,octetNum).Data());
	
	inc = 0;

	while ( infile >> binLow >> rate >> rateErr ) {
	  superSum_uk[3][inc] = rate;
	  superSumError_uk[3][inc] = rate!=0. ? rateErr : 0.;
	  inc++;
	}
	infile.close();

	for (int n=0; n<numBins;n++) {
	  superSumTotal_uk[3][n]+=superSumError_uk[3][n]>0.?(1./power(superSumError_uk[3][n],2))*superSum_uk[3][n]:0.;
	  superSumTotalError_uk[3][n]+=superSumError_uk[3][n]>0.?(1./power(superSumError_uk[3][n],2)):0.;
	  //if (enBins[n]>=enWinLow && enBins[n]<=enWinHigh) ukIntegral+=superSum_uk[3][n];
	}
	

	//Type 0 events
	infile.open(TString::Format("%s/Octet_%i/OctetAsymmetry/superSum_Octet%i_AnaChD.dat",getenv("ANALYSIS_RESULTS"),octetNum,octetNum).Data());
	
	inc = 0;

	while ( infile >> binLow >> rate >> rateErr ) {
	  superSum_uk[0][inc] = rate;
	  superSumError_uk[0][inc] = rate!=0. ? rateErr : 0.;
	  inc++;
	}
	infile.close();

	for (int n=0; n<numBins;n++) {
	  //std::cout << enBins[n] << " " << superSum[n] << " " << superSumError[n] << std::endl; 
	  superSumTotal_uk[0][n]+=superSumError_uk[0][n]>0.?(1./power(superSumError_uk[0][n],2))*superSum_uk[0][n]:0.;
	  superSumTotalError_uk[0][n]+=superSumError_uk[0][n]>0.?(1./power(superSumError_uk[0][n],2)):0.;
	  //if (enBins[n]>=enWinLow && enBins[n]<=enWinHigh) ukIntegral+=superSum_uk[0][n];
	}
       	


	//Type 1 events
	infile.open(TString::Format("%s/Octet_%i/OctetAsymmetry/superSum_Octet%i_AnaChF.dat",getenv("ANALYSIS_RESULTS"),octetNum,octetNum).Data());
	
	inc = 0;

	while ( infile >> binLow >> rate >> rateErr ) {
	  superSum_uk[1][inc] = rate;
	  superSumError_uk[1][inc] = rate!=0. ? rateErr : 0.;
	  inc++;
	}
	infile.close();

	for (int n=0; n<numBins;n++) {
	  superSumTotal_uk[1][n]+=superSumError_uk[1][n]>0.?(1./power(superSumError_uk[1][n],2))*superSum_uk[1][n]:0.;
	  superSumTotalError_uk[1][n]+=superSumError_uk[1][n]>0.?(1./power(superSumError_uk[1][n],2)):0.;
	}
	
		

	//Type 2/3 events
	infile.open(TString::Format("%s/Octet_%i/OctetAsymmetry/superSum_Octet%i_AnaChG.dat",getenv("ANALYSIS_RESULTS"),octetNum,octetNum).Data());
	
	inc = 0;

	while ( infile >> binLow >> rate >> rateErr ) {
	  superSum_uk[2][inc] = rate;
	  superSumError_uk[2][inc] = rate!=0. ? rateErr : 0.;
	  inc++;
	}
	infile.close();

	for (int n=0; n<numBins;n++) {
	  superSumTotal_uk[2][n]+=superSumError_uk[2][n]>0.?(1./power(superSumError_uk[2][n],2))*superSum_uk[2][n]:0.;
	  superSumTotalError_uk[2][n]+=superSumError_uk[2][n]>0.?(1./power(superSumError_uk[2][n],2)):0.;
	}
	
	

	
	
      }

      // SIMULATION

      { //Setting scope for second OctetAsymmetry so that it will be deleted to clear memory

	infile.open(TString::Format("%s/Octet_%i/OctetAsymmetry/superSum_Octet%i_AnaChA.dat",getenv("SIM_ANALYSIS_RESULTS"),octetNum,octetNum).Data());
	
	inc = 0;

	while ( infile >> binLow >> rate >> rateErr ) {
	  superSum_sim[3][inc] = rate;
	  superSumError_sim[3][inc] = rate!=0. ? rateErr : 0.;
	  inc++;
	}
	infile.close();


	//for (int n=0; n<numBins; n++) {
	// if (enBins[n]>=enWinLow && enBins[n]<=enWinHigh) simIntegral+=superSum_sim[3][n];
	//}
	//normFactor = ukIntegral/simIntegral;

	for (int n=0; n<numBins;n++) {
	  //superSum_sim[3][n] = normFactor*superSum_sim[3][n];
	  //superSumError_sim[3][n] = normFactor*superSumError_sim[3][n];
	  superSumTotal_sim[3][n]+=superSumError_sim[3][n]>0.?(1./power(superSumError_sim[3][n],2))*superSum_sim[3][n]:0.;
	  superSumTotalError_sim[3][n]+=superSumError_sim[3][n]>0.?(1./power(superSumError_sim[3][n],2)):0.;
	}
	

	//Type 0 events
	infile.open(TString::Format("%s/Octet_%i/OctetAsymmetry/superSum_Octet%i_AnaChD.dat",getenv("SIM_ANALYSIS_RESULTS"),octetNum,octetNum).Data());
	
	inc = 0;

	while ( infile >> binLow >> rate >> rateErr ) {
	  superSum_sim[0][inc] = rate;
	  superSumError_sim[0][inc] = rate!=0. ? rateErr : 0.;
	  inc++;
	}
	infile.close();

	/*for (int n=0; n<numBins; n++) {
	  if (enBins[n]>=enWinLow && enBins[n]<=enWinHigh) simIntegral+=superSum_sim[0][n];
	}
	normFactor = ukIntegral/simIntegral;*/

	for (int n=0; n<numBins;n++) {
	  //std::cout << enBins[n] << " " << superSum[n] << " " << superSumError[n] << std::endl; 
	  //superSum_sim[0][n] = normFactor*superSum_sim[0][n];
	  //superSumError_sim[0][n] = normFactor*superSumError_sim[0][n];
	  superSumTotal_sim[0][n]+=superSumError_sim[0][n]>0.?(1./power(superSumError_sim[0][n],2))*superSum_sim[0][n]:0.;
	  superSumTotalError_sim[0][n]+=superSumError_sim[0][n]>0.?(1./power(superSumError_sim[0][n],2)):0.;
	}
       	
	

	//Type 1 events
	infile.open(TString::Format("%s/Octet_%i/OctetAsymmetry/superSum_Octet%i_AnaChF.dat",getenv("SIM_ANALYSIS_RESULTS"),octetNum,octetNum).Data());
	
	inc = 0;

	while ( infile >> binLow >> rate >> rateErr ) {
	  superSum_sim[1][inc] = rate;
	  superSumError_sim[1][inc] = rate!=0. ? rateErr : 0.;
	  inc++;
	}
	infile.close();

	for (int n=0; n<numBins;n++) {
	  //superSum_sim[1][n] = normFactor*superSum_sim[1][n];
	  //superSumError_sim[1][n] = normFactor*superSumError_sim[1][n];
	  superSumTotal_sim[1][n]+=superSumError_sim[1][n]>0.?(1./power(superSumError_sim[1][n],2))*superSum_sim[1][n]:0.;
	  superSumTotalError_sim[1][n]+=superSumError_sim[1][n]>0.?(1./power(superSumError_sim[1][n],2)):0.;
	}
	
	

	//Type 2/3 events
	infile.open(TString::Format("%s/Octet_%i/OctetAsymmetry/superSum_Octet%i_AnaChG.dat",getenv("SIM_ANALYSIS_RESULTS"),octetNum,octetNum).Data());
	
	inc = 0;

	while ( infile >> binLow >> rate >> rateErr ) {
	  superSum_sim[2][inc] = rate;
	  superSumError_sim[2][inc] = rate!=0. ? rateErr : 0.;
	  inc++;
	}
	infile.close();

	for (int n=0; n<numBins;n++) {
	  //superSum_sim[2][n] = normFactor*superSum_sim[2][n];
	  //superSumError_sim[2][n] = normFactor*superSumError_sim[2][n];
	  superSumTotal_sim[2][n]+=superSumError_sim[2][n]>0.?(1./power(superSumError_sim[2][n],2))*superSum_sim[2][n]:0.;
	  superSumTotalError_sim[2][n]+=superSumError_sim[2][n]>0.?(1./power(superSumError_sim[2][n],2)):0.;
	}
	
	

	

	

      }
      
      std::string file1 = asymOn ? "superSumPlots/octet_"+itos(octetNum)+"_AsymOn.root" : "superSumPlots/octet_"+itos(octetNum)+"_AsymOff.root";
      TFile *f = new TFile(file1.c_str(),"RECREATE");
      
      TH1D Erecon0_uk("Erecon0_uk","Type 0 Super-Sum",numBins,0.,1200.);
      Erecon0_uk.SetLineColor(kBlue);
      //Erecon0_uk.SetMarkerStyle(21);
      Erecon0_uk.GetXaxis()->SetTitle("Energy (keV)");
      TH1D Erecon1_uk("Erecon1_uk","Type 1 Super-Sum",numBins,0.,1200.);
      Erecon1_uk.SetLineColor(kBlue);
      Erecon1_uk.GetXaxis()->SetTitle("Energy (keV)");
      TH1D Erecon23_uk("Erecon23_uk","Type 23 Super-Sum",numBins,0.,1200.);
      Erecon23_uk.SetLineColor(kBlue);
      Erecon23_uk.GetXaxis()->SetTitle("Energy (keV)");
      TH1D EreconALL_uk("EreconALL_uk","Type ALL Super-Sum",numBins,0.,1200.);
      EreconALL_uk.SetLineColor(kBlue);
      EreconALL_uk.GetXaxis()->SetTitle("Energy (keV)");

      TH1D Erecon0_sim("Erecon0_sim","Type 0 Super-Sum",numBins,0.,1200.);
      Erecon0_sim.SetLineColor(kBlue);
      //Erecon0_sim.SetMarkerStyle(21);
      Erecon0_sim.GetXaxis()->SetTitle("Energy (keV)");
      TH1D Erecon1_sim("Erecon1_sim","Type 1 Super-Sum",numBins,0.,1200.);
      Erecon1_sim.SetLineColor(kBlue);
      Erecon1_sim.GetXaxis()->SetTitle("Energy (keV)");
      TH1D Erecon23_sim("Erecon23_sim","Type 23 Super-Sum",numBins,0.,1200.);
      Erecon23_sim.SetLineColor(kBlue);
      Erecon23_sim.GetXaxis()->SetTitle("Energy (keV)");
      TH1D EreconALL_sim("EreconALL_sim","Type ALL Super-Sum",numBins,0.,1200.);
      EreconALL_sim.SetLineColor(kBlue);
      EreconALL_sim.GetXaxis()->SetTitle("Energy (keV)");
     
      for (int bin=0; bin<numBins; bin++) {
	Erecon0_uk.SetBinContent(bin+1,superSum_uk[0][bin]);
	Erecon1_uk.SetBinContent(bin+1,superSum_uk[1][bin]);
	Erecon23_uk.SetBinContent(bin+1,superSum_uk[2][bin]);
	EreconALL_uk.SetBinContent(bin+1,superSum_uk[3][bin]);

	Erecon0_uk.SetBinError(bin+1,superSumError_uk[0][bin]);
	Erecon1_uk.SetBinError(bin+1,superSumError_uk[1][bin]);
	Erecon23_uk.SetBinError(bin+1,superSumError_uk[2][bin]);
	EreconALL_uk.SetBinError(bin+1,superSumError_uk[3][bin]);

	Erecon0_sim.SetBinContent(bin+1,superSum_sim[0][bin]);
	Erecon1_sim.SetBinContent(bin+1,superSum_sim[1][bin]);
	Erecon23_sim.SetBinContent(bin+1,superSum_sim[2][bin]);
	EreconALL_sim.SetBinContent(bin+1,superSum_sim[3][bin]);

	Erecon0_sim.SetBinError(bin+1,superSumError_sim[0][bin]);
	Erecon1_sim.SetBinError(bin+1,superSumError_sim[1][bin]);
	Erecon23_sim.SetBinError(bin+1,superSumError_sim[2][bin]);
	EreconALL_sim.SetBinError(bin+1,superSumError_sim[3][bin]);
      }
      f->Write();
      f->Close();
    }

    catch(const char* ex){
      std::cerr << "Error: " << ex << std::endl;
    }
  }

  //Creating file which is summed over all octets in range
  if (octetNumStart!=octetNumEnd) {
    TString fileName = asymOn ? "superSumPlots/SuperSum_octets_"+itos(octetNumStart)+"-"+itos(octetNumEnd)+"_AsymOn.root" : "superSumPlots/SuperSum_octets_"+itos(octetNumStart)+"-"+itos(octetNumEnd)+"_AsymOff.root";
    //fileName+= octetNumStart;
    //fileName+= "-";
    //fileName+= octetNumEnd;
    //fileName+= "norm_"+ftos(enWinLow)+"-"+ftos(enWinHigh);
    //fileName+= ".root";
    TFile *f2 = new TFile(fileName,"RECREATE");
    
    for (int t=0; t<4; t++) {
      for (int bin=0; bin<numBins; bin++) {

	superSumTotal_uk[t][bin] = superSumTotalError_uk[t][bin]>0.? superSumTotal_uk[t][bin]/superSumTotalError_uk[t][bin] : 0.;
	superSumTotalError_uk[t][bin] = superSumTotalError_uk[t][bin]>0.? (1./TMath::Sqrt(superSumTotalError_uk[t][bin])) : 0.;
	superSumTotalError_uk[t][bin] = superSumTotalError_uk[t][bin] < (0.1*superSumTotal_uk[t][bin]) ? superSumTotalError_uk[t][bin] : (0.1*superSumTotal_uk[t][bin]);
	//std::cout << enBins[bin] << " " << superSumTotal_uk[0][bin] << " " << superSumTotalError_uk[0][bin] << std::endl;
	
	superSumTotal_sim[t][bin] = superSumTotalError_sim[t][bin]>0.? superSumTotal_sim[t][bin]/superSumTotalError_sim[t][bin] : 0.;
	superSumTotalError_sim[t][bin] = superSumTotalError_sim[t][bin]>0.? (1./TMath::Sqrt(superSumTotalError_sim[t][bin])) : 0.;
	superSumTotalError_sim[t][bin] = superSumTotalError_sim[t][bin] < (0.1*superSumTotal_sim[t][bin]) ? superSumTotalError_sim[t][bin] : (0.1*superSumTotal_sim[t][bin]);
      }
    }
    
    TH1D Erecon0_uk("Erecon0_uk","Type 0 Super-Sum",numBins,0.,1200.);
    Erecon0_uk.GetXaxis()->SetTitle("Energy (keV)");
    TH1D Erecon1_uk("Erecon1_uk","Type 1 Super-Sum",numBins,0.,1200.);
    Erecon1_uk.GetXaxis()->SetTitle("Energy (keV)");
    TH1D Erecon23_uk("Erecon23_uk","Type 23 Super-Sum",numBins,0.,1200.);
    Erecon23_uk.GetXaxis()->SetTitle("Energy (keV)");
    TH1D EreconALL_uk("EreconALL_uk","Type ALL Super-Sum",numBins,0.,1200.);
    EreconALL_uk.GetXaxis()->SetTitle("Energy (keV)");
    
    TH1D Erecon0_sim("Erecon0_sim","Type 0 Super-Sum",numBins,0.,1200.);
    Erecon0_sim.GetXaxis()->SetTitle("Energy (keV)");
    TH1D Erecon1_sim("Erecon1_sim","Type 1 Super-Sum",numBins,0.,1200.);
    Erecon1_sim.GetXaxis()->SetTitle("Energy (keV)");
    TH1D Erecon23_sim("Erecon23_sim","Type 23 Super-Sum",numBins,0.,1200.);
    Erecon23_sim.GetXaxis()->SetTitle("Energy (keV)");
    TH1D EreconALL_sim("EreconALL_sim","Type ALL Super-Sum",numBins,0.,1200.);
    EreconALL_sim.GetXaxis()->SetTitle("Energy (keV)");
    
    for (int bin=0; bin<numBins; bin++) {
      Erecon0_uk.SetBinContent(bin+1,superSumTotal_uk[0][bin]);
      Erecon1_uk.SetBinContent(bin+1,superSumTotal_uk[1][bin]);
      Erecon23_uk.SetBinContent(bin+1,superSumTotal_uk[2][bin]);
      EreconALL_uk.SetBinContent(bin+1,superSumTotal_uk[3][bin]);
      
      Erecon0_uk.SetBinError(bin+1,superSumTotalError_uk[0][bin]);
      Erecon1_uk.SetBinError(bin+1,superSumTotalError_uk[1][bin]);
      Erecon23_uk.SetBinError(bin+1,superSumTotalError_uk[2][bin]);
      EreconALL_uk.SetBinError(bin+1,superSumTotalError_uk[3][bin]);
      
      Erecon0_sim.SetBinContent(bin+1,superSumTotal_sim[0][bin]);
      Erecon1_sim.SetBinContent(bin+1,superSumTotal_sim[1][bin]);
      Erecon23_sim.SetBinContent(bin+1,superSumTotal_sim[2][bin]);
      EreconALL_sim.SetBinContent(bin+1,superSumTotal_sim[3][bin]);
      
      Erecon0_sim.SetBinError(bin+1,superSumTotalError_sim[0][bin]);
      Erecon1_sim.SetBinError(bin+1,superSumTotalError_sim[1][bin]);
      Erecon23_sim.SetBinError(bin+1,superSumTotalError_sim[2][bin]);
      EreconALL_sim.SetBinError(bin+1,superSumTotalError_sim[3][bin]);
    }
    
    f2->Write();
    f2->Close();
  }
  
  return 0;

}

    
    

