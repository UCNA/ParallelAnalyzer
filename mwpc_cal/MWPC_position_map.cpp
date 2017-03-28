//SWANKS EDIT
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <string>
#include <iomanip>

// ROOT libraries
#include <TH1D.h>

#include "DataTree.hh"
#include "MBUtils.hh"
#include "peaks.hh"
#include "positionMapHandler.hh"


using namespace std;


int main(int argc, char *argv[])
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  // Prompt for filename of run numbers
  int iXeRunPeriod;
  cout << "Enter Xenon run period: " << endl;
  cin  >> iXeRunPeriod;
  cout << endl;


  int nRuns = 0;
  int runList[500];

  char temp[500];
  sprintf(temp, "%s/run_lists/Xenon_Calibration_Run_Period_%i.dat", getenv("ANALYSIS_CODE"),iXeRunPeriod);
  ifstream fileRuns(temp);

  int ii = 0;
  while (fileRuns >> runList[ii]) {
    ii++;
    nRuns++;
  }

  cout << "... Number of runs: " << nRuns << endl;
  for (int i=0; i<nRuns; i++) {
    cout << runList[i] << endl;
  }

  
  double xyBinWidth = 5.; //2.5;
  MWPCPositionMap posmap(xyBinWidth,50.);
  Int_t nBinsXY = posmap.getNbinsXY();
  

  // Open output ntuple
  string tempOutBase;
  string tempOut;
  //sprintf(tempOut, "position_map_%s.root", argv[1]);
  tempOutBase = "MWPC_position_map_" + itos(iXeRunPeriod);
  
  tempOut =  getenv("MWPC_CALIBRATION")+std::string("/position_maps/")+tempOutBase+"_"+ftos(xyBinWidth)+"mm.root";
  TFile *fileOut = new TFile(tempOut.c_str(),"RECREATE");

  // Output histograms
  double EnergyBinWidth = 50.;
  int nEnergyBins = 14;
  double enBinStart = 100.;
  int nBinHist = 400;//1025;

  TH1D *hisxy_E[nEnergyBins][nBinsXY][nBinsXY];
  TH1D *hisxy_W[nEnergyBins][nBinsXY][nBinsXY];
  char *hisxyName = new char[10];

  for (int p=0; p<nEnergyBins; p++) {
    for (int i=0; i<posmap.getNbinsXY(); i++) {
      for (int j=0; j<posmap.getNbinsXY(); j++) {

	double elow = enBinStart + p*EnergyBinWidth;
	double ehigh = elow + EnergyBinWidth;
	
        hisxy_E[p][i][j] = new TH1D(TString::Format("East_%.0f-%.0f_x%i_y%i",elow,ehigh,i,j),
				    "", nBinHist,-0.5,3999.5);
	hisxy_W[p][i][j] = new TH1D(TString::Format("West_%.0f-%.0f_x%i_y%i",elow,ehigh,i,j),
				    "", nBinHist,-0.5,3999.5);

      }
    }
  }

  // Loop through input ntuples
  char tempIn[500];
  for (int i=0; i<nRuns; i++) {

    // Open input ntuple
    sprintf(tempIn, "%s/replay_pass3_%i.root",getenv("REPLAY_PASS3"), runList[i]);
    DataTree *t = new DataTree();
    t->setupInputTree(std::string(tempIn),"pass3");

    if ( !t->inputTreeIsGood() ) { 
      std::cout << "Skipping " << tempIn << "... Doesn't exist or couldn't be opened.\n";
      continue;
    }

    int nEvents = t->getEntries();
    cout << "Processing " << runList[i] << " ... " << endl;
    cout << "... nEvents = " << nEvents << endl;


    // Loop over events
    for (int i=0; i<nEvents; i++) {
      t->getEvent(i);

      if ( i%100000 == 0 ) std::cout << i <<  std::endl;
      
      // Select Type 0 electron events which are within our energy window
      if (t->PID != 1) continue;
      if (t->Type != 0) continue;
      if (t->Erecon < enBinStart || t->Erecon > (enBinStart+nEnergyBins*EnergyBinWidth) ) continue;
      
      //Cut out clipped events
      if ( t->Side==0 && ( t->xE.nClipped>0 || t->yE.nClipped>0 || t->xeRC<1 || t->xeRC>4 || t->yeRC<1 || t->yeRC>4 ) ) continue;
      else if ( t->Side==1 && ( t->xW.nClipped>0 || t->yW.nClipped>0 || t->xwRC<1 || t->xwRC>4 || t->ywRC<1 || t->ywRC>4) ) continue;

      // Check that the event has a reasonable position
      double r2E = t->xE.center*t->xE.center + t->yE.center*t->yE.center;
      double r2W = t->xW.center*t->xW.center + t->yW.center*t->yW.center;

      if ( r2E > 55.*55. || r2W > 55.*55. ) continue; 

		
      // Type 0 East Trigger
      int intBinX, intBinY;
      int intEnBin;
      
      if (t->Side == 0) {
	
	intBinX = posmap.getBinNumber(t->xE.center);
        intBinY = posmap.getBinNumber(t->yE.center);
	intEnBin = (int) ( ( t->Erecon - enBinStart ) / EnergyBinWidth );

	hisxy_E[intEnBin][intBinX][intBinY]->Fill(t->AnodeE);
      }

      // Type 0 West Trigger
      //moveOnX=moveOnY=true;
      else if (t->Side == 1) {

        intBinX = posmap.getBinNumber(t->xW.center);
        intBinY = posmap.getBinNumber(t->yW.center);
	intEnBin = (int) ( ( t->Erecon - enBinStart ) / EnergyBinWidth );

	hisxy_W[intEnBin][intBinX][intBinY]->Fill(t->AnodeW);
      }

    }

    // Close input ntuple
    delete t;

  }

  // Fit the histograms with a Landau peak and save the fitted values and then normalize to the center

  double mpv_E[nEnergyBins][nBinsXY][nBinsXY];
  double mpv_W[nEnergyBins][nBinsXY][nBinsXY];

  SinglePeakHist *peak;
  
  for (int p=0; p<nEnergyBins; p++) {
    for (int i=0; i<posmap.getNbinsXY(); i++) {
      for (int j=0; j<posmap.getNbinsXY(); j++) {

	peak = new SinglePeakHist(hisxy_E[p][i][j],0., 3999.5, true, 5, 2., 7., true);
	mpv_E[p][i][j] = peak->ReturnMean();
	if ( !peak->isGoodFit() ) mpv_E[p][i][j] = hisxy_E[p][i][j]->GetMean();
	delete peak;

	peak = new SinglePeakHist(hisxy_W[p][i][j],0., 3999.5, true, 5, 2., 7., true);
	mpv_W[p][i][j] = peak->ReturnMean();
	if ( !peak->isGoodFit() ) mpv_W[p][i][j] = hisxy_W[p][i][j]->GetMean();
	delete peak;
      }
    }
  }


  

  // Extract position maps
  double normE[nEnergyBins], normW[nEnergyBins];

  for (int p=0; p<nEnergyBins; p++) {
    normE[p] = mpv_E[p][nBinsXY/2][nBinsXY/2];
    normW[p] = mpv_W[p][nBinsXY/2][nBinsXY/2];
    cout << normE[p] << "\t" << normW[p] << endl;
  }

  //Checking for weird outliers
  for (int p=0; p<nEnergyBins; p++) {
    for (int i=0; i<nBinsXY; i++) {
      for (int j=0; j<nBinsXY; j++) {
	
	if ( mpv_E[p][i][j]<(0.25*normE[p]) || mpv_E[p][i][j]>(5.*normE[p]) ) 
	  mpv_E[p][i][j] = (0.25*normE[p]);
	if ( mpv_W[p][i][j]<(0.25*normW[p]) || mpv_W[p][i][j]>(5.*normW[p]) ) 
	  mpv_W[p][i][j] = (0.25*normW[p]);

      }
    }
  }

  double positionMap_E[nEnergyBins][nBinsXY][nBinsXY];
  double positionMap_W[nEnergyBins][nBinsXY][nBinsXY];
  
  for (int p=0; p<nEnergyBins; p++) {
    for (int i=0; i<nBinsXY; i++) {
      for (int j=0; j<nBinsXY; j++) {
        positionMap_E[p][i][j] = mpv_E[p][i][j] / normE[p];
        positionMap_W[p][i][j] = mpv_W[p][i][j] / normW[p];
      }
    }
  }

  // Write position maps to file

  for ( int p=0; p<nEnergyBins; ++p ) {
    
    string tempMap;
    tempMap = ( getenv("MWPC_CALIBRATION") + std::string("/position_maps/")+ tempOutBase + "_" +
		ftos( xyBinWidth )+ "mm_" + ftos( enBinStart + p*EnergyBinWidth ) +
		"-"+ftos( enBinStart + (p+1)*EnergyBinWidth )+".dat" );
    ofstream outMap(tempMap.c_str());
    outMap << std::setprecision(7);
    
    for (int i=0; i<nBinsXY; i++) {
      for (int j=0; j<nBinsXY; j++) {
	outMap << posmap.getBinCenter(i) << "  \t"
	       << posmap.getBinCenter(j) << "  \t"
	       << positionMap_E[p][i][j] << "  \t"
	       << positionMap_W[p][i][j] << endl;
      }
    }
    outMap.close();
  }

  /*// Write norms to file
  string tempNorm;
  tempNorm =  getenv("POSITION_MAPS") + std::string("norm_") + tempOutBase + "_" +ftos(xyBinWidth)+ "mm.dat";
  ofstream outNorm(tempNorm.c_str());

  for (int p=0; p<nPMT; p++) {
    outNorm << norm[p] << endl;
  }
  outNorm.close();*/

  // Write output ntuple
  fileOut->Write();
  delete fileOut;

  return 0;
}
