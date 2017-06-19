//SWANKS EDIT
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <string>

// ROOT libraries
#include "TRandom3.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TH1D.h>
#include <TSpectrum.h>

#include "fullTreeVariables.h"
#include "MWPCGeometry.h"
#include "pedestals.h"
#include "cuts.h"
#include "basic_reconstruction.h"
#include "DataTree.hh"
#include "MBUtils.hh"
#include "peaks.hh"
#include "BetaDecayTools.hh"


#include "positionMapHandler.hh"

#include "replay_pass2.h"

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

  /*bool allResponseClasses = true;
  int numResponseClasses = 0;
  vector <int> responseClasses;
  
  cout << "All Response Classes? (true=1/false=0): " << endl;
  cin  >> allResponseClasses;
  cout << endl;

  if (!allResponseClasses) {
    
    cout << "Enter number of response classes: " << endl;
    cin >> numResponseClasses;
    cout << endl;

    if (numResponseClasses<1 || numResponseClasses>9) {
      cout << "Bad number of response classes to include!!\n";
      exit(1);
    }
    responseClasses.resize(numResponseClasses,0);
    char quest[30];
    for (int i=0; i<numResponseClasses; i++) {
      sprintf(quest,"Enter class #%i: ",i+1);
      cout << quest;
      cin >> responseClasses[i];
      cout << endl;
      
      if (responseClasses[i]<0 || responseClasses[i]>8) {
	cout << "You entered a non-existent response class!\n";
	exit(1);
      }
    }
    }*/


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
  PositionMap posmap(xyBinWidth,50.);
  posmap.setRCflag(false); //telling the position map to not use the RC choices
  Int_t nBinsXY = posmap.getNbinsXY();
  

  // Open output ntuple
  string tempOutBase;
  string tempOut;
  //sprintf(tempOut, "position_map_%s.root", argv[1]);
  tempOutBase = "position_map_" + itos(iXeRunPeriod);
  /*if (!allResponseClasses) {
    tempOutBase+="_RC_";
    for (int i=0; i< numResponseClasses; i++) {
      tempOutBase+=itos(responseClasses[i]);
    }
    }*/
  tempOut =  getenv("POSITION_MAPS")+tempOutBase+"_"+ftos(xyBinWidth)+"mm.root";
  TFile *fileOut = new TFile(tempOut.c_str(),"RECREATE");

  // Output histograms
  int nPMT = 8;
  int nBinHist = 4100;//1025;

  TH1D *hisxy[nPMT][nBinsXY][nBinsXY];
  char *hisxyName = new char[10];

  for (int p=0; p<nPMT; p++) {
    for (int i=0; i<posmap.getNbinsXY(); i++) {
      for (int j=0; j<posmap.getNbinsXY(); j++) {
        if (p == 0)
          sprintf(hisxyName, "e0_%0.0f_%0.0f", posmap.getBinCenter(i), posmap.getBinCenter(j));
        if (p == 1)
          sprintf(hisxyName, "e1_%0.0f_%0.0f", posmap.getBinCenter(i), posmap.getBinCenter(j));
        if (p == 2)
          sprintf(hisxyName, "e2_%0.0f_%0.0f", posmap.getBinCenter(i), posmap.getBinCenter(j));
        if (p == 3)
          sprintf(hisxyName, "e3_%0.0f_%0.0f", posmap.getBinCenter(i), posmap.getBinCenter(j));
        if (p == 4)
          sprintf(hisxyName, "w0_%0.0f_%0.0f", posmap.getBinCenter(i), posmap.getBinCenter(j));
        if (p == 5)
          sprintf(hisxyName, "w1_%0.0f_%0.0f", posmap.getBinCenter(i), posmap.getBinCenter(j));
        if (p == 6)
          sprintf(hisxyName, "w2_%0.0f_%0.0f", posmap.getBinCenter(i), posmap.getBinCenter(j));
        if (p == 7)
          sprintf(hisxyName, "w3_%0.0f_%0.0f", posmap.getBinCenter(i), posmap.getBinCenter(j));

        hisxy[p][i][j] = new TH1D(hisxyName, "", nBinHist,-100.,4000.0);

      }
    }
  }

  // Loop through input ntuples
  char tempIn[500];
  for (int i=0; i<nRuns; i++) {

    // Open input ntuple
    sprintf(tempIn, "%s/replay_pass2_%i.root",getenv("REPLAY_PASS2"), runList[i]);
    DataTree *t = new DataTree();
    t->setupInputTree(std::string(tempIn),"pass2");

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

      // Select Type 0 events
      if (t->PID != 1) continue;
      if (t->Type != 0) continue;

      //Cut out clipped events
      if ( t->Side==0 && ( t->xE.nClipped>0 || t->yE.nClipped>0 || t->xeRC<1 || t->xeRC>4 || t->yeRC<1 || t->yeRC>4 ) ) continue;
      else if ( t->Side==1 && ( t->xW.nClipped>0 || t->yW.nClipped>0 || t->xwRC<1 || t->xwRC>4 || t->ywRC<1 || t->ywRC>4) ) continue;

		
      
      /*bool moveOnX = true, moveOnY=true; // Determining if the event is of the correct response class in x and y
     
	//Swank addition: Wire Chamber Response class. 
	for (int j=0; j<numResponseClasses; j++) {
	  if (t->xeRC == responseClasses[j]) {moveOnX=false;}
	  if (t->yeRC == responseClasses[j]) {moveOnY=false;}
	}
      
	if (moveOnX || moveOnY) continue;*/

      // Type 0 East Trigger
      int intBinX, intBinY; 
      
      if (t->Side == 0) {
	
	intBinX = posmap.getBinNumber(t->xE.center);
        intBinY = posmap.getBinNumber(t->yE.center);

        // Fill PMT histograms
        if (intBinX>-1 && intBinY>-1) hisxy[0][intBinX][intBinY]->Fill(t->ScintE.q1);
        if (intBinX>-1 && intBinY>-1) hisxy[1][intBinX][intBinY]->Fill(t->ScintE.q2);
        if (intBinX>-1 && intBinY>-1) hisxy[2][intBinX][intBinY]->Fill(t->ScintE.q3);
        if (intBinX>-1 && intBinY>-1) hisxy[3][intBinX][intBinY]->Fill(t->ScintE.q4);
      }

      // Type 0 West Trigger
      //moveOnX=moveOnY=true;
      else if (t->Side == 1) {

	//Swank Only Allow triangles!!!	  	
	//for (int j=0; j<numResponseClasses; j++) {
	// if (t->xwRC == responseClasses[j]) {moveOnX=false;}
	// if (t->ywRC == responseClasses[j]) {moveOnY=false;}
	//}
      
	//if (moveOnX || moveOnY) continue;
	
        intBinX = posmap.getBinNumber(t->xW.center);
        intBinY = posmap.getBinNumber(t->yW.center);

	// Fill PMT histograms 
        if (intBinX>-1 && intBinY>-1) hisxy[4][intBinX][intBinY]->Fill(t->ScintW.q1);
        if (intBinX>-1 && intBinY>-1) hisxy[5][intBinX][intBinY]->Fill(t->ScintW.q2);
        if (intBinX>-1 && intBinY>-1) hisxy[6][intBinX][intBinY]->Fill(t->ScintW.q3);
        if (intBinX>-1 && intBinY>-1) hisxy[7][intBinX][intBinY]->Fill(t->ScintW.q4);
      }


    }

    // Close input ntuple
    delete t;

  }


  //Rebinning the histograms based on the mean value...
  for (int p=0; p<nPMT; p++) {
    for (int i=0; i<nBinsXY; i++) {
      for (int j=0; j<nBinsXY; j++) {	

	hisxy[p][i][j]->GetXaxis()->SetRange(1,nBinHist);
	Double_t mean = hisxy[p][i][j]->GetMean();
	hisxy[p][i][j]->GetXaxis()->SetRange(0,nBinHist);
	double nGroup = 4.*mean/200.;
	nGroup = nGroup>1. ? nGroup : 1.;
	hisxy[p][i][j]->Rebin((int)nGroup);
	
      }
    }
  }

  // Extracting mean from 200 keV peak 

  // Define fit ranges
  double xLowBin[nPMT][nBinsXY][nBinsXY];
  double xHighBin[nPMT][nBinsXY][nBinsXY];
  double xLow[nPMT][nBinsXY][nBinsXY];
  double xHigh[nPMT][nBinsXY][nBinsXY];
  int maxBin[nPMT][nBinsXY][nBinsXY];
  double maxCounts[nPMT][nBinsXY][nBinsXY];
  double binCenterMax[nPMT][nBinsXY][nBinsXY];
  double meanVal[nPMT][nBinsXY][nBinsXY]; // Holds the mean of the distribution as defined by first fitting the peak
  double fitMean[nPMT][nBinsXY][nBinsXY]; // Holds the mean of the low energy peak
  double fitSigma[nPMT][nBinsXY][nBinsXY]; // Holds the sigma of the low energy peak
  double endpoint[nPMT][nBinsXY][nBinsXY]; // Holds the endpoint


  /////// First determine roughly where the low energy peak is
  TSpectrum *spec;

  for (int p=0; p<nPMT; p++) {
    for (int i=0; i<nBinsXY; i++) {
      for (int j=0; j<nBinsXY; j++) {	
	
	double r = sqrt(power(posmap.getBinCenter(j),2)+power(posmap.getBinCenter(i),2));
	
        // Find bin with maximum content
	hisxy[p][i][j]->GetXaxis()->SetRange(2,nBinHist);
        maxBin[p][i][j] = hisxy[p][i][j]->GetMaximumBin();
        maxCounts[p][i][j] = hisxy[p][i][j]->GetBinContent(maxBin[p][i][j]);
        binCenterMax[p][i][j] = hisxy[p][i][j]->GetBinCenter(maxBin[p][i][j]);
	
	  
	if (r<=(50.+2*xyBinWidth))
	      {
		spec = new TSpectrum(20);
		Int_t npeaks = spec->Search(hisxy[p][i][j],1.5,"",0.5);
		
		if (npeaks==0)
		  {
		    cout << "No peaks identified at PMT" << p << " position " << posmap.getBinCenter(i) << ", " << posmap.getBinCenter(j) << endl;
		  }
		else
		  {
		    Float_t *xpeaks = spec->GetPositionX(); // Note that newer versions of ROOT return a pointer to double...
		    TAxis *xaxis = (TAxis*)(hisxy[p][i][j]->GetXaxis());
		    Int_t peakBin=0;
		    Double_t BinSum=0.;
		    Double_t BinSumHold = 0.;
		    Int_t maxPeak=0.;
		    for (int pk=0;pk<npeaks;pk++) {
		      peakBin = xaxis->FindBin(xpeaks[pk]);
		      //Sum over 3 center bins of peak and compare to previos BinSum to see which peak is higher
		      BinSum = hisxy[p][i][j]->GetBinContent(peakBin) + hisxy[p][i][j]->GetBinContent(peakBin-1) + hisxy[p][i][j]->GetBinContent(peakBin+1);
		      if (BinSum>BinSumHold) {
			BinSumHold=BinSum;
			maxPeak=pk;
		      }
		    }
		    binCenterMax[p][i][j] = xpeaks[maxPeak];
		  }
		delete spec;
	      }
	xLow[p][i][j] = binCenterMax[p][i][j]*2./3.;
	xHigh[p][i][j] = 1.5*binCenterMax[p][i][j];

	
      }
    }
  }
  
  //////// Now fit the histograms for the peak

  for (int p=0; p<nPMT; p++) {
    for (int i=0; i<nBinsXY; i++) {
      for (int j=0; j<nBinsXY; j++) {	

	if ( hisxy[p][i][j]->Integral() > 500.) {// && r<=(50.+2*xBinWidth)) {

	  SinglePeakHist sing(hisxy[p][i][j], xLow[p][i][j], xHigh[p][i][j], true, 5, 0.8, 1.);

	  if (sing.isGoodFit() && sing.ReturnMean()>xLow[p][i][j] && sing.ReturnMean()<xHigh[p][i][j]) {
	    fitMean[p][i][j] = sing.ReturnMean();
	    fitSigma[p][i][j] = sing.ReturnSigma();
	  }

	  else  {
	    cout << "Can't converge on peak in PMT " << p << " at (" << posmap.getBinCenter(i) << ", " << posmap.getBinCenter(j) << "). Trying one more time......" << endl;
	    sing.SetRangeMin(xLow[p][i][j]);
	    sing.SetRangeMax(xHigh[p][i][j]);
	    sing.FitHist((double)maxBin[p][i][j], hisxy[p][i][j]->GetMean()/5., hisxy[p][i][j]->GetBinContent(maxBin[p][i][j]));

	    if (sing.isGoodFit() && sing.ReturnMean()>xLow[p][i][j] && sing.ReturnMean()<xHigh[p][i][j]) { 
	      fitMean[p][i][j] = sing.ReturnMean();
	      fitSigma[p][i][j] = sing.ReturnSigma();
	    }
	    else {
	      fitMean[p][i][j] = hisxy[p][i][j]->GetMean()/1.8;
	      int meanBin = hisxy[p][i][j]->GetXaxis()->FindBin(fitMean[p][i][j]);
	      int counts_check = hisxy[p][i][j]->GetBinContent(meanBin);
	      int counts = counts_check;
	      int bin=0;
	      if ( counts_check > 10 ) {
		while (counts>0.6*counts_check) {
		  bin++;
		  counts = hisxy[p][i][j]->GetBinContent(meanBin+bin);
		}
	      }
	  
	      xHighBin[p][i][j] = (meanBin+bin) < hisxy[p][i][j]->GetNbinsX()/3 ? (meanBin+bin): meanBin;
	      fitSigma[p][i][j] = hisxy[p][i][j]->GetBinCenter(xHighBin[p][i][j]) - fitMean[p][i][j];

	      cout << "Can't converge on peak in PMT " << p << " at bin (" << posmap.getBinCenter(i) << ", " << posmap.getBinCenter(j) << "). ";
	      cout << "**** replaced fit mean with hist_mean/1.8 " << fitMean[p][i][j] << endl;
	    }
	  }
	}
	else { 
	  fitMean[p][i][j] = hisxy[p][i][j]->GetMean()/1.8;
	  //if ( fitMean[p][i][j]>xLow[p][i][j] && fitMean[p][i][j]<xHigh[p][i][j] ) 
	  cout << "**** replaced fit mean with hist_mean/1.8 " << fitMean[p][i][j] << endl;
	  //else { 
	  //fitMean[p][i][j] = binCenterMax[p][i][j];
	  //cout << "**** replaced fit mean with binCenterMax " << fitMean[p][i][j] << endl;
	  //}
	  int meanBin = hisxy[p][i][j]->GetXaxis()->FindBin(fitMean[p][i][j]);
	  int counts_check = hisxy[p][i][j]->GetBinContent(meanBin);
	  int counts = counts_check;
	  int bin=0;
	  double frac = exp(-1/2.); // This should be the fraction of events for a gaussian at 1 sigma
	  if ( counts_check > 10 ) {
	    while (counts>frac*counts_check) {
	      bin++;
	      counts = hisxy[p][i][j]->GetBinContent(meanBin+bin);
	    }
	  }	  
	  xHighBin[p][i][j] = (meanBin+bin) < hisxy[p][i][j]->GetNbinsX()/3 ? (meanBin+bin): meanBin;
	  fitSigma[p][i][j] = hisxy[p][i][j]->GetBinCenter(xHighBin[p][i][j]) - fitMean[p][i][j];	 
	}

      }
    }
  }

  fileOut->Write(); // Writing the histograms with the peaks to file

  ////////// Now determine the mean of the Xe distribution in every position bin above the peak to avoid
  ////////// trigger effects
  
  double nSigmaFromMean = 1.5; // This is how far over from the peak we are starting the sum of the spectra

  for (int p=0; p<nPMT; p++) {
    for (int i=0; i<nBinsXY; i++) {
      for (int j=0; j<nBinsXY; j++) {

	hisxy[p][i][j]->GetXaxis()->SetRange(hisxy[p][i][j]->GetXaxis()->FindBin(fitMean[p][i][j]+nSigmaFromMean*fitSigma[p][i][j]), hisxy[p][i][j]->GetNbinsX()-1);
	meanVal[p][i][j] = hisxy[p][i][j]->GetMean();
	hisxy[p][i][j]->GetXaxis()->SetRange(0, hisxy[p][i][j]->GetNbinsX()); // Set the range back to the full range

      }
    }
  }

  ////////// Now determine the endpoint of the Xe distribution in every position bin

  double lowerBoundMult = 1.;
  double upperBoundMult = 2.;

  TFile *epfile = new TFile(TString::Format("%s/%s_%smm_endpoints.root",getenv("POSITION_MAPS"),tempOutBase.c_str(),ftos(xyBinWidth).c_str()),"RECREATE");

  TGraphErrors epgraph;
  
  KurieFitter kf;
  
  for (int p=0; p<nPMT; p++) {
    for (int i=0; i<nBinsXY; i++) {
      for (int j=0; j<nBinsXY; j++) {

	kf.FitSpectrum(hisxy[p][i][j], lowerBoundMult*meanVal[p][i][j], upperBoundMult*meanVal[p][i][j], 1.);
	endpoint[p][i][j] = kf.returnK0();

	epgraph = kf.returnKuriePlot();
	epgraph.SetName(TString::Format("pmt%i_x%0.1f_y%0.1f",p,posmap.getBinCenter(i), posmap.getBinCenter(j)));
	epgraph.Write();

      }
    }
  }

  delete epfile;
  delete fileOut;


  ///////////////////////// Extract position maps for meanVal ///////////////////////////////
  double norm[nPMT];
  for (int p=0; p<nPMT; p++) {
    norm[p] = meanVal[p][nBinsXY/2][nBinsXY/2];
    cout << norm[p] << endl;
  }

  //Checking for weird outliers
  for (int p=0; p<nPMT; p++) {
    for (int i=0; i<nBinsXY; i++) {
      for (int j=0; j<nBinsXY; j++) {
	
	if ( meanVal[p][i][j]<(0.1*norm[p]) || meanVal[p][i][j]>(8.*norm[p]) ) 
	  meanVal[p][i][j] = (0.1*norm[p]);

      }
    }
  }

  double positionMap[nPMT][nBinsXY][nBinsXY];
  for (int p=0; p<nPMT; p++) {
    for (int i=0; i<nBinsXY; i++) {
      for (int j=0; j<nBinsXY; j++) {
        positionMap[p][i][j] = meanVal[p][i][j] / norm[p];
      }
    }
  }

  // Write position maps to file
  TString mapFile = TString::Format("%s/%s_%0.1fmm.dat", getenv("POSITION_MAPS"),tempOutBase.c_str(),xyBinWidth);
  ofstream outMap(mapFile.Data());

  for (int i=0; i<nBinsXY; i++) {
    for (int j=0; j<nBinsXY; j++) {
      outMap << posmap.getBinCenter(i) << "  "
             << posmap.getBinCenter(j) << "  "
             << positionMap[0][i][j] << "  "
             << positionMap[1][i][j] << "  "
             << positionMap[2][i][j] << "  "
             << positionMap[3][i][j] << "  "
             << positionMap[4][i][j] << "  "
             << positionMap[5][i][j] << "  "
             << positionMap[6][i][j] << "  "
             << positionMap[7][i][j] << endl;
    }
  }
  outMap.close();

  // Write norms to file
  TString normFile = TString::Format("%s/norm_%s_%0.1fmm.dat", getenv("POSITION_MAPS"),tempOutBase.c_str(),xyBinWidth);
  ofstream outNorm(normFile.Data());

  for (int p=0; p<nPMT; p++) {
    outNorm << norm[p] << endl;
  }
  outNorm.close();

  ///////////////////////// Extract position maps for peaks ///////////////////////////////

  for (int p=0; p<nPMT; p++) {
    norm[p] = fitMean[p][nBinsXY/2][nBinsXY/2];
    cout << norm[p] << endl;
  }

  //Checking for weird outliers
  for (int p=0; p<nPMT; p++) {
    for (int i=0; i<nBinsXY; i++) {
      for (int j=0; j<nBinsXY; j++) {
	
	if ( fitMean[p][i][j]<(0.1*norm[p]) || fitMean[p][i][j]>(8.*norm[p]) ) 
	  fitMean[p][i][j] = (0.1*norm[p]);

      }
    }
  }

  for (int p=0; p<nPMT; p++) {
    for (int i=0; i<nBinsXY; i++) {
      for (int j=0; j<nBinsXY; j++) {
        positionMap[p][i][j] = fitMean[p][i][j] / norm[p];
      }
    }
  }

  // Write position maps to file
  mapFile = TString::Format("%s/%s_%0.1fmm_peakFits.dat", getenv("POSITION_MAPS"),tempOutBase.c_str(),xyBinWidth);
  outMap.open(mapFile.Data());
  

  for (int i=0; i<nBinsXY; i++) {
    for (int j=0; j<nBinsXY; j++) {
      outMap << posmap.getBinCenter(i) << "  "
             << posmap.getBinCenter(j) << "  "
             << positionMap[0][i][j] << "  "
             << positionMap[1][i][j] << "  "
             << positionMap[2][i][j] << "  "
             << positionMap[3][i][j] << "  "
             << positionMap[4][i][j] << "  "
             << positionMap[5][i][j] << "  "
             << positionMap[6][i][j] << "  "
             << positionMap[7][i][j] << endl;
    }
  }
  outMap.close();


  ///////////////////////// Extract position maps for endpoints ///////////////////////////////

  for (int p=0; p<nPMT; p++) {
    norm[p] = endpoint[p][nBinsXY/2][nBinsXY/2];
    cout << norm[p] << endl;
  }

  //Checking for weird outliers
  for (int p=0; p<nPMT; p++) {
    for (int i=0; i<nBinsXY; i++) {
      for (int j=0; j<nBinsXY; j++) {
	
	if ( endpoint[p][i][j]<(0.1*norm[p]) || endpoint[p][i][j]>(8.*norm[p]) ) 
	  endpoint[p][i][j] = (0.1*norm[p]);

      }
    }
  }

  for (int p=0; p<nPMT; p++) {
    for (int i=0; i<nBinsXY; i++) {
      for (int j=0; j<nBinsXY; j++) {
        positionMap[p][i][j] = endpoint[p][i][j] / norm[p];
      }
    }
  }

  // Write position maps to file
  mapFile = TString::Format("%s/%s_%0.1fmm_endpoints.dat", getenv("POSITION_MAPS"),tempOutBase.c_str(),xyBinWidth);
  outMap.open(mapFile.Data());

  for (int i=0; i<nBinsXY; i++) {
    for (int j=0; j<nBinsXY; j++) {
      outMap << posmap.getBinCenter(i) << "  "
             << posmap.getBinCenter(j) << "  "
             << positionMap[0][i][j] << "  "
             << positionMap[1][i][j] << "  "
             << positionMap[2][i][j] << "  "
             << positionMap[3][i][j] << "  "
             << positionMap[4][i][j] << "  "
             << positionMap[5][i][j] << "  "
             << positionMap[6][i][j] << "  "
             << positionMap[7][i][j] << endl;
    }
  }
  outMap.close();



  return 0;
}
