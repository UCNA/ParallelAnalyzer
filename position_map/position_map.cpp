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
#include "../source_peaks/peaks.hh"

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
  int nBinHist = 1025;

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
      if ( t->Side==0 && ( t->xE.nClipped>0 || t->yE.nClipped>0 ) ) continue;
      else if ( t->Side==1 && ( t->xW.nClipped>0 || t->yW.nClipped>0 ) ) continue;
	

	  
		

      
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

  // Extracting mean from 200 keV peak 

  // Define fit ranges
  double xLowBin[nPMT][nBinsXY][nBinsXY];
  double xHighBin[nPMT][nBinsXY][nBinsXY];
  double xLow[nPMT][nBinsXY][nBinsXY];
  double xHigh[nPMT][nBinsXY][nBinsXY];
  int maxBin[nPMT][nBinsXY][nBinsXY];
  double maxCounts[nPMT][nBinsXY][nBinsXY];
  double binCenterMax[nPMT][nBinsXY][nBinsXY];
  double meanVal[nPMT][nBinsXY][nBinsXY];

  /*for (int p=0; p<nPMT; p++) {
    for (int i=0; i<nBinsXY; i++) {
      for (int j=0; j<nBinsXY; j++) {	
	
	double r = sqrt(power(posmap.getBinCenter(j),2)+power(posmap.getBinCenter(i),2));
        // Find bin with maximum content
	hisxy[p][i][j]->GetXaxis()->SetRange(2,nBinHist);
        maxBin[p][i][j] = hisxy[p][i][j]->GetMaximumBin();
        maxCounts[p][i][j] = hisxy[p][i][j]->GetBinContent(maxBin[p][i][j]);
        binCenterMax[p][i][j] = hisxy[p][i][j]->GetBinCenter(maxBin[p][i][j]);
	
	//Determine Low edge bin
	int counts = maxCounts[p][i][j];
	int bin=0;
	while (counts>0.5*maxCounts[p][i][j]) {
	  bin++;
	  counts = hisxy[p][i][j]->GetBinContent(maxBin[p][i][j]-bin);
	}
	xLowBin[p][i][j] = (maxBin[p][i][j]-bin) > 0 ? (maxBin[p][i][j]-bin): 1;
	
      //Determine high edge bin
	counts = maxCounts[p][i][j];
	bin=0;
	while (counts>0.5*maxCounts[p][i][j]) {
	  bin++;
	  counts = hisxy[p][i][j]->GetBinContent(maxBin[p][i][j]+bin);
	}
	xHighBin[p][i][j] = (maxBin[p][i][j]+bin) < nBinHist ? (maxBin[p][i][j]+bin): nBinHist-1;

	//Determine mean of range
	//hisxy[p][i][j]->GetXaxis()->SetRange(xLowBin[p][i][j],xHighBin[p][i][j]);
	hisxy[p][i][j]->GetXaxis()->SetRange(xHighBin[p][i][j], nBinHist);
	meanVal[p][i][j] = hisxy[p][i][j]->GetMean();
	hisxy[p][i][j]->GetXaxis()->SetRange(0,nBinHist);
	
      }
    }
    }*/


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
		    Float_t *xpeaks = spec->GetPositionX();
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
  
  // Fit histograms
  //TF1 *gaussian_fit[nPMT][nPosBinsX][nPosBinsY];
  //double fitMean[nPMT][nPosBinsX][nPosBinsY];

  double sigmaMax[8] = {400.,250.,250.,250.,250.,300.,300.,500};

  for (int p=0; p<nPMT; p++) {
    for (int i=0; i<nBinsXY; i++) {
      for (int j=0; j<nBinsXY; j++) {

	double r = sqrt(power(posmap.getBinCenter(j),2)+power(posmap.getBinCenter(i),2));
	

	if ( hisxy[p][i][j]->Integral() > 500.) {// && r<=(50.+2*xBinWidth)) {

	  SinglePeakHist sing(hisxy[p][i][j], xLow[p][i][j], xHigh[p][i][j]);

	  if (sing.isGoodFit() && sing.ReturnMean()>xLow[p][i][j] && sing.ReturnMean()<xHigh[p][i][j]) {
	    meanVal[p][i][j] = sing.ReturnMean();
	  }

	  else  {
	    cout << "Can't converge on peak in PMT " << p << "at (" << posmap.getBinCenter(i) << ", " << posmap.getBinCenter(j) << "). Trying one more time......" << endl;
	    sing.SetRangeMin(xLow[p][i][j]);
	    sing.SetRangeMax(xHigh[p][i][j]);
	    sing.FitHist((double)maxBin[p][i][j], hisxy[p][i][j]->GetMean()/5., hisxy[p][i][j]->GetBinContent(maxBin[p][i][j]));

	    if (sing.isGoodFit() && sing.ReturnMean()>xLow[p][i][j] && sing.ReturnMean()<xHigh[p][i][j]) { 
	      meanVal[p][i][j] = sing.ReturnMean();
	    }
	    else {
	      meanVal[p][i][j] = hisxy[p][i][j]->GetMean()/1.8;
	      cout << "Can't converge on peak in PMT " << p << "at bin (" << posmap.getBinCenter(i) << ", " << posmap.getBinCenter(j) << "). ";
	      cout << "**** replaced fit mean with hist_mean/1.8 " << meanVal[p][i][j] << endl;
	    }
	  }
	}
	else { 
	  meanVal[p][i][j] = hisxy[p][i][j]->GetMean()/1.8;
	  if ( meanVal[p][i][j]>xLow[p][i][j] && meanVal[p][i][j]<xHigh[p][i][j] ) 
	    cout << "**** replaced fit mean with hist_mean/1.8 " << meanVal[p][i][j] << endl;
	  else { 
	    meanVal[p][i][j] = binCenterMax[p][i][j];
	    cout << "**** replaced fit mean with binCenterMax " << meanVal[p][i][j] << endl;
	  }
	}
      }
    }
  }

	  /*
        char fitName[500];
        sprintf(fitName, "gaussian_fit_%i_%i_%i.root", p, i, j);

        gaussian_fit[p][i][j] = new TF1(fitName,
       	     	                        "[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",
                                        xLow[p][i][j], xHigh[p][i][j]);
        gaussian_fit[p][i][j]->SetParameter(0,maxCounts[p][i][j]);
        gaussian_fit[p][i][j]->SetParameter(1,binCenterMax[p][i][j]);
        gaussian_fit[p][i][j]->SetParameter(2,100.);
	gaussian_fit[p][i][j]->SetParLimits(0,0.,5000.);
	gaussian_fit[p][i][j]->SetParLimits(1,0.,4000.);	
	gaussian_fit[p][i][j]->SetParLimits(2,0.,500.);

	

        if (maxCounts[p][i][j] > 0. && r<=(50.+2*xBinWidth)) {
          hisxy[p][i][j]->Fit(fitName, "LRQ");
          fitMean[p][i][j] = gaussian_fit[p][i][j]->GetParameter(1);
	  //attempt at manual fix of fitting problem. Constrain p1 and p2
	  if (fitMean[p][i][j]<xLow[p][i][j] || fitMean[p][i][j]>xHigh[p][i][j])
	    { 
	      string pmt;
	      if (p<4) pmt = "East";
	      else pmt = "West";
	      int pmtNum = p<4?p:(p-4);
	      cout << "Bad Fit in PMT " << pmt << " " << pmtNum << " at " << xBinCenter[i] << ", " << yBinCenter[j] <<  " -> " << fitMean[p][i][j] << endl;
	      Float_t lowerBoundPar1 = (float)xLow[p][i][j];
	      Float_t upperBoundPar1 = (float)xHigh[p][i][j];
	      Float_t lowerBoundPar2 = 0.;
	      Float_t upperBoundPar2 = sigmaMax[p];
	      gaussian_fit[p][i][j]->SetParameter(0,maxCounts[p][i][j]);
	      gaussian_fit[p][i][j]->SetParameter(1,binCenterMax[p][i][j]);
	      gaussian_fit[p][i][j]->SetParameter(2,100.);
	      gaussian_fit[p][i][j]->SetParLimits(1,lowerBoundPar1,upperBoundPar1);
	      gaussian_fit[p][i][j]->SetParLimits(2,lowerBoundPar2,upperBoundPar2);
	      hisxy[p][i][j]->Fit(fitName, "LRQ");
	      fitMean[p][i][j] = gaussian_fit[p][i][j]->GetParameter(1);
	      cout << "New Fit Mean -> " << fitMean[p][i][j] << endl;

	      //Check no. 2. Attempt to constrain p0 
	      if (fitMean[p][i][j]<(1.2*xLow[p][i][j]) || fitMean[p][i][j]>(0.8*xHigh[p][i][j]))
		{
		  Float_t lowerBoundPar0 = (float)(0.8*maxCounts[p][i][j]);
		  Float_t upperBoundPar0 = (float)(2.*maxCounts[p][i][j]);
		  gaussian_fit[p][i][j]->SetParameter(0,maxCounts[p][i][j]);
		  gaussian_fit[p][i][j]->SetParameter(1,binCenterMax[p][i][j]);
		  gaussian_fit[p][i][j]->SetParameter(2,100.);
		  gaussian_fit[p][i][j]->SetParLimits(0,lowerBoundPar0,upperBoundPar0);
		  gaussian_fit[p][i][j]->SetParLimits(1,lowerBoundPar1,upperBoundPar1);	
		  gaussian_fit[p][i][j]->SetParLimits(2,lowerBoundPar2,upperBoundPar2);
		  hisxy[p][i][j]->Fit(fitName, "LRQ");
		  fitMean[p][i][j] = gaussian_fit[p][i][j]->GetParameter(1);
		  cout << "New Fit Mean -> " << fitMean[p][i][j] << endl;
		
		  //Checking if fit mean is still out of range. If so, set value to max bin
		  if (fitMean[p][i][j]<xLow[p][i][j] || fitMean[p][i][j]>xHigh[p][i][j])
		    {
		      fitMean[p][i][j] = binCenterMax[p][i][j];
		      cout << "**** replaced fit mean with max bin " << fitMean[p][i][j] << endl;
		    }
		}
		}	      
        }
        else {
          meanVal[p][i][j] = 0.;
        }
	//cout << "PMT " << p << " at " << xBinCenter[i] << ", " << yBinCenter[j] << " fitMean = " << fitMean[p][i][j] << endl;

      }
    }
    }*/

  // Extract position maps
  double norm[nPMT];
  for (int p=0; p<nPMT; p++) {
    norm[p] = meanVal[p][nBinsXY/2][nBinsXY/2];
    cout << norm[p] << endl;
  }

  //Checking for weird outliers
  for (int p=0; p<nPMT; p++) {
    for (int i=0; i<nBinsXY; i++) {
      for (int j=0; j<nBinsXY; j++) {
	
	if ( meanVal[p][i][j]<(0.25*norm[p]) || meanVal[p][i][j]>(5.*norm[p]) ) 
	  meanVal[p][i][j] = (0.25*norm[p]);

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
  string tempMap;
  tempMap =  getenv("POSITION_MAPS") + tempOutBase + "_" +ftos(xyBinWidth)+ "mm.dat";
  ofstream outMap(tempMap.c_str());

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
  string tempNorm;
  tempNorm =  getenv("POSITION_MAPS") + std::string("norm_") + tempOutBase + "_" +ftos(xyBinWidth)+ "mm.dat";
  ofstream outNorm(tempNorm.c_str());

  for (int p=0; p<nPMT; p++) {
    outNorm << norm[p] << endl;
  }
  outNorm.close();

  // Write output ntuple
  fileOut->Write();
  delete fileOut;

  return 0;
}
