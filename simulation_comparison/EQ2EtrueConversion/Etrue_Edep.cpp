#include "Etrue_Edep.hh"
#include "cstdlib"
#include "fstream"
#include "iomanip"
#include "math.h"
#include "iostream"
#include "string"
#include "vector"
#include <TString.h>
#include <TRandom3.h>

#include "peaks.hh"

using std::ifstream;
using std::ofstream;
using namespace std;

DataTree::DataTree() : inputFile(NULL), inputTree(NULL) { };

DataTree::~DataTree() {


  if (inputTree) delete inputTree;
  if (inputFile) { inputFile->Close(); delete inputFile;}
};

void DataTree::setupInputTree(string inputFileName, string inputTreeName) {

  inputFile = new TFile(inputFileName.c_str(),"READ");
  inputTree = (TTree*)(inputFile->Get(inputTreeName.c_str()));
  
  inputTree->SetBranchAddress("primKE", &primKE);
  inputTree->SetBranchAddress("primTheta", &primTheta);
  inputTree->SetBranchAddress("Edep", &Edep);
  inputTree->SetBranchAddress("EdepQ", &EdepQ); 
  inputTree->SetBranchAddress("MWPCPos", &MWPCPos); 
  inputTree->SetBranchAddress("MWPCEnergy", &MWPCEnergy); 
  inputTree->SetBranchAddress("time", &time); 

  cout << "Prepared input tree " << inputTreeName << " in " << inputFileName << "...\n";
};


int main(int argc, char *argv[]) {

  if (argc!=2) {
    std::cout << "Usage: ./Etrue_Edep [geometry]\n";
    exit(0);
  }

  TString geometry = TString(argv[1]);

  if (geometry!=TString("2011-2012") && geometry!=TString("2012-2013") && geometry!=TString("2012-2013_isobutane")) {
    std::cout << "Bad choice in geometry!!\n";
    exit(0);
  }

  TFile *outputHists = new TFile(TString::Format("Hists_%s.root",geometry.Data()),"RECREATE");

  ofstream outdata;
  outdata.open (TString::Format("HistMeans_%s.dat",geometry.Data()).Data());
  outdata << setw(5) << "Nhist" << setw(7) << "Emin" << setw(7) << "Emax" 
                     << setw(11) << "EQType0E" << setw(11) << "meanErr" 
                     << setw(11) << "EQType0W" << setw(11) << "meanErr"  
                     << setw(11) << "EQType1E" << setw(11) << "meanErr" 
                     << setw(11) << "EQType1W" << setw(11) << "meanErr"  
                     << setw(11) << "EQType2E" << setw(11) << "meanErr" 
                     << setw(11) << "EQType2W" << setw(11) << "meanErr" 
                     << setw(11) << "EQType3E" << setw(11) << "meanErr" 
                     << setw(11) << "EQType3W" << setw(11) << "meanErr"
                     << setw(12) << "EQType23E" << setw(11) << "meanErr" 
                     << setw(12) << "EQType23W" << setw(11) << "meanErr" << endl; 

  char temp[100];
  int nFile = 1000;
  double MaxE = 800.0;
  double DeltaE = 10.0;

  

  double alpha = 0.4; // nPE/keV of roughly 400 PE per 1 GeV
  Double_t g_d = 16.;
  Double_t g_rest = 12500.;
  TRandom3 *rand = new TRandom3();
  TRandom3 *rand1 = new TRandom3(rand->Rndm()*1000);
  TRandom3 *rand2 = new TRandom3(rand->Rndm()*1000);

  #define pi M_PI

  vector<double> Emin;
  vector<double> Emax;

  vector<TH1D*> EDepQType0E;
  vector<TH1D*> EDepQType1E;
  vector<TH1D*> EDepQType2E;
  vector<TH1D*> EDepQType3E;
  vector<TH1D*> EDepQType23E;
  vector<TH1D*> EDepQType0W;
  vector<TH1D*> EDepQType1W;
  vector<TH1D*> EDepQType2W;
  vector<TH1D*> EDepQType3W;
  vector<TH1D*> EDepQType23W;

  int Nhist = 0;
  for (double E=DeltaE/2.; E<=MaxE; E+=DeltaE) Nhist++;
  int MaxNhist = Nhist;

  EDepQType0E.resize(MaxNhist);
  EDepQType1E.resize(MaxNhist);
  EDepQType2E.resize(MaxNhist);
  EDepQType3E.resize(MaxNhist);
  EDepQType23E.resize(MaxNhist);
  EDepQType0W.resize(MaxNhist);
  EDepQType1W.resize(MaxNhist);
  EDepQType2W.resize(MaxNhist);
  EDepQType3W.resize(MaxNhist);
  EDepQType23W.resize(MaxNhist);

  Nhist = 0;
  for (double E=DeltaE/2.; E<=MaxE; E+=DeltaE) 
      {
       Emin.push_back(E-DeltaE/2.);
       Emax.push_back(E+DeltaE/2.);

       sprintf(temp, "EDepQType0E_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType0E[Nhist] = new TH1D(temp, "Energy Deposition Quenched, TYPE 0 East", 200, 0. , 800.);

       sprintf(temp, "EDepQType0W_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType0W[Nhist] = new TH1D(temp, "Energy Deposition Quenched, TYPE 0 West", 200, 0. , 800.);

       sprintf(temp, "EDepQType1E_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType1E[Nhist] = new TH1D(temp, "Energy Deposition Quenched, TYPE 1 East", 200, 0. , 800.);

       sprintf(temp, "EDepQType1W_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType1W[Nhist] = new TH1D(temp, "Energy Deposition Quenched, TYPE 1 West", 200, 0. , 800.);

       sprintf(temp, "EDepQType2E_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType2E[Nhist] = new TH1D(temp, "Energy Deposition Quenched, TYPE 2 East", 200, 0. , 800.);

       sprintf(temp, "EDepQType2W_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType2W[Nhist] = new TH1D(temp, "Energy Deposition Quenched, TYPE 2 West", 200, 0. , 800.);

       sprintf(temp, "EDepQType3E_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType3E[Nhist] = new TH1D(temp, "Energy Deposition Quenched, TYPE 3 East", 200, 0. , 800.);

       sprintf(temp, "EDepQType3W_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType3W[Nhist] = new TH1D(temp, "Energy Deposition Quenched, TYPE 3 West", 200, 0. , 800.);

       sprintf(temp, "EDepQType23E_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType23E[Nhist] = new TH1D(temp, "Energy Deposition Quenched, TYPE 23 East", 200, 0. , 800.);

       sprintf(temp, "EDepQType23W_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType23W[Nhist] = new TH1D(temp, "Energy Deposition Quenched, TYPE 23 West", 200, 0. , 800.);

       Nhist++;
      }
 
  double FidRad = 45.0; // Set to 45 to eliminate any contamination from edge effects

  // We loop over both polarization states to mix the two...
  std::vector <TString> pols {"polE","polW"};

  for ( auto pol : pols ) {
  
    TString path = TString::Format("/extern/mabrow05/ucna/geant4work/output/official_%s/Beta_%s/",geometry.Data(),pol.Data());

    for (int nF = 0; nF<nFile; nF++) 
      {
	DataTree Tin;
	sprintf(temp, "%s/analyzed_%d.root", path.Data(), nF);
	Tin.setupInputTree (temp, "anaTree");
	int Nevt = Tin.getEntries();
	for (Int_t evt = 0; evt<Nevt; evt++) 
	  {
            Tin.getEvent(evt);
	    
	    double r2E = ( Tin.MWPCPos.MWPCPosE[0]*Tin.MWPCPos.MWPCPosE[0] + Tin.MWPCPos.MWPCPosE[1]*Tin.MWPCPos.MWPCPosE[1] ) * 0.6 * 100.;
	    double r2W = ( Tin.MWPCPos.MWPCPosW[0]*Tin.MWPCPos.MWPCPosW[0] + Tin.MWPCPos.MWPCPosW[1]*Tin.MWPCPos.MWPCPosW[1] ) * 0.6 * 100.;
	      
	    
	    if ( r2E > FidRad*FidRad || r2W > FidRad*FidRad ) continue; 

	    //Calculate smeared energies
	    //Tin.EdepQ.EdepQE = rand->Gaus(Tin.EdepQ.EdepQE, sqrt(Tin.EdepQ.EdepQE/alpha));
	    //Tin.EdepQ.EdepQW = rand->Gaus(Tin.EdepQ.EdepQW, sqrt(Tin.EdepQ.EdepQW/alpha));
	    Tin.EdepQ.EdepQE = (1./(alpha*g_d*g_rest)) * (rand->Poisson(g_rest*rand2->Poisson(g_d*rand1->Poisson(alpha*Tin.EdepQ.EdepQE))));
	    Tin.EdepQ.EdepQW = (1./(alpha*g_d*g_rest)) * (rand->Poisson(g_rest*rand2->Poisson(g_d*rand1->Poisson(alpha*Tin.EdepQ.EdepQW))));

	       
            for (Nhist=0; Nhist<MaxNhist; Nhist++) //Sets Nhist to the proper histogram based on the true energy of the event
	      if (Emin[Nhist]<=Tin.primKE && Tin.primKE<Emax[Nhist]) break; 
	    
            //TYPE 0E EVENTS
            if (Tin.MWPCEnergy.MWPCEnergyE > 0 && Tin.EdepQ.EdepQE > 0 &&
                Tin.MWPCEnergy.MWPCEnergyW == 0 && Tin.EdepQ.EdepQW == 0) 
	      EDepQType0E[Nhist]->Fill(Tin.EdepQ.EdepQE);
	    
            //TYPE 0W EVENTS
            else if (Tin.MWPCEnergy.MWPCEnergyE == 0 && Tin.EdepQ.EdepQE == 0 &&
                     Tin.MWPCEnergy.MWPCEnergyW > 0 && Tin.EdepQ.EdepQW > 0)
	      EDepQType0W[Nhist]->Fill(Tin.EdepQ.EdepQW);
	    
            //TYPE 1E & 1W EVENTS
            else if (Tin.MWPCEnergy.MWPCEnergyE > 0 && Tin.EdepQ.EdepQE > 0 &&
                     Tin.MWPCEnergy.MWPCEnergyW > 0 && Tin.EdepQ.EdepQW > 0)
	      {
		if (Tin.time.timeE < Tin.time.timeW) EDepQType1E[Nhist]->Fill(Tin.EdepQ.EdepQE+Tin.EdepQ.EdepQW);
		else EDepQType1W[Nhist]->Fill(Tin.EdepQ.EdepQW+Tin.EdepQ.EdepQE);
	      }
	    
            //TYPE 3E & 2E EVENTS
            else if (Tin.MWPCEnergy.MWPCEnergyE > 0 && Tin.EdepQ.EdepQE > 0 &&
                     Tin.MWPCEnergy.MWPCEnergyW > 0 && Tin.EdepQ.EdepQW == 0)
	      {
		EDepQType23E[Nhist]->Fill(Tin.EdepQ.EdepQE);
		if (Tin.primTheta>pi/2.) EDepQType3E[Nhist]->Fill(Tin.EdepQ.EdepQE);
		else EDepQType2E[Nhist]->Fill(Tin.EdepQ.EdepQE);
	      }
	    
            //TYPE 2W & 3W EVENTS
            else if (Tin.MWPCEnergy.MWPCEnergyE > 0 && Tin.EdepQ.EdepQE == 0 &&
                     Tin.MWPCEnergy.MWPCEnergyW > 0 && Tin.EdepQ.EdepQW > 0)
	      {
		EDepQType23W[Nhist]->Fill(Tin.EdepQ.EdepQW);
		if (Tin.primTheta<pi/2) EDepQType3W[Nhist]->Fill(Tin.EdepQ.EdepQW);
		else EDepQType2W[Nhist]->Fill(Tin.EdepQ.EdepQW);
	      }

	  }
	cout << "******* analyzed_" << nF << ".root done *******" << endl;
	cout << endl;
      }
    cout << "//////////////////" << pol.Data() << " finished //////////////////" << endl;
    cout << endl;
  }


  // NEED TO PUT IN FITS TO THESE SPECTRA
  vector<double> mean_EDepQType0E(MaxNhist,0.);
  vector<double> mean_EDepQType1E(MaxNhist,0.);
  vector<double> mean_EDepQType2E(MaxNhist,0.);
  vector<double> mean_EDepQType3E(MaxNhist,0.);
  vector<double> mean_EDepQType23E(MaxNhist,0.);
  vector<double> mean_EDepQType0W(MaxNhist,0.);
  vector<double> mean_EDepQType1W(MaxNhist,0.);
  vector<double> mean_EDepQType2W(MaxNhist,0.);
  vector<double> mean_EDepQType3W(MaxNhist,0.);
  vector<double> mean_EDepQType23W(MaxNhist,0.);

  vector<double> err_EDepQType0E(MaxNhist,0.);
  vector<double> err_EDepQType1E(MaxNhist,0.);
  vector<double> err_EDepQType2E(MaxNhist,0.);
  vector<double> err_EDepQType3E(MaxNhist,0.);
  vector<double> err_EDepQType23E(MaxNhist,0.);
  vector<double> err_EDepQType0W(MaxNhist,0.);
  vector<double> err_EDepQType1W(MaxNhist,0.);
  vector<double> err_EDepQType2W(MaxNhist,0.);
  vector<double> err_EDepQType3W(MaxNhist,0.);
  vector<double> err_EDepQType23W(MaxNhist,0.);
  
  SinglePeakHist *p;
  
  for (Nhist=0; Nhist<MaxNhist; Nhist++) {

    p = new SinglePeakHist(EDepQType0E[Nhist],0., 800., true, 5, 0.75,1.5);
    mean_EDepQType0E[Nhist] = p->ReturnMean();
    err_EDepQType0E[Nhist] = p->ReturnMeanError();
    delete p;

    p = new SinglePeakHist(EDepQType0W[Nhist],0., 800., true, 5, 0.75,1.5);
    mean_EDepQType0W[Nhist] = p->ReturnMean();
    err_EDepQType0W[Nhist] = p->ReturnMeanError();
    delete p;

    p = new SinglePeakHist(EDepQType1E[Nhist],0., 800., true, 5, 0.75,1.5);
    mean_EDepQType1E[Nhist] = p->ReturnMean();
    err_EDepQType1E[Nhist] = p->ReturnMeanError();
    delete p;

    p = new SinglePeakHist(EDepQType1W[Nhist],0., 800., true, 5, 0.75,1.5);
    mean_EDepQType1W[Nhist] = p->ReturnMean();
    err_EDepQType1W[Nhist] = p->ReturnMeanError();
    delete p;

    p = new SinglePeakHist(EDepQType2E[Nhist],0., 800., true, 5, 0.75,1.5);
    mean_EDepQType2E[Nhist] = p->ReturnMean();
    err_EDepQType2E[Nhist] = p->ReturnMeanError();
    delete p;

    p = new SinglePeakHist(EDepQType2W[Nhist],0., 800., true, 5, 0.75,1.5);
    mean_EDepQType2W[Nhist] = p->ReturnMean();
    err_EDepQType2W[Nhist] = p->ReturnMeanError();
    delete p;
    
    p = new SinglePeakHist(EDepQType3E[Nhist],0., 800., true, 5, 0.75,1.5);
    mean_EDepQType3E[Nhist] = p->ReturnMean();
    err_EDepQType3E[Nhist] = p->ReturnMeanError();
    delete p;

    p = new SinglePeakHist(EDepQType3W[Nhist],0., 800., true, 5, 0.75,1.5);
    mean_EDepQType3W[Nhist] = p->ReturnMean();
    err_EDepQType3W[Nhist] = p->ReturnMeanError();
    delete p;

    p = new SinglePeakHist(EDepQType23E[Nhist],0., 800., true, 5, 0.75,1.5);
    mean_EDepQType23E[Nhist] = p->ReturnMean();
    err_EDepQType23E[Nhist] = p->ReturnMeanError();
    delete p;

    p = new SinglePeakHist(EDepQType23W[Nhist],0., 800., true, 5, 0.75,1.5);
    mean_EDepQType23W[Nhist] = p->ReturnMean();
    err_EDepQType23W[Nhist] = p->ReturnMeanError();
    delete p;

  }
    

  for (Nhist=0; Nhist<MaxNhist; Nhist++)
      {
       outdata << setw(5)  << Nhist << setw(7) << Emin[Nhist] << setw(7) << Emax[Nhist]
               << setw(11) << mean_EDepQType0E[Nhist] << setw(11) << err_EDepQType0E[Nhist] 
               << setw(11) << mean_EDepQType0W[Nhist] << setw(11) << err_EDepQType0W[Nhist]
               << setw(11) << mean_EDepQType1E[Nhist] << setw(11) << err_EDepQType1E[Nhist] 
               << setw(11) << mean_EDepQType1W[Nhist] << setw(11) << err_EDepQType1W[Nhist]
               << setw(11) << mean_EDepQType2E[Nhist] << setw(11) << err_EDepQType2E[Nhist] 
               << setw(11) << mean_EDepQType2W[Nhist] << setw(11) << err_EDepQType2W[Nhist]
               << setw(11) << mean_EDepQType3E[Nhist] << setw(11) << err_EDepQType3E[Nhist] 
               << setw(11) << mean_EDepQType3W[Nhist] << setw(11) << err_EDepQType3W[Nhist]  
               << setw(12) << mean_EDepQType23E[Nhist] << setw(11) << err_EDepQType23E[Nhist]
               << setw(12) << mean_EDepQType23W[Nhist] << setw(11) << err_EDepQType23W[Nhist] << endl;  
       }
  
  outputHists->Write();
  outputHists->Close();
}
  
