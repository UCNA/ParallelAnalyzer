/* Code to take a run number, retrieve it's runperiod, and construct the 
weighted spectra which would be seen as a reconstructed energy on one side
of the detector. Also applies the trigger functions */

#include "revCalSim.h"
#include "posMapReader.h"
#include "sourcePeaks.h"
#include "runInfo.h"

#include "../Asymmetry/SQLinterface.hh"

#include <vector>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include <TRandom3.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH1D.h>
#include <TChain.h>
#include <TMath.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>


using namespace std;


int getPolarization(int run) {
 
  std::string dbAddress = std::string(getenv("UCNADBADDRESS"));
  std::string dbname = std::string(getenv("UCNADB"));
  std::string dbUser = std::string(getenv("UCNADBUSER"));
  std::string dbPass = std::string(getenv("UCNADBPASS"));
  
  char cmd[200];
  sprintf(cmd,"SELECT flipper FROM run WHERE run_number=%i;",run);

  SQLdatabase *db = new SQLdatabase(dbname, dbAddress, dbUser, dbPass);
  db->fetchQuery(cmd);
  std::string flipperStatus = db->returnQueryEntry();
  delete db;

  std::cout << flipperStatus << std::endl;
  if (flipperStatus=="On") return 1;
  else if (flipperStatus=="Off") return -1;
  else {
    std::cout <<  "Polarization isn't applicaple or you chose a Depol Run";
    return 0;
  }
};

vector < vector < Double_t > > getTriggerFunctionParams(Int_t XeRunPeriod, Int_t nParams=8) {
  Char_t file[500];
  sprintf(file,"%s/trigger_functions_XePeriod_%i.dat",getenv("TRIGGER_FUNC"),XeRunPeriod);
  ifstream infile(file);
  vector < vector <Double_t> > func;
  func.resize(2,vector <Double_t> (nParams,0.));
  //cout << "made it here\n";
  for (Int_t side = 0; side<2; side++) {
    Int_t param = 0;
    while (param<nParams) {
      infile >> func[side][param];
      cout << func[side][param] << " ";
      param++;
    }
    cout << endl;
  }
  infile.close();
  return func;
}

Double_t triggerProbability(vector <Double_t> params, Double_t En) {
  //Double_t prob = params[0]+params[1]*TMath::Erf((En-params[2])/params[3])
  //+ params[4]*TMath::Gaus(En,params[5],params[6]);
  if (En<100.) return (params[0]+params[1]*TMath::Erf((En-params[2])/params[3]))*(0.5-0.5*TMath::TanH((En-params[2])/params[4])) + 
    (0.5+0.5*TMath::TanH((En-params[2])/params[4]))*(params[5]+params[6]*TMath::TanH((En-params[2])/params[7]));
  else return 1.;
}

vector <vector <double> > returnSourcePosition (Int_t runNumber, string src) {
  Char_t temp[500];
  sprintf(temp,"%s/source_list_%i.dat",getenv("SOURCE_LIST"),runNumber);
  ifstream file(temp);
  cout << src << endl;
  int num = 0;
  file >> num;
  cout << num << endl;
  int srcNum = 0;
  string src_check;
  for (int i=0; i<num;srcNum++,i++) {
    file >> src_check;
    cout << src_check << endl;
    if (src_check==src) break;   
  }
  cout << "The source Number is: " << srcNum << endl;
  if (srcNum==num) {
    cout << "Didn't find source in that run\n"; exit(0);
  }
  file.close();
  
  sprintf(temp,"%s/source_positions_%i.dat",getenv("SOURCE_POSITIONS"),runNumber);
  file.open(temp);
  
  vector < vector < double > > srcPos;
  srcPos.resize(2,vector <double> (3,0.));
  
  for (int i=0; i<srcNum+1; i++) {
    for (int j=0; j<2; j++) {
      for (int jj=0; jj<3; jj++) {
	file >> srcPos[j][jj];
      }
    }
  }
  return srcPos;
}

vector <Int_t> getPMTQuality(Int_t runNumber) {
  //Read in PMT quality file
  cout << "Reading in PMT Quality file ...\n";
  vector <Int_t>  pmtQuality (8,0);
  Char_t temp[200];
  sprintf(temp,"%s/residuals/PMT_runQuality_master.dat",getenv("ANALYSIS_CODE")); 
  ifstream pmt;
  std::cout << temp << std::endl;
  pmt.open(temp);
  Int_t run_hold;
  while (pmt >> run_hold >> pmtQuality[0] >> pmtQuality[1] >> pmtQuality[2]
	 >> pmtQuality[3] >> pmtQuality[4] >> pmtQuality[5]
	 >> pmtQuality[6] >> pmtQuality[7]) {
    if (run_hold==runNumber) break;
    if (pmt.fail()) break;
  }
  pmt.close();
  if (run_hold!=runNumber) {
    cout << "Run not found in PMT quality file!" << endl;
    exit(0);
  }
  return pmtQuality;
}

vector < Double_t > GetAlphaValues(Int_t runPeriod)
{
  Char_t temp[500];
  vector < Double_t > alphas (8,0.);
  sprintf(temp,"%s/simulation_comparison/nPE_per_keV/nPE_per_keV_%i.dat",getenv("ANALYSIS_CODE"),runPeriod);
  ifstream infile;
  infile.open(temp);
  Int_t i = 0;

  while (infile >> alphas[i]) { std::cout << alphas[i] << std::endl; i++; }
  return alphas;
}

//Returns the average eta value from the position map for the Sn source used to calculate alpha
vector < Double_t > getMeanEtaForAlpha(Int_t run) 
{
  Char_t temp[500];
  vector < Double_t > eta (8,0.);
  sprintf(temp,"%s/nPE_meanEtaVal_%i.dat",getenv("NPE_WEIGHTS"),run);
  ifstream infile;
  infile.open(temp);
  Int_t i = 0;

  while (infile >> eta[i]) i++;
  return eta;
}
  
//Get the conversion from EQ2Etrue
std::vector < std::vector < std::vector <double> > > getEQ2EtrueParams(int runNumber) {
  ifstream infile;
  std::string basePath = getenv("ANALYSIS_CODE"); 
  if (runNumber<16000) basePath+=std::string("/simulation_comparison/EQ2EtrueConversion/2011-2012_EQ2EtrueFitParams.dat");
  else if (runNumber<20000) basePath+=std::string("/simulation_comparison/EQ2EtrueConversion/2011-2012_EQ2EtrueFitParams.dat");
  else if (runNumber<24000) basePath+=std::string("/simulation_comparison/EQ2EtrueConversion/2011-2012_EQ2EtrueFitParams.dat");
  else {
    std::cout << "Bad runNumber passed to getEQ2EtrueParams\n";
    exit(0);
  }
  infile.open(basePath.c_str());
  std::vector < std::vector < std::vector < double > > > params;
  params.resize(2,std::vector < std::vector < double > > (3, std::vector < double > (6,0.)));

  char holdType[10];
  int side=0, type=0;
  while (infile >> holdType >> params[side][type][0] >> params[side][type][1] >> params[side][type][2] >> params[side][type][3] >> params[side][type][4] >> params[side][type][5]) { 
    std::cout << holdType << " " << params[side][type][0] << " " << params[side][type][1] << " " << params[side][type][2] << " " << params[side][type][3] << " " << params[side][type][4] << " " << params[side][type][5] << std::endl;
    type+=1;
    if (type==3) {type=0; side=1;}
  }
  return params;
}


void SetUpTree(TTree *tree) {
  tree->Branch("PID", &PID, "PID/I");
  tree->Branch("side", &side, "side/I");
  tree->Branch("type", &type, "type/I");
  tree->Branch("Erecon", &Erecon,"Erecon/D");
  
  tree->Branch("Evis",&evis,"EvisE/D:EvisW");
  tree->Branch("Edep",&edep,"EdepE/D:EdepW");
  tree->Branch("EdepQ",&edepQ,"EdepQE/D:EdepQW");
  tree->Branch("Eprim",&Eprim,"Eprim/D");
  tree->Branch("AsymWeight",&AsymWeight,"AsymWeight/D");
  
  tree->Branch("time",&Time,"timeE/D:timeW");
  tree->Branch("MWPCEnergy",&mwpcE,"MWPCEnergyE/D:MWPCEnergyW");
  tree->Branch("MWPCPos",&mwpc_pos,"MWPCPosE[3]/D:MWPCPosW[3]");
  tree->Branch("ScintPos",&scint_pos,"ScintPosE[3]/D:ScintPosW[3]");
  tree->Branch("ScintPosAdjusted",&scint_pos_adj,"ScintPosAdjE[3]/D:ScintPosAdjW[3]");
  tree->Branch("PMT",&pmt,"Evis0/D:Evis1:Evis2:Evis3:Evis4:Evis5:Evis6:Evis7:etaEvis0/D:etaEvis1:etaEvis2:etaEvis3:etaEvis4:etaEvis5:etaEvis6:etaEvis7:nPE0/D:nPE1:nPE2:nPE3:nPE4:nPE5:nPE6:nPE7");
  
}
  

void revCalSimulation (Int_t runNumber, string source) 
{
  cout << "Running reverse calibration for run " << runNumber << " and source " << source << endl;

  //Start by getting the source position if not a Beta decay run
  vector < vector <double> > srcPos;

  if (source!="Beta") {
    string srcShort = source;
    srcShort.erase(2);
    srcPos = returnSourcePosition(runNumber, srcShort);
  }

  //If the source is In114, checking which side the Indium was facing
  if (source=="In114") {
    string side = getIndiumSide(runNumber);
    if (side=="East") source+="E";
    else if (side=="West") source+="W";
    else {cout << "can't find Indium in run requested\n"; exit(0);}
  }

  //First get the number of total electron events around the source position from the data file and also the length of the run
  Char_t temp[500],tempE[500],tempW[500];
  //sprintf(temp,"%s/replay_pass4_%i.root",getenv("REPLAY_PASS4"),runNumber);
  sprintf(temp,"%s/replay_pass2_%i.root",getenv("REPLAY_PASS2"),runNumber);
  TFile *dataFile = new TFile(temp,"READ");
  TTree *data = (TTree*)(dataFile->Get("pass2"));

  string outputBase;
  if (source!="Beta") {
    sprintf(tempE,"Type<4 && PID==1 && Side==0 && (xE.center>(%f-2.*%f) && xE.center<(%f+2.*%f) && yE.center>(%f-2.*%f) && yE.center<(%f+2.*%f))",srcPos[0][0],fabs(srcPos[0][2]),srcPos[0][0],fabs(srcPos[0][2]),srcPos[0][1],fabs(srcPos[0][2]),srcPos[0][1],fabs(srcPos[0][2]));
    sprintf(tempW,"Type<4 && PID==1 && Side==1 && (xW.center>(%f-2.*%f) && xW.center<(%f+2.*%f) && yW.center>(%f-2.*%f) && yW.center<(%f+2.*%f))",srcPos[1][0],fabs(srcPos[1][2]),srcPos[1][0],fabs(srcPos[1][2]),srcPos[1][1],fabs(srcPos[1][2]),srcPos[1][1],fabs(srcPos[1][2]));
    outputBase = string(getenv("REVCALSIM")) + "sources/";
  }
  else if (source=="Beta") {
    sprintf(tempE,"Type<4 && PID==1 && Side==0 && (xE.center*xE.center+yE.center*yE.center)<2500.");
    sprintf(tempW,"Type<4 && PID==1 && Side==1 && (xW.center*xW.center+yW.center*yW.center)<2500.");
    outputBase = string(getenv("REVCALSIM")) + "beta/";
  }
  UInt_t BetaEvents = data->GetEntries(tempE) + data->GetEntries(tempW);
  cout << "Electron Events in Data file: " << BetaEvents << endl;
  cout << "East: " << data->GetEntries(tempE) << "\tWest: " << data->GetEntries(tempW) << endl;
  delete data;
  dataFile->Close();

  //Create simulation output file
  Char_t outputfile[500];
  sprintf(outputfile,"%s/revCalSim_%i_%s.root",outputBase.c_str(),runNumber,source.c_str());
  
  TFile *outfile = new TFile(outputfile, "RECREATE");


  ///////////////////////// SETTING GAIN OF FIRST DYNYODE
  Double_t g_d = 16.;

  /////// Loading other run dependent quantities
  vector <Int_t> pmtQuality = getPMTQuality(runNumber); // Get the quality of the PMTs for that run
  UInt_t calibrationPeriod = getSrcRunPeriod(runNumber); // retrieve the calibration period for this run
  UInt_t XePeriod = getXeRunPeriod(runNumber); // Get the proper Xe run period for the Trigger functions
  GetPositionMap(XePeriod);
  vector <Double_t> alpha = GetAlphaValues(calibrationPeriod); // fill vector with the alpha (nPE/keV) values for this run period
  vector <Double_t> meanEta = getMeanEtaForAlpha(runNumber); // fill vector with the calculated mean position correction for the source used to calculate alpha
  vector < vector <Double_t> > triggerFunc = getTriggerFunctionParams(XePeriod,8); // 2D vector with trigger function for East side and West side in that order
  std::vector < std::vector < std::vector <double> > > EQ2Etrue = getEQ2EtrueParams(runNumber);
  Int_t pol = getPolarization(runNumber);

  //Setup the output tree
  TTree *tree = new TTree("revCalSim", "revCalSim");
  SetUpTree(tree); //Setup the output tree and branches

  //Histograms of event types for quick checks
  vector <TH1D*> finalEn (6,NULL);
  finalEn[0] = new TH1D("finalE0", "Simulated Weighted Sum East Type 0", 400, 0., 1200.);
  finalEn[1] = new TH1D("finalW0", "Simulated Weighted Sum West Type 0", 400, 0., 1200.);
  finalEn[2] = new TH1D("finalE1", "Simulated Weighted Sum East Type 1", 400, 0., 1200.);
  finalEn[3] = new TH1D("finalW1", "Simulated Weighted Sum West Type 1", 400, 0., 1200.);
  finalEn[4] = new TH1D("finalE23", "Simulated Weighted Sum East Type 2/3", 400, 0., 1200.);
  finalEn[5] = new TH1D("finalW23", "Simulated Weighted Sum West Type 2/3", 400, 0., 1200.);

  //Decide which simulation to use...
  std::string simLocation;
  TChain *chain = new TChain("anaTree");
  
  if (runNumber<20000) simLocation = string(getenv("SIM_2011_2012"));
  else if (runNumber>21087 && runNumber<21679) simLocation = string(getenv("SIM_2012_2013_ISOBUTANE"));
  else simLocation = string(getenv("SIM_2012_2013"));

  std::cout << "Using simulation from " << simLocation << "...\n";

  //Read in simulated data and put in a TChain
  TRandom3 *randFile = new TRandom3(0);
  int numFiles = source=="Beta"?400:250;
  int fileNum = source=="Beta" ? (int)(randFile->Rndm()*400) : 0;
  delete randFile;
  for (int i=0; i<numFiles; i++) {
    if (fileNum==numFiles) fileNum=0;
    sprintf(temp,"%s/%s/analyzed_%i.root",simLocation.c_str(),source.c_str(),fileNum);
    //sprintf(temp,"/extern/mabrow05/ucna/XuanSim/%s/xuan_analyzed.root",source.c_str());
    chain->AddFile(temp);
    fileNum++;
  }

  
  // Set the addresses of the information read in from the simulation file
  chain->SetBranchAddress("MWPCEnergy",&mwpcE);
  chain->SetBranchAddress("time",&Time);
  chain->SetBranchAddress("Edep",&edep);
  chain->SetBranchAddress("EdepQ",&edepQ);
  chain->SetBranchAddress("MWPCPos",&mwpc_pos);
  chain->SetBranchAddress("ScintPos",&scint_pos);
  chain->SetBranchAddress("primKE",&Eprim);
  chain->SetBranchAddress("primTheta",&primTheta);
  //chain->GetBranch("EdepQ")->GetLeaf("EdepQE")->SetAddress(&EdepQE);
  //chain->GetBranch("EdepQ")->GetLeaf("EdepQW")->SetAddress(&EdepQW);
  //chain->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyE")->SetAddress(&MWPCEnergyE);
  //chain->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyW")->SetAddress(&MWPCEnergyW);

  //These are for feeding in Xuan's simulations... this needs to be updated so that I can pass a flag and change these on the fly
  //chain->SetBranchAddress("PrimaryParticleSpecies",&primaryID);
  //chain->SetBranchAddress("PrimaryParticleSpecies",&primaryID);
  //chain->SetBranchAddress("mwpcEnergy",&mwpcE);
  //chain->SetBranchAddress("scintTimeToHit",&Time);
  //chain->SetBranchAddress("scintillatorEdep",&edep);
  //chain->SetBranchAddress("scintillatorEdepQuenched",&edepQ);
  //chain->SetBranchAddress("MWPCPos",&mwpc_pos);
  //chain->SetBranchAddress("ScintPos",&scint_pos);
  //chain->SetBranchAddress("primaryKE",&Eprim);


  

  //Trigger booleans
  bool EastScintTrigger, WestScintTrigger, EMWPCTrigger, WMWPCTrigger;
  Double_t MWPCThreshold=0.001;

  //Set random number generator
  TRandom3 *rand = new TRandom3(0);
  TRandom3 *rand2 = new TRandom3(0);
  
  //Get total number of events in TChain
  UInt_t nevents = chain->GetEntries();
  cout << "events = " << nevents << endl;
 
  //Start from random position in evt sequence
  UInt_t evtStart = rand->Rndm()*nevents;
  UInt_t evtTally = 0; //To keep track of the number of events 
  UInt_t evt = evtStart; //current event number
  vector < vector <Int_t> > gridPoint;
  

  //Read in events and determine evt type based on triggers
  while (evtTally<=BetaEvents) {
    if (evt>=nevents) evt=0; //Wrapping the events back to the beginning
    EastScintTrigger = WestScintTrigger = EMWPCTrigger = WMWPCTrigger = false; //Resetting triggers each event

    chain->GetEvent(evt);

    //Checking that the event occurs within the fiducial volume in the simulation to minimize
    // contamination from edge effects and interactions with detector walls
    Double_t fidCut = 45.;
    if (source=="Beta") fidCut=50.; //Since I cut at 50 when determining the number of beta events
    else if (sqrt(scint_pos.ScintPosE[0]*scint_pos.ScintPosE[0]+scint_pos.ScintPosE[1]+scint_pos.ScintPosE[1])*sqrt(0.6)*10.>fidCut
	     || sqrt(scint_pos.ScintPosW[0]*scint_pos.ScintPosW[0]+scint_pos.ScintPosW[1]+scint_pos.ScintPosW[1])*sqrt(0.6)*10.>fidCut) {evt++; continue;} //For source events, 
    // We don't want edge contamination, and I use a cut of 45 mm when selecting what sources are present in calibration runs.


    //Calculate event weight
    if (source=="Beta") AsymWeight = 1+(-0.12)*pol*sqrt(1-(1/((Eprim/511.+1.)*(Eprim/511.+1.))))*cos(primTheta);
    else AsymWeight = 1.;
    
    scint_pos_adj.ScintPosAdjE[0] = source=="Beta"?scint_pos.ScintPosE[0]*sqrt(0.6)*10.:rand2->Gaus(srcPos[0][0], fabs(srcPos[0][2]));
    scint_pos_adj.ScintPosAdjE[1] = source=="Beta"?scint_pos.ScintPosE[1]*sqrt(0.6)*10.:rand2->Gaus(srcPos[0][1], fabs(srcPos[0][2]));
    scint_pos_adj.ScintPosAdjW[0] = source=="Beta"?scint_pos.ScintPosW[0]*sqrt(0.6)*10.:rand2->Gaus(srcPos[1][0], fabs(srcPos[1][2]));//sqrt(0.6)*10.*scint_pos.ScintPosW[0]+displacementX;
    scint_pos_adj.ScintPosAdjW[1] = source=="Beta"?scint_pos.ScintPosW[1]*sqrt(0.6)*10.:rand2->Gaus(srcPos[1][1], fabs(srcPos[1][2]));//sqrt(0.6)*10.*scint_pos.ScintPosW[1]+displacementY;
    scint_pos_adj.ScintPosAdjE[2] = scint_pos.ScintPosE[2]*10.;
    scint_pos_adj.ScintPosAdjW[2] = scint_pos.ScintPosW[2]*10.;

          
    //retrieve point on grid for each side of detector [E/W][x/y]
    gridPoint = getGridPoint(scint_pos_adj.ScintPosAdjE[0],scint_pos_adj.ScintPosAdjE[1],scint_pos_adj.ScintPosAdjW[0],scint_pos_adj.ScintPosAdjW[1]);

    Int_t intEastBinX = gridPoint[0][0];
    Int_t intEastBinY = gridPoint[0][1];
    Int_t intWestBinX = gridPoint[1][0];
    Int_t intWestBinY = gridPoint[1][1];
      
    //MWPC triggers
    if (mwpcE.MWPCEnergyE>MWPCThreshold) EMWPCTrigger=true;
    if (mwpcE.MWPCEnergyW>MWPCThreshold) WMWPCTrigger=true;

    Double_t pmtEnergyLowerLimit = 1.; //To put a hard cut on the weight
    

    //East Side smeared PMT energies

    for (UInt_t p=0; p<4; p++) {
      if (pmtQuality[p]) { //Check to make sure PMT was functioning
	//cout << p << " " << positionMap[p][intEastBinX][intEastBinY] << endl;
	
	if (positionMap[p][intEastBinX][intEastBinY]>0.) {
	  pmt.etaEvis[p] = (1./(alpha[p]*g_d)) * rand->Poisson(g_d*rand2->Poisson(alpha[p]*positionMap[p][intEastBinX][intEastBinY]*edepQ.EdepQE));
	  pmt.Evis[p] = pmt.etaEvis[p]/positionMap[p][intEastBinX][intEastBinY];
	}
	else { //To avoid dividing by zero.. these events won't be used in analysis since they are outside the fiducial cut
	  pmt.etaEvis[p] = (1./(alpha[p]*g_d)) * rand->Poisson(g_d*rand2->Poisson(alpha[p]*edepQ.EdepQE));
	  pmt.Evis[p] = pmt.etaEvis[p];
	}

	pmt.nPE[p] = alpha[p]*pmt.etaEvis[p];
      }	
      // If PMT quality failed, set the evis and the weight to zero for this PMT
      else {
	pmt.Evis[p] = 0.;
	pmt.Evis[p] = 0.;
	pmt.nPE[p] = 0.;
      }
    }
      
    //Calculate the weighted energy on a side
    Double_t numer=0., denom=0.;
    for (UInt_t p=0;p<4;p++) {
      numer += pmt.nPE[p];
      denom += pmt.nPE[p]>0. ? positionMap[p][intEastBinX][intEastBinY]*alpha[p] : 0.;
    }

    //Now we apply the trigger probability
    Double_t totalEnE = numer>0. ? numer/denom : 0.;
    evis.EvisE = totalEnE;
    Double_t triggProb = triggerProbability(triggerFunc[0],totalEnE);
      
    //Set East Scint Trigger to true if event passes triggProb
    if (rand->Rndm(0)<triggProb) { // && evis.EvisE>0.) {
      EastScintTrigger=true;
    }
      
      
    //West Side
    for (UInt_t p=4; p<8; p++) {
      if (pmtQuality[p]) { //Check to make sure PMT was functioning
	
	if (positionMap[p][intWestBinX][intWestBinY]>0.) {
	  pmt.etaEvis[p] = (1./(alpha[p]*g_d)) * rand->Poisson(g_d*rand2->Poisson(alpha[p]*positionMap[p][intWestBinX][intWestBinY]*edepQ.EdepQW));
	  pmt.Evis[p] = pmt.etaEvis[p]/positionMap[p][intWestBinX][intWestBinY];
	}
	else { //To avoid dividing by zero.. these events won't be used in analysis since they are outside the fiducial cut
	  pmt.etaEvis[p] = (1./(alpha[p]*g_d)) * rand->Poisson(g_d*rand2->Poisson(alpha[p]*edepQ.EdepQW));
	  pmt.Evis[p] = pmt.etaEvis[p];
	}

	pmt.nPE[p] = alpha[p]*pmt.etaEvis[p];
      }	
      // If PMT quality failed, set the evis and the weight to zero for this PMT
      else {
	pmt.Evis[p] = 0.;
	pmt.Evis[p] = 0.;
	pmt.nPE[p] = 0.;
      }
    }
      
    //Calculate the total weighted energy
    numer=denom=0.;
    for (UInt_t p=4;p<8;p++) {
      numer+=pmt.nPE[p];
      denom += pmt.nPE[p]>0. ? positionMap[p][intWestBinX][intWestBinY]*alpha[p] : 0.;
    }
    //Now we apply the trigger probability
    Double_t totalEnW = numer>0. ? numer/denom : 0.;
    evis.EvisW = totalEnW;
    triggProb = triggerProbability(triggerFunc[1],totalEnW);
    //Fill histograms if event passes trigger function
    if (rand->Rndm(0)<triggProb) { // && evis.EvisW>0.) {
      WestScintTrigger = true;      
    }
            
    //Fill proper total event histogram based on event type
    PID=-1;
    //Type 0 East
    if (EastScintTrigger && EMWPCTrigger && !WestScintTrigger && !WMWPCTrigger) {
      PID=1;
      type=0;
      side=0;
      //finalEn[0]->Fill(evis.EvisE);
      //cout << "Type 0 East E = " << totalEnE << endl;
    }
    //Type 0 West
    else if (WestScintTrigger && WMWPCTrigger && !EastScintTrigger && !EMWPCTrigger) {
      PID=1;
      type=0;
      side=1;
      //finalEn[1]->Fill(totalEnW);
    }
    //Type 1 
    else if (EastScintTrigger && EMWPCTrigger && WestScintTrigger && WMWPCTrigger) {
      PID=1;
      type=1;
      //East
      if (Time.timeE<Time.timeW) {
	//finalEn[2]->Fill(totalEnE);
	side=0;
      }
      //West
      else if (Time.timeE>Time.timeW) {
	//finalEn[3]->Fill(totalEnW);
	side=1;
      }
    }
    //Type 2/3 East
    else if (EastScintTrigger && EMWPCTrigger && !WestScintTrigger && WMWPCTrigger) {
      PID=1;
      type=2;
      side=0;
      //finalEn[4]->Fill(totalEnE);
      //cout << "Type 2/3 East E = " << totalEnE << endl;
    }
    //Type 2/3 West
    else if (!EastScintTrigger && EMWPCTrigger && WestScintTrigger && WMWPCTrigger) {
      PID=1;
      type=2;
      side=1;
      //finalEn[5]->Fill(totalEnW);
      //cout << "Type 2/3 East W = " << totalEnW << endl;
    }   
    //Gamma events and missed events (Type 4)
    else {
      if (!WMWPCTrigger && !EMWPCTrigger) {
	if (EastScintTrigger && !WestScintTrigger) {
	  PID=0;
	  type=0;
	  side=0;
	}
	else if (!EastScintTrigger && WestScintTrigger) {
	  PID=0;
	  type=0;
	  side=1;
	}
	else if (EastScintTrigger && WestScintTrigger) {
	  PID=0;
	  type=1;
	  if (Time.timeE<Time.timeW) {
	    side=0;
	  }
	  else {
	    side=1;
	  }
	}
	else {
	  PID=6;
	  type=4;
	  side=2;
	}
      }
      else {
	PID=1;
	type=4;
	side = (WMWPCTrigger && EMWPCTrigger) ? 2 : (WMWPCTrigger ? 1 : 0);
      }
    }
    
    //Calculate Erecon
    Int_t typeIndex = (type==0 || type==4) ? 0:(type==1 ? 1:2); //for retrieving the parameters from EQ2Etrue
    if (side==0) {
      Double_t totalEvis = type==1 ? (evis.EvisE+evis.EvisW):evis.EvisE;
      if (totalEvis>0.) {
	Erecon = EQ2Etrue[0][typeIndex][0]+EQ2Etrue[0][typeIndex][1]*totalEvis+EQ2Etrue[0][typeIndex][2]/(totalEvis+EQ2Etrue[0][typeIndex][3])+EQ2Etrue[0][typeIndex][4]/((totalEvis+EQ2Etrue[0][typeIndex][5])*(totalEvis+EQ2Etrue[0][typeIndex][5]));
	if (type==0) finalEn[0]->Fill(Erecon); 
	else if (type==1) finalEn[2]->Fill(Erecon); 
	else if(type==2 ||type==3) finalEn[4]->Fill(Erecon);
      }
      else Erecon=-1.;
    }
    if (side==1) {
      Double_t totalEvis = type==1 ? (evis.EvisE+evis.EvisW):evis.EvisW;
      if (totalEvis>0.) {
	Erecon = EQ2Etrue[1][typeIndex][0]+EQ2Etrue[1][typeIndex][1]*totalEvis+EQ2Etrue[1][typeIndex][2]/(totalEvis+EQ2Etrue[1][typeIndex][3])+EQ2Etrue[1][typeIndex][4]/((totalEvis+EQ2Etrue[1][typeIndex][5])*(totalEvis+EQ2Etrue[1][typeIndex][5]));
	
	if (type==0) finalEn[1]->Fill(Erecon); 
	else if (type==1) finalEn[3]->Fill(Erecon); 
	else if(type==2 ||type==3) finalEn[5]->Fill(Erecon);
      }
      else Erecon=-1.;
    }    
  
    // Increment the event tally if the event was PID = 1 (electron) and the event was inside the fiducial radius used to determine num of events in data file
    if (PID==1 && sqrt(scint_pos.ScintPosE[0]*scint_pos.ScintPosE[0]+scint_pos.ScintPosE[1]+scint_pos.ScintPosE[1])*sqrt(0.6)*10.<fidCut
	&& sqrt(scint_pos.ScintPosW[0]*scint_pos.ScintPosW[0]+scint_pos.ScintPosW[1]+scint_pos.ScintPosW[1])*sqrt(0.6)*10.<fidCut) evtTally++;
    evt++;
    if (PID>=0) tree->Fill();
    //cout << evtTally << endl;
    //if (evtTally%1000==0) {cout << evtTally << endl;}//cout << "filled event " << evt << endl;
  }
  cout << endl;
  delete chain;
  delete rand;
  delete rand2;
  outfile->Write();
  outfile->Close();
  
}

int main(int argc, char *argv[]) {
  string source = string(argv[2]);
  revCalSimulation(atoi(argv[1]),source);

  //tests
  /*UInt_t XePeriod = getXeRunPeriod(atoi(argv[1]));
  vector < vector <Double_t> > triggerFunc = getTriggerFunctionParams(XePeriod,7);
  Double_t triggProbE = triggerProbability(triggerFunc[0],25.);
  Double_t triggProbW = triggerProbability(triggerFunc[1],25.);
  cout << triggProbE << " " << triggProbW << endl;*/


}
  
