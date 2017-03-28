#ifndef EVTRATEHANDLER_HH
#define EVTRATEHANDLER_HH

#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TH1.h>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <string>



class EvtRateHandler { // base clase for reading in runs and calculating rates and errors
public:
  EvtRateHandler(std::vector<int> rn, bool fg, std::string anaCh, double enBinWidth=10., double fidCut=100., bool ukdata=true, bool unblind=false);
  virtual ~EvtRateHandler();

  double getFiducialCut() { return fiducialCut; };
  bool isUKdata() { return UKdata; }; 

  int polarization(int run);       // determines the polarization. Assigns value to pol and returns polarization
                                   // This comes from the database. Flipper on -> 1, Flipper Off -> -1

  double returnRunLength(int side) { return side==0 ? totalRunLengthE : totalRunLengthW; }       // Return the length of the run (s)

  virtual void CalcRates();       // Calculate the rates and use the reference spectra to fill in errors if necessary

  
  
  std::vector <double> getRateVectors(int side) { 
    return side==0 ? rateEvec : ( side==1 ? rateWvec : throw "BAD SIDE GIVEN TO getRateVectors"); 
  } ;

  std::vector <double> getRateErrors(int side) {
    return side==0 ? rateEerr : ( side==1 ? rateWerr : throw "BAD SIDE GIVEN TO getRateErrors"); 
  } ;

protected:
  virtual void dataReader();       //Read in data and fill histograms
  void loadReferenceSpectra(); //Loads the reference rates from the proper files
  double referenceError(int side, int bin);       // goes to the proper reference spectra and returns the error

  std::vector <int> runs;       // runs
  bool FG;        // holds whether the run is a Foreground or a background run
  std::string analysisChoice; 
  double fiducialCut;       //definition of a fiducial volume
  bool UKdata;       //This is true if the tree format is UK style, False if it's official analyzer style
  bool unblinded;       //Holds whether the result should be unblinded
  int pol;       // polarization of the foreground run if polarization has been run;
  unsigned int numEnergyBins;
  std::vector < std::vector <double> > runLength;       // E/W for each run in run vector
  double totalRunLengthE;       // Holds the sum of the run lengths of the applicable runs
  double totalRunLengthW;
  double totalCountsE;       // Holds the total number of events in the runs of interest
  double totalCountsW;
  
  std::vector <double>  UCNMonIntegral;
  
  TH1D *hisCounts[2];       // histogram for counts in an energy bin
  
  std::vector <double> rateWvec;       //Stores the rate for each bin
  std::vector <double> rateEvec;
  
  std::vector <double> rateEerr;       //Stores the statistical error for each Bin
  std::vector <double> rateWerr;

  std::vector <double> refSpectraE;       //Stores the reference runs event rate per bin
  std::vector <double> refSpectraW;       //Stores the reference runs event rate per bin
  double totalRefTime;       //UNBLINDED total reference time
  double totalRefTimeE;       // BLINDED time
  double totalRefTimeW;       // Blinded time
  double totalRefCountsE;
  double totalRefCountsW;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
 
class SimEvtRateHandler: public EvtRateHandler {
public:
  SimEvtRateHandler(std::vector <int> run, bool fg, std::string anaCh, double enBinWidth=10., double fidCut=100., bool unblind=false): EvtRateHandler(run, fg, anaCh, enBinWidth, fidCut, true, unblind) {}

protected:
  void dataReader();       //Different set of variables for reverse calibrated simulated data
  
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class BGSubtractedRate {
public:
  BGSubtractedRate(std::vector <int> fgruns, std::vector <int> bgruns, std::string anaCh, double enBin, double fidCut, bool ukdata=true, bool sim=false, bool unblind=false);
  ~BGSubtractedRate() {}

  std::vector<int> getCurrentForegroundRuns() { return FGruns; };       // vector of foreground runs (possibly more than one when runs spilled over
  std::vector<int> getCurrentBackgroundRuns() { return BGruns; };       // vector of background runs (possibly more than one when runs spilled over
  
  std::vector<double> returnRunLengths(bool beta=true);       // for BG do beta=false
  void calcBGSubtRates();       //Loads the rates and calculates BG subtr rates and errors
  std::vector<double> returnBGSubtRate(int side);       //returns a vector holding all the BG subtracted rates.
  std::vector<double> returnBGSubtRateError(int side);       //Returns background subtracted rate statistical error

private:

  std::vector <int> FGruns;      // vector of foreground runs
  std::vector <int> BGruns;      // vector of background runs
  std::string analysisChoice;       // string holding the analysis choice
  double EnergyBinWidth;       //Width of energy bins used
  double fiducialCut;
  bool UKdata;       //Whether the replay was done using UK code or MPM code
  bool Simulation;       //Whether the rate is from simulation or not, in which case there is no background run
  bool UNBLIND;       //SHOULD BE FALSE UNTIL UNBLINDING

  std::vector <double> BetaRateE;       //Save the rates here for beta run
  std::vector <double> BetaRateErrorE;       //Save the error here for beta run
  std::vector <double> BGRateE;         // Save the background rates here
  std::vector <double> BGRateErrorE;         // Save the background rates here
  std::vector <double> FinalRateE;       //This is the difference in the rates
  std::vector <double> FinalRateErrorE;       //This is the statistical error in the difference in the rates
  

  std::vector <double> BetaRateW;       //Save the rates here for beta run
  std::vector <double> BetaRateErrorW;       //Save the error here for beta run
  std::vector <double> BGRateW;         // Save the background rates here
  std::vector <double> BGRateErrorW;         // Save the background rates here
  std::vector <double> FinalRateW;       //This is the difference in the rates
  std::vector <double> FinalRateErrorW;       //This is the statistical error in the difference in the rates

  std::vector <double> runLengthBG;       // E/W
  std::vector <double> runLengthBeta;
  
  void LoadRatesByBin();       // Load rates and save them to vectors
  void CalcFinalRate();
  
  
};
#endif
