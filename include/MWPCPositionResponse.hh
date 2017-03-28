#ifndef MWPCPOSRESPONSE_HH
#define MWPCPOSRESPONSE_HH

#include <vector>



class MWPCCathodeHandler {

public:

  MWPCCathodeHandler(); // In case you just need access to wire positions
  
  MWPCCathodeHandler(double *ex,double *ey,double *wx,double *wy,
		     double *pedex,double *pedey,double *pedwx,double *pedwy);

  MWPCCathodeHandler(double *ex,double *ey,double *wx,double *wy); // If they are already pedestal
                                                                   // subtracted
  ~MWPCCathodeHandler() {};

  //void setCathodeSignals(double *ex,double *ey,double *wx,double *wy,
  //			 double *pedex,double *pedey,
  //			 double *pedwx,double *pedwy);
    
  void setCathodeThreshold(double *ex,double *ey,double *wx,double *wy);
  void setClippingThreshold(double *ex,double *ey,double *wx,double *wy); 
  void setSigma(double s) { _sigma = s; };
  void loadGainFactors(int run); // Loads the gain factors if you want gain
                                 // gain corrected cathode signals
  void purgeGainFactors(); //Erases previously set gain factors

  void loadCathodeModelParams(int run); // When doing simulation, you can load the model params for clipping and triggering


  void PrintSignals();
  void findAllPositions(bool gaus=true, bool gausAllEvents=false); // You only need to call this to return all position.
                                         // If gaus is false, it just uses the weighted average for all events
                                         // as was originally done. If gaus is true and gausAllEvents is true,
                                         // it will look for three wires for all events, not just clipped events
  
  std::vector<double> getPosEX() { return posEX; }
  std::vector<double> getPosEY() { return posEY; }
  std::vector<double> getPosWX() { return posWX; }
  std::vector<double> getPosWY() { return posWY; }

  int getnClippedEX() { return clippedEX.size(); }
  int getnClippedEY() { return clippedEY.size(); }
  int getnClippedWX() { return clippedWX.size(); }
  int getnClippedWY() { return clippedWY.size(); }

  int getMultEX() { return signalEX.size(); }
  int getMultEY() { return signalEY.size(); }
  int getMultWX() { return signalWX.size(); }
  int getMultWY() { return signalWY.size(); }
  
  int getMaxWireEX() { return getMaxWire(signalEX,pedSubtrEX); }
  int getMaxWireEY() { return getMaxWire(signalEY,pedSubtrEY); }
  int getMaxWireWX() { return getMaxWire(signalWX,pedSubtrWX); }
  int getMaxWireWY() { return getMaxWire(signalWX,pedSubtrWY); }

  double getMaxSignalEX() { return getMaxSignal(pedSubtrEX) ; } // These are pedestal subtracted
  double getMaxSignalEY() { return getMaxSignal(pedSubtrEY) ; }
  double getMaxSignalWX() { return getMaxSignal(pedSubtrWX) ; }
  double getMaxSignalWY() { return getMaxSignal(pedSubtrWY) ; }

  double getCathSumEX() { return getCathSum(pedSubtrEX,16) ; } // These are pedestal subtracted
  double getCathSumEY() { return getCathSum(pedSubtrEY,16) ; }
  double getCathSumWX() { return getCathSum(pedSubtrWX,16) ; }
  double getCathSumWY() { return getCathSum(pedSubtrWY,16) ; }
  
  double getWirePosEX(int i) { return wireposEX[i]; }
  double getWirePosEY(int i) { return wireposEY[i]; }
  double getWirePosWX(int i) { return wireposWX[i]; }
  double getWirePosWY(int i) { return wireposWY[i]; }
  




private:
  
  double cathodeThresholdEX[16], cathodeThresholdEY[16], cathodeThresholdWX[16], cathodeThresholdWY[16]; // Value where trigger occurs
  double clipThresholdEX[16], clipThresholdEY[16], clipThresholdWX[16], clipThresholdWY[16]; // value where clipping occurs
  double cathEX[16], cathEY[16], cathWX[16], cathWY[16];
  double pedEX[16], pedEY[16], pedWX[16], pedWY[16];
  double pedSubtrEX[16],pedSubtrEY[16],pedSubtrWX[16],pedSubtrWY[16];
  std::vector<double> posEX, posEY, posWX, posWY; // Final position of the event (mean,width,height)
  //int nClippedEX, nClippedEY, nClippedWX, nClippedWY;
  std::vector <int> clippedEX, clippedEY, clippedWX, clippedWY; // Vector of clipped wires
  std::vector <int> signalEX, signalEY, signalWX, signalWY; // Holds wires with signal above threshold
  std::vector<double> gainEX, gainEY, gainWX, gainWY;

  bool _bGaus, _bGausAllEvents;
  double _sigma; // This is the characteristic sigma of an event
  
  void doPedestalSubtraction(); // Subtracts pedestals and fills pedSubtr** 
  void doThresholdCheck(); // Sets all of the wires that were above threshold in signal**
  std::vector <int>  doClipping(std::vector<int> wires, double *sig, double *thresh); //Checks for clipped wires
  std::vector<double> fitCathResponse(std::vector <int> wires, std::vector<int> clip, double *sig, const double *pos); // // returns gaussian mean and width and height of a signal in a three element vector
  int getMaxWire(std::vector <int> wires, double *sig); // returns max wire
  double getMaxSignal(double *sig); // returns max signal
  std::vector <int> getNonClippedSorted( const std::vector<int>& wires, const std::vector<int>& clipWires, double *sig);
  std::vector <int> getNonClippedSequential( const std::vector<int>& wires, const std::vector<int>& clipWires, double *sig);
  double getCathSum(double *sig, int length);

  std::vector<double> fitGaus(std::vector<double> x, std::vector<double> y); // returns gaussian mean and width and height of data arrays in a 3 element vector
  std::vector<double> fitParabola(std::vector<double> x, std::vector<double> y); // returns gaussian mean and width and height of data arrays in a 3 element vector
  std::vector<double> fitGaus2Points(std::vector<double> x, std::vector<double> y); // If there are only two points, this is used

  
  bool boolPedSubtr{false}; // Whether or not the pedestal subtraction has been calculated
  bool boolThresholdCheck{false}; // Whether or not the threshold has been checked has been calculated
 
  const double wireposEX[16] {76.20, 66.04, 55.88, 45.72, 35.56, 25.40, 15.24, 5.08, 
      -5.08, -15.24, -25.40, -35.56, -45.72, -55.88, -66.04, -76.20};
  const double wireposEY[16] {76.20, 66.04, 55.88, 45.72, 35.56, 25.40, 15.24, 5.08, 
      -5.08, -15.24, -25.40, -35.56, -45.72, -55.88, -66.04, -76.20};
  const double wireposWY[16] {76.20, 66.04, 55.88, 45.72, 35.56, 25.40, 15.24, 5.08, 
      -5.08, -15.24, -25.40, -35.56, -45.72, -55.88, -66.04, -76.20};
  const double wireposWX[16] {-76.20, -66.04, -55.88, -45.72, -35.56, -25.40, -15.24, -5.08,
      5.08, 15.24, 25.40, 35.56, 45.72, 55.88, 66.04, 76.20};
  
};


#endif


  
