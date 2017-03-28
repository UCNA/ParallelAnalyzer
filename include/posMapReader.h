#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>

using namespace std;

// Position bins
const int nPMT = 8;
const int nPosBinsX = 23; //43;
const int nPosBinsY = 23; //43;  
const double xBinWidth = 5.;//2.5;
const double yBinWidth = 5.;//2.5;
double xBinLower[nPosBinsX];
double xBinUpper[nPosBinsX];
double xBinCenter[nPosBinsX];
double yBinLower[nPosBinsY];
double yBinUpper[nPosBinsY];
double yBinCenter[nPosBinsY];
int intXBinCenter[nPosBinsX];
int intYBinCenter[nPosBinsY];

double positionMap[nPMT][nPosBinsX][nPosBinsY];

void GetPositionMap(int XePeriod) {
  for (int k=0; k<nPosBinsX; k++) {
    xBinLower[k]     = -(double)nPosBinsX*xBinWidth/2. + ((double) k)*xBinWidth;
    xBinUpper[k]     = -(double)nPosBinsX*xBinWidth/2. + ((double) k)*xBinWidth + xBinWidth;
    xBinCenter[k]    = (xBinLower[k] + xBinUpper[k])/2.;
    intXBinCenter[k] = (int) xBinCenter[k];
    //cout << xBinLower[k] << " " << intXBinCenter[k] << " " << xBinUpper[k] << endl;
  }
  
  for (int k=0; k<nPosBinsY; k++) {
    yBinLower[k]     = -(double)nPosBinsY*yBinWidth/2. + ((double) k)*yBinWidth;
    yBinUpper[k]     = -(double)nPosBinsY*yBinWidth/2. + ((double) k)*yBinWidth + yBinWidth;
    yBinCenter[k]    = (yBinLower[k] + yBinUpper[k])/2.;
    intYBinCenter[k] = (int) yBinCenter[k];
    //cout << yBinLower[k] << " " << intYBinCenter[k] << " " << yBinUpper[k] << endl;
  }
  
  // Determine position map to use
  char tempFileXePositionMap[500];
  sprintf(tempFileXePositionMap, "%s/position_map_%i_RC_123_%0.1fmm.dat",getenv("POSITION_MAPS"),XePeriod, xBinWidth);
  
  cout << "... Reading: " << tempFileXePositionMap << endl;

  // Read position map
  double x, y;
  
  ifstream fileXePositionMap(tempFileXePositionMap);
  for (int i=0; i<nPosBinsX; i++) {
    for (int j=0; j<nPosBinsY; j++) {
      fileXePositionMap >> x >> y
                        >> positionMap[0][i][j]
			>> positionMap[1][i][j]
			>> positionMap[2][i][j]
			>> positionMap[3][i][j]
			>> positionMap[4][i][j]
			>> positionMap[5][i][j]
			>> positionMap[6][i][j]
			>> positionMap[7][i][j];
    }
  }
}

vector < vector < int > > getGridPoint(double xE, double yE, double xW, double yW) {
  // Determine (x,y) bin
  vector <vector <int> > gridPoint; 
  gridPoint.resize(2,vector <int> (2,-1));

  for (int m=0; m<nPosBinsX; m++) {
    if ( (xE >= xBinLower[m]) && (xE < xBinUpper[m]) )  gridPoint[0][0] = m;
    if ( (xW >= xBinLower[m]) && (xW < xBinUpper[m]) )  gridPoint[1][0] = m;
  }
  
  for (int m=0; m<nPosBinsY; m++) {
    if ( (yE >= yBinLower[m]) && (yE < yBinUpper[m]) ) gridPoint[0][1] = m;
    if ( (yW >= yBinLower[m]) && (yW < yBinUpper[m]) ) gridPoint[1][1] = m;
  }
  return gridPoint;
}
