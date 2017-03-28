#include "positionMapHandler.hh"

#include <string>
#include <fstream>
#include <cstdlib>
#include "MBUtils.hh"

#include <TMath.h>



PositionMap::PositionMap(Double_t bin_width, Double_t r_max) : binWidth(bin_width), rmax(r_max), bUseRC(useRC) {
  
  Int_t nBinsHold = (int)(110./binWidth);
  nBinsXY = ((int)nBinsHold%2)==0 ? (int)nBinsHold+1 : (int)nBinsHold;
  nBinsTotal = nBinsXY*nBinsXY;

  xyBinLower.resize(nBinsXY,0.);
  xyBinUpper.resize(nBinsXY,0.);
  xyBinCenter.resize(nBinsXY,0.);
  
  posMap.resize(nBinsXY, std::vector < std::vector < Double_t > > (nBinsXY, std::vector <Double_t> (8,0.)));
  
  setBinValues();
  
  funcs.push_back(new TF1("func0", "-0.5*(1.-x)*(1.-x)*x", 0., 1.));
  funcs.push_back(new TF1("func1", "(1.+x-3./2.*x*x)*(1.-x)", 0., 1.));
  funcs.push_back(new TF1("func2", "(1+(1-x)-3./2.*(1.-x)*(1.-x))*x", 0., 1.));
  funcs.push_back(new TF1("func3", "-0.5*(1.-x)*x*x", 0., 1.));
};

PositionMap::~PositionMap() {
  for (UInt_t i=0; i<funcs.size(); i++) {
    if (funcs[i]) delete funcs[i];
  }
};

void PositionMap::setBinValues() {
  for (Int_t k=0; k<nBinsXY; k++) {
    xyBinLower[k]     = -(Double_t)nBinsXY*binWidth/2. + ((Double_t) k)*binWidth;
    xyBinUpper[k]     = -(Double_t)nBinsXY*binWidth/2. + ((Double_t) k)*binWidth + binWidth;
    xyBinCenter[k]    = (xyBinLower[k] + xyBinUpper[k])/2.;
    
    //cout << xBinLower[k] << " " << intXBinCenter[k] << " " << xBinUpper[k] << endl;
  }
};

void PositionMap::readPositionMap(Int_t XeRunPeriod) {
  XePeriod = XeRunPeriod;
  std::string file = "";
  if ( bUseRC ) file = std::string(getenv("POSITION_MAPS")) + "/position_map_" + itos(XePeriod) + "_" + "RC_123_" + ftos(binWidth) +  "mm.dat";
  else file = std::string(getenv("POSITION_MAPS")) + "/position_map_" + itos(XePeriod) + "_" + ftos(binWidth) +  "mm.dat";
  std::cout << "Reading position map from " << file << std::endl;
  std::ifstream infile(file.c_str());
  
  Int_t binValx, binValy;
  std::vector <Double_t> p(8,0.);
  
  while (infile >> binValx >> binValy >> p[0] >> p[1] >> p[2] >> p[3] >> p[4] >> p[5] >> p[6] >> p[7])
    {
      setPositionMapPoint(getBinNumber(binValx), getBinNumber(binValy), p);
    }

  std::cout << "Read in Position Map from Xe Period " << XePeriod << std::endl;
};

Int_t PositionMap::getBinNumber(Double_t pos) {
  for (Int_t m=0; m<nBinsXY; m++) {
    if ( (pos >= xyBinLower[m]) && (pos < xyBinUpper[m]) )  return m;
  }
  return -1;
};

void PositionMap::setPositionMapPoint(Int_t xBin, Int_t yBin, std::vector <Double_t> params) {
  posMap[xBin][yBin] = params;
};




std::vector <Double_t> PositionMap::getInterpolatedEta(Double_t xE, Double_t yE, Double_t xW, Double_t yW) {

  std::vector <Double_t> val(8, 0.);

  if ( (xE*xE+yE*yE) >= 65.*65. || (xW*xW+yW*yW) >= 65.*65. ) { return std::vector<Double_t>(8,1.); }
 
  
  if ( (xE*xE+yE*yE) >= rmax*rmax) {

    //First we go to a point at the same angle on the position map, only with radius = 50. mm
    Double_t tanTheta = TMath::Abs(xE)>1.e-5 ? TMath::Abs(yE)/TMath::Abs(xE) : 1.;
    
    Double_t x = TMath::Abs(xE)>1.e-5 ? (xE>0. ? rmax/TMath::Sqrt(tanTheta*tanTheta+1.) : -rmax/TMath::Sqrt(tanTheta*tanTheta+1.)) : 0.;
    Double_t y = TMath::Abs(xE)>1.e-5 ? (yE>0. ? tanTheta*rmax/TMath::Sqrt(tanTheta*tanTheta+1.) : -tanTheta*rmax/TMath::Sqrt(tanTheta*tanTheta+1.)) : (yE>0. ? rmax : -rmax);

    //std::cout << "x = " << x << "  y = " << y << std::endl;
    
    Int_t xBin = x > getBinCenter(getBinNumber(x)) ? getBinNumber(x) : getBinNumber(x)-1;
    Int_t yBin = y > getBinCenter(getBinNumber(y)) ? getBinNumber(y) : getBinNumber(y)-1;

    Double_t xp = (x-getBinCenter(xBin))/binWidth;
    Double_t yp = (y-getBinCenter(yBin))/binWidth; 

    // if we lie directly on the center of a bin, use only one bin to the left/down if val is neg. in interpolation (Shift bins)
    if (xp==1. && x<0.) {xp=0.; xBin+=1;}
    if (yp==1. && y<0.) {yp=0.; yBin+=1;}

    //std::cout << "xBin = " << xBin << "  yBin = " << yBin << std::endl;
    //std::cout << "xp = " << xp << "  yp = " << yp << std::endl;

    for (int p=0; p<4; p++) {
      Int_t xx=0;
      for (int xb=xBin-1; xb<=xBin+2; xb++) {
	Int_t yy=0;
	for (int yb=yBin-1; yb<=yBin+2; yb++) {
	  val[p]+=(funcs[xx]->Eval(xp))*(funcs[yy]->Eval(yp))*posMap[xb][yb][p];
	  yy++;
	}
	xx++;
      }
    }
  }
  else {
    Int_t xBin = xE > getBinCenter(getBinNumber(xE)) ? getBinNumber(xE) : getBinNumber(xE)-1;
    Int_t yBin = yE > getBinCenter(getBinNumber(yE)) ? getBinNumber(yE) : getBinNumber(yE)-1;
    
    Double_t xp = (xE-getBinCenter(xBin))/binWidth;
    Double_t yp = (yE-getBinCenter(yBin))/binWidth; 
    
    for (Int_t p=0; p<4; p++) {
      Int_t xx=0;
      for (int xb=xBin-1; xb<=xBin+2; xb++) {
	Int_t yy=0;
	for (int yb=yBin-1; yb<=yBin+2; yb++) {
	  val[p]+=(funcs[xx]->Eval(xp))*(funcs[yy]->Eval(yp))*posMap[xb][yb][p];
	  yy++;
	}
	xx++;
      }
      //std::cout << p << " " << val[p] << std::endl;
    }
  }
  
  

  if ( (xW*xW+yW*yW) >= rmax*rmax) {

    Double_t tanTheta = TMath::Abs(xW)>1.e-5 ? TMath::Abs(yW)/TMath::Abs(xW) : 1.;
    Double_t x = TMath::Abs(xW)>1.e-5 ? (xW>0. ? rmax/TMath::Sqrt(tanTheta*tanTheta+1.) : -rmax/TMath::Sqrt(tanTheta*tanTheta+1.)) : 0.;
    Double_t y = TMath::Abs(xW)>1.e-5 ? (yW>0. ? tanTheta*rmax/TMath::Sqrt(tanTheta*tanTheta+1.) : -tanTheta*rmax/TMath::Sqrt(tanTheta*tanTheta+1.)) : (yW>0. ? rmax : -rmax);
    Int_t xBin = x > getBinCenter(getBinNumber(x)) ? getBinNumber(x) : getBinNumber(x)-1;
    Int_t yBin = y > getBinCenter(getBinNumber(y)) ? getBinNumber(y) : getBinNumber(y)-1;

    Double_t xp = (x-getBinCenter(xBin))/binWidth;
    Double_t yp = (y-getBinCenter(yBin))/binWidth; 

    if (xp==1. && x<0.) {xp=0.; xBin+=1;}
    if (yp==1. && y<0.) {yp=0.; yBin+=1;}

    for (int p=4; p<8; p++) {
      Int_t xx=0;
      for (int xb=xBin-1; xb<=xBin+2; xb++) {
	Int_t yy=0;
	for (int yb=yBin-1; yb<=yBin+2; yb++) {
	  val[p]+=(funcs[xx]->Eval(xp))*(funcs[yy]->Eval(yp))*posMap[xb][yb][p];
	  yy++;
	}
	xx++;
      }
    }
  }
  
  else {
    
    Int_t xBin = xW > getBinCenter(getBinNumber(xW)) ? getBinNumber(xW) : getBinNumber(xW)-1;
    Int_t yBin = yW > getBinCenter(getBinNumber(yW)) ? getBinNumber(yW) : getBinNumber(yW)-1;
    
    Double_t xp = (xW-getBinCenter(xBin))/binWidth;
    Double_t yp = (yW-getBinCenter(yBin))/binWidth; 
    
    for (Int_t p=4; p<8; p++) {
      Int_t xx=0;
      for (int xb=xBin-1; xb<=xBin+2; xb++) {
	Int_t yy=0;
	for (int yb=yBin-1; yb<=yBin+2; yb++) {
	  val[p]+=(funcs[xx]->Eval(xp))*(funcs[yy]->Eval(yp))*posMap[xb][yb][p];
	  yy++;
	}
	xx++;
      }
      // std::cout << p << " " << val[p] << std::endl;
    }
  }
  return val;
  
};



//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

MWPCPositionMap::MWPCPositionMap(Double_t bin_width, Double_t r_max) : binWidth(bin_width), rmax(r_max) {
  
  Int_t nBinsHold = (int)(110./binWidth);
  nBinsXY = ((int)nBinsHold%2)==0 ? (int)nBinsHold+1 : (int)nBinsHold;
  nBinsTotal = nBinsXY*nBinsXY;

  xyBinLower.resize(nBinsXY,0.);
  xyBinUpper.resize(nBinsXY,0.);
  xyBinCenter.resize(nBinsXY,0.);
  
  posMap.resize(nBinsXY, std::vector < std::vector < Double_t > > (nBinsXY, std::vector <Double_t> (2,0.)));
  
  setBinValues();
  
  funcs.push_back(new TF1("func0", "-0.5*(1.-x)*(1.-x)*x", 0., 1.));
  funcs.push_back(new TF1("func1", "(1.+x-3./2.*x*x)*(1.-x)", 0., 1.));
  funcs.push_back(new TF1("func2", "(1+(1-x)-3./2.*(1.-x)*(1.-x))*x", 0., 1.));
  funcs.push_back(new TF1("func3", "-0.5*(1.-x)*x*x", 0., 1.));
};

MWPCPositionMap::~MWPCPositionMap() {
  for (UInt_t i=0; i<funcs.size(); i++) {
    if (funcs[i]) delete funcs[i];
  }
};

void MWPCPositionMap::setBinValues() {
  for (Int_t k=0; k<nBinsXY; k++) {
    xyBinLower[k]     = -(Double_t)nBinsXY*binWidth/2. + ((Double_t) k)*binWidth;
    xyBinUpper[k]     = -(Double_t)nBinsXY*binWidth/2. + ((Double_t) k)*binWidth + binWidth;
    xyBinCenter[k]    = (xyBinLower[k] + xyBinUpper[k])/2.;
    
    //cout << xBinLower[k] << " " << intXBinCenter[k] << " " << xBinUpper[k] << endl;
  }
};

void MWPCPositionMap::readMWPCPositionMap(Int_t XeRunPeriod, Double_t elow, Double_t ehigh) {

  XePeriod = XeRunPeriod;
  
  std::string file = ( std::string(getenv("MWPC_CALIBRATION")) + "/position_maps/MWPC_position_map_" + itos(XePeriod) + "_" + ftos(binWidth) +  "mm_"
		       + ftos( elow )+"-"+ftos( ehigh )+".dat" );
  std::cout << "Reading position map from " << file << std::endl;
  std::ifstream infile(file.c_str());
  
  Int_t binValx, binValy;
  std::vector <Double_t> p(2,0.);
  
  while ( infile >> binValx >> binValy >> p[0] >> p[1] )
    {
      setMWPCPositionMapPoint(getBinNumber(binValx), getBinNumber(binValy), p);
    }

  std::cout << "Read in Position Map from Xe Period " << XePeriod << std::endl;
};

Int_t MWPCPositionMap::getBinNumber(Double_t pos) {
  for (Int_t m=0; m<nBinsXY; m++) {
    if ( (pos >= xyBinLower[m]) && (pos < xyBinUpper[m]) )  return m;
  }
  return -1;
};

void MWPCPositionMap::setMWPCPositionMapPoint(Int_t xBin, Int_t yBin, std::vector <Double_t> params) {
  posMap[xBin][yBin] = params;
};




std::vector <Double_t> MWPCPositionMap::getInterpolatedEta(Double_t xE, Double_t yE, Double_t xW, Double_t yW) {

  std::vector <Double_t> val(2, 0.);

  if ( (xE*xE+yE*yE) >= 65.*65. || (xW*xW+yW*yW) >= 65.*65. ) { return std::vector<Double_t>(2,1.); }
  
  if (xE*xE+yE*yE >= rmax*rmax) {

    //First we go to a point at the same angle on the position map, only with radius = 50. mm
    Double_t tanTheta = TMath::Abs(xE)>1.e-5 ? TMath::Abs(yE)/TMath::Abs(xE) : 1.;
    
    Double_t x = TMath::Abs(xE)>1.e-5 ? (xE>0. ? rmax/TMath::Sqrt(tanTheta*tanTheta+1.) : -rmax/TMath::Sqrt(tanTheta*tanTheta+1.)) : 0.;
    Double_t y = TMath::Abs(xE)>1.e-5 ? (yE>0. ? tanTheta*rmax/TMath::Sqrt(tanTheta*tanTheta+1.) : -tanTheta*rmax/TMath::Sqrt(tanTheta*tanTheta+1.)) : (yE>0. ? rmax : -rmax);

    //std::cout << "x = " << x << "  y = " << y << std::endl;
    
    Int_t xBin = x > getBinCenter(getBinNumber(x)) ? getBinNumber(x) : getBinNumber(x)-1;
    Int_t yBin = y > getBinCenter(getBinNumber(y)) ? getBinNumber(y) : getBinNumber(y)-1;

    Double_t xp = (x-getBinCenter(xBin))/binWidth;
    Double_t yp = (y-getBinCenter(yBin))/binWidth; 

    // if we lie directly on the center of a bin, use only one bin to the left/down if val is neg. in interpolation (Shift bins)
    if (xp==1. && x<0.) {xp=0.; xBin+=1;}
    if (yp==1. && y<0.) {yp=0.; yBin+=1;}

    //std::cout << "xBin = " << xBin << "  yBin = " << yBin << std::endl;
    //std::cout << "xp = " << xp << "  yp = " << yp << std::endl;

    Int_t xx=0;
    for (int xb=xBin-1; xb<=xBin+2; xb++) {
      Int_t yy=0;
      for (int yb=yBin-1; yb<=yBin+2; yb++) {
	val[0]+=(funcs[xx]->Eval(xp))*(funcs[yy]->Eval(yp))*posMap[xb][yb][0];
	yy++;
      }
      xx++;
    }
  }
  
  else {
    Int_t xBin = xE > getBinCenter(getBinNumber(xE)) ? getBinNumber(xE) : getBinNumber(xE)-1;
    Int_t yBin = yE > getBinCenter(getBinNumber(yE)) ? getBinNumber(yE) : getBinNumber(yE)-1;
    
    Double_t xp = (xE-getBinCenter(xBin))/binWidth;
    Double_t yp = (yE-getBinCenter(yBin))/binWidth; 
    
    Int_t xx=0;
    for (int xb=xBin-1; xb<=xBin+2; xb++) {
      Int_t yy=0;
      for (int yb=yBin-1; yb<=yBin+2; yb++) {
	val[0]+=(funcs[xx]->Eval(xp))*(funcs[yy]->Eval(yp))*posMap[xb][yb][0];
	  yy++;
      }
      xx++;
    }
    //std::cout << p << " " << val[p] << std::endl;
  }
  
    
  if (xW*xW+yW*yW >= rmax*rmax) {

    Double_t tanTheta = TMath::Abs(xW)>1.e-5 ? TMath::Abs(yW)/TMath::Abs(xW) : 1.;
    Double_t x = TMath::Abs(xW)>1.e-5 ? (xW>0. ? rmax/TMath::Sqrt(tanTheta*tanTheta+1.) : -rmax/TMath::Sqrt(tanTheta*tanTheta+1.)) : 0.;
    Double_t y = TMath::Abs(xW)>1.e-5 ? (yW>0. ? tanTheta*rmax/TMath::Sqrt(tanTheta*tanTheta+1.) : -tanTheta*rmax/TMath::Sqrt(tanTheta*tanTheta+1.)) : (yW>0. ? rmax : -rmax);
    Int_t xBin = x > getBinCenter(getBinNumber(x)) ? getBinNumber(x) : getBinNumber(x)-1;
    Int_t yBin = y > getBinCenter(getBinNumber(y)) ? getBinNumber(y) : getBinNumber(y)-1;
    
    Double_t xp = (x-getBinCenter(xBin))/binWidth;
    Double_t yp = (y-getBinCenter(yBin))/binWidth; 
    
    if (xp==1. && x<0.) {xp=0.; xBin+=1;}
    if (yp==1. && y<0.) {yp=0.; yBin+=1;}
    
    
    Int_t xx=0;
    for (int xb=xBin-1; xb<=xBin+2; xb++) {
      Int_t yy=0;
      for (int yb=yBin-1; yb<=yBin+2; yb++) {
	val[1]+=(funcs[xx]->Eval(xp))*(funcs[yy]->Eval(yp))*posMap[xb][yb][1];
	yy++;
      }
      xx++;
    }
  }
  
  
  else {
    
    Int_t xBin = xW > getBinCenter(getBinNumber(xW)) ? getBinNumber(xW) : getBinNumber(xW)-1;
    Int_t yBin = yW > getBinCenter(getBinNumber(yW)) ? getBinNumber(yW) : getBinNumber(yW)-1;
    
    Double_t xp = (xW-getBinCenter(xBin))/binWidth;
    Double_t yp = (yW-getBinCenter(yBin))/binWidth; 
    
    Int_t xx=0;
    for (int xb=xBin-1; xb<=xBin+2; xb++) {
      Int_t yy=0;
      for (int yb=yBin-1; yb<=yBin+2; yb++) {
	val[1]+=(funcs[xx]->Eval(xp))*(funcs[yy]->Eval(yp))*posMap[xb][yb][1];
	yy++;
      }
      xx++;
    }
    // std::cout << p << " " << val[p] << std::endl;
  }
  
  return val;
  
};

//////////////////////////////////////////////////////////////////////////////////

/*int main(int argc, char* argv[]) {
  
  int XePeriod = 3; 
  PositionMap map(5.);
  map.readPositionMap(XePeriod);

  Double_t xE = atof(argv[1]);
  Double_t yE = atof(argv[2]);
  Double_t xW = atoi(argv[3]);
  Double_t yW = atoi(argv[4]);

  std::vector <Double_t> eta = map.getInterpolatedEta(xE, yE, xW, yW);
  
  for (int i=0; i<8; i++) std::cout << eta[i] << std::endl;
  }*/
