#include "calibrationTools.hh"
#include <TString.h>


TString getGeometry(int run) {

  if ( run < 20000 ) return TString("2011-2012");
  else if ( run >= 21087 && run < 21679) return TString("2012-2013_isobutane");
  else return TString("2012-2013");

};



LinearityCurve::LinearityCurve(Int_t period, bool useTanhSmear) : sourceCalPeriod(period), useTanh(useTanhSmear) {
  
  if (useTanh) {
    nParams = 6;
    linCurve = new TF1("linCurve","([0] + [1]*x + [2]*x*x)*(0.5+0.5*TMath::TanH((x-[4])/[5]))+([3]*x)*(0.5-0.5*TMath::TanH((x-[4])/[5]))", -200.5, 4096.5);
  }
  else {
    nParams = 3;
    linCurve = new TF1("linCurve","([0] + [1]*x + [2]*x*x)", -200.5,4096.5);
  }
  
  pmtParams.resize(8, std::vector <Double_t> (nParams,0.));
  
  
  //
  readLinearityCurveParams(sourceCalPeriod);

}

LinearityCurve::~LinearityCurve() {
  delete linCurve;
}


void LinearityCurve::readLinearityCurveParams(Int_t srcPeriod) {
  // Read linearity curve                                                                                                                      
  char tempFileLinearityCurve[500];

  sprintf(tempFileLinearityCurve, "%s/lin_curves_srcCal_Period_%i.dat",getenv("LINEARITY_CURVES"),srcPeriod);
  std::cout << "Reading in linearity curve from " << tempFileLinearityCurve << ":\n";

  if (useTanh) std::cout << "p0\tp1\tp2\tp3\tp4\tp5\n";
  else std::cout << "p0\tp1\tp2\n";

  std::ifstream fileLinearityCurve(tempFileLinearityCurve);
  //std::vector <Double_t> p(nParams);                                                           
  Double_t p;

  Int_t i = 0;
  Int_t ii=0;
  while (fileLinearityCurve >> p) {
    pmtParams[i][ii] = p;
    std::cout << p << " ";
    ii++;
    if (ii==nParams) { ii = 0; i++; std::cout << "\n"; }

    if (fileLinearityCurve.fail()) break;
  }
                                                          
  sourceCalPeriod = srcPeriod;
}

Double_t LinearityCurve::applyLinCurve(Int_t pmt, Double_t x) {
  return linCurve->EvalPar(&x,&pmtParams[pmt][0]);
}

Double_t LinearityCurve::applyInverseLinCurve(Int_t pmt, Double_t y) {
  if (useTanh) linCurve->SetParameters(pmtParams[pmt][0], pmtParams[pmt][1], pmtParams[pmt][2], pmtParams[pmt][3], pmtParams[pmt][4], pmtParams[pmt][5]);
  else linCurve->SetParameters(pmtParams[pmt][0], pmtParams[pmt][1], pmtParams[pmt][2]);
  return linCurve->GetX(y, -50., 3000.);
}

////////////////////////////////////////////////////////////////////////


WirechamberCal::WirechamberCal(Int_t run) : _run(run) {
  
  _nParams = 5;

  _calFunc = new TF1("calFunc","([0] + [1]*x)", 0.,4096.5);
  _extrap = new TF1("extrap","[0]*TMath::Power(x,[1])", -10.,4096.5);

  _params.resize(2, std::vector <Double_t> (_nParams,0.));
  
  readWirechamberCalParams(_run);

}

WirechamberCal::~WirechamberCal() {
  delete _calFunc;
  delete _extrap;
}


void WirechamberCal::readWirechamberCalParams(Int_t run) {                                                                                                                     
  char tempFileCal[500];

  sprintf(tempFileCal, "%s/runs/MWPC_cal_%i.dat",getenv("MWPC_CALIBRATION"),run);
  std::cout << "Reading in MWPC calibration from " << tempFileCal << ":\n";

  std::cout << "p0\tp1\n";

  std::ifstream fileCal(tempFileCal);
 
  std::string label;

  if ( fileCal.is_open() ) {
    
    fileCal >> label;
    for ( int i=0; i<_nParams; ++i ) { fileCal >> _params[0][i]; std::cout << _params[0][i] << "\t"; }
    std::cout << "\n";
    fileCal >> label;
    for ( int i=0; i<_nParams; ++i ) { fileCal >> _params[1][i]; std::cout << _params[1][i] << "\t"; }
    std::cout << "\n";
    fileCal.close();
  }
  
  else  std::cout << TString::Format("Couldn't read wirechamber calibration for run %i",run).Data() << std::endl;
    
                                                       
  _run = run;
}

Double_t WirechamberCal::applyCal(Int_t side, Double_t adc) {

  if ( adc > _params[side][2] ) { 
    Double_t par[] {_params[side][0],_params[side][1]};
    return _calFunc->EvalPar(&adc,par);
  }
  else if ( adc > 0. ) {
    Double_t par[] {_params[side][3],_params[side][4]};
    return _extrap->EvalPar(&adc,par);
  }
  else return 0.;
    
}


///////////////////////////////////////////////////////

BackscatterSeparator::BackscatterSeparator() : _run(0),_geometry("") {

  params.resize(3,0.);
  _func = new TF1("_func","[0]+[1]*TMath::Exp(-x/[2])",0., 1200.);
}

void BackscatterSeparator::LoadCutCurve(int run) {

  _run = run;
  TString _geometry = getGeometry(_run);

  TString filename = TString::Format("%s/backscSepParameters_%s.dat",
				     getenv("MWPC_CALIBRATION"),_geometry.Data());
  std::ifstream infile(filename.Data());

  std::string hold;
  infile >> hold >> hold >> hold;
  infile >> params[0] >> params[1] >> params[2];
  
  infile.close();
}

Int_t BackscatterSeparator::separate23(Double_t en) {

  Double_t cut = _func->EvalPar(&en,&params[0]);

  Int_t type = ( en > cut ) ? 3 : 2;  
  
  return type;
};







///////////////////////////////////////////////////////////////////////

TriggerFunctions::TriggerFunctions(Int_t run) : currentRun(run) {

  triggerFuncs.resize(8,std::vector<Double_t>(8,0.));
  readTriggerFunctions(currentRun);
  
  func = new TF1("func","([0]+[1]*TMath::Erf((x-[2])/[3]))*(0.5-.5*TMath::TanH((x-[2])/[4]))+(0.5+.5*TMath::TanH((x-[2])/[4]))*([5]+[6]*TMath::TanH((x-[2])/[7]))",-100.,4096.);

}

TriggerFunctions::~TriggerFunctions() {
  
  delete func;

}


void TriggerFunctions::readTriggerFunctions(Int_t run) {

  std::ifstream file(TString::Format("%s/runs/%i_thresholds.dat",getenv("TRIGGER_FUNC"),run).Data());

  std::vector <Double_t> params(8,0.);

  if (file.is_open()) {
    Int_t pmt = 0;
    while (file >> params[0] >> params[1] >> params[2] >> params[3] >> params[4] >>
	   params[5] >> params[6] >> params[7]) {
      
      for (Int_t param=0; param<8; param++) {
	triggerFuncs[pmt][param] = params[param];
	std::cout << triggerFuncs[pmt][param] << " ";
      }
      std::cout << std::endl;
      pmt++;
    }
  }
  else { 
    std::cout << "Couldn't open trigger function file for run " << run << std::endl;
    exit(0);
  }
  file.close();

  //Making sure the trigger probability is zero for WPMT2 when in the bad run range...
  if (run<=17249 || run>=16983) {
    triggerFuncs[5][0] = triggerFuncs[5][1] = triggerFuncs[5][5] = triggerFuncs[5][6] = 0.;
  }

}

bool TriggerFunctions::decideEastTrigger(std::vector<Double_t> ADC, std::vector<Double_t> rands) {

  Int_t numTriggs = 0;
  Double_t adc = 0.;
  
  adc = ADC[0];
  numTriggs = rands[0] < func->EvalPar(&adc,&triggerFuncs[0][0]) ? numTriggs+1 : numTriggs;
  adc = ADC[1];
  numTriggs = rands[1] < func->EvalPar(&adc,&triggerFuncs[1][0]) ? numTriggs+1 : numTriggs;
  
  if (numTriggs<2) {
    adc = ADC[2];
    numTriggs = rands[2] < func->EvalPar(&adc,&triggerFuncs[2][0]) ? numTriggs+1 : numTriggs;
  }
  else return true;
  
  if (numTriggs<2) {
    adc = ADC[3];
    numTriggs = rands[3] < func->EvalPar(&adc,&triggerFuncs[3][0]) ? numTriggs+1 : numTriggs;
  }
  else return true;

  return numTriggs>=2 ? true : false;

}

bool TriggerFunctions::decideWestTrigger(std::vector<Double_t> ADC, std::vector<Double_t> rands) {

  Int_t numTriggs = 0;
  Double_t adc = 0.;
  
  adc = ADC[0];
  numTriggs = rands[0] < func->EvalPar(&adc,&triggerFuncs[4][0]) ? numTriggs+1 : numTriggs;

  //std::cout << rands[0] << " " << func->EvalPar(&adc,&triggerFuncs[4][0]) << std::endl;

  adc = ADC[1];
  numTriggs = rands[1] < func->EvalPar(&adc,&triggerFuncs[5][0]) ? numTriggs+1 : numTriggs;
  
  if (numTriggs<2) {
    adc = ADC[2];
    numTriggs = rands[2] < func->EvalPar(&adc,&triggerFuncs[6][0]) ? numTriggs+1 : numTriggs;
  }
  else return true;

  if (numTriggs<2) {
    adc = ADC[3];
    numTriggs = rands[3] < func->EvalPar(&adc,&triggerFuncs[7][0]) ? numTriggs+1 : numTriggs;
  }
  else return true;

  return numTriggs>=2 ? true : false;

};




EreconParameterization::EreconParameterization(Int_t runNumber) {

  if (runNumber<16000) geometry = TString("2011-2012");
  else if (runNumber<20000) geometry = TString("2011-2012");
  else if (runNumber<=21628 && runNumber>21087) geometry = TString("2012-2013_isobutane");
  else if (runNumber<24000) geometry = TString("2012-2013");
  else {
    std::cout << "Bad runNumber passed to getEQ2EtrueParams\n";
    exit(0);
  }
  params.resize(2,std::vector < std::vector < double > > (3, std::vector < double > (7,0.)));
  readParams();

};

void EreconParameterization::readParams() {

  std::ifstream infile;
  TString basePath = TString::Format("%s/simulation_comparison/EQ2EtrueConversion/",getenv("ANALYSIS_CODE")); 
  basePath += TString::Format("%s_EQ2EtrueFitParams.dat",geometry.Data());
  
  infile.open(basePath.Data());
  
  char holdType[10];
  int side=0, type=0;
  while (infile >> holdType >> params[side][type][0] 
	 >> params[side][type][1] >> params[side][type][2] 
	 >> params[side][type][3] >> params[side][type][4] 
	 >> params[side][type][5] >> params[side][type][6]) { 

    std::cout << holdType << " " << params[side][type][0] << " " 
	      << params[side][type][1] << " " << params[side][type][2]
	      << " " << params[side][type][3] << " " << params[side][type][4] 
	      << " " << params[side][type][5] << std::endl;

    type+=1;
    if (type==3) {type=0; side=1;}
  }
};

Double_t EreconParameterization::getErecon(Int_t side, Int_t type, Double_t Evis) {

  Double_t EnCutoff = params[side][type][4]; // Below the EQ, we extrapolate to zero

  // Above EnCutoff: p0 + p1*eQ + p2/eQ + p3/eQ*eQ
  // Below : C*eQ^alpha where p4 is the cutoff, p5 is alpha, and p6 is C

  if ( Evis > EnCutoff ) {

    return ( params[side][type][0] + params[side][type][1]*Evis + 
	     params[side][type][2]/Evis + params[side][type][3]/(Evis*Evis) );
  }

  else {
    
    return params[side][type][6] * TMath::Power(Evis,params[side][type][5]);

  }

};
