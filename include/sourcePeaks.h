//Holds the smeared source peaks which can be changed in single location

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <vector>
#include <cstdlib>
#include <string>

#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

using namespace std;

// These are the weighted true peaks from the KLM CE sources
const static double peakCe = 130.3;// 131.5956;//80.5;
const static double peakIn = 174.35;
const static double peakSn = 368.4938;//317.8;
const static double peakBiLow = 498.;//501.420;//448.8;
const static double peakBiHigh = 993.789;//926.;

const static double peakCe_EQ = 80.5;
const static double peakIn_EQ = 120.;
const static double peakSn_EQ = 317.8;
const static double peakBiLow_EQ = 448.8;
const static double peakBiHigh_EQ = 926.;

//Below "type" should be "EQ" or "Etrue"
vector < vector <double> > returnPeaks(int srcPeriod, string type) {
  vector < vector <double> > peaks;
  peaks.resize( 8, vector <double> (4,0.));
  ifstream peakfile;
  char filename[500];
  if (type=="EQ") {
    sprintf(filename, "../smeared_peaks/weighted_smeared_peaks/fits/weightedSimPeaks_PMTbyPMT_runPeriod_%i.dat",srcPeriod);}
  else {
    sprintf(filename, "../smeared_peaks/%s_smeared_%i.dat",type.c_str(),srcPeriod);}
  peakfile.open(filename);
  
  for (int pmt = 0; pmt<8; pmt++) {
    peakfile >> peaks[pmt][0] >> peaks[pmt][1] >> peaks[pmt][2] >> peaks[pmt][3];
  }
  return peaks;
}

vector < vector <double> > EQ2EtrueFit(int srcPeriod) {
  vector < vector <double> > fits;
  fits.resize( 8, vector <double> (3,0.));
  ifstream sourcefile;
  char filename[500];
  sprintf(filename, "../smeared_peaks/EQ2Etrue_srcPeriod_%i.dat",srcPeriod);
  sourcefile.open(filename);
  
  for (int pmt = 0; pmt<8; pmt++) {
    sourcefile >> fits[pmt][0] >> fits[pmt][1] >> fits[pmt][2];
  }
  return fits;
}

string getIndiumSide(int runNumber) {
  std::string dbAddress = std::string(getenv("UCNADBADDRESS"));
  std::string dbname = std::string(getenv("UCNADB"));
  std::string dbUser = std::string(getenv("UCNADBUSER"));
  std::string dbPass = std::string(getenv("UCNADBPASS"));
  std::string port = "3306";
  std::string dbAddressFull = "mysql://"+dbAddress+":"+port+"/"+dbname;
  
  std::string qresult;
  std::string Side="";
  bool passFlag = false;

  char cmd[200];
  sprintf(cmd,"SELECT sourcetype FROM sources WHERE run_number=%i;",runNumber);

  TSQLServer *db = TSQLServer::Connect(dbAddressFull.c_str(), dbUser.c_str(), dbPass.c_str());
  if (!db) cout << "Couldn't connect to database\n";
  else {
    //cout << "Connected to DB Server\n";
    TSQLResult *res = db->Query(cmd);
    int rows = res->GetRowCount();
    //cout << rows << endl;
    TSQLRow *row = res->Next();
   
    if (!row) {
      cout << "This run wasn't a source run of any type\n";
    }
    else {
      for (int i=0; i<rows; i++) {
	qresult = std::string(row->GetField(0));
	if (qresult=="In114W") {Side="West"; passFlag = true; continue;}
	else if (qresult=="In114E") {Side="East"; passFlag = true; continue;}
	row = res->Next();
      }
      if (passFlag) std::cout << "Run " << runNumber << " Indium Facing " << Side << std::endl;
      else cout << "This run is not an Indium Run\n";
    }
  }
  //delete (row);
  //delete (res); 
  db->Close();
  return Side;
}

