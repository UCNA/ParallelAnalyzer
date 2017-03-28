#include <map>
using namespace std;


// Returns the calibration period to be used for this run
unsigned int getSrcRunPeriod(int runNumber) {
  unsigned int calibrationPeriod=0;
  if (runNumber <= 17297) {
    calibrationPeriod=1;
  }
  else if (runNumber <= 17439) {
    calibrationPeriod=2;
  }
  else if (runNumber <= 17734) {
    calibrationPeriod=3;
  }
  else if (runNumber <= 17955) {
    calibrationPeriod=4;
  }
  else if (runNumber <= 18386) {
    calibrationPeriod=5;
  }
  else if (runNumber <= 18683) {
    calibrationPeriod=6;
  }
  else if (runNumber <= 18994) {
    calibrationPeriod=7;
  }
  else if (runNumber <= 19239) {
    calibrationPeriod=8;
  }
  else if (runNumber <= 19544) {
    calibrationPeriod=9;
  }
  else if (runNumber <= 20000) {
    calibrationPeriod=11;
  }
  else if (runNumber <= 20741) { // Octets 60-65
    calibrationPeriod=13;
  }
  else if (runNumber <= 20837) { // Octets 66
    calibrationPeriod=14;
  }
  else if (runNumber <= 21237) { // Octets 67-71
    calibrationPeriod=16;
  }
  else if (runNumber <= 21605) { // Octets 72-79
    calibrationPeriod=17;
  }
  else if (runNumber <= 21863) { // Octets 80-85
    calibrationPeriod=18;
  }
  else if (runNumber <= 22118) { // Octets 86-91
    calibrationPeriod=19;
  }
  else if (runNumber <= 22238) { // Octets 92-95
    calibrationPeriod=20;
  }
  else if (runNumber <= 22630) { // Octets 96-105
    calibrationPeriod=22;
  }
  else if (runNumber <= 23173) { // Octets 106-121
    calibrationPeriod=23;
  }
  return calibrationPeriod;
}

unsigned int getXeRunPeriod(int runNumber) {
  unsigned int calibrationPeriod=0;
  if (runNumber <= 18080) {
    calibrationPeriod=2;
  }
  else if (runNumber <= 18389) {
    calibrationPeriod=3;
  }
  else if (runNumber <= 18711) {
    calibrationPeriod=4;
  }
  else if (runNumber <= 19238) {
    calibrationPeriod=5;
  }
  //else if (runNumber <= 19872) {
  //calibrationPeriod=6;
  //}
  else if (runNumber <= 20000) {
    calibrationPeriod=7;
  }
  else if (runNumber <= 21605) {
    calibrationPeriod=8;
  }
  else if (runNumber <= 22238) {
    calibrationPeriod=9;
  }
  else if (runNumber <= 23173) {
    calibrationPeriod=10;
  }
  else {
    cout << "Bad run number\n";
    exit(0);}

  return calibrationPeriod;
}

unsigned int getXeRunPeriodForMWPCmap(int runNumber) {
  unsigned int calibrationPeriod=0;
  if (runNumber <= 18080) {
    calibrationPeriod=2;
  }
  else if (runNumber <= 18389) {
    calibrationPeriod=3;
  }
  else if (runNumber <= 18711) {
    calibrationPeriod=4;
  }
  else if (runNumber <= 19544) { // This is different due to sudden change in 
                                 // West anode spectra (gain cranked way up)
    calibrationPeriod=5;
  }
  //else if (runNumber <= 19872) {
  //calibrationPeriod=6;
  //}
  else if (runNumber <= 20000) {
    calibrationPeriod=7;
  }
  else if (runNumber <= 21605) {
    calibrationPeriod=8;
  }
  else if (runNumber <= 22238) {
    calibrationPeriod=9;
  }
  else if (runNumber <= 23173) {
    calibrationPeriod=10;
  }
  else {
    cout << "Bad run number\n";
    exit(0);}

  return calibrationPeriod;
}

unsigned int getGainReferenceRun(int runNumber) {
  map <unsigned int,unsigned int> refRun;
  refRun[1]=17238;
  refRun[2]=17370;
  refRun[3]=17521;
  refRun[4]=17892;
  refRun[5]=18361;
  refRun[6]=18621;
  refRun[7]=18749;
  refRun[8]=19232;
  refRun[9]=19359;
  refRun[10]=19511;
  refRun[11]=19857;
  refRun[12]=19899;
  refRun[13]=20519;
  refRun[14]=20820;
  refRun[15]=20905;
  refRun[16]=21091;
  refRun[17]=21315; // start weird EPMT4 Bi Pulser
  refRun[18]=21683;
  refRun[19]=21918;
  refRun[20]=22219;
  refRun[21]=22298;
  refRun[22]=22441;
  refRun[23]=22771;
  refRun[24]=22925;
  unsigned int calPeriod = getSrcRunPeriod(runNumber);
  map <unsigned int, unsigned int>::iterator it = refRun.find(calPeriod);
  return it->second;
}
