#include "MBUtils.hh"

std::string itos(int val) {
  char temp[32];
  sprintf(temp,"%i",val);
  return std::string(temp);
};

std::string ftos(double val) {
  char temp[32];
  sprintf(temp,"%.1f",val);
  return std::string(temp);
};

double power(double num, int pow) {
  double result = 1.;
  if (pow==0) return result;
  else if (pow>0) {
    for (int i=0; i<pow; i++) result = result*num;
    return result;
  }
  else if (pow<0) {
    for (int i=0; i>pow; i--) result = result/num;   
    return result;
  }
  else throw "Error in power somewhere";
};

std::vector <int> sortVecInt(std::vector<int> vec,
			     bool hightolow) {

  unsigned int veclen = vec.size();
  
  for ( unsigned int i=1; i<veclen; ++i ) {
    
    for ( unsigned int j=0; j<i; ++j ) {

      if ( (hightolow && vec[i]>vec[j]) || (!hightolow && vec[i]<vec[j]) ) {
	vec.insert(vec.begin()+j,1,vec[i]);
	std::cout << i << " " << j << std::endl;
	vec.erase(vec.begin()+i+1);
	continue;
      }
    }
  }

  return vec;
};

std::vector <double> sortVecDouble(std::vector<double> vec,
			     bool hightolow) {

  unsigned int veclen = vec.size();
  
  for ( unsigned int i=1; i<veclen; ++i ) {
    
    for ( unsigned int j=0; j<i; ++j ) {

      if ( (hightolow && vec[i]>vec[j]) || (!hightolow && vec[i]<vec[j]) ) {
	vec.insert(vec.begin()+j,1,vec[i]);
	std::cout << i << " " << j << std::endl;
	vec.erase(vec.begin()+i+1);
	continue;
      }
    }
  }

  return vec;
};
    
