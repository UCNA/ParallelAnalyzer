#ifndef MBUTILS_HH
#define MBUTILS_HH

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <string>


std::string itos(int val); //Take integer to a string

std::string ftos(double val); //Take float or double to a string

double power(double num, int pow); //calculate the integer power of a number

std::vector <int> sortVecInt(std::vector<int> vec,
			     bool hightolow); //Sort vector of ints

std::vector <double> sortVecDouble(std::vector<double> vec,
			     bool hightolow); //Sort vector of ints



#endif
