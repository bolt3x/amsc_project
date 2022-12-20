#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <iostream>
#include <vector>
#include <complex>
#include <stdlib.h>
#include <chrono>
#include <functional>
#include <fstream>
#include <iterator>


//------------------------------UTILITIES----------------------------------


// print function
void PrintIt(const std::vector<std::complex<double>> &v,const std::string &msg = "", std::ostream &out = std::cout);

// reverse bit function

size_t ReverseBit(size_t num, const unsigned n);

// read signal from file function
void ReadIt(const std::string& file,std::vector<std::complex<double>> &v,size_t &n);

//random or specific generator of signals
std::vector<std::complex<double>> RandomGen(const unsigned long &N, const bool &specific = true);

#endif
