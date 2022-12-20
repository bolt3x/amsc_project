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

// timing function for class methods
template<class C>
void TimeIt(void (C::*func)(),
            C& obj,
            std::ostream &out = std::cout);

// print function
void PrintIt(const std::vector<std::complex<double>> &v,const std::string &msg = "", std::ostream &out = std::cout);

// reverse bit function

size_t ReverseBit(size_t num, const unsigned n);

// read signal from file function
void ReadIt(const std::string& file,std::vector<std::complex<double>> &v,size_t &n);

#endif
