#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <iostream>
#include <vector>
#include <complex>
#include <stdlib.h>
#include <chrono>
#include <functional>
#include <iomanip>


//------------------------------UTILITIES----------------------------------

// timing function for class methods
template<class C, typename T>
void TimeIt(void (C::*func)(T),
            C& obj,
            const T arg,
            std::ostream &out = std::cout);


// print function
void PrintIt(const std::vector<std::complex<double>> &v,std::string msg,std::ostream &out = std::cout);

// reverse bit function

size_t ReverseBit(size_t num, const unsigned n);

#endif