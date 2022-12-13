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
template <typename T> 
void TimeIt(void (*func)(T),
            T arg,
            std::ostream &out);


// print function

void PrintIt(const std::vector<std::complex<double>> &v,std::string msg,std::ostream &out = std::cout);

// reverse bit function

size_t ReverseBit(size_t num, const unsigned n);

// random signal generator

std::vector<std::complex<double>> RandomGen(unsigned long N, bool specific = true);

#endif