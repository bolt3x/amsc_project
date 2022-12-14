#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <iostream>
#include <vector>
#include <complex>
#include <stdlib.h>
#include <chrono>
#include <functional>


//------------------------------UTILITIES----------------------------------

// timing function for class methods
template<class C>
void TimeIt(void (C::*func)(),
            C& obj,
            std::ostream &out = std::cout);

// print function
void PrintIt(const std::vector<std::complex<double>> &v,std::string msg = "",std::ostream &out = std::cout);

// reverse bit function

size_t ReverseBit(size_t num, const unsigned n);

// Template definitions
/*-----------------------------------------------------
TIMEIT function: to time the execution time of a specified parameter-given routine. Specify:
                * function to be used to compute FFT;
                * input signal and size;
                * output vector for transform
                * ostream: where to print the results
-----------------------------------------------------*/
template<class C>
void TimeIt(void (C::*func)(),
            C& obj,
            std::ostream &out)
{
    auto start = std::chrono::high_resolution_clock::now();
    (obj.*func)();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    out << "Time required to run: "<<duration.count() << " microseconds."<< std::endl;

}



#endif
