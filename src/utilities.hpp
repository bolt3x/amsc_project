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
            std::ostream &out = std::cout)
{
    auto start = std::chrono::high_resolution_clock::now();
    (obj.*func)();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    out << "Time required to run: "<<duration.count() << " microseconds."<< std::endl;

}

template <typename C >
void TimeIt(void (C::*func)( const unsigned long&),
            C& obj,
            const unsigned long& threads,
            std::ostream &out = std::cout)
{
    auto start = std::chrono::high_resolution_clock::now();

    (obj.*func)(threads);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    out << "Time required to run per process#: "<<duration.count() << " microseconds."<< std::endl;

    return;
}




// print function
void PrintIt(const std::vector<std::complex<double>> &v,const std::string &msg = "", std::ostream &out = std::cout);

// reverse bit function

size_t ReverseBit(size_t num, const unsigned n);

// read signal from file function
void ReadIt(const std::string& file,std::vector<std::complex<double>> &v,size_t &n);

//random or specific generator of signals
std::vector<std::complex<double>> RandomGen(const unsigned long &N, const bool &specific = true);

#endif
