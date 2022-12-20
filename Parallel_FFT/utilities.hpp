#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <iostream>
#include <vector>
#include <complex>
#include <stdlib.h>
#include <chrono>
#include <functional>
#include <iomanip>
#include <mpi.h>


//------------------------------UTILITIES----------------------------------

// timing function for class methods

template <typename T>
void TimeIt(void (*func)(T&),
            T& arg, 
            std::ostream &out = std::cout)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    auto start = std::chrono::high_resolution_clock::now();

    (*func)(arg);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    if(rank ==0)
        out << "Time required to run per process#" << rank <<": "<<duration.count() << " microseconds."<< std::endl;

    return;
}


// print function

void PrintIt(const std::vector<std::complex<double>> &v,const std::string &msg,std::ostream &out = std::cout);

// reverse bit function

size_t ReverseBit(size_t num, const unsigned n);

// random signal generator

std::vector<std::complex<double>> RandomGen(unsigned long N, bool specific = true);

#endif