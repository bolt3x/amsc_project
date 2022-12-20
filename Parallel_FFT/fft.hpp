#ifndef FFT_HPP
#define FFT_HPP

#include <iostream>
#include <vector>
#include <numbers>
#include <complex>
#include <stdlib.h>
#include <functional>
#include <random>
#include <algorithm>


#define N_STANDARD 150


//Out of class declaration of Parallel FFT function (MPI)
//The parallel functio has to be defined outside the class, otherwise each process would have to define
// a class.
//To avoid that possibility, we define the parallel vesion as a free function

void Parallel_FFT(std::vector<std::complex<double>> &x);



#endif