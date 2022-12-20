#ifndef FFT_HPP
#define FFT_HPP

#include <iostream>
#include <vector>
#include <complex>
#include <stdlib.h>
#include <functional>
#include <random>
#include <algorithm>
#include <omp.h>
#include <mpi.h>



//---------------------------------FFT CLASS------------------------------

//Out of class declaration of Parallel FFT function (MPI)

void MPI_FFT(std::vector<std::complex<double>> &x);



#endif
