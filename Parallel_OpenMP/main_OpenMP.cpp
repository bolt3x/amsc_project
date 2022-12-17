#include <iostream>

#include "fft_OpenMP.hpp"
#include "utilities_OpenMP.hpp"



int main(int argc, char ** argv){


    if(argc != 3){
        
        std::cout << "Wrong number of parameters to start the Parallel FFT"<< std::endl;

        return 0;
    }

    unsigned int n = std::stoi(argv[1]); /*power of 2*/;
    // the number of elements chosen is already known to all the processors

    unsigned long num_threads = std::stoi(argv[2]); /* power of 2, less than n*/
    // number of threads that process the elements of the signal

    if(!((n & (n-1)) == 0) && !((num_threads & (num_threads-1)) == 0) && num_threads > n){

        std::cout << "Either the number of elements or the number of process is incorrect"<< std::endl;
        return 0;
    }

    
    // after the construction is done, each process calls the parallel function...

    std::vector<std::complex<double>> x = RandomGen(n);
    
    PrintIt(x, "Input signal specifically generated");

    TimeIt(&Parallel_FFT, x, num_threads);


    return 0;
}
