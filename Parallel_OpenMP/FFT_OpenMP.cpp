#include "fft_OpenMP.hpp"
#include "utilities_OpenMP.hpp"



/*-----------------------------------------------------
PARALLEL FFT function: 
-----------------------------------------------------*/

// Implememtation of parallel function:
// First of all, each avaiable processor handles n/p elements.
// We identify 3 different phases:
// ---->  I PHASE  :: swap the elements of a into y with bit reverse. To do so, each processor handles only a section of the vectors
// a and y, so each processos will have to scatter specific elemnts among specific processors. The algorithm to do
// so is reported below.


void Parallel_FFT(std::vector<std::complex<double>> &x, const unsigned long &num_threads ){


    unsigned long N (x.size());

    
    #pragma omp parallel num_threads(num_threads) shared (x)

    // once we have created the threads with the specified number, we let OpenMP distribute the work among them
    #pragma omp for 
        // each thread swaps some elements of the input signal
        for (unsigned int i = 1; i < N/2; i++) {
            std::swap(x[ReverseBit(i, N)],x[i]);
        }

    // Once the first phase is complete, start with the second one

    #pragma omp for
        for(size_t j =1; j <= static_cast<size_t>(std::log2(N)); ++j){
            unsigned int d = 2<<(j-1); // size
            unsigned int d2 = d >> 1; // m2 = m/2
            // principle root of nth complex root of unity.
            std::complex<double> w_d(std::polar(1.0, static_cast<double>(- 2)*(std::numbers::pi/static_cast<double>(d))));
            std::complex<double> w(1.0,0.0);

            // k :: number of iterations representing the number of "consecutive combining pattern"
            
            for(size_t k = 0; k < d2; ++k) {

                // m :: number of iteration representing the number of total of "combining patterns" of a certain type 
                for(size_t m = k; m < N; m += d) {
                    std::complex<double> t = w * x[m + d2];
                    std::complex<double> u = x[m];

                    // similar calculating y[m]
                    x[m] = (u + t);

                    // similar calculating ym+n/2]
                    x[m + d2] = (u - t);
                }
                w *= w_d;


                }
        }

    //-------------DEBUGGING SECTION------------------

    PrintIt(x, "SwappedVec: ");
    


    //-------------DEBUGGING SECTION------------------

   
    

    return;
}
