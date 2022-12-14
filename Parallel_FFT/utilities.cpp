#include "utilities.hpp"
#include "fft.hpp"

using namespace std::chrono;


/*-----------------------------------------------------
TIMEIT function: to time the execution time of a specified parameter-given routine. Specify:
                * function to be used to compute FFT;
                * input signal and size;
                * output vector for transform
                * ostream: where to print the results
-----------------------------------------------------*/

template <typename T>
void TimeIt(void (*func)(T&),
            T& arg, 
            std::ostream &out)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    auto start = high_resolution_clock::now();

    (*func)(arg);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    out << "Time required to run per process#" << rank <<": "<<duration.count() << " microseconds."<< std::endl;

}


template void TimeIt <std::vector<std::complex<double>>>(void (*)(std::vector<std::complex<double>>&), std::vector<std::complex<double>>&, std::ostream&);

/*-----------------------------------------------------
PRINTIT function: print the result of a specified vector
-----------------------------------------------------*/

void PrintIt(const std::vector<std::complex<double>> &v,std::string msg,std::ostream &out){

    for(int i = 0; i < 10; i++){
        out << "-";
    }
    out << std::endl;
    
    out << msg << std::endl;
    for(auto &i : v){
        out <<std::setprecision(2) << i << " ";
    }
    out << std::endl;
    
    for(int i = 0; i < 10; i++){
        out << "-";

    }
    out << std::endl;
    return;
}

/*----------------------------------------------------
REVERSE BIT : given a number num (the index representig a certain element of the original), 
              and the size of the orignal signal n, generate the bit reversed number
------------------------------------------------------*/

size_t ReverseBit(size_t num, const unsigned n)
{

    // * num is the number to reverse;

    unsigned int bits = static_cast<unsigned int>(std::log2(n));
    // * bits expresses the number of valid bits of num, for a correct bit reversal,
    //   dependent on the size of the vector considered

	// num = (((num & 0xAAAAAAAA) >> 1) | ((num & 0x55555555) << 1));
	// num = (((num & 0xCCCCCCCC) >> 2) | ((num & 0x33333333) << 2));
	// num = (((num & 0xF0F0F0F0) >> 4) | ((num & 0x0F0F0F0F) << 4));
	// num = (((num & 0xFF00FF00) >> 8) | ((num & 0x00FF00FF) << 8));
	
    // return((num >> 16) | (num << 16));

        
    size_t nrev;
    unsigned int count, N;  
    N = 1<<bits;
    count = bits-1;   // initialize the count variable
    nrev = num;
    for(num>>=1; num; num>>=1)
    {
        nrev <<= 1;
        nrev |= num & 1;
        count--;
    }

    nrev <<= count;
    nrev &= N - 1;

    return nrev;


}


/*-----------------------------------------------------
RANDOM GEN function: give a specified dimension, generate a signal of random positive numbers
-----------------------------------------------------*/


std::vector<std::complex<double>> RandomGen(unsigned long N, bool specific)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<std::complex<double>> signal;
    if(rank == 0){
        signal.resize(N);
        if(!specific){
            std::mt19937 gen32; // MersenneTWISTER
            for(size_t i = 0; i< N; ++i){

                signal[i]=(static_cast<double>(gen32()%10)/*module*/*std::polar(1.0,2*std::numbers::pi / (gen32()%gen32())))/*random complex number of module 1*/;

            }
        }
        else{

            for(size_t i = 0; i < N; i++){

                signal[i] = std::complex<double>(static_cast<double>(i),/*std::pow(-1.0, i +1)*/1);
            
            }
        }
        
    }

    return signal;

}