#include "utilities.hpp"
#include "fft.hpp"


/*-----------------------------------------------------
TIMEIT function: to time the execution time of a specified parameter-given routine. Specify:
                * function to be used to compute FFT;
                * input signal;
                * eventually the number of threads
                * ostream: where to print the results
-----------------------------------------------------*/

template void TimeIt<FFTGenerator>(void (FFTGenerator::*)(), FFTGenerator&, std::ostream&);

template void TimeIt<FFTGenerator>(void (FFTGenerator::*)(const unsigned long&), FFTGenerator&, const unsigned long&, std::ostream&);

/*-----------------------------------------------------
PRINTIT function: print the result of a specified vector
-----------------------------------------------------*/

void PrintIt(const std::vector<std::complex<double>> &v,const std::string &msg, std::ostream &out){

    if(v.size() == 0){
        return;
    }
    
    for(int i = 0; i < 10; i++){
        out << "-";
    }
    out << std::endl;
    
    out << msg << std::endl;
    for(auto &i : v){
        out << i << " ";
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

/*----------------------------------------------------
READIT : given a string it will try to read a vector
            from the file with that name
------------------------------------------------------*/

void ReadIt(const std::string& file,std::vector<std::complex<double>> &v,size_t &n){
    
    std::ifstream input(file);
    std::vector<std::complex<double>> buffer{
            std::istream_iterator<std::complex<double>>(input), 
            std::istream_iterator<std::complex<double>>() };
    
    v = buffer;
    n = v.size();
    
    return;
}


/*-----------------------------------------------------
RANDOM GEN function: give a specified dimension, generate a signal of random positive numbers
-----------------------------------------------------*/

std::vector<std::complex<double>> RandomGen(const unsigned long &N, const bool &specific){

    std::vector<std::complex<double>> signal;
    signal.resize(N);
    if(!specific){
        std::mt19937 gen32; // MersenneTWISTER
        for(size_t i = 0; i< N; ++i){

            signal[i]=(static_cast<double>(gen32()%10)/*module*/*std::polar(1.0,2*M_PI/(gen32()%10)))/*random complex number of module 1*/;

        }
    }
    else{

        for(size_t i = 0; i< N; ++i){

            signal[i] = std::complex<double>(i<N/2?1.0:0.0, 0.0);

        }
    }
    

    return signal;

}