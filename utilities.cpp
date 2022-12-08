#include "utilities.hpp"
#include "fft.hpp"




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


template void TimeIt<FFTGenerator>(void (FFTGenerator::*)(), FFTGenerator&, std::ostream&);

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