#include <iostream>
#include <vector>
#include <complex>
#include "fft.hpp"
#include "utilities.hpp"

int main(){

    size_t n = 8;
    std::vector<std::complex<double>> input{1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 };

    FFTGenerator fft(input,n);
    
    TimeIt(&FFTGenerator::Recursive_FFT,fft);

    fft.Iterative_FFT();

    PrintIt(input,"Input");
    PrintIt(fft.getRecT(),"Recursive");
    PrintIt(fft.getIterT(),"Iterative");

    input = fft.getIterT();

    FFTGenerator ifft(input,n);

    ifft.Inverse_FFT();

    PrintIt(input,"Inverse Input");
    PrintIt(ifft.getInvT(),"Inverse");

    return 0;

}