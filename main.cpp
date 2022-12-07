#include <iostream>

#include "fft.hpp"
#include "utilities.hpp"


int main(){
    unsigned int n =4/*power of 2*/;
    std::vector<std::complex<double>> input(n);
    for(size_t i = 0; i < n; i++){
        //input[i] = std::complex<double>(static_cast<double>(rand())/RAND_MAX,static_cast<double>(rand())/RAND_MAX);
        input[i] = std::complex<double>(static_cast<double>(i),std::pow(-1.0, i +1));
        //std::cout << input[i] << "\n";
        std::cout << ReverseBit(i,n) << std::endl;
    }
    std::cout << "\n" <<std::endl;

    FFTGenerator FFT(input,n);

    std::cout<< std::endl;


    PrintIt(FFT.getSignal());

    std::cout << std::endl;


    std::cout << "Starting construction of Recursive_FFT:" <<std::endl;
    TimeIt<FFTGenerator>(&FFTGenerator::Recursive_FFT(), FFT );
    PrintIt(FFT.getRecT());
    std::cout << std::endl;
    std::cout << "Starting construction of Iterative_FFT:" <<std::endl;
    FFT.Iterative_FFT();
    PrintIt(FFT.getIterT());



    return 0;
}
