#include <iostream>

#include "fft.hpp"
#include "utilities.hpp"


int main(int argc, char ** argv){

    if(argc != 2)
    {
        std::cout << "Wrong number of parameters to start the Parallel FFT"<< std::endl;
        return 0;
    }
    unsigned int n = std::stoi(argv[1])/*power of 2*/;
    std::vector<std::complex<double>> input(n);
    std::cout<< "Input signal is:" << std::endl;
    for(size_t i = 0; i < n; i++){
        //input[i] = std::complex<double>(static_cast<double>(rand())/RAND_MAX,static_cast<double>(rand())/RAND_MAX);
        input[i] = std::complex<double>(static_cast<double>(i),/*std::pow(-1.0, i +1)*/1);
        //std::cout << input[i] << "\n";
    }

    FFTGenerator FFT(input,n);

    std::cout<< std::endl;


    PrintIt(input,"Input");

    std::cout << std::endl;


    std::cout << "Starting construction of Recursive_FFT:" <<std::endl;
    TimeIt<FFTGenerator>(&FFTGenerator::Recursive_FFT, FFT , false);
    PrintIt(FFT.getRecT(),"Recursive Direct");

    std::cout << std::endl;
    std::cout << "Starting construction of Recursive_IFFT:" <<std::endl;
    TimeIt<FFTGenerator>(&FFTGenerator::Recursive_FFT, FFT, true );
    PrintIt(FFT.getRecT(),"Recursive Inverse");

    std::cout << std::endl;
    std::cout << "Starting construction of Iterative_FFT:" <<std::endl;
    TimeIt<FFTGenerator>(&FFTGenerator::Iterative_FFT, FFT , false);
    PrintIt(FFT.getIterT(),"Iterative Direct");

    std::cout << std::endl;
    std::cout << "Starting construction of Iterative_IFFT:" <<std::endl;
    TimeIt<FFTGenerator>(&FFTGenerator::Iterative_FFT, FFT, true);
    PrintIt(FFT.getRecT(),"Iterative Inverse");



    return 0;
}
