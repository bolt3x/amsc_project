#include <iostream>

#include "fft.hpp"
#include "utilities.hpp"



int main(int argc, char ** argv){

    // initializing MPI
    MPI_Init(nullptr, nullptr);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(argc != 2){
        if(rank == 0)
        {
            std::cout << "Wrong number of parameters to start the Parallel FFT"<< std::endl;
        }

        MPI_Finalize();

        return 0;
    }

    unsigned int n = std::stoi(argv[1]); /*power of 2*/;
    // the number of elements chosen is already known to all the processors

    
    // after the construction is done, each process calls the parallel function...

    std::vector<std::complex<double>> x = RandomGen(n);
    
    if(rank == 0){
        std::cout << "Process "<< rank << " says:"<< std::endl;
        PrintIt(x, "Input signal specifically generated");
    }

    Parallel_FFT(x);

    if(rank == 0)
        PrintIt(x,"FFT computed using MPI");

    //TimeIt(&Parallel_FFT, x);
    
    MPI_Finalize();


    return 0;
}
