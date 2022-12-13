#include "fft.hpp"
#include "utilities.hpp"



/*-----------------------------------------------------
PARALLEL FFT function: 
-----------------------------------------------------*/

// Implememtation of parallel function:
// First of all, each avaiable processor handles n/p elements.
// We identify 3 different phases:
// ---->  I PHASE  :: swap the elements of a into y with bit reverse. To do so, each processor handles only a section of the vectors
// a and y, so each processos will have to scatter specific elemnts among specific processors. The algorithm to do
// so is reported below.


void Parallel_FFT(std::vector<std::complex<double>> &x){
    // each process calls the function and does its work...

    //first, each process must determine the his id and total number of processes

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // We first obsedrve that the numeber of pocesses to be used for the parallel vrersion of the Cooley-Tukey
    //algorithm has to be a power of 2. As such, we assess the number of processes currently involved and
    // we return in the case this number is not satisfactory...

    

    if(!((size & (size-1)) == 0)){
        if(rank == 0){
            std::cout << "Number of called processes is not correct: for the algorithm to work properly\
            the number of processes has to be a power of 2, exiting..."<< std::endl;
        }
        return;
    }

    unsigned long N(x.size()); //variable to store number of total elements, assumed to be only known to rank 0 
    
    if(!((N & (N-1)) == 0)){
        if(rank == 0){
            std::cout << "Number elements of the input signal is not a power of 2: for the algorithm to work properly\
            the number of elements has to be a power of 2, exiting.."<< std::endl;
        }
        return;
    }
    // Now, as a first step, process of rank 0 has to broadcast the number of total elments to the others

    MPI_Bcast(&N, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD); 

    // Each processor calculates the elements it has to deal with...                                                                 

    unsigned long width(N/size); // number of elements locally owned

    std::vector<std::complex<double>> x_local(width /* + (N%size > rank)*/); // segment of input singal locally owned
    std::vector<std::complex<double>> y_local(width /* + (N%size > rank)*/); // segment of output signal locally handled
                                                                         // the dimensions are calculated on the assumpion the number of elements
                                                                         // are multiple of p. This assumpion simpifies the work

    // Note that rank 0 still has to share the elements of x between all the other processes

    // Now, we present two different approaches:
    //      1) : we can let process 0 alone handle the permutation of the elements of the signal (phase I)
    //          and then scatter this elements among the different processes. By doing so, rank 0 is doing
    //          more work as it is compelled to compute all the changes in position by itself, but, on the other
    //          hand the overall communication cost is realtively small...
    //      2) : we can let process 0 broadcast the elements of the signal to all the different processes,
    //           then each one of then broadcast in the right position the element of the singal into the
    //           ouput vector, by bit reversal. This leads to a reduced computational job for rank 0
    //           but the overall communication cost is increased instead
    //
    // We choose the first approach

    
    if(rank == 0){

        for (unsigned int i = 1, j = 0; i < N; i++) {
            unsigned bit = N >> 1;
            for (; j & bit; bit >>= 1)
                j ^= bit;
            j ^= bit;

            if (i < j)
                swap(x[i], x[j]);
        }

    }

    MPI_Scatter(x.data(),width , MPI_DOUBLE_COMPLEX,
               x_local.data(), width , MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    



    // Once the elements are redistributed among the processes, everyone of them can start the second phase
    //  ----> 2 PHASE : log2(n) - log2(p) iterations of the FFT algorithm without communication: simple additions
    //                  and subtruction on complex phases.

    for(size_t j =1; j <= std::log2(n)-std::log2(size); ++j){
        unsigned int d = 2<<(j-1); // size
        unsigned int d2 = d >> 1; // m2 = m/2
        // principle root of nth complex root of unity.
        std::complex<double> w_d(std::polar(1.0, (inverse*4 - 2)*std::numbers::pi/d));
        std::complex<double> w(1.0,0.0);

        // k :: number of iterations representing the number of "consecutive combining pattern"
        
        for(size_t k = 0; k < d2; ++k) {
            // m :: number of iteration representing the number of total of "combining patterns" of a certain type 
            for(size_t m = k; m < n; m += d) {
                std::complex<double> t = w * y_local[m + d2];
                std::complex<double> u = y_local[m];

                // similar calculating y[m]
                y_local[m] = (u + t);

                // similar calculating ym+n/2]
                y_local[m + d2] = (u - t);
            }
            w *= w_d;
        }

    }

    //Then the third phase starts:
    //  ----> 3 PHASE : log2(p) iterations durig which every process swaps copies of his values
    //                  with a process adiacient across some dimension of the hypercube
    //                  Note that in this representation, then number of processes p would represent
    //                  the dimension of the hypercube ==> p = 2^dim

    for(size_t j =std::log2(n) - std::log2(size) + 1 ; j <= log2(size); ++j){

        //the dimension across which the copies are swapped between processes are dependent
        //form the iterations of the first loop

        std::vector<std::complex<double>> y_adiacent(y_local);
        unsigned int ad(rank);
        for(int i = 0; i< size; i++){
            //Process sends the elements of y_local to the other process by first coping them in y_adiacent
            MPI_Send(y_adiacent.data(), size, MPI_COMPLEX_DOUBLE, i + ad, 1, MPI_COMM_WORLD);
        }

        for(int i = 0; i< size; i++){
            //Then receives the elements for the same process in y_adiacent
            MPI_Recv(y_adiacent.data(), size, MPI_COMPLEX_DOUBLE, i + ad, 1, MPI_COMM_WORLD);
        }

        unsigned int d = 2<<(j-1); // size
        unsigned int d2 = d >> 1; // m2 = m/2
        // principle root of nth complex root of unity.
        std::complex<double> w_d(std::polar(1.0, (inverse*4 - 2)*std::numbers::pi/d));
        std::complex<double> w(1.0,0.0);

        // k :: number of iterations representing the number of "consecutive combining pattern"
        
        for(size_t k = 0; k < d2; ++k) {
            // m :: number of iteration representing the number of total of "combining patterns" of a certain type 
            for(size_t m = k; m < n; m += d) {
                

                std::complex<double> t = w * y_local[m + d2];
                std::complex<double> u = y_local[m];

                // similar calculating y[m]
                y_local[m] = (u + t);

                // similar calculating ym+n/2]
                y_local[m + d2] = (u - t);
            }
            w *= w_d;
        }

    }



    MPI_Barrier(MPI_COMM_WORLD);

    int mpi_rank = 0;
    
    while(mpi_rank < size){
        if(mpi_rank == rank){
            std::cout << "----------------" << std::endl;
            std::cout<<"process #"<< rank << " elements" << std::endl;
            for(auto i : x_local){
                std::cout << " " << i << " " <<std::endl;
            }
        }

        mpi_rank++;
        MPI_Barrier(MPI_COMM_WORLD);
    }

    return;
}
