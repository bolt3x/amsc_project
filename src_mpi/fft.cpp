#include "fft.hpp"
#include "utilities.hpp"





/*-----------------------------------------------------
MPI_FFT function: a second possible parallel implementation of the FFT done with MPI
-----------------------------------------------------*/

// Implememtation of parallel function:
// First of all, each avaiable processor handles n/p elements.
// We identify 3 different phases, and for each of those phases we calculate an estimate of the time required to 
// run them, as to assess the possible bottlenecks
// ---->  I PHASE  :: swap the elements of x into y with bit reverse.


void MPI_FFT(std::vector<std::complex<double>> &x){
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
            std::cout << "Number of called processes is not correct: for the algorithm to work properly the number of processes has to be a power of 2. Exiting..."<< std::endl;
        }
        return;
    }

    // We measure the time required for the boradcast of n and the sorting
    auto start_tot = std::chrono::high_resolution_clock::now();

    unsigned long N(x.size()); //variable to store number of total elements, assumed to be only known to rank 0 
    
    
    // Now, as a first step, process of rank 0 has to broadcast the number of total elments to the others

    MPI_Bcast(&N, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD); 


    if(!((N & (N-1)) == 0)){
        if(rank == 0)
        {
            std::cout << "Number elements of the input signal is not a power of 2: for the algorithm to work properly the number of elements has to be a power of 2. Eexiting.."<< std::endl;
        }
        return;
    }

    if((int)N < size){
        if(rank == 0)
        {
            std::cout << "Number signal elements have to be bigger than size. Exiting "<< std::endl;
        }
        return;
    }

    // Each processor calculates the elements it has to deal with...                                                                 

    unsigned long width(N/size); // number of elements locally owned

    //std::vector<std::complex<double>> x_local(width /* + (N%size > rank)*/); // segment of input singal locally owned
    std::vector<std::complex<double>> y_local(width /* + (N%size > rank)*/); // segment of output signal locally handled
                                                                         // the dimensions are calculated on the assumpion the number of elements
                                                                         // are multiple of p. This assumpion simpifies the work

    // Note that rank 0 still has to share the elements of x between all the other processes

    // Now, we present different approaches:
    //      1) : we can let process 0 alone handle the permutation of the elements of the signal (phase I)
    //          and then scatter this elements among the different processes. By doing so, rank 0 is doing
    //          more work as it is compelled to compute all the changes in position by itself, but, on the other
    //          hand the overall communication cost is realtively small...
    //      2) : we can let process 0 broadcast the elements of the signal to all the different processes,
    //           then each one of them broadcast in the right position the element of the singal into the
    //           ouput vector, by bit reversal. This leads to a reduced computational job for rank 0
    //           but the overall communication cost is increased instead
    //      3) : with point to point comunication rank 0 sends the element required by a process to the process itself
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

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start_tot);

    if(rank ==0)
        std::cout << "Time required to broadcast the dimension, create temporary variables and sort the elements by process 0 as per process " << rank <<": "<<duration.count() << " microseconds."<< std::endl;


    // Now we measure the time required for the scatter
    auto start = std::chrono::high_resolution_clock::now();


    MPI_Scatter(x.data(),width , MPI_DOUBLE_COMPLEX,
               y_local.data(), width , MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);



    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    if(rank ==0)
        std::cout << "Time required to scatter the elements as per process " << rank <<": "<<duration.count() << " microseconds."<< std::endl;

    

    // Once the elements are redistributed among the processes, everyone of them can start the second phase
    //  ----> 2 PHASE : log2(n) - log2(p) iterations of the FFT algorithm without communication: simple additions
    //                  and subtruction on complex phases.

    //Now time phase 2

    start = std::chrono::high_resolution_clock::now();

    for(size_t j =1; j <= std::log2(N)-std::log2(size); ++j){
        unsigned int d = 2<<(j-1); // size
        unsigned int d2 = d >> 1; // m2 = m/2
        // principle root of nth complex root of unity.
        std::complex<double> w_d(std::polar(1.0, - 2*M_PI/d));
        std::complex<double> w(1.0,0.0);

        // k :: number of iterations representing the number of "consecutive combining pattern"
        
        for(size_t k = 0; k < d2; ++k) {
            // m :: number of iteration representing the number of total of "combining patterns" of a certain type 
            for(size_t m = k; m < width; m += d) {
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

    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    if(rank ==0)
        std::cout << "Time required to prerform the first log2(n)-log2(p) iterations of butterfly pattern as per process " << rank <<": "<<duration.count() << " microseconds."<< std::endl;


    

    //Then the third phase starts:
    //  ----> 3 PHASE : log2(p) iterations durig which every process swaps copies of his values
    //                  with a process adiacient across some dimension of the hypercube
    //                  Note that in this representation, the number of processes p would represent
    //                  the dimension of the hypercube ==> p = 2^dim



    // Finally, we measure the time necessary to perfom the last phase
    start = std::chrono::high_resolution_clock::now();


    for(size_t iter =1 ; iter <= std::log2(size); ++iter){

        //the dimension across which the copies are swapped between processes are dependent
        //form the iterations of the first loop

        std::vector<std::complex<double>> y_adiacent(y_local);
        // the receiver of the data is determined by the rank and the iteration of the
        //outer loop

        // ad is the rank of the adiacent process
        const int ad = rank + (((rank/iter)%2) *(-2) +1)*iter;
        unsigned int d = 2<<(iter + static_cast<unsigned int>(std::log2(N)- std::log2(size))-1);
        
        
        //Process sends the elements of y_local to the adiacent process by first coping them in y_adiacent
        MPI_Send(y_adiacent.data(), width, MPI_DOUBLE_COMPLEX, ad, 0, MPI_COMM_WORLD);

        //Then receives the elements for the same process in y_adiacent
        MPI_Recv(y_adiacent.data(), width, MPI_DOUBLE_COMPLEX, ad, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // // principle root of nth complex root of unity.
        std::complex<double> w_d(std::polar(1.0,-2*M_PI/d));
        std::complex<double> w(1.0,0.0);

        // k :: number of iterations representing the number of "consecutive combining pattern"
        
        for(size_t k = 0; k < width; ++k) {
            // m :: number of iteration representing the number of total of "combining patterns" of a certain type 
                

                std::complex<double> t = ((1.0 - w) * static_cast<double>(rank < ad) + w) * y_local[k];
                std::complex<double> u = ((1.0 - w) * static_cast<double>(rank > ad) + w) * y_adiacent[k];

                // Wether to add or to subtract t and u dipends on the specific
                // processor and itereration of the outer loop
                // similar calculating y[m]
                y_local[k] = (u + (2.0*static_cast<double>(rank < ad) - 1.0)*t);
                w *= w_d;

        }

        

        
        

    }

    // Once finished the third phase, all processes send their local data to rank 0 using a Gather

    MPI_Gather(y_local.data(), width, MPI_DOUBLE_COMPLEX, x.data(),width,MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    auto stop_tot = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop_tot - start);
    auto duration_tot = std::chrono::duration_cast<std::chrono::microseconds>(stop_tot - start_tot);

    if(rank ==0){
        std::cout << "Time required to perform the last log2(p) iterations of the butterfly pattern as per process" << rank <<": "<<duration.count() << " microseconds.\n"<< std::endl;
        std::cout << "Total time required as per process "<< rank<<" is :"<<duration_tot.count() << " microseconds."<< std::endl;
    }
        


    //--------DEBUGGING SECTION-------

    /*
    MPI_Barrier(MPI_COMM_WORLD);

    int mpi_rank = 0;
    
    while(mpi_rank < size){
        if(mpi_rank == rank){
            std::cout << "----------------" << std::endl;
            std::cout<<"process #"<< rank << " elements" << std::endl;
            for(auto i : y_local){
                std::cout << " " << i << " " <<std::endl;
            }
        }

        mpi_rank++;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    */



    //--------DEBUGGING SECTION-------

    

    return;
}
