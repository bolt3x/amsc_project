#include <iostream>
#include <vector>
#include <complex>
#include <string>

#include <unistd.h>

#include "fft.hpp"
#include "utilities.hpp"

int main(int argc, char **argv){


  char *signal_file = NULL;
	
	int c;
	
	while ((c = getopt (argc, argv, "rtips:")) != -1)
    switch (c)
	{
      case 's':
        signal_file = optarg;
        break;
      case '?':
        if (optopt == 's')
          fprintf (stderr, "Option -%s requires an argument.\n",optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%s'.\n", optopt);
        else
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        return 1;
      default:
        abort ();
	}

  unsigned int print(0);
	
	
	// initializing MPI
    MPI_Init(nullptr, nullptr);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<std::complex<double>> input(0);
    size_t n(0);

    if(rank == 0){

      if(!signal_file){
        n = 256;
        input.resize(n);
        input = RandomGen(n);
      }

      else {
        std::string s(signal_file);
        ReadIt(s,input,n);
      }
      
      if(print)
        PrintIt(input,"Input");

    }

    MPI_FFT(input);

    if(print)
      PrintIt(input,"MPI -- Please note the implementation is running on processors specified when lauched");

    MPI_Finalize();



    return 0;

}
