#include <iostream>
#include <vector>
#include <complex>

#include <unistd.h>

#include "fft.hpp"
#include "utilities.hpp"

int main(int argc, char **argv){
	
	
	int rec_flag = 0;
	int iter_flag = 0;
	int inv_flag = 0;
	int par_flag = 0;
	char *signal_value = NULL;
	
	int c;
	
	while ((c = getopt (argc, argv, "rtips:")) != -1)
    switch (c)
	{
      case 'r':
        rec_flag = 1;
        break;
      case 't':
        iter_flag = 1;
        break;
      case 'i':
        inv_flag = 1;
        break;
      case 'p':
        par_flag = 1;
        break;
      case 's':
        signal_value = optarg;
        break;
      case '?':
        if (optopt == 's')
          fprintf (stderr, "Option -%s requires an argument.\n", optopt);
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
	
	size_t n = 8;
	std::vector<std::complex<double>> input{1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 };

	PrintIt(input,"Input");
	FFTGenerator fft(input,n);
	
	if(rec_flag){

		fft.Recursive_FFT();

		PrintIt(fft.getRecT(),"Recursive");
	}
	if(iter_flag){

		fft.Iterative_FFT();

		PrintIt(fft.getIterT(),"Iterative");
	}
	if(par_flag){
		std::cout << "yet to be implemented..." << std::endl;
	}
	if(inv_flag){
		
		fft.Inverse_FFT();

		PrintIt(fft.getInvT(),"Inverse");
	}
	



    return 0;

}
