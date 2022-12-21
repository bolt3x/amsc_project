#include <iostream>
#include <vector>
#include <complex>
#include <string>

#include <unistd.h>

#include "fft.hpp"
#include "utilities.hpp"

int main(int argc, char **argv){
	
	
	int rec_flag = 0;
	int iter_flag = 0;
	int inv_flag = 0;
	int OMP_flag = 0;
	char *signal_file = NULL;
	
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
        OMP_flag = 1;
        break;
      case 's':
        signal_file = optarg;
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
	


  std::vector<std::complex<double>> input;
  size_t n;
  //--------------------

    if(!signal_file){
      n = 32;
      input.resize(n);
      input = RandomGen(n);
    }

    else {
      std::string s(signal_file);
      ReadIt(s,input,n);
    }

    PrintIt(input,"Input");
    FFTGenerator fft(input,n);

  //-------------------- 
  
	
	
	if(rec_flag){ 

  

    TimeIt<FFTGenerator>(&FFTGenerator::Recursive_FFT, fft);
    
		//fft.Recursive_FFT();

		PrintIt(fft.getRecT(),"Recursive");
	}
	else if(iter_flag){

    TimeIt<FFTGenerator>(&FFTGenerator::Iterative_FFT, fft);

		//fft.Iterative_FFT();

		PrintIt(fft.getIterT(),"Iterative");
	}
	else if(OMP_flag){

   

    std::cout << "Select number of threads:"<< std::endl;

    int num_threads;

    std::cin >>num_threads;

    TimeIt<FFTGenerator>(&FFTGenerator::OPENMP_FFT, fft, num_threads);
		
    //fft.OPENMP_FFT(num_threads);

    PrintIt(fft.getOpenMPT(),"OpenMp -- Please note the number of default threads is 2");
    
	}
	else if(inv_flag){

    TimeIt<FFTGenerator>(&FFTGenerator::Inverse_FFT, fft);
		
		//fft.Inverse_FFT();

		PrintIt(fft.getInvT(),"Inverse");
	}
  return 0;

}
