#include "fft.hpp"
#include "utilities.hpp"

/*----------------------------------------------------
NAIVE RECURSIVE FFT :  basic, less-performant recursive version of the FFT.
                        * creates 4 new vectors of decreasing sizes for each recursion
                        * serial version
                        * performs n operations for every recursion and n is halved each recursion
                        * time complexity: T(n) = 2T(n/2) * O(n) = O(n*log(n))
------------------------------------------------------*/


void FFTGenerator::Recursive_FFT(){
    unsigned long N = this->getN();
    std::vector<std::complex<double>> x(this->getSignal());
    
    if(N==1){
        this->setRecT(x);
        return;
    }

    // create vectors with even and odd coefficients
    std::vector<complex<double>> x_even(N/2), x_odd(N/2);

    for(size_t it=0; it*2<N; ++it){
        x_even[it] = x[it*2];
        x_odd[it] = x[it*2+1];
    }

    // create class istances for recursive calls
    FFTGenerator FFT_even(x_even, N/2);
    FFTGenerator FFT_odd(x_odd, N/2);

    // recursive calls
    FFT_even.Recursive_FFT();
    FFT_odd.Recursive_FFT();

    // merging results from recursive calls
    for(size_t k=0; k*2<N; ++k){
      std::complex<double> Wn = std::polar(1.0, -2 * std::numbers::pi * k/N);
      this->setRecT(k, FFT_even.getRecT(k) + (Wn * FFT_odd.getRecT(k)));
      this->setRecT(k + N/2, FFT_even.getRecT(k) - (Wn * FFT_odd.getRecT(k)));
    }

    return;
 }

// OUT-OF-CLASS DEFINITION...

/*

void FFTGenerator::Recursive_FFT(const std::vector<std::complex<double>> &x //signal values, 
                   const unsigned &n//signal dimension,
                   std::vector<std::complex<double>> &y//discrete Fourier transform value
                   //WARNING:: we hearby implement the function without the reference operator as to use
                   //the method in synbiosis with the reference getter method){
    
    NEED TO IMPLEMENT A CHECK ON THE DIMENSION OF N : METHOD IS VALID ONLY FOR DIMENSIONS MULTIPLE OF 2

    if(n == 1){
        y[0] = x[0];
        return;
    }

    std::complex<double> w (1.0,0.0);
    const std::complex<double> w_n(std::polar(1.0,2*std::numbers::pi/n));

    std::vector<std::complex<double>> even (n/2);
    std::vector<std::complex<double>> odd  (n/2);

    for(size_t i = 0; i < n/2; ++i){
        even[i] = x[i*2];
        odd[i] = x[i*2 +1];
    }

    std::vector<std::complex<double>>y_0(n/2);
    std::vector<std::complex<double>>y_1(n/2);
    Recursive_FFT(even, n/2, y_0);
    Recursive_FFT(odd, n/2, y_1);

    for(size_t k = 0; k < y_0.size(); k++){
        y[k] = y_0[k] + w * y_1[k];
        y[k + n/2] = y_0[k] - w * y_1[k];
        w = w * w_n;

    }
    return y;
}

*/

/*----------------------------------------------------
ITERATIVE FFT :  basic iterative version of the FFT.
                        * doesn't create new vectores, simply updates the values of the existing ones
                        * serial version
                        * performs n operattions, log(n) times
                        * time complexity: T(n) = O(n*log(n))
                        * need a reverse bit function
------------------------------------------------------*/



void FFTGenerator::Iterative_FFT(){

    // temp variables...

    unsigned long n = this->getN();
    std::vector<std::complex<double>> x(this->getSignal());
    std::vector<std::complex<double>> y(n);

    //NEED TO IMPLEMENT A CHECK ON THE DIMENSION OF N : METHOD IS VALID ONLY FOR DIMENSIONS MULTIPLE OF 2

    for(size_t i = 0; i < n; ++i)
        y[ReverseBit(i,n)] = x[i];

    for(size_t j = 1; j <= std::log2(n); ++j){
        unsigned int d = 1<<j;
        std::complex<double> w_d(std::polar(1.0, 2*std::numbers::pi/d));
        std::complex<double> w(1.0,0.0);
        for(size_t k = 0; k < d/2-1; ++k){
            for(size_t m = k; m < n; m += d){
                std::complex<double> t = w*y[m + d/2];
                std::complex<double> a = y[k];
                y[k] = a + t;
                y[k + d/2] = a - t;
            }
        }
        w *= w_d;
    }

    this->setIterT(y);

    return;

}

// OUT-OF-CLASS DEFINITION...

/*

void FFTGenerator::Iterative_FFT(const std::vector<std::complex<double>> &x //signal values, 
                   const unsigned &n//signal dimension,
                   std::vector<std::complex<double>> &y//discrete Fourier transform value
                   //WARNING:: we hearby implement the function without the reference operator as to use
                   //the method in synbiosis with the reference getter method){

    
    //NEED TO IMPLEMENT A CHECK ON THE DIMENSION OF N : METHOD IS VALID ONLY FOR DIMENSIONS MULTIPLE OF 2

    for(size_t i = 0; i < x.size(); ++i)
        y[//bit_reverse of i
          ReverseBit(i, //number of bits required to represent the max index
          n) ] = x[i];

    for(size_t j = 1; j <= std::log2(n); ++j){
        unsigned int d = 1<<j;
        std::complex<double> w_d(std::polar(1.0, 2*std::numbers::pi/d));
        std::complex<double> w(1.0,0.0);
        for(size_t k = 0; k < d/2-1; ++k){
            for(size_t m = k; m < n; m += d){
                std::complex<double> t = w*y[m + d/2];
                std::complex<double> a = y[k];
                y[k] = a + t;
                y[k + d/2] = a - t;
            }
        }
        w *= w_d;
    }

    return;

}

*/



//-------------------------------------- SUPPORTING FUNCTIONS --------------------------------------

/*-----------------------------------------------------
GET ERROR funciton: to compare two diffennt FFT implementations and 
return discrepancy between the two:
                    * first we determine the difference vector between the two
                    compolex-valued FFTs, compared at a given k index
                    * then we compute the norm for finite dimensional colpex spaces
                     of such quantities: ||z|| = sqrt(abs(z_1)^2 + ... +abs(z_n)^2)
-----------------------------------------------------*/

void FFTGenerator::GetError(std::ostream &output)const{

    // by construction the to FFT to be compared have the same dimension as
    // n is a class member (dimension parameter) common to all vectors of the classe
    double norm_err(0);
    output << "Norm_2 of error -- RECURSIVE vs. ITERATIVE:  ";
    for(size_t i = 0; i< this->getN(); ++i){
         norm_err += std::pow(std::abs(this->getRecT(i) - this->getIterT(i)),2);
    }
    output << std::sqrt(norm_err) << std::endl;
}


// OUT-OF-CLASS DEFINITION...

/*
void FFTGenerator::GetError(const std::vector<std::complex<double>> &ExactFFT,
              const std::vector<std::complex<double>> &EmpiricalFFT,
              std::ostream &output)const{

    if(ExactFFT.size() != EmpiricalFFT.size()){
        output << "Error: mismatching size of the two FFTs"<< std::endl;
        return;
    }
    double norm_err(0);
    output << "Norm_2 of error:  ";
    for(size_t i = 0; i< ExactFFT.size(); ++i){
         norm_err += std::pow(std::abs(ExactFFT[i] - EmpiricalFFT[i]),2);
    }
    output << std::sqrt(norm_err) << std::endl;
}

*/


/*-----------------------------------------------------
RANDOM GEN function: give a specified dimension, generate a signal of random positive numbers
-----------------------------------------------------*/

 std::vector<std::complex<double>> FFTGenerator::RandomGen()
{
    unsigned long n = this->getN();
    std::mt19937 gen32; // MersenneTWISTER
    std::vector<std::complex<double>> signal(n);
    for(size_t i = 0; i< n; ++i){
        signal[i]=(static_cast<double>(gen32()%10)/*module*/*std::polar(1.0,2*std::numbers::pi / (gen32()%gen32())))/*random complex number of module 1*/;
    }

    return signal;

}


/*----------------------------------------------------
IS MULTIPLE funcition : to check whether a number is a multiple of 2
------------------------------------------------------*/

bool FFTGenerator::IsMultiple()const
{

    unsigned long n = this->getN();
    while(n){
        if(n%2 == 1) return false;
        n <<= 2;
    }
    return true;
}


//----------------------------STANDARD GETTERS----------------------------

// dimesion getter...

unsigned int FFTGenerator::getN() const{
    return m_n;
}

// input signal getter...

std::vector<std::complex<double>> FFTGenerator::getSignal() const{
    return m_signal;
}

// element-of-input getter...

std::complex<double> FFTGenerator::getSignal(unsigned int const pos) const{
    return m_signal[pos];
}

//rec. tranform getter...

std::vector<std::complex<double>> FFTGenerator::getRecT()const{
    return m_rec;
}

// element-of-RecT getter...

std::complex<double> FFTGenerator::getRecT(unsigned int const pos) const{
    return m_rec[pos];
}


//iter. tranform getter...

std::vector<std::complex<double>> FFTGenerator::getIterT()const{
    return m_iter;
}

// element-of-IterT getter...

std::complex<double> FFTGenerator::getIterT(unsigned int const pos) const{
    return m_iter[pos];
}

//----------------------------STANDARD GETTERS----------------------------

// by-element Recursive setter

void FFTGenerator::setRecT(unsigned long const pos, std::complex<double> const val){
    m_rec[pos] = val;
}

// global Recursive setter

void FFTGenerator::setRecT(std::vector<std::complex<double>> y){
    m_rec = y;
}

// by-element Iterative setter

void FFTGenerator::setIterT(unsigned long const pos, std::complex<double> const val){
    m_iter[pos] = val;
}

//global Iterative setter

void FFTGenerator::setIterT(std::vector<std::complex<double>> y){
    m_iter = y;
}