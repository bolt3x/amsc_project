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
    
    //temp variables...

    unsigned long n = this->getN();
    std::vector<std::complex<double>> x(this->getSignal());

    if(n == 1){
        this->setRecT(x);
        return;
    }

    std::complex<double> w (1.0,0.0);
    const std::complex<double> w_n(std::polar(1.0,-2*M_PI/(this)->getN()));

    std::vector<std::complex<double>> even (n/2);
    std::vector<std::complex<double>> odd  (n/2);

    for(size_t i = 0; i < n/2 ; ++i){
        even[i] = x[i*2];
        odd[i] = x[i*2 +1];
    }

    FFTGenerator FFT_0(even,n/2);
    FFTGenerator FFT_1(odd,n/2);

    FFT_0.Recursive_FFT();
    FFT_1.Recursive_FFT();

    /*std::vector<std::complex<double>>y_0(n/2);
    std::vector<std::complex<double>>y_1(n/2);
    Recursive_FFT(even, n/2, y_0);
    Recursive_FFT(odd, n/2, y_1);*/

    for(unsigned long k = 0; k < n/2; k++){
        //y[k] = FFT_0.getRecT(k) + w * FFT_1.getRecT(k);
        this->setRecT(k, FFT_0.getRecT(k) + w * FFT_1.getRecT(k));
        //y[k + n/2] = FFT_0.getRecT(k) - w * FFT_0.getRecT(k);
        this->setRecT(k + n/2, FFT_0.getRecT(k) - w * FFT_1.getRecT(k));
        w = w * w_n;
    }

    //this->setRecT(y);

    return;
}

/*----------------------------------------------------
ITERATIVE FFT :  basic iterative version of the FFT.
                        * doesn't create new vectres, simply updates the values of the existing ones
                        * serial version
                        * performs n operattions, log(n) times
                        * time complexity: T(n) = O(n*log(n))
                        * needs a reverse bit function
------------------------------------------------------*/



void FFTGenerator::Iterative_FFT(){

    // temp variables...

    unsigned long n = this->getN();
    std::vector<std::complex<double>> x(this->getSignal());
    std::vector<std::complex<double>> y(n);
    

    //CHECK ON THE DIMENSION OF N : METHOD IS VALID ONLY FOR DIMENSIONS MULTIPLE OF 2

    if(!((n & (n-1)) == 0)){
        //the dimension is not power of 2
        // return and notify the user
        std::cout << "Dimension of the input signal is not a power of 2:\n as such the Iterative_FFT is not adequate to calculate the FFT "<< std::endl;
        return;
    }

    for(size_t i = 0; i < n; ++i)
        y[ReverseBit(i,n)] = x[i];

    for(size_t j =1; j <= std::log2(n); ++j){
        unsigned int d = 2<<(j-1); // size
        unsigned int d2 = d >> 1; // m2 = m/2
        // principle root of nth complex root of unity.
        std::complex<double> w_d(std::polar(1.0, -2*M_PI/d));
        std::complex<double> w(1.0,0.0);


        
        for(size_t k = 0; k < d2; ++k) {
            for(size_t m = k; m < n; m += d) {
                std::complex<double> t = w * y[m + d2];
                std::complex<double> u = y[m];

                // similar calculating y[m]
                y[m] = u + t;

                // similar calculating ym+n/2]
                y[m + d2] = u - t;
            }
            w *= w_d;
        }
    }

    this->setIterT(y);

    return;

}

/*----------------------------------------------------
INVERSE FFT :  basic iterative version of the IFFT.
                        * uses previous iterative FFT
                        * serial version
------------------------------------------------------*/

void FFTGenerator::Inverse_FFT(){

   // temp variables...

    unsigned long n = this->getN();
    std::vector<std::complex<double>> x(this->getSignal());
    std::vector<std::complex<double>> y(n);

    //CHECK ON THE DIMENSION OF N : METHOD IS VALID ONLY FOR DIMENSIONS MULTIPLE OF 2

    if(!((n & (n-1)) == 0)){
        //the dimension is not power of 2
        // return and notify the user
        std::cout << "Dimension of the input signal is not a power of 2:\n as such the Iterative_FFT is not adequate to calculate the FFT "<< std::endl;
        return;
    }

    //compute the conjugate of the values
    std::for_each(x.begin(),x.end(),[n](auto &elem){elem = std::conj(elem);});

    for(size_t i = 0; i < n; ++i)
        y[ReverseBit(i,n)] = x[i];

    for(size_t j =1; j <= std::log2(n); ++j){
        unsigned int d = 2<<(j-1); // size
        unsigned int d2 = d >> 1; // m2 = m/2
        // principle root of nth complex root of unity.
        std::complex<double> w_d(std::polar(1.0, -2*M_PI/d));
        std::complex<double> w(1.0,0.0);


        
        for(size_t k = 0; k < d2; ++k) {
            for(size_t m = k; m < n; m += d) {
                std::complex<double> t = w * y[m + d2];
                std::complex<double> u = y[m];

                // similar calculating y[m]
                y[m] = u + t;

                // similar calculating ym+n/2]
                y[m + d2] = u - t;
            }
            w *= w_d;
        }
    }

    //compute the conjugate of the values and scale them
    std::for_each(y.begin(),y.end(),[n](auto &elem){elem = std::conj(elem) / std::complex<double>{n,0};});

    this->setInvT(y);

    return;

}

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


/*-----------------------------------------------------
RANDOM GEN function: give a specified dimension, generate a signal of random positive numbers
-----------------------------------------------------*/

 std::vector<std::complex<double>> FFTGenerator::RandomGen()
{
    unsigned long n = this->getN();
    std::mt19937 gen32; // MersenneTWISTER
    std::vector<std::complex<double>> signal(n);
    for(size_t i = 0; i< n; ++i){
        signal[i]=(static_cast<double>(gen32()%10)/*module*/*std::polar(1.0,2*M_PI / (gen32()%gen32())))/*random complex number of module 1*/;
    }

    return signal;

}


/*----------------------------------------------------
IS Power funcition : to check whether a number is a Power of 2
------------------------------------------------------*/

bool FFTGenerator::IsPower()const
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

// element-of-InvT getter...

std::complex<double> FFTGenerator::getIterT(unsigned int const pos) const{
    return m_iter[pos];
}

std::vector<std::complex<double>> FFTGenerator::getInvT()const{
    return m_inv;
}

// element-of-InvT getter...

std::complex<double> FFTGenerator::getInvT(unsigned int const pos) const{
    return m_inv[pos];
}


//----------------------------STANDARD SETTERS----------------------------

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

// by-element Inverse setter

void FFTGenerator::setInvT(unsigned long const pos, std::complex<double> const val){
    m_inv[pos] = val;
}

//global Inverse setter

void FFTGenerator::setInvT(std::vector<std::complex<double>> y){
    m_inv = y;
}