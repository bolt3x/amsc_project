#ifndef FFT_HPP
#define FFT_HPP

#include <iostream>
#include <vector>
#include <complex>
#include <stdlib.h>
#include <functional>
#include <random>
#include <algorithm>
#include <omp.h>

inline constexpr unsigned N_STANDARD=150;

//---------------------------------FFT CLASS------------------------------

class FFTGenerator{
    public:


    //-------------------------------CONSTRUCTORS---------------------------------

    // Default constructor: neither dimension or signal vector given
    FFTGenerator():
        m_n(N_STANDARD),
        m_signal(this->RandomGen()),
        m_rec(m_n,0.),
        m_iter(m_n,0.)
        {};

    // Partial constructor: just the dimension is given
    FFTGenerator(const unsigned &n):
        m_n(n),
        m_signal(this->RandomGen()),
        m_rec(n,0.),
        m_iter(n,0.)
        {};

    // Full constructor: both the input signal and the dimensions are given
    FFTGenerator(const std::vector<std::complex<double>> &signal, const unsigned &n):
        m_n(n),
        m_signal(signal),
        m_rec(n,0.),
        m_iter(n,0.)
        {};

    // Destructor

    ~FFTGenerator(){};

    // Declaration of FFT Recursive function

    void Recursive_FFT();

    // Declaration of FFT Iterative funcion

    void Iterative_FFT();

    // Declaration of Inverse FFT

    void Inverse_FFT();

    //-------------------------------------- SUPPORTING FUNCTIONS --------------------------------------
    // SEE SOURCE FILE FOR DETAILED EXPLANATION...

    // error function
    void GetError(std::ostream &output = std::cout)const;
    
    // random signal generator
    std::vector<std::complex<double>> RandomGen();

    // check for multiple of 2
    bool IsPower()const;

    


    // standard getters for class protected variables

    unsigned int getN() const;

    std::vector<std::complex<double>> getSignal() const;

    std::complex<double> getSignal(unsigned int const pos) const;

    std::vector<std::complex<double>> getRecT()const;

    std::complex<double> getRecT(unsigned int const pos) const;

    std::vector<std::complex<double>> getIterT()const;

    std::complex<double> getIterT(unsigned int const pos) const;

    std::vector<std::complex<double>> getInvT()const;

    std::complex<double> getInvT(unsigned int const pos) const;

    // standard setters to write the solution

    void setRecT(unsigned long const pos, std::complex<double> const val);

    void setRecT(std::vector<std::complex<double>> y);

    void setIterT(unsigned long const pos, std::complex<double> const val);

    void setIterT(std::vector<std::complex<double>> y);

    void setInvT(unsigned long const pos, std::complex<double> const val);

    void setInvT(std::vector<std::complex<double>> y);


    protected:

        // signal dimensions
        const unsigned int m_n;

        // signal values
        std::vector<std::complex<double>> m_signal;

        //FFT rec values
        std::vector<std::complex<double>> m_rec;

        //FFT iter values
        std::vector<std::complex<double>> m_iter;

        //FFT inv values
        std::vector<std::complex<double>> m_inv;
        

};



#endif
