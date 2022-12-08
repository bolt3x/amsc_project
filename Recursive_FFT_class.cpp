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