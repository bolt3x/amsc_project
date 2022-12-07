#include <iostream>
#include <vector>
#include <complex>
#include <numbers>

using namespace std;

void Recursive_FFT(vector<complex<double>>& A, int N)
{
    if(N==1) return;

    const double pi = acos(-1);
    const complex<double> i(0, 1);

    // create vectors with even and odd coefficients
    vector<complex<double>> A_0(N/2), A_1(N/2);

    for(size_t it=0; it*2<N; ++it){
        A_0[it]=A[it*2];
        A_1[it]=A[it*2+1];
    }

    // recursive call
    Recursive_FFT(A_0, N/2);
    Recursive_FFT(A_1, N/2);

    for(size_t k=0; k*2<N; ++k){
      complex<double> Wn = std::polar(1.0, 2 * pi * k/N);
      A[k]=A_0[k]+(Wn*A_1[k]);
      A[k+N/2]=A_0[k]-(Wn*A_1[k]);
    }

    return;
}

int main(){
  size_t n = 4;
  vector<complex<double>> input(n);
  for(size_t i = 1; i <= n; i++){
    input[i-1] = (complex<double>(i, 1));
    cout << input[i-1] << " ";
  }
  cout << endl;

  
  Recursive_FFT(input,n);

  for(size_t i = 0; i < n; i++){
    cout << input[i] << " ";
  }
  cout << endl;
  return 0;
}
