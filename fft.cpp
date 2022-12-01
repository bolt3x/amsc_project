#include <iostream>
#include <vector>
#include <complex>

using namespace std;

void Recursive_FFT(vector<complex<double>> a, size_t n,vector<complex<double>>& y){
    if(n == 1){
        y = a;
        cout << y[0] << std::endl;
        return;
    }

    complex<double> w = 1;
    const complex<double> w_n (std::cos(2 * 3.14 / n),std::sin(2* 3.14/ n));

    vector<complex<double>> even(n / 2);
    vector<complex<double>> odd(n / 2);
    for(size_t i = 0; i < n/2; i++){
        even[i] = a[i*2];
        odd[i] = a[i*2+1];
    } 
    vector<complex<double>> y_o(n/2);
    Recursive_FFT(even,n/2,y_o);
    vector<complex<double>> y_e(n/2);
    Recursive_FFT(odd,n/2,y_e);

    for(size_t k = 0; k < n / 2 - 1; k++){
        y[k] = y_o[k] + w * y_e[k];
        y[k + n/2] = y_o[k] - w * y_e[k];
        w = w * w_n;
    }

    return;
}

int main(){
    size_t n = 10;
    vector<complex<double>> input(n);
    for(size_t i = 0; i < n; i++){
        input[i] = (complex<double>(1, i));
        cout << input[i] << " ";
    }
    cout << endl;

    vector<complex<double>> y(n); 
    Recursive_FFT(input,n,y);

    for(size_t i = 0; i < n; i++){
        cout << y[i] << " ";
    }
    cout << endl;
    return 0;
}
