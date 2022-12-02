#include <iostream>
#include <vector>
#include <complex>

using namespace std;

vector<complex<double>> Recursive_FFT(vector<complex<double>> x, size_t n){
    if(n == 1){

        return x;
    }

    vector<complex<double>> y(n);

    complex<double> w(1,0);
    const complex<double> w_n(std::cos(2 * 3.14 / n),std::sin(2* 3.14/ n));

    vector<complex<double>> even(n / 2);
    vector<complex<double>> odd(n / 2);
    
    for(int i = 0; i < n/2; i++){
        even[i] = x[i*2];
        odd[i] = x[i*2+1];
    } 
    vector<complex<double>> y_o = Recursive_FFT(even,n/2);
    vector<complex<double>> y_e = Recursive_FFT(odd,n/2);

    for(int k = 0; k < n / 2; k++){
        
        y[k] = y_o[k] + w * y_e[k];
        y[k + n/2] = y_o[k] - w * y_e[k];
    
        w = w * w_n;

    }
    return y;
}

int main(){
    size_t n = 4;
    vector<complex<double>> input(n);
    for(size_t i = 0; i < n; i++){
        input[i] = (complex<double>(i, 0));
        cout << input[i] << " ";
    }
    cout << endl;

    vector<complex<double>> y = Recursive_FFT(input,n);

    for(size_t i = 0; i < n; i++){
        cout << y[i] << " ";
    }
    cout << endl;
    return 0;
}
