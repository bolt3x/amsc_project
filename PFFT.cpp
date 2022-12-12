#include <iostream>
#include <vector>
#include <complex>
#include <numbers>
#include <tuple>
#include <mpi.h>

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

int mod(int a , int b){
  int res=a%b;
  if(res<0) res+=b;
  return res;
}

bool isprime(int x)
{
  if (x == 0 || x == 1 || (x % 2 == 0 && x > 2)) return false;
  if(x<=3) return true;
  int n=floor(sqrt(x));
  for(size_t i=5; i<n; ++i){
    if(x%i==0) return false;
  }
  return true;
}

int euclidean_extended(int a, int b, int& u, int& v)
{
  // given a and b, returns u and v such that ua+vb=1 (final a1 is the gcd(a,b))
  u = 1, v = 0;
  int x = 0, y = 1, a1 = a, b1 = b;
  while (b1) {
    int q = a1 / b1;
    tie(u, x) = make_tuple(x, u - q * x);
    tie(v, y) = make_tuple(y, v - q * y);
    tie(a1, b1) = make_tuple(b1, a1 - q * b1);
  }
  return a1;
}

void coprime_factorization(int N, int& R, int& S, int& Sr, int& Rs)
{
  R=1;
  S=N;
  int n=floor(sqrt(N));
  if(isprime(N) || isprime(n)) return; // errore o warning? "Algorithm will run serially"

  for(; R<n; ++R){
    if(N%R==0 && euclidean_extended(R, N/R, Sr, Rs)==1){
      S=N/R;
      return;
    }
  }
  return;
}

vector<complex<double>> Parallel_FFT(const vector<complex<double>>& A, size_t N)
{
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const double pi = acos(-1);
  const complex<double> i(0, 1);

  vector<complex<double>> A_vect(N);
  int R, S, Rs, Sr;

  if(!rank){
    // compute R and S as relatively prime, as well as Rs and Sr
    coprime_factorization(N, R, S, Sr, Rs); 
  
    // convert A(k) in A(k_0,k_1)
    size_t k_0, k_1; //redefinition! fuori da ogni if
    for(size_t k=0; k<N; ++k){
      k_0=mod(k*Rs,S);
      k_1=mod(k*Sr,R);
      A_vect[k_1*S+k_0]=A[k]; 
    }
  }
  
  // send data to each processor
  MPI_Bcast(A_vect.data(), N, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(&R, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
  MPI_Bcast(&S, 1, MPI_INTEGER, 0, MPI_COMM_WORLD); 

  vector<complex<double>> P(S);
  
  if(rank<R){ 
    complex<double> w;

    // compute partial results (S R-transforms)
    for(size_t k_1=0; k_1<R; ++k_1){
      w = std::polar(1.0, 2 * pi * rank*k_1/R); // j_0=rank
      for(size_t k_0=0; k_0<S; ++k_0){
        P[k_0]+=A_vect[k_1*S+k_0]*w; // P contains A1(rank,i) for i=0,...,S-1
      }
    }

    // each processor compute a S-points transform of vector P (using a FFT algorithm)
    vector<complex<double>> X_loc(S);
    
    Recursive_FFT(P,S);
  }
    // send results to rank 0, gathering R vectors of S elements X(rank,i) i=0,...,S
    vector<complex<double>> X_vect(N);

    MPI_Gather(P.data(), S, MPI_C_DOUBLE_COMPLEX, X_vect.data(), S, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

  // convert X(j_0,j_1) in X(j)
  vector<complex<double>> X(N);
  if(rank==0){
    size_t j_0, j_1; //redefinition! fuori da ogni if
    for(size_t j=0; j<N; ++j){
      j_0=mod(j,R);
      j_1=mod(j,S);
      X[j]=X_vect[j_1*R+j_0]; 
    }
  }

  return X;
} 


int main(){
  size_t n = 15;
  vector<complex<double>> input(n);
  for(size_t i = 0; i < n; i++){
    input[i] = (complex<double>(i, 0));
    cout << input[i] << " ";
  }
  cout << endl;

  vector<complex<double>> y = Parallel_FFT(input,n);

  for(size_t i = 0; i < n; i++){
    cout << y[i] << " ";
  }
  cout << endl;
  return 0;
}



