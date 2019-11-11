#include <iostream>
#include <vector>
#include <complex>
#include "omp.h"
#include "Util.hpp"
#include <set>

#define INTYPE double
#define RTYPE std::complex<INTYPE>

using namespace std;

void check(RTYPE refpsi, RTYPE psi) {
  if (abs(psi - refpsi) > 0.01f)
  {
    cout << "fail" << endl;
    cout << "Ref: " << refpsi << endl;
    cout << "Got: " << psi << endl;
  }
  else
    cout << "pass" << endl;
};
int pick_vec(std::vector<int> & v) {
  const int n = v.size();
  int r = rand() % n; // random index in vector
  int ret = v[r]; // get the value
  v.erase(v.begin() + r); // remove returned element
  return ret;
}
void fill_index(const int n, int* c, const int cache_line, const float fill) {
  int i, j;
  int num_cache_lines = n/cache_line; // full cache lines
  int remainder = n % cache_line;  // remaining elements
  int r = remainder > 0 ? 1 : 0; // add a cache line if there's a remainder
  int num_stride = cache_line * fill; // filled with stride 1
  int num_rand = cache_line - num_stride; // fillled randomly

  // create a vector of leftover indices
  // leftovers from each cache line + randoms left in last cache line
  int k = 0;
  int num_leftovers = num_cache_lines * num_rand + (remainder % num_stride);
  std::vector<int> leftover(num_leftovers);
  for (i = 0; i < num_cache_lines + r; i++) {
    for (j = 0; j < num_rand; j++) {
      int fill = i * cache_line + num_rand + j;
      if (fill < n) {
        leftover[j + num_rand * i] = fill;
        cout << "rem[" << j + num_rand * i << "]=" << fill << endl;
      }
    }
  }
  
  

  for (i = 0; i < num_cache_lines + r; i++) {
    for (j = 0; j < num_stride; j++) {
      if (i * cache_line + j < n)
        c[i * cache_line + j] = i * cache_line + j;
    }
    for (j = 0; j < num_rand; j++) {
      if (i * cache_line + num_stride + j < n)
        c[i * cache_line + num_stride + j] = pick_vec(leftover);
    }
  }

  for (i = 0; i < n; i++)
    std::cout << "c[" << i << "] = " << c[i] << std::endl;
}
int main(int argc, char** argv) {
  srand(1234);

  const int N = atoi(argv[1]);
  cout << "Using N = " << N << endl;
  std::vector<int>det0(N, 0);
  std::vector<int>det1(N, 0);
  std::vector<RTYPE>detValues0(N);
  std::vector<RTYPE>detValues1(N);
  ComplexSoA mydetValues0(N);
  ComplexSoA mydetValues1(N);
  INTYPE realdetValues0[N];
  INTYPE imagdetValues0[N];
  INTYPE realdetValues1[N];
  INTYPE imagdetValues1[N];
  RTYPE psi;

  INTYPE t;

  cout << "Initializing... det0\n";
  t = omp_get_wtime();
  fill_index(N, det0.data(), 6 , 0.5);
  t = omp_get_wtime() - t;
  cout << t << " sec" << endl;

//  cout << "Initializing...";
//  t = omp_get_wtime();
//  for (int i = 0; i < N; i++)
//  {
//    int ii=0, kk=0;
//    detValues0[i] = complex<INTYPE>(ii, kk);
//    detValues1[i] = complex<INTYPE>(kk, ii);
//  }
//  t = omp_get_wtime() - t;
//  cout << t << " sec" << endl;
//  cout << "\n\n\n\n";
//
//  psi = 0;
//  cout << "Compute...";
//  t = omp_get_wtime();
//  for (int i = 0; i < N; i++)
//    psi += detValues0[det0[i]] * detValues1[det1[i]];
//  t = omp_get_wtime() - t;
//  cout << t << " sec" << endl;
//  RTYPE psiref = psi;
//  
//
//
//  INTYPE psi_r = 0, psi_i = 0;
//  psi = 0;
//  cout << "Compute...";
//  t = omp_get_wtime();
//  for (int i = 0; i < N; i++) 
//  {
//    psi_r += detValues0[det0[i]].real() * detValues1[det1[i]].real();
//    psi_i += detValues0[det0[i]].imag() * detValues1[det1[i]].imag();
//  }
//  t = omp_get_wtime() - t;
//  psi = complex<INTYPE>(psi_r, psi_i);
//  cout << t << " sec" << endl;
//  check(psiref, psi);
//
//
//  psi = 0;
//  psi_r = 0;
//  psi_i = 0;
//  cout << "Compute...";
//  t = omp_get_wtime();
//#pragma omp parallel for simd reduction(+:psi_r, psi_i) schedule(static, 1)
//  for (int i = 0; i < N; i++) 
//  {
//    psi_r += mydetValues0.real(det0[i]) * mydetValues1.real(det1[i]);
//    psi_i += mydetValues0.imag(det0[i]) * mydetValues1.imag(det1[i]);
//  }
//  t = omp_get_wtime() - t;
//  psi = complex<INTYPE>(psi_r, psi_i);
//  cout << t << " sec" << endl;
//  check(psiref, psi);
//
//
//  psi = 0;
//  psi_r = 0;
//  psi_i = 0;
//  cout << "Compute...";
//  t = omp_get_wtime();
//#pragma omp parallel for simd reduction(+:psi_r, psi_i)
//  for (int i = 0; i < N; i++) 
//  {
//    psi_r += realdetValues0[det0[i]] * realdetValues1[det1[i]];
//    psi_i += imagdetValues0[det0[i]] * imagdetValues1[det1[i]];
//  }
//  t = omp_get_wtime() - t;
//  psi = complex<INTYPE>(psi_r, psi_i);
//  cout << t << " sec" << endl;
//  check(psiref, psi);
//
  return 0;
}
