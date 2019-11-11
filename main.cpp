#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include "omp.h"
#include "Util.hpp"
#include <set>
#include <random>
#include <cmath>
#define INTYPE double
#define RTYPE std::complex<INTYPE>
#define CACHELINE 64
#define MEGA 1000000.f
#define GIGA 1000000000.f

using namespace std;

void check(RTYPE refpsi, RTYPE psi) {
  if (abs(psi - refpsi) > 0.01f)
  {
    cout << "fail" << endl;
    cout << "Ref: " << refpsi << endl;
    cout << "Got: " << psi << endl;
  }
  else
    return;
    //cout << "pass" << endl;
};
int pick_vec(std::mt19937 & mt, std::vector<int> & v) {
  const int n = v.size();
  if (n == 0) {
    cout << "Crap" << endl;
    exit(1);
  }
  std::uniform_real_distribution<double> dist(0, n-1);
  int r = dist(mt); // random index in vector
  int ret = v[r]; // get the value
  v.erase(v.begin() + r); // remove returned element
  return ret;
}
void fill_index(const int n, int* c, const int cache_line_bytes, const float ratio) {
  std::random_device rd;
  std::mt19937 mt(rd());

  int cache_line = cache_line_bytes / sizeof(int);
  int i, j;
  int num_cache_lines = n/cache_line; // full cache lines
  int remainder = n % cache_line;  // remaining elements
  int r = remainder > 0 ? 1 : 0; // add a cache line if there's a remainder
  int num_stride = cache_line * ratio; // filled with stride 1
  if (num_stride < 1) num_stride++; // why cna't I use min?
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
//        cout << "rem[" << j + num_rand * i << "]=" << fill << endl;
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
        c[i * cache_line + num_stride + j] = pick_vec(mt, leftover);
    }
  }

//  for (i = 0; i < n; i++)
//    std::cout << "c[" << i << "] = " << c[i] << std::endl;
}
RTYPE calc0(const int N, vector<RTYPE> & detValues0, vector<RTYPE> & detValues1, vector<int> & det0, vector<int> & det1) {
  RTYPE psi = 0;
  for (int i = 0; i < N; i++)
    psi += detValues0[det0[i]] * detValues1[det1[i]];
  return psi;
}
RTYPE calc1(const int N, const vector<RTYPE> & detValues0, const vector<RTYPE> & detValues1, const vector<int> & det0, const vector<int> & det1) {
  INTYPE psi_r = 0, psi_i = 0;
  RTYPE psi = 0;
  for (int i = 0; i < N; i++) 
  {
    psi_r += detValues0[det0[i]].real() * detValues1[det1[i]].real() - detValues0[det0[i]].imag() * detValues1[det1[i]].imag();
    psi_i += detValues0[det0[i]].real() * detValues1[det1[i]].real() + detValues0[det0[i]].imag() * detValues1[det1[i]].imag();
  }
  psi = complex<INTYPE>(psi_r, psi_i);
  return psi;
}
RTYPE calc1_simdht(const int N, const vector<RTYPE> & detValues0, const vector<RTYPE> & detValues1, const vector<int> & det0, const vector<int> & det1) {
  INTYPE psi_r = 0, psi_i = 0;
  RTYPE psi = 0;
#pragma omp parallel for simd reduction(+:psi_r, psi_i) schedule(static, 1)
  for (int i = 0; i < N; i++) 
  {
    psi_r += detValues0[det0[i]].real() * detValues1[det1[i]].real() - detValues0[det0[i]].imag() * detValues1[det1[i]].imag();
    psi_i += detValues0[det0[i]].real() * detValues1[det1[i]].real() + detValues0[det0[i]].imag() * detValues1[det1[i]].imag();
  }
  psi = complex<INTYPE>(psi_r, psi_i);
  return psi;
}
RTYPE calc2(const int N, ComplexSoA & mydetValues0, ComplexSoA & mydetValues1, vector<int> & det0, vector<int> & det1) {
  RTYPE psi = 0;
  INTYPE psi_r = 0;
  INTYPE psi_i = 0;
  auto t = omp_get_wtime();
#pragma omp parallel for simd reduction(+:psi_r, psi_i) schedule(static, 1)
  for (int i = 0; i < N; i++) 
  {
    psi_r += mydetValues0.real(det0[i]) * mydetValues1.real(det1[i]) - mydetValues0.imag(det0[i]) * mydetValues1.imag(det1[i]);
    psi_i += mydetValues0.real(det0[i]) * mydetValues1.real(det1[i]) + mydetValues0.imag(det0[i]) * mydetValues1.imag(det1[i]);
  }
  psi = complex<INTYPE>(psi_r, psi_i);
  return psi;
}
RTYPE calc3(const int N, INTYPE* realdetValues0, INTYPE* realdetValues1, INTYPE* imagdetValues0, INTYPE* imagdetValues1, vector<int> & det0, vector<int> & det1) {
  RTYPE psi = 0;
  INTYPE psi_r = 0;
  INTYPE psi_i = 0;
#pragma omp parallel for simd reduction(+:psi_r, psi_i) schedule(static, 1)
  for (int i = 0; i < N; i++) 
  {
    psi_r += realdetValues0[det0[i]] * realdetValues1[det1[i]] - imagdetValues0[det0[i]] * imagdetValues1[det1[i]];
    psi_i += realdetValues0[det0[i]] * realdetValues1[det1[i]] + imagdetValues0[det0[i]] * imagdetValues1[det1[i]];
  }
  psi = complex<INTYPE>(psi_r, psi_i);
  return psi;
}

int main(int argc, char** argv) {

  const int N = atoi(argv[1]);
  const int M = atoi(argv[2]);
  const float R = atof(argv[3]);
  float size = (2 * 2 * sizeof(INTYPE)) + (2 * sizeof(int)) * N / MEGA;
  cout << "Using N = " << N << endl;
  cout << "Using M(outer loop) = " << M << endl;
  cout << "Footprint = " << size << " MB" << endl;
  cout << "Using fill = " << R << endl;


  std::vector<int>det0(N);
  std::vector<int>det1(N);
  std::vector<RTYPE>detValues0(N, complex<INTYPE>(1, 1));
  std::vector<RTYPE>detValues1(N, complex<INTYPE>(1, 1));
  ComplexSoA mydetValues0(N);
  ComplexSoA mydetValues1(N);
  INTYPE realdetValues0[N];
  INTYPE imagdetValues0[N];
  INTYPE realdetValues1[N];
  INTYPE imagdetValues1[N];
  RTYPE psi, psiref;

  for (int i = 0; i < N; i++) {
    mydetValues0._real[i] = detValues0[i].real();
    mydetValues1._real[i] = detValues1[i].real();
    mydetValues0._imag[i] = detValues0[i].imag();
    mydetValues1._imag[i] = detValues1[i].imag();
  }

  cout << "Initialize... det0, det1\n";
  double t = omp_get_wtime();
  fill_index(N, det0.data(), CACHELINE , R);
  fill_index(N, det1.data(), CACHELINE , R);
  t = omp_get_wtime() - t;
  cout << t << " sec" << endl;


  double t0 = omp_get_wtime();
  for (int i = 0; i < M; i++)
#pragma noinline
    psiref = calc0(N, detValues0, detValues1, det0, det1);
  t0 = omp_get_wtime() - t0;

  double t1ht = omp_get_wtime();
  for (int i = 0; i < M; i++)
#pragma noinline
    psi = calc1_simdht(N, detValues0, detValues1, det0, det1);
  t1ht = omp_get_wtime() - t1ht;
  check(psiref, psi);


  double t1 = omp_get_wtime();
  for (int i = 0; i < M; i++)
#pragma noinline
    psi = calc1(N, detValues0, detValues1, det0, det1);
  t1 = omp_get_wtime() - t1;
  check(psiref, psi);

  double t2 = omp_get_wtime();
  for (int i = 0; i < M; i++)
#pragma noinline
    psi = calc2(N, mydetValues0, mydetValues1, det0, det1);
  t2 = omp_get_wtime() - t2;
  check(psiref, psi);

  double t3 = omp_get_wtime();
  for (int i = 0; i < M; i++)
#pragma noinline
    psi = calc3(N, realdetValues0, realdetValues1, imagdetValues0, imagdetValues1, det0, det1);
  t3 = omp_get_wtime() - t3;
  check(psiref, psi);

  cout << "-------------- RESULT -------------------" << endl;
  cout << std::left << std::setw(8) << t0/t0 << std::endl;
  cout << std::left << std::setw(8) << t0/t1 << std::endl;
  cout << std::left << std::setw(8) << t0/t1ht << std::endl;
  cout << std::left << std::setw(8) << t0/t2 << std::endl;
  cout << std::left << std::setw(8) << t0/t3 << std::endl;

  return 0;
}
