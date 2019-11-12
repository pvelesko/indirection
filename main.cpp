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
#define CRINTPTR const int* __restrict__
#define CRINTYPEPTR const INTYPE* __restrict__
#ifndef DEBUG
#define DEBUG 0
#endif
using namespace std;

std::random_device rd;
std::mt19937 mt(rd());

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
int pick_vec(std::vector<int> & v) {
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
      // TODO bug here
      if (leftover.size() > 0 and i * cache_line + num_stride + j < n)
        c[i * cache_line + num_stride + j] = pick_vec(leftover);
    }
  }

//  for (i = 0; i < n; i++)
//    std::cout << "c[" << i << "] = " << c[i] << std::endl;
}
RTYPE calc0(const int N, vector<RTYPE> & detValues0, vector<RTYPE> & detValues1, CRINTPTR det0, CRINTPTR det1) {
  RTYPE psi = 0;
  for (int i = 0; i < N; i++)
    psi += detValues0[det0[i]] * detValues1[det1[i]];
  return psi;
}
RTYPE calc1(const int N, const vector<RTYPE> & detValues0, const vector<RTYPE> & detValues1, CRINTPTR det0, CRINTPTR det1) {
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
RTYPE calc2(const int N, const vector<RTYPE> & detValues0, const vector<RTYPE> & detValues1, CRINTPTR det0, CRINTPTR det1) {
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
RTYPE calc3(const int N, ComplexSoA & mydetValues0, ComplexSoA & mydetValues1, CRINTPTR det0, CRINTPTR det1) {
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
RTYPE calc4(const int N, CRINTYPEPTR realdetValues0, CRINTYPEPTR realdetValues1, CRINTYPEPTR imagdetValues0, CRINTYPEPTR imagdetValues1, CRINTPTR det0, CRINTPTR det1) {
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
  int num_t;

//  std::mt19937 mt(rd());
  #pragma omp parallel
  {
    #pragma omp master
    num_t = omp_get_num_threads();
  }
  const int N = atoi(argv[1]);
  const int M = atoi(argv[2]);
  const float R = atof(argv[3]);
  float size = ((2 * 2 * sizeof(INTYPE)) + (2 * sizeof(int))) * N / MEGA;
  cout << "Using N = " << N << endl;
  cout << "Using M(outer loop) = " << M << endl;
  cout << "Footprint = " << size << " MB" << endl;
  cout << "Using fill = " << R << endl;


  int* det0 = static_cast<int*>(malloc(sizeof(int) * N));
  int* det1 = static_cast<int*>(malloc(sizeof(int) * N));
  std::vector<RTYPE>detValues0(N, complex<INTYPE>(1, 1));
  std::vector<RTYPE>detValues1(N, complex<INTYPE>(1, 1));
  ComplexSoA mydetValues0(N);
  ComplexSoA mydetValues1(N);
  INTYPE* realdetValues0 = static_cast<INTYPE*>(malloc(sizeof(INTYPE) * N));
  INTYPE* imagdetValues0 = static_cast<INTYPE*>(malloc(sizeof(INTYPE) * N));
  INTYPE* realdetValues1 = static_cast<INTYPE*>(malloc(sizeof(INTYPE) * N));
  INTYPE* imagdetValues1 = static_cast<INTYPE*>(malloc(sizeof(INTYPE) * N));
  RTYPE psi, psiref;

  for (int i = 0; i < N; i++) {
    mydetValues0._real[i] = detValues0[i].real();
    mydetValues1._real[i] = detValues1[i].real();
    mydetValues0._imag[i] = detValues0[i].imag();
    mydetValues1._imag[i] = detValues1[i].imag();
    realdetValues0[i] = detValues0[i].real();
    realdetValues1[i] = detValues1[i].real();
    imagdetValues0[i] = detValues0[i].imag();
    imagdetValues1[i] = detValues1[i].imag();
  }

  cout << "Initialize... det0, det1\n";
  double t = omp_get_wtime();
  fill_index(N, det0, CACHELINE , R);
  fill_index(N, det1, CACHELINE , R);
  t = omp_get_wtime() - t;
  cout << t << " sec" << endl;


  cout << "calc0\n";
  double t0 = omp_get_wtime();
  for (int i = 0; i < M; i++)
#pragma noinline
    psiref = calc0(N, detValues0, detValues1, det0, det1);
  t0 = omp_get_wtime() - t0;

  cout << "calc1\n";
  double t1 = omp_get_wtime();
  for (int i = 0; i < M; i++)
#pragma noinline
    psi = calc1(N, detValues0, detValues1, det0, det1);
  t1 = omp_get_wtime() - t1;
  check(psiref, psi);

  cout << "calc2\n";
  double t2 = omp_get_wtime();
  for (int i = 0; i < M; i++)
#pragma noinline
    psi = calc2(N, detValues0, detValues1, det0, det1);
  t2 = omp_get_wtime() - t2;
  check(psiref, psi);

  cout << "calc3\n";
  double t3 = omp_get_wtime();
  for (int i = 0; i < M; i++)
#pragma noinline
    psi = calc3(N, mydetValues0, mydetValues1, det0, det1);
  t3 = omp_get_wtime() - t3;
  check(psiref, psi);

  cout << "calc4\n";
  double t4 = omp_get_wtime();
  for (int i = 0; i < M; i++)
#pragma noinline
    psi = calc4(N, realdetValues0, realdetValues1, imagdetValues0, imagdetValues1, det0, det1);
  t4 = omp_get_wtime() - t4;
  check(psiref, psi);

  cout << "-------------- RESULT -------------------" << endl;
  cout << "OpenMP Threads: " << num_t << endl;
  cout << std::left << std::setprecision(3) << std::setw(10) << t0    << " Runtime" <<  std::endl;
  cout << std::left << std::setprecision(3) << std::setw(10) << t0/t0 << " Test0 Complex" <<  std::endl;
  cout << std::left << std::setprecision(3) << std::setw(10) << t0/t1 << " Test1 Real/Imag" << std::endl;
  cout << std::left << std::setprecision(3) << std::setw(10) << t0/t2 << " Test2 Real/Imag SIMD HT" <<  std::endl;
  cout << std::left << std::setprecision(3) << std::setw(10) << t0/t3 << " Test3 ComplexSoA" <<  std::endl;
  cout << std::left << std::setprecision(3) << std::setw(10) << t0/t4 << " Test4 Complex Arrays" <<  std::endl;

  //free(det0);
  //free(det1);
  ////free(detValues0);
  ////free(detValues1);
  ////~mydetValues0();
  ////~mydetValues1();
  //free(realdetValues0);
  //free(imagdetValues0);
  //free(realdetValues1);
  //free(imagdetValues1);


  return 0;
}
