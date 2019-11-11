#include <cstdlib>
class ComplexSoA
{
  public:
  double* __restrict__  _real;
  double* __restrict__  _imag;

  public:
  ComplexSoA(int n);

  inline double real(int i) {return _real[i];};
  inline double imag(int i) {return _imag[i];};
  //operator[](const int i);
  //ComplexSoA operator*(const ComplexSoA);
  //ComplexSoA operator=(const ComplexSoA);
};

ComplexSoA::ComplexSoA(int n)
{
  _real = (double*)malloc(sizeof(double) * n);
  _imag = (double*)malloc(sizeof(double) * n);
};


