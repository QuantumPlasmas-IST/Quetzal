#ifndef __ALGOS__
#define __ALGOS__

#include "stdlib.h"
#include "math.h"
#include <complex.h>
#include "compile.h"

void LUdecomp(REAL** A, REAL** L, REAL** U, int N);
void LUinverse(REAL** A, REAL** L, REAL** U, int N);
void LUsolve(REAL** L, REAL** U, REAL* b, REAL* x, int N);

void fastFourierAuxiliary(int N, double complex* ft, double complex* gw, int s);
void fastInverseFourierAuxiliary(int N, double complex* ft, double complex* gw, int s);
void fastFourier(int N, double complex* ft, double complex* gw, int s);
void fastInverseFourier(int N, double complex* ft, double complex* gw, int s);
void NdFourier(int n_dims, int* N, double complex* ft, double complex* gw);
void NdInverseFourier(int n_dims, int* N, double complex* ft, double complex* gw);
void fourierTransform(int n_dims, int* N, double complex* ft, double complex* gw, double* dt, double* t0);
void inverseFourierTransform(int n_dims, int* N, double complex* ft, double complex* gw, double* dt, double* t0);

void FFTConvolution(int n_dims, int* N, double complex* f1, double complex* f2, double complex* g0);

__host__ __device__ void hadamardProduct(int n_dims, int* N, double complex* f1, double complex* f2, double complex* g0);
__host__ __device__ void kroneckerProduct(int n_dims, int* N1, int* N2, double complex* f1, double complex* f2, double complex* g0);
__host__ __device__ void tensorProduct(int n_dims, int* N, double complex* f1, double complex* f2, double complex* g0);

#endif