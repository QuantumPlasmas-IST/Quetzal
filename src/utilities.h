#ifndef __UTILS__
#define __UTILS__

#include "aux.h"
#include "compile.h"

__host__ __device__ int binomial_coeff(int n, int m);

__host__ __device__ int factorial(int n);

__host__ __device__ int ipow(int base, int exp);

__host__ __device__ int* list_combinations(int n, int m, int* list, int k);
__host__ __device__ int* list_combinations_all(int n, int m, int* list, int k);

__host__ __device__ int list_index(int n_dims, int* Ns, int* is);
__host__ __device__ void axis_index(int n_dims, int* Ns, int* is, int index);

__host__ __device__ REAL minmod3(REAL x1, REAL x2, REAL x3);
__host__ __device__ REAL maxmod3(REAL x1, REAL x2, REAL x3);
__host__ __device__ REAL minmod2(REAL x1, REAL x2);
__host__ __device__ REAL maxmod2(REAL x1, REAL x2);
__host__ __device__ REAL min2(REAL x1, REAL x2);
__host__ __device__ REAL max2(REAL x1, REAL x2);
__host__ __device__ REAL min3(REAL x1, REAL x2, REAL x3);
__host__ __device__ REAL max3(REAL x1, REAL x2, REAL x3);
__host__ __device__ REAL minmod_slopeLimiter(REAL u0, REAL u1, REAL u2, REAL Dx);
__host__ __device__ REAL superbee_slopeLimiter(REAL u0, REAL u1, REAL u2, REAL Dx);
__host__ __device__ REAL MC_slopeLimiter(REAL u0, REAL u1, REAL u2, REAL Dx);
__host__ __device__ REAL vanLeer_slopeLimiter(REAL u0, REAL u1, REAL u2, REAL Dx);
__host__ __device__ REAL minmod_fluxLimiter(REAL u0, REAL u1, REAL u2);
__host__ __device__ REAL superbee_fluxLimiter(REAL u0, REAL u1, REAL u2);
__host__ __device__ REAL MC_fluxLimiter(REAL u0, REAL u1, REAL u2);
__host__ __device__ REAL vanLeer_fluxLimiter(REAL u0, REAL u1, REAL u2);

#endif