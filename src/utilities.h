#ifndef __UTILS__
#define __UTILS__

#include "aux.h"
#include "compile.h"

int binomial_coeff(int n, int m);

int factorial(int n);

int ipow(int base, int exp);

int* list_combinations(int n, int m, int* list, int k);
int* list_combinations_all(int n, int m, int* list, int k);

int list_index(int n_dims, int* Ns, int* is);
void axis_index(int n_dims, int* Ns, int* is, int index);

REAL minmod3(REAL x1, REAL x2, REAL x3);
REAL maxmod3(REAL x1, REAL x2, REAL x3);
REAL minmod2(REAL x1, REAL x2);
REAL maxmod2(REAL x1, REAL x2);
REAL min2(REAL x1, REAL x2);
REAL max2(REAL x1, REAL x2);
REAL min3(REAL x1, REAL x2, REAL x3);
REAL max3(REAL x1, REAL x2, REAL x3);
REAL minmod_slopeLimiter(REAL u0, REAL u1, REAL u2, REAL Dx);
REAL superbee_slopeLimiter(REAL u0, REAL u1, REAL u2, REAL Dx);
REAL MC_slopeLimiter(REAL u0, REAL u1, REAL u2, REAL Dx);
REAL vanLeer_slopeLimiter(REAL u0, REAL u1, REAL u2, REAL Dx);
REAL minmod_fluxLimiter(REAL u0, REAL u1, REAL u2);
REAL superbee_fluxLimiter(REAL u0, REAL u1, REAL u2);
REAL MC_fluxLimiter(REAL u0, REAL u1, REAL u2);
REAL vanLeer_fluxLimiter(REAL u0, REAL u1, REAL u2);

#endif