#ifndef __DISP__
#define __DISP__

#include "aux.h"
#include "compile.h"

void linear_dispersion(int pos_dims, int mom_dims, REAL* p);

void quadratic_dispersion(int pos_dims, int mom_dims, REAL* p);

void anisotropic_quadratic_dispersion(int pos_dims, int mom_dims, REAL* p);

#endif