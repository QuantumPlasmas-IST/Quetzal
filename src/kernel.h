#ifndef __KERNEL__
#define __KERNEL__

#include "aux.h"
#include "compile.h"

COMPLEX density(int n_dims, COMPLEX density, REAL* wavevector);

COMPLEX poisson_2D(int n_dims, COMPLEX density, REAL* wavevector);

COMPLEX poisson_3D(int n_dims, COMPLEX density, REAL* wavevector);

COMPLEX softcore_2D(int n_dims, COMPLEX density, REAL* wavevector);

COMPLEX softcore_3D(int n_dims, COMPLEX density, REAL* wavevector);

COMPLEX maxwell_2D(int n_dims, COMPLEX density, REAL* wavevector);

COMPLEX maxwell_3D(int n_dims, COMPLEX density, REAL* wavevector);

COMPLEX null_kernel(int n_dims, COMPLEX density, REAL* wavevector);

#endif