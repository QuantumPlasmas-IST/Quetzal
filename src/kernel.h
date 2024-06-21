#ifndef __KERNEL__
#define __KERNEL__

#include "aux.h"
#include "compile.h"

COMPLEX poisson_2D(int n_dims, COMPLEX density, REAL* wavevector);

COMPLEX poisson_3D(int n_dims, COMPLEX density, REAL* wavevector);

#endif