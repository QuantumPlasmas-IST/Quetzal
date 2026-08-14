#ifndef __FORCE__
#define __FORCE__

#include "aux.h"
#include "compile.h"

__host__ __device__ REAL gradient_force(input_t* input, COMPLEX** fields, int fld, int j, int k, int dim);

__host__ __device__ REAL electric0_force(input_t* input, COMPLEX** fields, int fld, int j, int k, int dim);

__host__ __device__ REAL electric1_force(input_t* input, COMPLEX** fields, int fld, int j, int k, int dim);

__host__ __device__ REAL magnetic_force(input_t* input, COMPLEX** fields, int fld, int j, int k, int dim);

#endif