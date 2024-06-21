#ifndef __FORCE__
#define __FORCE__

#include "aux.h"
#include "compile.h"

REAL gradient_force(input_t* input, COMPLEX** fields, int fld, int j, int k, int dim);

REAL electric0_force(input_t* input, COMPLEX** fields, int fld, int j, int k, int dim);

REAL electric1_force(input_t* input, COMPLEX** fields, int fld, int j, int k, int dim);

REAL magnetic_force(input_t* input, COMPLEX** fields, int fld, int j, int k, int dim);

#endif