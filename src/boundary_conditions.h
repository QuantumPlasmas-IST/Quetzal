#ifndef __BOUNDCOND__
#define __BOUNDCOND__

#include "compile.h"
#include "aux.h"
#include "utilities.h"

void apply_bound_cond(input_t* input, REAL** species);

void dirichelet_bound(int total_dims, int* Ns, int dim, int padding, REAL* species);
void periodic_bound(int total_dims, int* Ns, int dim, int padding, REAL* species);

void boundary_shift(input_t* input, COMPLEX* start, COMPLEX* final, int* aux);
void inverse_boundary_shift(input_t* input, COMPLEX* start, COMPLEX* final, int* aux);

#endif