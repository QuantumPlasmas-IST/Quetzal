#ifndef __BOUNDCOND__
#define __BOUNDCOND__

#include "compile.h"
#include "aux.h"
#include "utilities.h"
#include "mpi.h"

void apply_bound_cond(input_t* input, REAL** species, REAL*** left_buffer, REAL*** right_buffer);
void field_bound_cond(input_t* input, COMPLEX** fields, COMPLEX*** left_buffer, COMPLEX*** right_buffer, int fld);

void dirichelet_bound(int total_dims, int* Ns, int dim, int padding, REAL* species);
void periodic_bound(int total_dims, int* Ns, int dim, int padding, REAL* species);
void forced_periodic_bound(int total_dims, int* Ns, int dim, int padding, REAL* species);
void forced_field_periodic_bound(int total_dims, int* Ns, int dim, int padding, COMPLEX* species);

void boundary_shift(input_t* input, COMPLEX* start, COMPLEX* final, int* aux);
void inverse_boundary_shift(input_t* input, COMPLEX* start, COMPLEX* final, int* aux);
void species_boundary_shift(input_t* input, REAL* start, REAL* final, int* aux);
void species_inverse_boundary_shift(input_t* input, REAL* start, REAL* final, int* aux);

void ghost_cell_transfer(input_t* input, REAL* species, REAL*** left_buffer, REAL*** right_buffer);
void field_cell_transfer(input_t* input, COMPLEX* field, COMPLEX*** left_buffer, COMPLEX*** right_buffer);

#endif