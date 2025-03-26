#ifndef __BOUNDCOND__
#define __BOUNDCOND__

#include "compile.h"
#include "aux.h"
#include "utilities.h"
#include "mpi.h"

void apply_bound_cond(input_t* input, REAL** species, REAL*** left_buffer, REAL*** right_buffer, int* aux_is);
void field_bound_cond(input_t* input, COMPLEX** fields, COMPLEX*** left_buffer, COMPLEX*** right_buffer, int fld);

void dirichelet_bound(input_t* input, int dim, int* aux_is, REAL* species, int side);
void periodic_bound(input_t* input, int dim, int* aux_is, REAL* species, int side);
void wall_bound(input_t* input, int dim, int* aux_is, REAL* species, int side);

void boundary_shift(input_t* input, COMPLEX* start, COMPLEX* final, int* aux);
void inverse_boundary_shift(input_t* input, COMPLEX* start, COMPLEX* final, int* aux);
void species_boundary_shift(input_t* input, REAL* start, REAL* final, int* aux);
void species_inverse_boundary_shift(input_t* input, REAL* start, REAL* final, int* aux);

void ghost_cell_transfer(input_t* input, REAL* species, REAL*** left_buffer, REAL*** right_buffer);
void field_cell_transfer(input_t* input, COMPLEX* field, COMPLEX*** left_buffer, COMPLEX*** right_buffer);

#endif