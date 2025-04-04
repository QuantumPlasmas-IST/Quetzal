#ifndef __PUSHER__
#define __PUSHER__

#include "stdlib.h"
#include "aux.h"
#include "initial_conditions.h"
#include "boundary_conditions.h"
#include "output.h"
#include "algos.h"
#include "math.h"

REAL** solve(input_t* input);

void integrate_source(input_t* input, REAL** species, COMPLEX** sources, int* aux_is, REAL* aux_momentum);
void convolve_source(input_t* input, COMPLEX** sources, COMPLEX** fields_fft, COMPLEX** fields_deriv, int* aux_is, REAL* aux_wavevector, int fld);
void invert_field(input_t* input, COMPLEX** fields, COMPLEX** fields_fft, int* aux_is, int fld);
void finite_volume1_fields(input_t* input, COMPLEX** fields, COMPLEX** fields_deriv, COMPLEX** fields_fft, int fld1);
void finite_volume2_fields(input_t* input, COMPLEX** fields, COMPLEX** fields_deriv, COMPLEX** fields_fft, int fld1);

// Spatial Operators
REAL central_diff(input_t* input, REAL** species, COMPLEX** fields, int i_sp, int index, int* aux_is, REAL* momentum);
void finite_volume2(input_t* input, REAL** species, COMPLEX** fields, int* aux_is, REAL* aux_momentum, REAL*** flows);

// Time Pushers
void leapfrog(input_t* input, REAL** species, REAL** aux1, REAL** aux2, COMPLEX**fields, int* aux_is);
void rungeKutta2(input_t* input, REAL** species, REAL** species_aux, COMPLEX** sources, COMPLEX** fields, COMPLEX** fields_aux, COMPLEX** fields_fft, COMPLEX** fields_deriv, REAL*** flows, int* aux_is, REAL* aux_momentum, REAL*** left_buffer, REAL*** right_buffer, COMPLEX*** left_field_buffer, COMPLEX*** right_field_buffer);
#endif