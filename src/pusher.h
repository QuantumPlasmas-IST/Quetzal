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

void integrate_source(input_t* input, REAL** species, COMPLEX** sources, COMPLEX** sources_aux);
void convolute_field(input_t* input, COMPLEX** sources, COMPLEX** fields, COMPLEX*** left_buffer, COMPLEX*** right_buffer);

// Spatial Operators
REAL central_diff(input_t* input, REAL** species, COMPLEX** fields, int i_sp, int index, int* aux_is, REAL* momentum);
void finite_volume2(input_t* input, REAL** species, COMPLEX** fields, int* aux_is, REAL* aux_momentum, REAL*** flows);
void finite_volumeNL2(input_t* input, REAL** species, COMPLEX** fields, int* aux_is, REAL* aux_momentum, REAL*** flows);

// Time Pushers
void leapfrog(input_t* input, REAL** species, REAL** aux1, REAL** aux2, COMPLEX**fields, int* aux_is);
void rungeKutta2(input_t* input, REAL** species, REAL** aux1, REAL** aux2, COMPLEX** sources, COMPLEX** sources_aux, COMPLEX** fields, REAL*** flows, int* aux_is, REAL* aux_momentum, REAL*** left_buffer, REAL*** right_buffer, COMPLEX*** left_field_buffer, COMPLEX*** right_field_buffer);
#endif