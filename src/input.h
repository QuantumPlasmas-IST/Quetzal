#ifndef __INPUT__
#define __INPUT__

#include "stdlib.h"
#include "aux.h"
#include <time.h>
#include <ctype.h>
#include <string.h>
#include "stdio.h"
#include "initial_conditions.h"
#include "kernel.h"
#include "boundary_conditions.h"
#include "dispersion.h"
#include "force.h"

input_t* read_input(const char* name);

void initialize_input(input_t* input, REAL* pos_parameters, REAL* mom_parameters, REAL* field_parameters);

void print_input(input_t* input);

void free_input(input_t* input);

#endif