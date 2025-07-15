#ifndef __INITCOND__
#define __INITCOND__

#include "aux.h"
#include "compile.h"
#include "utilities.h"
#include <gsl/gsl_sf_fermi_dirac.h>

void apply_init_cond(input_t* input, REAL** species);
void field_init_cond(input_t* input, COMPLEX** fields);

REAL gaussian(int n, REAL* x, REAL* p);
REAL linear_gaussian(int n, REAL* x, REAL* p);
REAL anisotropic_gaussian(int n, REAL* x, REAL* p);
REAL uniform(int n, REAL* x, REAL* p);
REAL double_gaussian(int n, REAL* x, REAL* p);
REAL noise(int n, REAL* x, REAL* p);
REAL linear_fermi(int n, REAL* x, REAL* p);
REAL quad_fermi(int n, REAL* x, REAL* p);
REAL linear_bose(int n, REAL* x, REAL* p);
REAL quad_bose(int n, REAL* x, REAL* p);
REAL plasmon_bose(int n, REAL* x, REAL* p);
REAL linear_maxwell(int n, REAL* x, REAL* p);
REAL quad_maxwell(int n, REAL* x, REAL* p);
REAL double_quad_maxwell(int n, REAL* x, REAL* p);
REAL circle(int n, REAL* x, REAL* p);
REAL sphere(int n, REAL* x, REAL* p);
REAL flattened_quad_fermi(int n, REAL* x, REAL* p);
REAL dirac_delta(int n, REAL* x, REAL* p);
REAL sine(int n, REAL* x, REAL* p);
REAL anisotropic_quad_maxwell(int n, REAL* x, REAL* p);
REAL anisotropic_quad_fermi(int n, REAL* x, REAL* p);
REAL anisotropic_linear_fermi(int n, REAL* x, REAL* p);


#endif
