#ifndef __AUX__
#define __AUX__

#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include <time.h>
#include <ctype.h>
#include <string.h>
#include <complex.h>
#include "compile.h"
#include "mpi.h"
#include "fft3d_wrap.h"
#include "fft2d_wrap.h"
#include "remap3d_wrap.h"
#include "remap2d_wrap.h"

typedef void(*disp_t)(int, int, REAL*); // type: pointer to dispersion relation (band) function

typedef REAL(*init_t)(int, REAL*, REAL*); // type: pointer to initial condition function

typedef void(*bound_t)(int, int*, int, int, REAL*); //type: pointer to boundary condition function

typedef COMPLEX(*kernel_t)(int, COMPLEX, REAL*); // type: pointer to convolution kernel function

typedef COMPLEX(*dynamics_t)(int, int, COMPLEX**, REAL* ); // type: pointer to field dynamics function

typedef struct input{
   
    char* filename;
    int pos_diag_freq;
    int mom_diag_freq;

    int pos_dims;
    int mom_dims;
    
    int* pos_points;
    int* fft_points;
    int* mom_points;
    int pos_total;
    int mom_total;

    REAL* pos_delta;
    REAL* mom_delta;
    REAL* lambda;
    REAL* dk;
    
    REAL* pos_min;
    REAL* pos_max;
    REAL* mom_min;
    REAL* mom_max;
    
    int n_species;
    int n_fields;

    REAL deltaT;
    int n_timesteps;
    
    char** dispersion_names;
    disp_t* dispersions;

    char* pusher;
    char* operator;
    
    char** kernel_names;
    kernel_t* kernels;

    char** dynamics_names;
    dynamics_t* dynamics;
    
    char** force_names;
    REAL(**forces)(struct input*, COMPLEX**, int, int, int, int);

    char** field_type;

    char** pos_init_names;
    char** mom_init_names;
    char** field_init_names;
    init_t* pos_inits;
    init_t* mom_inits;
    init_t* field_inits;
    REAL** pos_init_params;
    REAL** mom_init_params;
    REAL** field_init_params;

    char** pos_bound_names;
    char** mom_bound_names;
    bound_t* pos_bounds;
    bound_t* mom_bounds;

    REAL** source_charges;
    REAL** force_charges;
    int** source_moments;

    int padding;

    int* grid_factor;
    int* space_factor;

    int rank;
    int size;
    int* procs;
    int* parallel_pos;
    int* parallel_factor;
    void* fft;

} input_t;

typedef REAL(*force_t)(input_t*, COMPLEX**, int, int, int, int); //type: pointer to force function

#endif