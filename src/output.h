#ifndef __OUTPUT__
#define __OUTPUT__

#include "aux.h"
#include "hdf5.h"
#include "mpi.h"
#include "H5FDmpio.h"
#include "H5FDmpi.h"

typedef struct hdf5_complex{
    REAL re;
    REAL im;
} hdf5_complex_t;

void open_file(input_t* input);

void write_solution(input_t* input, REAL** results, int timestep);
void write_fields(input_t* input, COMPLEX** fields, int timestep);
void write_sources(input_t* input, COMPLEX** sources, int timestep);

#endif