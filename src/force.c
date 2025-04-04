#include "force.h"

REAL gradient_force(input_t* input, COMPLEX** fields, int fld, int j, int k, int dim){

    if(dim>=input->pos_dims) return 0;
    
    REAL forw = creal(fields[fld][j + input->space_factor[dim]]);
    REAL back = creal(fields[fld][j - input->space_factor[dim]]);
    return - (forw - back) / (2 * input->pos_delta[dim]);
}

REAL electric0_force(input_t* input, COMPLEX** fields, int fld, int j, int k, int dim){
    
    return (1-dim) * creal(fields[fld][j]);

}

REAL electric1_force(input_t* input, COMPLEX** fields, int fld, int j, int k, int dim){
    
    return dim * creal(fields[fld][j]);

}

REAL magnetic_force(input_t* input, COMPLEX** fields, int fld, int j, int k, int dim){

    //TODO Check this carefully
    int mom_index = dim * (k/input->mom_points[0]) + (1-dim) * (k%input->mom_points[0]);
    return creal((-2*dim+1) * fields[fld][j] * (input->mom_min[1-dim] + input->mom_delta[1-dim] *(mom_index-input->padding)));

}