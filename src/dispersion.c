#include "dispersion.h"

void linear_dispersion(int pos_dims, int mom_dims, REAL* p){

    REAL total_p = 0;

    for(int i = 0; i < mom_dims; ++i){
        total_p += p[i]*p[i];
    }
    total_p = sqrt(total_p);

    for(int i = 0; i < pos_dims; ++i){
        p[i] = p[i] / total_p;
    }
}

// Yes, it's supposed to do nothing
void quadratic_dispersion(int pos_dims, int mom_dims, REAL* p){}


// TODO add user-specified parameters to this
void anisotropic_quadratic_dispersion(int pos_dims, int mom_dims, REAL* p){
    
    p[1] /= 3.;

}