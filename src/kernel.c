#include "kernel.h"

COMPLEX poisson_2D(int n_dims, COMPLEX density, REAL* wavevector){
    REAL q = 0;
    for(int i = 0; i < n_dims; ++i){
        q += wavevector[i]*wavevector[i];
    }
    q = sqrt(q);

    if (q < 1e-8) return 0;
    return density/q;
}

COMPLEX poisson_3D(int n_dims, COMPLEX density, REAL* wavevector){
    REAL q = 0;
    for(int i = 0; i < n_dims; ++i){
        q += wavevector[i]*wavevector[i];
    }

    if (q < 1e-8) return 0;
    return density/q;
}