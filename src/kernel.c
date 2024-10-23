#include "kernel.h"

COMPLEX density(int n_dims, COMPLEX density, REAL* wavevector){
    return density;
}

COMPLEX poisson_2D(int n_dims, COMPLEX density, REAL* wavevector){
    REAL q = 0;
    for(int i = 0; i < n_dims; ++i){
        q += wavevector[i]*wavevector[i];
    }
    q = sqrt(q);

    if (q < 1e-8) return 0;
    return density/(2*q);
}

COMPLEX poisson_3D(int n_dims, COMPLEX density, REAL* wavevector){
    REAL q = 0;
    for(int i = 0; i < n_dims; ++i){
        q += wavevector[i]*wavevector[i];
    }

    if (q < 1e-8) return 0;
    return density/q;
}

COMPLEX softcore_2D(int n_dims, COMPLEX density, REAL* wavevector){
    REAL q = 0;
    for(int i = 0; i < n_dims; ++i){
        q += wavevector[i]*wavevector[i];
    }
    q = sqrt(q);

    if (q < 1e-8) return 1/2;

    return density * j1(q) / q;
}

COMPLEX softcore_3D(int n_dims, COMPLEX density, REAL* wavevector){
    REAL q = 0;
    for(int i = 0; i < n_dims; ++i){
        q += wavevector[i]*wavevector[i];
    }
    q = sqrt(q);

    if (q < 1e-8) return 1/3;
    return density * (sin(q) - q * cos(q))/(q*q*q);
}

COMPLEX maxwell_2D(int n_dims, COMPLEX density, REAL* wavevector){

    REAL q = 0;
    for(int i = 0; i < n_dims; ++i){
        q += wavevector[i]*wavevector[i];
    }
    q = sqrt(q);

    return - q * density / 2;
}

COMPLEX maxwell_3D(int n_dims, COMPLEX density, REAL* wavevector){
    
    return -density;
}