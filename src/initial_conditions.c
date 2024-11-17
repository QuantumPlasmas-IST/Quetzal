#include "initial_conditions.h"

void apply_init_cond(input_t* input, REAL** species){

    int* indices = malloc((input->pos_dims + input->mom_dims) * sizeof(int));
    REAL* positions = malloc(input->pos_dims * sizeof(REAL));
    REAL* momenta = malloc(input->mom_dims * sizeof(REAL));

    for (int i = 0; i < input->n_species; ++i){
        for (int j = 0; j < input->pos_total * input->mom_total; ++j){
            species[i][j] = 1;

            axis_index(input->pos_dims+input->mom_dims, input->pos_points, indices, j);

            for(int pdim = 0; pdim < input->pos_dims; ++pdim){
                positions[pdim] = (indices[pdim]-input->padding) * input->pos_delta[pdim] + input->pos_min[pdim];
            }
            for(int mdim = 0; mdim < input->mom_dims; ++mdim){
                momenta[mdim] = (indices[input->pos_dims+mdim]-input->padding) * input->mom_delta[mdim] + input->mom_min[mdim];
            }

            species[i][j] *= input->pos_inits[i](input->pos_dims, positions, input->pos_init_params[i]);

            species[i][j] *= input->mom_inits[i](input->mom_dims, momenta, input->mom_init_params[i]);
        }
    }

    free(indices);
    free(positions);
    free(momenta);
}

void field_init_cond(input_t* input, COMPLEX** fields){

    int* indices = malloc(input->pos_dims * sizeof(int));
    REAL* positions = malloc(input->pos_dims * sizeof(REAL));

    for (int i = 0; i < input->n_fields; ++i){
        for (int j = 0; j < input->pos_total; ++j){
            fields[i][j] = 1;

            axis_index(input->pos_dims, input->pos_points, indices, j);

            for(int pdim = 0; pdim < input->pos_dims; ++pdim){
                positions[pdim] = (indices[pdim]-input->padding) * input->pos_delta[pdim] + input->pos_min[pdim];
            }
            
            fields[i][j] *= input->field_inits[i](input->pos_dims, positions, input->field_init_params[i]);
        
        }
    }

    free(indices);
    free(positions);
}

REAL gaussian(int n, REAL* x, REAL* p){
    x[0] = (x[0] - p[2])*(x[0] - p[2]);
    for(int i = 1; i < n; ++i){
        x[0] += (x[i] - p[2+i])*(x[i] - p[2+i]);
    }
    //printf("%f\n", p[0] * exp( -x[0] / (p[1]*p[1])) / (sqrt(M_PI) * p[1]));
    return p[0] * exp( -x[0] / (p[1]*p[1]));
}

REAL uniform(int n, REAL* x, REAL* p){
    for (int i = 0; i < n; ++i){
        if (x[i] < p[1+2*i]) return 0;
        if (x[i] > p[2+2*i]) return 0;
    }
    return p[0];
}

REAL double_quad_maxwell(int n, REAL* x, REAL* p){
    REAL aux = x[0];
    REAL m1 = quad_maxwell(n, x, p);
    x[0] = aux;
    return quad_maxwell(n, x, p+2+n) + m1;
}

REAL double_gaussian(int n, REAL* x, REAL* p){
    REAL aux = x[0];
    REAL m1 = gaussian(n, x, p);
    x[0] = aux;
    return gaussian(n, x, p+2+n) + m1;
}

REAL noise(int n, REAL* x, REAL* p){

    return p[0] + rand() * p[1]/RAND_MAX;
}

REAL linear_fermi(int n, REAL* x, REAL* p){
    x[0] = (x[0] - p[2])*(x[0] - p[2]);
    for(int i = 1; i < n; ++i){
        x[0] += (x[i] - p[2+i])*(x[i] - p[2+i]);
    }
    x[0] = sqrt(x[0]);
    
    return 1.0 / (1 + exp((x[0] - p[0]) / p[1]));
}

REAL quad_fermi(int n, REAL* x, REAL* p){
    x[0] = (x[0] - p[2])*(x[0] - p[2]);
    for(int i = 1; i < n; ++i){
        x[0] += (x[i] - p[2+i])*(x[i] - p[2+i]);
    }

    return 1 / (1 + exp((x[0]/2 - p[0]) / p[1]));
}

REAL linear_bose(int n, REAL* x, REAL* p){
    x[0] = (x[0] - p[2])*(x[0] - p[2]);
    for(int i = 1; i < n; ++i){
        x[0] += (x[i] - p[2+i])*(x[i] - p[2+i]);
    }
    x[0] = sqrt(x[0]);
    
    return 1.0 / (-1 + exp((x[0] - p[0]) / p[1]));
}

REAL quad_bose(int n, REAL* x, REAL* p){
    x[0] = (x[0] - p[2])*(x[0] - p[2]);
    for(int i = 1; i < n; ++i){
        x[0] += (x[i] - p[2+i])*(x[i] - p[2+i]);
    }

    return 1 / (-1 + exp((x[0]/2 - p[0]) / p[1]));
}

REAL linear_maxwell(int n, REAL* x, REAL* p){
    x[0] = (x[0] - p[2])*(x[0] - p[2]);
    for(int i = 1; i < n; ++i){
        x[0] += (x[i] - p[2+i])*(x[i] - p[2+i]);
    }
    x[0] = sqrt(x[0]);

    return p[0] * exp(- x[0] / p[1]) / (2 * pow(sqrt(M_PI) * p[1],n)*tgamma(n)/tgamma(((double)(n))/2.));
}

REAL quad_maxwell(int n, REAL* x, REAL* p){
    x[0] = (x[0] - p[2])*(x[0] - p[2]);
    for(int i = 1; i < n; ++i){
        x[0] += (x[i] - p[2+i])*(x[i] - p[2+i]);
    }
    
    return p[0] * exp(- x[0]/2 / p[1]) / pow(2 * M_PI * p[1], ((double)(n))/2.);
}

REAL flattened_quad_fermi(int n, REAL* x, REAL* p){
    x[0] = (x[0] - p[2])*(x[0] - p[2]);
    for(int i = 1; i < n; ++i){
        x[0] += (x[i] - p[2+i])*(x[i] - p[2+i]);
    }

    return sqrt(2*M_PI*p[1]) * gsl_sf_fermi_dirac_mhalf(-(x[0]/2 - p[0]) / p[1]);
} 

REAL dirac_delta(int n, REAL* x, REAL* p){
    int check = 1;
    for (int i = 0; i < n; ++i){
        if (fabs(x[i]-p[1+i])>1e-5) check = 0;
    }
    if (check) return p[0];
    return 0;
}

REAL sine(int n, REAL* x, REAL* p){
    x[0] = sin(p[2] * x[0]);
    for(int i = 1; i < n; ++i){
        x[0] *= sin(p[2+i] * x[i]);
    }
    return p[0] + p[1] * x[0];
}

REAL anisotropic_quad_maxwell(int n, REAL* x, REAL* p){
    x[0] = (x[0] - p[n+1])*(x[0] - p[n+1])/(2 * p[1]);
    REAL norm = sqrt(2 * M_PI * p[1]);
    for(int i = 1; i < n; ++i){
        x[0] += (x[i] - p[n+1+i])*(x[i] - p[n+1+i])/(2 * p[1+i]);
        norm *= sqrt(2 * M_PI * p[1+i]);
    }
    
    return p[0] * exp(- x[0] ) / norm;
}