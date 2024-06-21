#include "boundary_conditions.h"

void apply_bound_cond(input_t* input, REAL** species){
    
    for (int i = 0; i < input->n_species; ++i){
        for (int j = 0; j < input->pos_dims; ++j){
            input->pos_bounds[j](input->pos_dims+input->mom_dims, input->pos_points, j, input->padding, species[i]);
        }
        for (int j = 0; j < input->mom_dims; ++j){
            input->mom_bounds[j](input->pos_dims+input->mom_dims, input->pos_points, input->pos_dims+j, input->padding, species[i]);
        }
    }
}

void dirichelet_bound(int total_dims, int* Ns, int dim, int padding, REAL* species){

    int above = 1;
    int below = 1;
    int Nj = Ns[dim];
    for(int j = 0; j < dim; ++j){
        above *= Ns[j];
    }
    for(int j = dim+1; j < total_dims; ++j){
        below *= Ns[j];
    }

    for(int i = 0; i < above; ++i){
        for(int k = 0; k < below; ++k){
            for (int p = 0; p < padding; ++p){
                species[(i * Nj + p) * below + k] = 0;
                species[(Nj * (i + 1) - p - 1) * below + k] = 0;
            }
        }
    }
}

void periodic_bound(int total_dims, int* Ns, int dim, int padding, REAL* species){

    int above = 1;
    int below = 1;
    int Nj = Ns[dim];
    for(int j = 0; j < dim; ++j){
        above *= Ns[j];
    }
    for(int j = dim+1; j < total_dims; ++j){
        below *= Ns[j];
    }

    for(int i = 0; i < above; ++i){
        for(int k = 0; k < below; ++k){
            for (int p = 0; p < padding; ++p){
                species[(i * Nj + p) * below + k] = species[(Nj * (i + 1) - 2* padding + p) * below + k];
                species[(Nj * (i + 1) - padding + p) * below + k] = species[(i * Nj + padding + p) * below + k];
            }
        }
    }
}

void boundary_shift(input_t* input, COMPLEX* start, COMPLEX* final, int* aux){

    int count = 0;
    int inside = 0;

    for(int i = 0; i < input->pos_total; ++i){

        inside = 1;

        axis_index(input->pos_dims, input->pos_points, aux, i);

        for (int dim = 0; dim < input->pos_dims; ++dim){
            if(aux[dim] < input->padding || aux[dim] >= input->pos_points[dim] - input->padding){
                inside = 0;
            }
        }

        if (inside){
            final[count] = start[i];
            count++;
        }
    }

    for(int i = 0; i < input->pos_dims; ++i){
        input->pos_points[i] -= 2 * input->padding;
    }

}

void inverse_boundary_shift(input_t* input, COMPLEX* start, COMPLEX* final, int* aux){

    for(int i = 0; i < input->pos_dims; ++i){
        input->pos_points[i] += 2 * input->padding;
    }

    int count = 0;
    int inside = 0;

    for(int i = 0; i < input->pos_total; ++i){

        inside = 1;

        axis_index(input->pos_dims, input->pos_points, aux, i);

        for (int dim = 0; dim < input->pos_dims; ++dim){
            for(int p = 0; p < input->padding; ++p){
                if(aux[dim] == p){
                    inside = 0;
                    aux[dim] = input->pos_points[dim]- 2 * input-> padding + p;
                }
                if(aux[dim] == input->pos_points[dim] - input->padding + p){
                    inside = 0;
                    aux[dim] = input->padding + p;
                }
            }
            aux[dim] -= input->padding;
        }

        if (inside){
            final[i] = start[count];
            count++;
        }

        else{
            for(int j = 0; j < input->pos_dims; ++j){
                input->pos_points[j] -= 2 * input->padding;
            }
            final[i] = start[list_index(input->pos_dims, input->pos_points, aux)];
            for(int j = 0; j < input->pos_dims; ++j){
                input->pos_points[j] += 2 * input->padding;
            }
        }
    }

}