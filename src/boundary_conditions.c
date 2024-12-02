#include "boundary_conditions.h"

void apply_bound_cond(input_t* input, REAL** species, REAL*** left_buffer, REAL*** right_buffer){
    
    for (int i = 0; i < input->n_species; ++i){

        ghost_cell_transfer(input, species[i], left_buffer, right_buffer);

        for (int j = 0; j < input->pos_dims; ++j){
            input->pos_bounds[j](input->pos_dims+input->mom_dims, input->pos_points, j, input->padding, species[i]);
        }


        for (int j = 0; j < input->mom_dims; ++j){
            input->mom_bounds[j](input->pos_dims+input->mom_dims, input->pos_points, input->pos_dims+j, input->padding, species[i]);
        }
    }
}

void field_bound_cond(input_t* input, COMPLEX** fields, COMPLEX*** left_buffer, COMPLEX*** right_buffer, int fld){

    field_cell_transfer(input, fields[fld], left_buffer, right_buffer);

    for (int ps = 0; ps < input->pos_total; ++ps){
        //fields[fld][ps] = creal(fields[fld][ps]);
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


    /*int above = 1;
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
    }*/

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

    //for(int i = 0; i < input->pos_dims; ++i){
    //    input->pos_points[i] -= 2 * input->padding;
    //}

}

void inverse_boundary_shift(input_t* input, COMPLEX* start, COMPLEX* final, int* aux){

    //for(int i = 0; i < input->pos_dims; ++i){
    //    input->pos_points[i] += 2 * input->padding;
    //}

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

void species_boundary_shift(input_t* input, REAL* start, REAL* final, int* aux){

    int count = 0;
    int inside = 0;

    for(int i = 0; i < input->pos_total * input->mom_total; ++i){

        inside = 1;

        axis_index(input->pos_dims+input->mom_dims, input->pos_points, aux, i);

        for (int dim = 0; dim < input->pos_dims+input->mom_dims; ++dim){
            if(aux[dim] < input->padding || aux[dim] >= input->pos_points[dim] - input->padding){
                inside = 0;
            }
        }

        if (inside){
            final[count] = start[i];
            count++;
        }
    }

    for(int i = 0; i < input->pos_dims+input->mom_dims; ++i){
        input->pos_points[i] -= 2 * input->padding;
    }

}

void species_inverse_boundary_shift(input_t* input, REAL* start, REAL* final, int* aux){

    for(int i = 0; i < input->pos_dims+input->mom_dims; ++i){
        input->pos_points[i] += 2 * input->padding;
    }

    int count = 0;
    int inside = 0;

    for(int i = 0; i < input->pos_total * input->mom_total; ++i){

        inside = 1;

        axis_index(input->pos_dims+input->mom_dims, input->pos_points, aux, i);

        for (int dim = 0; dim < input->pos_dims+input->mom_dims; ++dim){
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
            for(int j = 0; j < input->pos_dims+input->mom_dims; ++j){
                input->pos_points[j] -= 2 * input->padding;
            }
            final[i] = start[list_index(input->pos_dims+input->mom_dims, input->pos_points, aux)];
            for(int j = 0; j < input->pos_dims+input->mom_dims; ++j){
                input->pos_points[j] += 2 * input->padding;
            }
        }
    }
}

void ghost_cell_transfer(input_t* input, REAL* species, REAL*** left_buffer, REAL*** right_buffer){

    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Request Req[2*input->pos_dims*input->padding];
    MPI_Status  Stat[2*input->pos_dims*input->padding];

    for(int i = 0; i < input->pos_dims; ++i){

        // Calculate Number of points
        int above = 1;
        int below = 1;
        int Nj = input->pos_points[i];
        for(int j = 0; j < i; ++j){
            above *= input->pos_points[j];
        }
        for(int j = i+1; j < input->pos_dims+input->mom_dims; ++j){
            below *= input->pos_points[j];
        }

        // Find neighboring processes
        int left_proc = input->rank - input->parallel_factor[i];
        int right_proc = input->rank + input->parallel_factor[i];
        if(input->parallel_pos[i]%input->procs[i] == 0){
            left_proc = input->rank + (input->procs[i]-1) * input->parallel_factor[i];
        }
        if(input->parallel_pos[i]%input->procs[i] == input->procs[i]-1){
            right_proc = input->rank - (input->procs[i]-1) * input->parallel_factor[i];
        }

        //Fill left-facing buffer
        for (int p = 0; p < input->padding; ++p){
            for(int j = 0; j < above; ++j){
                for(int k = 0; k < below; ++k){
                    left_buffer[i][p][j*below + k] = species[(Nj * (j + 1) - 2* input->padding + p) * below + k];
                }
            }
        }

        // Send left-facing data to right-facing buffer of left neighbor
        for (int p = 0; p < input->padding; ++p){
            MPI_Send(left_buffer[i][p], above*below, QTZ_MPI_REAL, left_proc, input->rank, MPI_COMM_WORLD);
            MPI_Recv(right_buffer[i][p], above*below, QTZ_MPI_REAL, right_proc, right_proc, MPI_COMM_WORLD, &Stat[2*(i*input->padding+p)+1]);
        }

        // Fill ghost-cells with right-facing buffer
        for (int p = 0; p < input->padding; ++p){
            for(int j = 0; j < above; ++j){
                for(int k = 0; k < below; ++k){
                    species[(Nj * (j + 1) - input->padding + p) * below + k] = right_buffer[i][p][j*below + k];
                }
            }
        }

        // Fill right-facing buffer
        for (int p = 0; p < input->padding; ++p){
            for(int j = 0; j < above; ++j){
                for(int k = 0; k < below; ++k){
                    right_buffer[i][p][j*below + k] = species[(j * Nj + input->padding + p) * below + k];
                }
            }
        }

        // Send right-facing data to left-facing buffer of right neighbor
        for (int p = 0; p < input->padding; ++p){
            MPI_Send(right_buffer[i][p], above*below, QTZ_MPI_REAL, right_proc, input->rank, MPI_COMM_WORLD);
            MPI_Recv(left_buffer[i][p], above*below, QTZ_MPI_REAL, left_proc, left_proc, MPI_COMM_WORLD, &Stat[2*(i*input->padding+p)+1]);
        }

        // Fill ghost-cells with left-facing buffer
        for (int p = 0; p < input->padding; ++p){
            for(int j = 0; j < above; ++j){
                for(int k = 0; k < below; ++k){
                    species[(j * Nj + p) * below + k] = left_buffer[i][p][j*below + k];
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

void field_cell_transfer(input_t* input, COMPLEX* field, COMPLEX*** left_buffer, COMPLEX*** right_buffer){

    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Request Req[2*input->pos_dims*input->padding];
    MPI_Status  Stat[2*input->pos_dims*input->padding];

    for(int i = 0; i < input->pos_dims; ++i){

        int above = 1;
        int below = 1;
        int Nj = input->pos_points[i];
        for(int j = 0; j < i; ++j){
            above *= input->pos_points[j];
        }
        for(int j = i+1; j < input->pos_dims; ++j){
            below *= input->pos_points[j];
        }

        for (int p = 0; p < input->padding; ++p){
            for(int j = 0; j < above; ++j){
                for(int k = 0; k < below; ++k){

                    left_buffer[i][p][j*below + k] = field[(Nj * (j + 1) - 2* input->padding + p) * below + k];
                    right_buffer[i][p][j*below + k] = field[(j * Nj + input->padding + p) * below + k];
                }
            }
        }

        int left_proc = input->rank - input->parallel_factor[i];
        int right_proc = input->rank + input->parallel_factor[i];
        if(input->parallel_pos[i]%input->procs[i] == 0){
            left_proc = input->rank + (input->procs[i]-1) * input->parallel_factor[i];
        }
        if(input->parallel_pos[i]%input->procs[i] == input->procs[i]-1){
            right_proc = input->rank - (input->procs[i]-1) * input->parallel_factor[i];
        }

        for (int p = 0; p < input->padding; ++p){
            
            //MPI_Isendrecv_replace(left_buffer[i][p], above*below, MPI_DOUBLE, left_proc, input->rank, right_proc, right_proc, MPI_COMM_WORLD, &Req[2*(i*input->padding+p)]);
            //MPI_Isendrecv_replace(right_buffer[i][p], above*below, MPI_DOUBLE, right_proc, input->rank, left_proc, left_proc, MPI_COMM_WORLD, &Req[2*(i*input->padding+p)+1]);

            MPI_Sendrecv_replace(left_buffer[i][p], above*below, QTZ_MPI_COMPLEX, left_proc, input->rank, right_proc, right_proc, MPI_COMM_WORLD, &Stat[2*(i*input->padding+p)]);
            MPI_Sendrecv_replace(right_buffer[i][p], above*below, QTZ_MPI_COMPLEX, right_proc, input->rank, left_proc, left_proc, MPI_COMM_WORLD, &Stat[2*(i*input->padding+p)+1]);

        }

        //MPI_Waitall(2*input->padding, &Req[2*i*input->padding], &Stat[2*i*input->padding]);

        for (int p = 0; p < input->padding; ++p){
            for(int j = 0; j < above; ++j){
                for(int k = 0; k < below; ++k){

                    field[(j * Nj + p) * below + k] = left_buffer[i][p][j*below + k];
                    field[(Nj * (j + 1) - input->padding + p) * below + k] = right_buffer[i][p][j*below + k];
                }
            }
        }

    }

    MPI_Barrier(MPI_COMM_WORLD);
}