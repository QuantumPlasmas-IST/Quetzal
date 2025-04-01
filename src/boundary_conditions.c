#include "boundary_conditions.h"

void apply_bound_cond(input_t* input, REAL** species, REAL*** left_buffer, REAL*** right_buffer, int* aux_is){
    
    for (int i = 0; i < input->n_species; ++i){

        ghost_cell_transfer(input, species[i], left_buffer, right_buffer);

        for (int j = 0; j < input->pos_dims; ++j){
            if(input->parallel_pos[j]==0){
                input->pos_bounds[j](input, j, aux_is, species[i],0);
            }
            if(input->parallel_pos[j]==input->procs[j]-1){
                input->pos_bounds[j](input, j, aux_is, species[i],1);
            }
        }

        
        for (int j = 0; j < input->mom_dims; ++j){
            input->mom_bounds[j](input, input->pos_dims+j, aux_is, species[i], 0);
            input->mom_bounds[j](input, input->pos_dims+j, aux_is, species[i], 1);
        }
    }
}

void field_bound_cond(input_t* input, COMPLEX** fields, COMPLEX*** left_buffer, COMPLEX*** right_buffer, int* aux_is, int fld){

    field_cell_transfer(input, fields[fld], left_buffer, right_buffer);

    for (int j = 0; j < input->pos_dims; ++j){
        if(input->parallel_pos[j]==0){
            input->field_bounds[j](input, j, aux_is, fields[fld],0);
        }
        if(input->parallel_pos[j]==input->procs[j]-1){
            input->field_bounds[j](input, j, aux_is, fields[fld],1);
        }
    }
}

void dirichelet_bound(input_t* input, int dim, int* aux_is, REAL* species, int side){

    int above = 1;
    int below = 1;
    int Nj = input->pos_points[dim];
    for(int j = 0; j < dim; ++j){
        above *= input->pos_points[j];
    }
    for(int j = dim+1; j < input->pos_dims+input->mom_dims; ++j){
        below *= input->pos_points[j];
    }

    if(side){
        for(int i = 0; i < above; ++i){
            for(int k = 0; k < below; ++k){
                for (int p = 0; p < input->padding; ++p){
                    species[(Nj * (i + 1) - p - 1) * below + k] = 0;
                }
            }
        }
    }
    else{
        for(int i = 0; i < above; ++i){
            for(int k = 0; k < below; ++k){
                for (int p = 0; p < input->padding; ++p){
                    species[(i * Nj + p) * below + k] = 0;
                }
            }
        }
    }
}

void periodic_bound(input_t* input, int dim, int* aux_is, REAL* species, int side){

}

void wall_bound(input_t* input, int dim, int* aux_is, REAL* species, int side){

    if(dim >= input->mom_dims) return;

    int above = 1;
    int below = 1;
    int Nj = input->pos_points[dim];
    for(int j = 0; j < dim; ++j){
        above *= input->pos_points[j];
    }
    for(int j = dim+1; j < input->pos_dims; ++j){
        below *= input->pos_points[j];
    }

    int new_mp=0;

    if(side){
        for(int mp = 0; mp < input->mom_total; ++mp){
            //Calculate new momentum point for velocity reflection along 'dim' axis
            axis_index(input->mom_dims, input->mom_points, aux_is, mp);
            //printf("%d\n",aux_is[dim]);
            aux_is[dim] = -2 * (int)(round(input->mom_min[dim]/input->mom_delta[dim]))-aux_is[dim] + 2*input->padding;
            //printf("%d\n\n",aux_is[dim]);
            new_mp = list_index(input->mom_dims, input->mom_points, aux_is);
            //Fill ghost cells
            for(int ab = 0; ab < above; ++ab){
                for(int bl = 0; bl < below; ++bl){
                    for(int p = 0; p < input->padding; ++p){
                        // Positive boundary
                        species[(((ab+1) * Nj - input->padding + p) * below + bl) * input->mom_total + mp] = species[(((ab+1) * Nj - 2*input->padding + p) * below + bl) * input->mom_total + new_mp];
                    }
                }
            }
        }
    }
    else{
        for(int mp = 0; mp < input->mom_total; ++mp){
            //Calculate new momentum point for velocity reflection along 'dim' axis
            axis_index(input->mom_dims, input->mom_points, aux_is, mp);
            //printf("%d\n",aux_is[dim]);
            aux_is[dim] = -2 * (int)(round(input->mom_min[dim]/input->mom_delta[dim]))-aux_is[dim] + 2*input->padding;
            //printf("%d\n\n",aux_is[dim]);
            new_mp = list_index(input->mom_dims, input->mom_points, aux_is);
            //Fill ghost cells
            for(int ab = 0; ab < above; ++ab){
                for(int bl = 0; bl < below; ++bl){
                    for(int p = 0; p < input->padding; ++p){
                        // Negative Boundary
                        species[((ab * Nj + p) * below + bl) * input->mom_total + mp] = species[((ab * Nj + 2*input->padding-p-1) * below + bl) * input->mom_total + new_mp];
                    }
                }
            }
        }
    }
}

void field_dirichelet_bound(input_t* input, int dim, int* aux_is, COMPLEX* fields, int side){

    int above = 1;
    int below = 1;
    int Nj = input->pos_points[dim];
    for(int j = 0; j < dim; ++j){
        above *= input->pos_points[j];
    }
    for(int j = dim+1; j < input->pos_dims; ++j){
        below *= input->pos_points[j];
    }

    if(side){
        for(int i = 0; i < above; ++i){
            for(int k = 0; k < below; ++k){
                for (int p = 0; p < input->padding; ++p){
                    fields[(Nj * (i + 1) - p - 1) * below + k] = 0;
                }
            }
        }
    }
    else{
        for(int i = 0; i < above; ++i){
            for(int k = 0; k < below; ++k){
                for (int p = 0; p < input->padding; ++p){
                    fields[(i * Nj + p) * below + k] = 0;
                }
            }
        }
    }
}

void field_periodic_bound(input_t* input, int dim, int* aux_is, COMPLEX* fields, int side){

}

void field_neumann_bound(input_t* input, int dim, int* aux_is, COMPLEX* fields, int side){

    int above = 1;
    int below = 1;
    int Nj = input->pos_points[dim];
    for(int j = 0; j < dim; ++j){
        above *= input->pos_points[j];
    }
    for(int j = dim+1; j < input->pos_dims; ++j){
        below *= input->pos_points[j];
    }

    if(side){
        for(int i = 0; i < above; ++i){
            for(int k = 0; k < below; ++k){
                for (int p = 0; p < input->padding; ++p){
                    fields[(Nj * (i + 1) - p - 1) * below + k] = fields[(Nj * (i + 1) - input->padding) * below + k];
                }
            }
        }
    }
    else{
        for(int i = 0; i < above; ++i){
            for(int k = 0; k < below; ++k){
                for (int p = 0; p < input->padding; ++p){
                    fields[(i * Nj + p) * below + k] = fields[(i * Nj + input->padding) * below + k];
                }
            }
        }
    }
}

void field_second_bound(input_t* input, int dim, int* aux_is, COMPLEX* fields, int side){

    int above = 1;
    int below = 1;
    int Nj = input->pos_points[dim];
    for(int j = 0; j < dim; ++j){
        above *= input->pos_points[j];
    }
    for(int j = dim+1; j < input->pos_dims; ++j){
        below *= input->pos_points[j];
    }

    if(side){
        for(int i = 0; i < above; ++i){
            for(int k = 0; k < below; ++k){
                for (int p = 0; p < input->padding; ++p){
                    fields[(Nj * (i + 1) - p - 1) * below + k] = 2 * fields[(Nj * (i + 1) - input->padding) * below + k] - fields[(Nj * (i + 1) - input->padding - 1) * below + k];
                }
            }
        }
    }
    else{
        for(int i = 0; i < above; ++i){
            for(int k = 0; k < below; ++k){
                for (int p = 0; p < input->padding; ++p){
                    fields[(i * Nj + p) * below + k] = 2 * fields[(i * Nj + input->padding) * below + k] - fields[(i * Nj + input->padding + 1) * below + k];
                }
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

        //Fill buffers
        for (int p = 0; p < input->padding; ++p){
            for(int j = 0; j < above; ++j){
                for(int k = 0; k < below; ++k){
                    left_buffer[i][p][j*below + k] = species[(j * Nj + input->padding + p) * below + k];
                    right_buffer[i][p][j*below + k] = species[(Nj * (j + 1) - 2* input->padding + p) * below + k];
                }
            }
        }

        // Send data to buffer of neighbors
        MPI_Barrier(MPI_COMM_WORLD);
        for (int p = 0; p < input->padding; ++p){
            MPI_Sendrecv_replace(left_buffer[i][p], above*below, QTZ_MPI_REAL, left_proc, input->rank, right_proc, right_proc, MPI_COMM_WORLD, &Stat[2*(i*input->padding+p)]);
            MPI_Sendrecv_replace(right_buffer[i][p], above*below, QTZ_MPI_REAL, right_proc, input->rank, left_proc, left_proc, MPI_COMM_WORLD, &Stat[2*(i*input->padding+p)+1]);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Fill ghost-cells with buffers
        for (int p = 0; p < input->padding; ++p){
            for(int j = 0; j < above; ++j){
                for(int k = 0; k < below; ++k){
                    species[(Nj * (j + 1) - input->padding + p) * below + k] = left_buffer[i][p][j*below + k];
                    species[(j * Nj + p) * below + k] = right_buffer[i][p][j*below + k];
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

        // Calculate Number of points
        int above = 1;
        int below = 1;
        int Nj = input->pos_points[i];
        for(int j = 0; j < i; ++j){
            above *= input->pos_points[j];
        }
        for(int j = i+1; j < input->pos_dims; ++j){
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

        //Fill buffers
        for (int p = 0; p < input->padding; ++p){
            for(int j = 0; j < above; ++j){
                for(int k = 0; k < below; ++k){
                    left_buffer[i][p][j*below + k] = field[(j * Nj + input->padding + p) * below + k];
                    right_buffer[i][p][j*below + k] = field[(Nj * (j + 1) - 2* input->padding + p) * below + k];
                }
            }
        }

        // Send data to buffer of neighbors
        MPI_Barrier(MPI_COMM_WORLD);
        for (int p = 0; p < input->padding; ++p){
            MPI_Sendrecv_replace(left_buffer[i][p], above*below, QTZ_MPI_COMPLEX, left_proc, input->rank, right_proc, right_proc, MPI_COMM_WORLD, &Stat[2*(i*input->padding+p)]);
            MPI_Sendrecv_replace(right_buffer[i][p], above*below, QTZ_MPI_COMPLEX, right_proc, input->rank, left_proc, left_proc, MPI_COMM_WORLD, &Stat[2*(i*input->padding+p)+1]);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Fill ghost-cells with buffers
        for (int p = 0; p < input->padding; ++p){
            for(int j = 0; j < above; ++j){
                for(int k = 0; k < below; ++k){
                    field[(Nj * (j + 1) - input->padding + p) * below + k] = left_buffer[i][p][j*below + k];
                    field[(j * Nj + p) * below + k] = right_buffer[i][p][j*below + k];
                }
            }
        }  
    }

    MPI_Barrier(MPI_COMM_WORLD);
}