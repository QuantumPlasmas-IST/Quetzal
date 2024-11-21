#include "pusher.h"

REAL** solve(input_t* input){

    // Setting up MPI-FFT //

    // Sizes and ranges of axes for each process
    int* fft_ns = malloc(sizeof(int) * input->pos_dims);
    int* fft_hi = malloc(sizeof(int) * input->pos_dims);
    int* fft_lo = malloc(sizeof(int) * input->pos_dims);
    
    int fftsize,sendsize,recvsize;
    int pos_size = input->pos_total;

    for (int i = 0; i < input->pos_dims; ++i){
        fft_lo[i] = input->parallel_pos[i] * input->fft_points[i];
        fft_hi[i] = (input->parallel_pos[i]+1) * input->fft_points[i] - 1;
        fft_ns[i] = input->procs[i] * input->fft_points[i];
    }
    
    if(input->pos_dims == 2){
        fft2d_create(MPI_COMM_WORLD, 2, &input->fft);
        fft2d_setup(input->fft, fft_ns[1], fft_ns[0], fft_lo[1], fft_hi[1], fft_lo[0], fft_hi[0], fft_lo[1], fft_hi[1], fft_lo[0], fft_hi[0], 0, &fftsize, &sendsize, &recvsize);
        if(fftsize > input->pos_total){
            pos_size = fftsize;
        }
    }

    if(input->pos_dims == 3){
        fft3d_create(MPI_COMM_WORLD, 2, &input->fft);
        fft3d_setup(input->fft, fft_ns[2], fft_ns[1], fft_ns[0], fft_lo[2], fft_hi[2], fft_lo[1], fft_hi[1], fft_lo[0], fft_hi[0], fft_lo[2], fft_hi[2], fft_lo[1], fft_hi[1], fft_lo[0], fft_hi[0], 0, &fftsize, &sendsize, &recvsize);
        if(fftsize > input->pos_total){
            pos_size = fftsize;
        }
    }

    free(fft_ns);
    free(fft_hi);
    free(fft_lo);

    // Allocating all the memory //

    REAL** species = malloc(input->n_species * sizeof(REAL*));
    REAL** species_aux = malloc(input->n_species * sizeof(REAL*));

    COMPLEX** fields = malloc(input->n_fields * sizeof(COMPLEX*));
    COMPLEX** fields_aux = malloc(input->n_fields * sizeof(COMPLEX*));
    COMPLEX** fields_deriv = malloc(input->n_fields * sizeof(COMPLEX*));
    COMPLEX** fields_fft = malloc(input->n_fields * sizeof(COMPLEX*));

    COMPLEX** sources = malloc(input->n_species * (1+input->mom_dims) * sizeof(COMPLEX*));
    COMPLEX** sources_aux = malloc(input->n_species * (1+input->mom_dims) * sizeof(COMPLEX*));

    REAL*** flows = malloc(input->n_species * sizeof(REAL**));

    REAL*** left_ghost_buffer = malloc(input->pos_dims * sizeof(REAL**));
    REAL*** right_ghost_buffer = malloc(input->pos_dims * sizeof(REAL**));

    COMPLEX*** left_field_buffer = malloc(input->pos_dims * sizeof(COMPLEX**));
    COMPLEX*** right_field_buffer = malloc(input->pos_dims * sizeof(COMPLEX**));
        
    for(int i = 0; i < input->pos_dims; ++i){

        int buffer_size = input->pos_total*input->mom_total / input->pos_points[i];
        int buffer_size2 = input->pos_total / input->pos_points[i];

        left_ghost_buffer[i]   =  malloc(input->padding * sizeof(REAL*));
        right_ghost_buffer[i]  =  malloc(input->padding * sizeof(REAL*));
        left_field_buffer[i]   =  malloc(input->padding * sizeof(COMPLEX*));
        right_field_buffer[i]  =  malloc(input->padding * sizeof(COMPLEX*));

        for(int p = 0; p < input->padding; ++p){
            left_ghost_buffer[i][p]   =  malloc(buffer_size * sizeof(REAL));
            right_ghost_buffer[i][p]  =  malloc(buffer_size * sizeof(REAL));
            left_field_buffer[i][p]   =  malloc(buffer_size * sizeof(COMPLEX));
            right_field_buffer[i][p]  =  malloc(buffer_size * sizeof(COMPLEX));
        }
    }

    species[0] = malloc(input->n_species * input->pos_total * input->mom_total * sizeof(REAL));
    species_aux[0] = malloc(input->n_species * input->pos_total * input->mom_total * sizeof(REAL));

    sources[0] = malloc(input->n_species * (1+input->mom_dims) * pos_size * sizeof(COMPLEX));
    sources_aux[0] = malloc(input->n_species * (1+input->mom_dims) * pos_size * sizeof(COMPLEX));
    
    fields[0] = malloc(input->n_fields * pos_size * sizeof(COMPLEX));
    fields_aux[0] = malloc(input->n_fields * pos_size * sizeof(COMPLEX));
    fields_deriv[0] = malloc(input->n_fields * pos_size * sizeof(COMPLEX));
    fields_fft[0] = malloc(input->n_fields * pos_size * sizeof(COMPLEX));

    for (int i = 0; i < input->n_species; ++i){
        species[i] = species[0] + i * input->pos_total * input->mom_total;
        species_aux[i] = species_aux[0] + i * input->pos_total * input->mom_total;

        for(int d = 0; d < input->mom_dims+1; ++d){
            sources[i * (1+input->mom_dims) + d] = sources[0] + (i * (1+input->mom_dims) + d) * pos_size;
            sources_aux[i * (1+input->mom_dims) + d] = sources_aux[0] + (i * (1+input->mom_dims) + d) * pos_size;
        }
    
        flows[i] = malloc((input->pos_dims + input->mom_dims) * sizeof(REAL*));
        for(int dim = 0; dim < input->pos_dims+input->mom_dims; ++dim){
            flows[i][dim] = malloc((input->pos_total * input->mom_total) * sizeof(REAL));
        }
    }

    for (int i = 0; i < input->n_fields; ++i){
        fields[i] = fields[0] + i * pos_size;
        fields_aux[i] = fields_aux[0] + i * pos_size;
        fields_deriv[i] = fields_deriv[0] + i * pos_size;
        fields_fft[i] = fields_fft[0] + i * pos_size;
    }

    int* is = malloc(sizeof(int) * (input->pos_dims+input->mom_dims));
    REAL* aux_momentum = malloc((input->pos_dims+input->mom_dims) * sizeof(REAL));

    // Applying initial and boundary conditions conditions //

    apply_init_cond(input, species);
    apply_init_cond(input, species_aux);
    field_init_cond(input, fields);
    if(!input->rank) printf("Initial Conditions Applied\n");
    apply_bound_cond(input, species, left_ghost_buffer, right_ghost_buffer);
    apply_bound_cond(input, species_aux, left_ghost_buffer, right_ghost_buffer);
    for (int fld = 0; fld < input->n_fields; ++fld) field_bound_cond(input, fields, left_field_buffer, right_field_buffer, fld);
    if(!input->rank) printf("Boundary Conditions Applied\n");
    integrate_source(input, species, sources, is, aux_momentum);

    if(!input->rank) printf("Initialization Complete\n");

    // Algorithm //

    if(!strcmp(input->pusher,"runge kutta")){
        for(int i = 0; i < input->n_timesteps; ++i){

            if(!(i%input->pos_diag_freq)){
                write_sources(input, sources, sources_aux, is, i);
                write_fields(input, fields, fields_fft, is, i);
                if(!input->rank){
                    printf("Timestep %d\n",i);
                    fflush(stdout);
                }
            }

            if(!(i%input->mom_diag_freq)){
                write_solution(input, species, species_aux, is, i);
                apply_bound_cond(input, species, left_ghost_buffer, right_ghost_buffer);
            }

            rungeKutta2(input, species, species_aux, sources, fields, fields_aux, fields_fft, fields_deriv, flows, is, aux_momentum, left_ghost_buffer, right_ghost_buffer, left_field_buffer, right_field_buffer);

        }
    }

    // Free all the memory

    if(input->pos_dims==2) fft2d_destroy(input->fft);
    if(input->pos_dims==3) fft3d_destroy(input->fft);

    free(is);
    free(aux_momentum);

    for (int i = 0; i < input->n_species; ++i){
        for(int dim = 0; dim < input->pos_dims+input->mom_dims; ++dim){
            free(flows[i][dim]);
        }
        free(flows[i]);
    }
    free(flows);

    for (int dim = 0; dim < input->pos_dims; ++dim){
        for (int p = 0; p < input->padding; ++p){
            free(left_field_buffer[dim][p]);
            free(right_field_buffer[dim][p]);
            free(left_ghost_buffer[dim][p]);
            free(right_ghost_buffer[dim][p]);
        }
        free(left_field_buffer[dim]);
        free(right_field_buffer[dim]);
        free(left_ghost_buffer[dim]);
        free(right_ghost_buffer[dim]);
    }
    free(left_field_buffer);
    free(right_field_buffer);
    free(left_ghost_buffer);
    free(right_ghost_buffer);

    free(species_aux[0]);
    free(sources[0]);
    free(sources_aux[0]);
    free(fields[0]);
    free(fields_aux[0]);
    free(species_aux);
    free(sources);
    free(sources_aux);
    free(fields);
    free(fields_aux);

    return species;
}

void leapfrog(input_t* input, REAL** species, REAL** aux1, REAL** aux2, COMPLEX**fields, int* aux_is){

    REAL* momentum = malloc(input->mom_dims * sizeof(REAL));

    for(int i = 0; i < input->n_species; ++i){
        for(int k = 0; k < input->mom_total; ++k){

            axis_index(input->pos_dims + input->mom_dims, input->pos_points, aux_is, k);
            for(int j = 0; j < input->mom_dims; ++j){
                momentum[j] = input->mom_min[j] + aux_is[j+input->pos_dims] * input->mom_delta[j];
            }           
            input->dispersions[i](input->pos_dims, input->mom_dims, momentum); 

            for(int j = 0; j < input->pos_total; ++j){
                aux2[i][j*input->mom_total + k] = aux1[i][j*input->mom_total + k] + 2*input->deltaT*central_diff(input, species, fields, i, j*input->mom_total + k, aux_is, momentum);
            }
        }
    }
}

void integrate_source(input_t* input, REAL** species, COMPLEX** sources, int* aux_is, REAL* aux_momentum){
    
    for(int i = 0; i < input->n_species * (1+input->mom_dims); ++i){
        for(int j = 0; j < input->pos_total; ++j){
            sources[i][j] = 0;
        }
    }

    for(int i = 0; i < input->n_species; ++i){

        for(int k = 0; k < input->mom_total; ++k){

            // Calculate velocity
            axis_index(input->mom_dims, input->pos_points + input->pos_dims, aux_is + input->pos_dims, k);
            for(int mdim = 0; mdim < input->mom_dims; ++mdim){
                aux_momentum[mdim] = input->mom_min[mdim] + (aux_is[mdim+input->pos_dims] - input->padding) * input->mom_delta[mdim];
            }           

            input->dispersions[i](input->pos_dims, input->mom_dims, aux_momentum);

            for(int j = 0; j < input->pos_total; ++j){
                sources[i*(1+input->mom_dims)][j] += species[i][j*input->mom_total + k];
                for(int d = 0; d < input->mom_dims; ++d){
                    sources[i*(1+input->mom_dims)+1+d][j] += aux_momentum[d] * species[i][j*input->mom_total + k];
                }
            }
        }
    }

    for(int i = 0; i < input->n_species * (1+input->mom_dims); ++i){
        for(int j = 0; j < input->pos_total; ++j){
            for(int d = 0; d < input->mom_dims; ++d){
                sources[i][j] *= input->mom_delta[d];
            }
        }
    }
}

void convolve_source(input_t* input, COMPLEX** sources, COMPLEX** fields_fft, COMPLEX** fields_deriv, int* aux_is, REAL* aux_wavevector, int fld){

    // Combine appropriate moments of each species for each field
    for(int ps = 0; ps < input->pos_total; ++ps){

        fields_fft[fld][ps] = 0;

        for(int spc = 0; spc < input->n_species; ++spc){
            fields_fft[fld][ps] += input->source_charges[spc][fld] * sources[spc*(1+input->mom_dims)+input->source_moments[spc][fld]][ps];
        }
    }

    //Calculate FFT
    boundary_shift(input, fields_fft[fld], fields_deriv[fld], aux_is);

    switch (input->pos_dims){
        case 1:
            fastFourier(input->fft_points[0], fields_deriv[fld], fields_fft[fld],1);
            break;
        case 2:
            fft2d_compute(input->fft, (double*) fields_deriv[fld], (double*) fields_fft[fld], 1);
            break;
        case 3:
            fft3d_compute(input->fft, (double*) fields_deriv[fld], (double*) fields_fft[fld], 1);
            break;
    }

    for(int ps = 0; ps < input->pos_total; ++ps){
        
        // Calculate wavevectors
        axis_index(input->pos_dims, input->fft_points, aux_is, ps);

        for(int dim = 0; dim < input->pos_dims+input->mom_dims; ++dim){
            aux_wavevector[dim] = 0;
        }

        for(int pdim = 0; pdim < input->pos_dims; ++pdim){

            aux_is[pdim] += input->parallel_pos[pdim] * input->fft_points[pdim];

            if (aux_is[pdim] < (input->fft_points[pdim]*input->procs[pdim])/2){
                aux_wavevector[pdim] = input->dk[pdim] * aux_is[pdim];
            }
            else{
                aux_wavevector[pdim] = input->dk[pdim] * (-input->fft_points[pdim]*input->procs[pdim] + aux_is[pdim]);
            }
        }

        //Convolve the source with the kernel
        fields_deriv[fld][ps] = input->kernels[fld](input->pos_dims, fields_fft[fld][ps], aux_wavevector);
    }

}

void transform_field(input_t* input, COMPLEX** fields, COMPLEX** fields_fft, int* aux_is, int fld){

    boundary_shift(input, fields[fld], fields_fft[fld], aux_is);

    switch (input->pos_dims){
        case 1:
            //NdFourier(input->pos_dims, input->pos_points, fields[k], sources[k]);
            fastFourier(input->fft_points[0], fields_fft[fld], fields[fld],1);
            break;
        case 2:
            fft2d_compute(input->fft, (double*) fields_fft[fld], (double*) fields[fld], 1);
            break;
        case 3:
            fft3d_compute(input->fft, (double*) fields_fft[fld], (double*) fields[fld], 1);
            break;
    }

}

void invert_field(input_t* input, COMPLEX** fields, COMPLEX** fields_fft, int* aux_is, int fld){

    switch (input->pos_dims){
        case 1:
            //NdFourier(input->pos_dims, input->pos_points, fields[k], sources[k]);
            fastInverseFourier(input->fft_points[0], fields[fld], fields_fft[fld],1);
            break;
        case 2:
            fft2d_compute(input->fft, (double*) fields[fld], (double*) fields_fft[fld], -1);
            break;
        case 3:
            fft3d_compute(input->fft, (double*) fields[fld], (double*) fields_fft[fld], -1);
            break;
    }

    inverse_boundary_shift(input, fields_fft[fld], fields[fld], aux_is);

}

void field_dynamics(input_t* input, COMPLEX** fields, COMPLEX** fields_deriv, int* aux_is, REAL* aux_wavevector, int fld){
        
    for(int ps = 0; ps < input->pos_total; ++ps){
        
        // Calculate wavevectors
        axis_index(input->pos_dims, input->fft_points, aux_is, ps);

        for(int dim = 0; dim < input->pos_dims+input->mom_dims; ++dim){
            aux_wavevector[dim] = 0;
        }

        for(int pdim = 0; pdim < input->pos_dims; ++pdim){

            aux_is[pdim] += input->parallel_pos[pdim] * input->fft_points[pdim];

            if (aux_is[pdim] < (input->fft_points[pdim]*input->procs[pdim])/2){
                aux_wavevector[pdim] = input->dk[pdim] * aux_is[pdim];
            }
            else{
                aux_wavevector[pdim] = input->dk[pdim] * (-input->fft_points[pdim]*input->procs[pdim] + aux_is[pdim]);
            }
        }

        // Calculate the field derivative
        fields_deriv[fld][ps] += input->dynamics[fld](fld, ps, fields, aux_wavevector);
    }
}

REAL central_diff(input_t* input, REAL** species, COMPLEX** fields, int i_sp, int index, int* aux_is, REAL* velocity){

    axis_index(input->pos_dims + input->mom_dims, input->pos_points, aux_is, index);

    for (int dim = 0; dim < input->pos_dims+input->mom_dims; ++dim){
        if(aux_is[dim] == 0 || aux_is[dim] == input->pos_points[dim]-1) return 0;
    }

    REAL res = 0;

    for(int j = 0; j < input->pos_dims; ++j){
        aux_is[j] += 1;
        REAL forw = species[i_sp][list_index(input->pos_dims+input->mom_dims,input->pos_points,aux_is)];
        aux_is[j] -= 2;
        REAL back = species[i_sp][list_index(input->pos_dims+input->mom_dims,input->pos_points,aux_is)];
        aux_is[j] += 1;
        res += velocity[j] * (forw - back) / (2 * input->pos_delta[j]); 
    }   

    for(int j = 0; j < input->mom_dims; ++j){
        aux_is[input->pos_dims + j] += 1;
        REAL forw = species[i_sp][list_index(input->pos_dims+input->mom_dims,input->pos_points,aux_is)];
        aux_is[input->pos_dims + j] -= 2;
        REAL back = species[i_sp][list_index(input->pos_dims+input->mom_dims,input->pos_points,aux_is)];
        aux_is[input->pos_dims + j] += 1;

        REAL force = 0;
        for(int k = 0; k < input->n_fields; ++k){
            aux_is[j] += 1;
            REAL e_forw = creal(fields[k][list_index(input->pos_dims,input->pos_points,aux_is)]);
            aux_is[j] -= 2;
            REAL e_back = creal(fields[k][list_index(input->pos_dims,input->pos_points,aux_is)]);
            aux_is[j] += 1;
            force -= (e_forw - e_back) / (2 * input->pos_delta[j]);
        }

        res += force * (forw - back) / (2 * input->mom_delta[j]); 
    }   

    return -res;    

}

void finite_volume2(input_t* input, REAL** species, COMPLEX** fields, int* aux_is, REAL* aux_momentum, REAL*** flows){

    REAL cent, forw, back, delta, force;
    
    // Loop through different species
    for(int i = 0; i < input->n_species; ++i){
        //Loop through different momenta
        for(int k = 1; k < input->mom_total-2; ++k){
            
            // Calculate velocity of this momentum point
            axis_index(input->mom_dims, input->pos_points + input->pos_dims, aux_is + input->pos_dims, k);
            for(int mdim = 0; mdim < input->mom_dims; ++mdim){
                aux_momentum[mdim] = input->mom_min[mdim] + (aux_is[mdim+input->pos_dims] - input->padding) * input->mom_delta[mdim];
            }           
            input->dispersions[i](input->pos_dims, input->mom_dims, aux_momentum); 

            //Loop through different positions
            for(int j = input->space_factor[0]; j < input->pos_total-2*input->space_factor[0]; ++j){
                
                for(int pdim = 0; pdim < input->pos_dims; ++pdim){

                    if(pdim == input->mom_dims){
                        flows[i][pdim][j*input->mom_total + k] = 0;
                        break;
                    }

                    if(aux_momentum[pdim]>=0){
                        
                        cent = species[i][j*input->mom_total+k];
                        forw = species[i][j*input->mom_total+k+input->grid_factor[pdim]];
                        back = species[i][j*input->mom_total+k-input->grid_factor[pdim]];

                        delta = (forw - cent) * MC_fluxLimiter(back, cent, forw);

                        flows[i][pdim][j*input->mom_total + k] = aux_momentum[pdim] * (cent + delta/2);
                    }

                    else{

                        back = species[i][j*input->mom_total+k];
                        cent = species[i][j*input->mom_total+k+input->grid_factor[pdim]];
                        forw = species[i][j*input->mom_total+k+2*input->grid_factor[pdim]];

                        delta = (forw - cent) * MC_fluxLimiter(back, cent, forw);

                        flows[i][pdim][j*input->mom_total + k] = aux_momentum[pdim] * (cent - delta/2);
                    }
                }   

                for(int mdim = 0; mdim < input->mom_dims; ++mdim){

                    force = 0;

                    for(int fld = 0; fld < input->n_fields; ++fld){
                        force += input->force_charges[i][fld] * input->forces[fld](input, fields, fld, j, k, mdim);
                    }

                    if(force>=0){
                        
                        cent = species[i][j*input->mom_total+k];
                        forw = species[i][j*input->mom_total+k+input->grid_factor[input->pos_dims+mdim]];
                        back = species[i][j*input->mom_total+k-input->grid_factor[input->pos_dims+mdim]];

                        delta = (forw - cent) * MC_fluxLimiter(back, cent, forw);

                        flows[i][input->pos_dims + mdim][j*input->mom_total + k] = force * (cent + delta/2);
                    }

                    else{

                        back = species[i][j*input->mom_total+k];
                        cent = species[i][j*input->mom_total+k+input->grid_factor[input->pos_dims+mdim]];
                        forw = species[i][j*input->mom_total+k+2*input->grid_factor[input->pos_dims+mdim]];

                        delta = (forw - cent) * MC_fluxLimiter(back, cent, forw);

                        flows[i][input->pos_dims + mdim][j*input->mom_total + k] = force * (cent - delta/2);
                    }
                }  
            }
        }
    }
}

void finite_volumeNL2(input_t* input, REAL** species, COMPLEX** fields, int* aux_is, REAL* aux_momentum, REAL*** flows){

    REAL cent, forw, back, delta, force;
    
    // Loop through different species
    for(int i = 0; i < input->n_species; ++i){
        //Loop through different momenta
        for(int k = 1; k < input->mom_total-2; ++k){
            
            // Calculate velocity of this momentum point
            axis_index(input->mom_dims, input->pos_points + input->pos_dims, aux_is + input->pos_dims, k);
            for(int mdim = 0; mdim < input->mom_dims; ++mdim){
                aux_momentum[mdim] = input->mom_min[mdim] + aux_is[mdim+input->pos_dims] * input->mom_delta[mdim];
            }           
            input->dispersions[i](input->pos_dims, input->mom_dims, aux_momentum); 

            //Loop through different positions
            for(int j = 1; j < input->pos_total-2; ++j){
                
                for(int pdim = 0; pdim < input->pos_dims; ++pdim){

                    if(pdim == input->mom_dims){
                        flows[i][pdim][j*input->mom_total + k] = 0;
                        break;
                    }

                    if(aux_momentum[pdim]>=0){
                        
                        cent = species[i][j*input->mom_total+k];
                        forw = species[i][j*input->mom_total+k+input->grid_factor[pdim]];
                        back = species[i][j*input->mom_total+k-input->grid_factor[pdim]];

                        delta = (forw - cent) * MC_fluxLimiter(back, cent, forw);

                        flows[i][pdim][j*input->mom_total + k] = aux_momentum[pdim] * (cent + (1 - aux_momentum[pdim] * input->deltaT/input->pos_delta[pdim]) * delta/2);
                    }

                    else{

                        back = species[i][j*input->mom_total+k];
                        cent = species[i][j*input->mom_total+k+input->grid_factor[pdim]];
                        forw = species[i][j*input->mom_total+k+2*input->grid_factor[pdim]];

                        delta = (forw - cent) * MC_fluxLimiter(back, cent, forw);

                        flows[i][pdim][j*input->mom_total + k] = aux_momentum[pdim] * (cent - (1 - aux_momentum[pdim] * input->deltaT/input->pos_delta[pdim]) * delta/2);
                    }
                }   

                for(int mdim = 0; mdim < input->mom_dims; ++mdim){

                    force = 0;

                    for(int fld = 0; fld < input->n_fields; ++fld){
                        force += input->force_charges[i][fld] * input->forces[fld](input, fields, fld, j, k, mdim);
                    }

                    if(force>=0){
                        
                        cent = species[i][j*input->mom_total+k];
                        forw = species[i][j*input->mom_total+k+input->grid_factor[input->pos_dims+mdim]];
                        back = species[i][j*input->mom_total+k-input->grid_factor[input->pos_dims+mdim]];

                        delta = (forw - cent) * MC_fluxLimiter(back, cent, forw);

                        flows[i][input->pos_dims + mdim][j*input->mom_total + k] = force * (cent + (1 - force * input->deltaT/input->mom_delta[mdim]) * delta/2);
                    }

                    else{

                        back = species[i][j*input->mom_total+k];
                        cent = species[i][j*input->mom_total+k+input->grid_factor[input->pos_dims+mdim]];
                        forw = species[i][j*input->mom_total+k+2*input->grid_factor[input->pos_dims+mdim]];

                        delta = (forw - cent) * MC_fluxLimiter(back, cent, forw);

                        flows[i][input->pos_dims + mdim][j*input->mom_total + k] = force * (cent - (1 - force * input->deltaT/input->mom_delta[mdim]) * delta/2);
                    }
                }  
            }
        }
    }
}

void rungeKutta2(input_t* input, REAL** species, REAL** species_aux, COMPLEX** sources, COMPLEX** fields, COMPLEX** fields_aux, COMPLEX** fields_fft, COMPLEX** fields_deriv, REAL*** flows, int* aux_is, REAL* aux_momentum, REAL*** left_buffer, REAL*** right_buffer, COMPLEX*** left_field_buffer, COMPLEX*** right_field_buffer){

    // First Step //
    integrate_source(input, species, sources, aux_is, aux_momentum);
    
    for(int fld = 0; fld < input->n_fields; ++fld){
        if(!strcmp(input->field_type[fld],"static")){
            convolve_source(input, sources, fields_fft, fields, aux_is, aux_momentum, fld);
            invert_field(input, fields, fields_fft, aux_is, fld);
            field_bound_cond(input, fields, left_field_buffer, right_field_buffer, fld);
        }
    }
    
    finite_volume2(input, species, fields, aux_is, aux_momentum, flows);

    for(int fld = 0; fld < input->n_fields; ++fld){
        if(!strcmp(input->field_type[fld],"dynamic")){
            convolve_source(input, sources, fields_fft, fields_deriv, aux_is, aux_momentum, fld);
            transform_field(input, fields, fields_fft, aux_is, fld);
        }
    }
    for(int fld = 0; fld < input->n_fields; ++fld){
        if(!strcmp(input->field_type[fld],"dynamic")){
            field_dynamics(input, fields, fields_deriv, aux_is, aux_momentum, fld);
        }
    }
    for(int fld = 0; fld < input->n_fields; ++fld){
        if(!strcmp(input->field_type[fld],"dynamic")){
            for(int j = 0; j < input->pos_total; ++j){
                fields_aux[fld][j] = fields[fld][j] + input->deltaT * fields_deriv[fld][j];
            }
            invert_field(input, fields, fields_fft, aux_is, fld);
            field_bound_cond(input, fields, left_field_buffer, right_field_buffer, fld);
            invert_field(input, fields_aux, fields_fft, aux_is, fld);
            field_bound_cond(input, fields_aux, left_field_buffer, right_field_buffer, fld);
        }
    }

    for(int i = 0; i < input->n_species; ++i){
        for(int j = 2*input->space_factor[0]; j < input->pos_total-2*input->space_factor[0]; ++j){
            for(int k = 1; k < input->mom_total; ++k){
                species_aux[i][j * input->mom_total + k] = species[i][j * input->mom_total + k];
                for(int dim = 0; dim < input->pos_dims+input->mom_dims; ++dim){
                    species_aux[i][j * input->mom_total + k] -= input->lambda[dim] * (flows[i][dim][j * input->mom_total + k] - flows[i][dim][j * input->mom_total + k - input->grid_factor[dim]]);
                }
            }
        }
    }

    apply_bound_cond(input, species_aux, left_buffer, right_buffer);    

    // Second Step //
    integrate_source(input, species_aux, sources, aux_is, aux_momentum);
    
    for(int fld = 0; fld < input->n_fields; ++fld){
        if(!strcmp(input->field_type[fld],"static")){
            convolve_source(input, sources, fields_fft, fields_aux, aux_is, aux_momentum, fld);
            invert_field(input, fields_aux, fields_fft, aux_is, fld);
            field_bound_cond(input, fields_aux, left_field_buffer, right_field_buffer, fld);
        }
    }

    finite_volume2(input, species_aux, fields_aux, aux_is, aux_momentum, flows);

    for(int fld = 0; fld < input->n_fields; ++fld){
        if(!strcmp(input->field_type[fld],"dynamic")){
            convolve_source(input, sources, fields_fft, fields_deriv, aux_is, aux_momentum, fld);
            transform_field(input, fields_aux, fields_fft, aux_is, fld);
        }
    }
    for(int fld = 0; fld < input->n_fields; ++fld){
        if(!strcmp(input->field_type[fld],"dynamic")){
            field_dynamics(input, fields_aux, fields_deriv, aux_is, aux_momentum, fld);
        }
    }
    for(int fld = 0; fld < input->n_fields; ++fld){
        if(!strcmp(input->field_type[fld],"dynamic")){
            for(int j = 0; j < input->pos_total; ++j){
                fields_aux[fld][j] += input->deltaT * fields_deriv[fld][j];
            }
            invert_field(input, fields_aux, fields_fft, aux_is, fld);
            field_bound_cond(input, fields_aux, left_field_buffer, right_field_buffer, fld);
        }
    }

    for(int i = 0; i < input->n_species; ++i){
        for(int j = 2*input->space_factor[0]; j < input->pos_total-2*input->space_factor[0]; ++j){
            for(int k = 1; k < input->mom_total; ++k){
                for(int dim = 0; dim < input->pos_dims+input->mom_dims; ++dim){
                    species_aux[i][j * input->mom_total + k] -= input->lambda[dim] * (flows[i][dim][j * input->mom_total + k] - flows[i][dim][j * input->mom_total + k - input->grid_factor[dim]]);
                }
            }
        }
    }

    apply_bound_cond(input, species_aux, left_buffer, right_buffer);

    // Average Steps //
    for(int i = 0; i < input->n_species; ++i){
        for(int j = 2*input->space_factor[0]; j < input->pos_total-2*input->space_factor[0]; ++j){
            for(int k = 1; k < input->mom_total; ++k){
                species[i][j * input->mom_total + k] = 0.5*(species[i][j * input->mom_total + k] + species_aux[i][j * input->mom_total + k]);
            }
        }
    }
    for(int fld = 0; fld < input->n_fields; ++fld){
        if(!strcmp(input->field_type[fld],"dynamic")){
            for(int j = 0; j < input->pos_total; ++j){
                fields[fld][j] = 0.5*(fields[fld][j] + fields_aux[fld][j]);
            }
            field_bound_cond(input, fields, left_field_buffer, right_field_buffer, fld);
        }
    }
    
    apply_bound_cond(input, species, left_buffer, right_buffer);
    
}