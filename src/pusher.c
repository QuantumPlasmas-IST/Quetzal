#include "pusher.h"

REAL** solve(input_t* input){

    // Allocating AAALLLLLL the memory //

    REAL** species = malloc(input->n_species * sizeof(REAL*));
    REAL** species_aux1 = malloc(input->n_species * sizeof(REAL*));
    REAL** species_aux2 = malloc(input->n_species * sizeof(REAL*));

    COMPLEX** fields = malloc(input->n_fields * sizeof(COMPLEX*));

    COMPLEX** sources = malloc(input->n_fields * sizeof(COMPLEX*));
    COMPLEX** sources_aux = malloc(input->n_species * sizeof(COMPLEX*));

    REAL*** flows = malloc(input->n_species * sizeof(REAL**));

    species[0] = malloc(input->n_species * input->pos_total * input->mom_total * sizeof(REAL));
    species_aux1[0] = malloc(input->n_species * input->pos_total * input->mom_total * sizeof(REAL));
    species_aux2[0] = malloc(input->n_species * input->pos_total * input->mom_total * sizeof(REAL));
    sources_aux[0] = malloc(input->n_species * input->pos_total * sizeof(COMPLEX));
    sources[0] = malloc(input->n_fields * input->pos_total * sizeof(COMPLEX));
    fields[0] = malloc(input->n_fields * input->pos_total * sizeof(COMPLEX));

    for (int i = 0; i < input->n_species; ++i){
        species[i] = species[0] + i * input->pos_total * input->mom_total;
        species_aux1[i] = species_aux1[0] + i * input->pos_total * input->mom_total;
        species_aux2[i] = species_aux2[0] + i * input->pos_total * input->mom_total;
        sources_aux[i] = sources_aux[0] + i * input->pos_total;

        flows[i] = malloc((input->pos_dims + input->mom_dims) * sizeof(REAL*));
        for(int dim = 0; dim < input->pos_dims+input->mom_dims; ++dim){
            flows[i][dim] = malloc((input->pos_total * input->mom_total) * sizeof(REAL));
        }
    }

    for (int i = 0; i < input->n_fields; ++i){
        fields[i] = fields[0] + i * input->pos_total;
        sources[i] = sources[0] + i * input->pos_total;
    }

    int* is = malloc(sizeof(int) * (input->pos_dims+input->mom_dims));
    REAL* aux_momentum = malloc(input->mom_dims * sizeof(REAL));

    printf("Memory Allocation Complete\n");
    // Memory successfully allocated *thumbs up emoji* //

    // Applying initial conditions //

    apply_init_cond(input, species);
    apply_init_cond(input, species_aux1);
    apply_init_cond(input, species_aux2);
    field_init_cond(input, fields);
    apply_bound_cond(input, species);
    apply_bound_cond(input, species_aux1);
    apply_bound_cond(input, species_aux2);
    integrate_source(input, species, sources, sources_aux);

    // Algorithm //

    printf("Initialization Complete\n");

    if(!strcmp(input->pusher,"leapfrog")){
        for(int i = 0; i < input->n_timesteps; ++i){
            convolute_field(input, sources, fields);
            leapfrog(input, species, species_aux1, species_aux2, fields, is);
            REAL** swap = species;
            species = species_aux2;
            species_aux2 = species_aux1;
            species_aux1 = swap;
            apply_bound_cond(input, species);
            integrate_source(input, species, sources, sources_aux);
            if(!(i%input->diag_freq)){
                write_fields(input, fields, i);
                write_solution(input, species, i);
                write_sources(input, sources, i);
                printf("Timestep %d\r",i);
                fflush(stdout);
            }
        }
    }

    if(!strcmp(input->pusher,"runge kutta")){
        for(int i = 0; i < input->n_timesteps; ++i){
            if(!(i%input->diag_freq)){
                write_fields(input, fields, i);
                write_solution(input, species, i);
                write_sources(input, sources_aux, i);
                printf("Timestep %d\r",i);
                fflush(stdout);
            }
            rungeKutta2(input, species, species_aux1, species_aux2, sources, sources_aux, fields, flows, is, aux_momentum);
            apply_bound_cond(input, species);
        }
    }

    // Free all the memory

    free(is);
    free(aux_momentum);

    for (int i = 0; i < input->n_species; ++i){
        for(int dim = 0; dim < input->pos_dims+input->mom_dims; ++dim){
            free(flows[i][dim]);
        }
        free(flows[i]);
    }
    free(flows);

    free(species_aux1[0]);
    free(species_aux2[0]);
    free(sources[0]);
    free(sources_aux[0]);
    free(fields[0]);
    free(species_aux1);
    free(species_aux2);
    free(sources);
    free(sources_aux);
    free(fields);

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

void integrate_source(input_t* input, REAL** species, COMPLEX** sources, COMPLEX** sources_aux){
    
    for(int i = 0; i < input->n_species; ++i){
        for(int j = 0; j < input->pos_total; ++j){
            sources_aux[i][j] = 0;
            for(int k = 0; k < input->mom_total; ++k){
                sources_aux[i][j] += species[i][j*input->mom_total + k];
            }
            for(int d = 0; d < input->mom_dims; ++d){
                sources_aux[i][j] *= input->mom_delta[d];
            }
        }
    }

    for(int k = 0; k < input->n_fields; ++k){
        if(strcmp(input->kernel_names[k],"fixed")){
            for(int j = 0; j < input->pos_total; ++j){
                sources[k][j] = 0;
                for(int i = 0; i < input->n_species; ++i){
                    sources[k][j] += sources_aux[i][j] * input->charges[i][k];
                }
            }
        }
    }

}

void convolute_field(input_t* input, COMPLEX** sources, COMPLEX** fields){

    REAL* wavevector = malloc(input->pos_dims * sizeof(REAL));
    int* indices = malloc(input->pos_dims * sizeof(int));

    for(int k = 0; k < input->n_fields; ++k){

        if(!strcmp(input->kernel_names[k],"density")){
            COMPLEX* aux_ccPtr = fields[k];
            fields[k] = sources[k];
            sources[k] = aux_ccPtr;
        }

        else if(strcmp(input->kernel_names[k],"fixed")){

            boundary_shift(input, sources[k], fields[k], indices);

            NdFourier(input->pos_dims, input->pos_points, fields[k], sources[k]);

            for(int i = 0; i < input->pos_total; ++i){

                axis_index(input->pos_dims, input->pos_points, indices, i);

                for(int j = 0; j < input->pos_dims; ++j){
                    if (indices[j]<input->pos_points[j]/2){
                        wavevector[j] = 2 * M_PI * indices[j] / (input->pos_max[j]-input->pos_min[j]);
                    }
                    else{
                        wavevector[j] = 2 * M_PI * (-input->pos_points[j] + indices[j]) / (input->pos_max[j]-input->pos_min[j]);
                    }
                }
                
                fields[k][i] = input->kernels[k](input->pos_dims, sources[k][i], wavevector);
            } 
            
            NdInverseFourier(input->pos_dims, input->pos_points, fields[k], sources[k]);

            inverse_boundary_shift(input, sources[k], fields[k], indices);

        }

    }

    free(wavevector);
    free(indices);
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
                aux_momentum[mdim] = input->mom_min[mdim] + aux_is[mdim+input->pos_dims] * input->mom_delta[mdim];
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
                        force += input->charges[i][fld] * input->forces[fld](input, fields, fld, j, k, mdim);
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
                        force += input->charges[i][fld] * input->forces[fld](input, fields, fld, j, k, mdim);
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

void rungeKutta2(input_t* input, REAL** species, REAL** aux1, REAL** aux2, COMPLEX** sources, COMPLEX** sources_aux, COMPLEX** fields, REAL*** flows, int* aux_is, REAL* aux_momentum){
    
    integrate_source(input, species, sources, sources_aux);
    convolute_field(input, sources, fields);
    finite_volume2(input, species, fields, aux_is, aux_momentum, flows);

    for(int i = 0; i < input->n_species; ++i){
        for(int j = 2*input->space_factor[0]; j < input->pos_total-2*input->space_factor[0]; ++j){
            for(int k = 1; k < input->mom_total; ++k){
                aux1[i][j * input->mom_total + k] = species[i][j * input->mom_total + k];
                for(int dim = 0; dim < input->pos_dims+input->mom_dims; ++dim){
                    aux1[i][j * input->mom_total + k] -= input->lambda[dim] * (flows[i][dim][j * input->mom_total + k] - flows[i][dim][j * input->mom_total + k - input->grid_factor[dim]]);
                }
            }
        }
    }

    apply_bound_cond(input, aux1);
    integrate_source(input, aux1, sources, sources_aux);
    convolute_field(input, sources, fields);
    finite_volume2(input, aux1, fields, aux_is, aux_momentum, flows);

    for(int i = 0; i < input->n_species; ++i){
        for(int j = 2*input->space_factor[0]; j < input->pos_total-2*input->space_factor[0]; ++j){
            for(int k = 1; k < input->mom_total; ++k){
                aux2[i][j * input->mom_total + k] = aux1[i][j * input->mom_total + k];
                for(int dim = 0; dim < input->pos_dims+input->mom_dims; ++dim){
                    aux2[i][j * input->mom_total + k] -= input->lambda[dim] * (flows[i][dim][j * input->mom_total + k] - flows[i][dim][j * input->mom_total + k - input->grid_factor[dim]]);
                }
            }
        }
    }

    for(int i = 0; i < input->n_species; ++i){
        for(int j = 2*input->space_factor[0]; j < input->pos_total-2*input->space_factor[0]; ++j){
            for(int k = 1; k < input->mom_total; ++k){
                species[i][j * input->mom_total + k] = 0.5*(species[i][j * input->mom_total + k] + aux2[i][j * input->mom_total + k]);
            }
        }
    }
}