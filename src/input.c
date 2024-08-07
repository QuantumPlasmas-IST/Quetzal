#include "input.h"
#define NPARAMS 20
#define STR_SIZE 100

void initialize_input(input_t* input, REAL* pos_parameters, REAL* mom_parameters, REAL* field_parameters){

    int processes = 1;
    for(int i = 0; i < input->pos_dims; ++i){
        processes *= input->procs[i];
    }
    if(!(input->rank || processes==input->size)){
        printf("Wrong Number of Processes Given!\n");
        return;
    }
    for(int i = 0; i < input->pos_dims; ++i){
        processes /= input->procs[i];
        input->parallel_factor[i] = processes;
    }

    axis_index(input->pos_dims, input->procs, input->parallel_pos, input->rank);

    input->pos_delta = malloc(input->pos_dims * sizeof(REAL));
    input->mom_delta = malloc(input->mom_dims * sizeof(REAL));
    input->dk = malloc(input->pos_dims * sizeof(REAL));
    for(int i = 0; i < input->pos_dims; ++i){
        input->pos_delta[i] = (input->pos_max[i]-input->pos_min[i])/input->pos_points[i];
        input->dk[i] = 2 * M_PI / (input->pos_max[i]-input->pos_min[i]);
        input->pos_points[i] = input->pos_points[i]/input->procs[i];
        input->pos_min[i] += input->parallel_pos[i] * input->pos_points[i] * input->pos_delta[i];
        input->pos_max[i] = input->pos_min[i] + input->pos_points[i]*input->pos_delta[i];
        input->pos_points[i] += 2 * input->padding;
    }
    for(int i = 0; i < input->mom_dims; ++i){
        input->mom_delta[i] = (input->mom_max[i]-input->mom_min[i])/input->mom_points[i];
        input->mom_points[i] += 2 * input->padding;
    }

    int* aux_points = malloc((input->pos_dims + input->mom_dims) * sizeof(int));
    for(int i = 0; i < input->pos_dims; ++i){
        aux_points[i] = input->pos_points[i];
    }
    for(int i = 0; i < input->mom_dims; ++i){
        aux_points[input->pos_dims+i] = input->mom_points[i];
    }

    REAL* aux_delta = malloc((input->pos_dims + input->mom_dims) * sizeof(REAL));
    input->lambda = malloc((input->pos_dims + input->mom_dims) * sizeof(REAL));
    for(int i = 0; i < input->pos_dims; ++i){
        aux_delta[i] = input->pos_delta[i];
        input->lambda[i] = input->deltaT/aux_delta[i];
    }
    for(int i = 0; i < input->mom_dims; ++i){
        aux_delta[input->pos_dims+i] = input->mom_delta[i];
        input->lambda[input->pos_dims+i] = input->deltaT/aux_delta[input->pos_dims+i];
    }

    free(input->pos_points);
    free(input->mom_points);
    free(input->pos_delta);
    free(input->mom_delta);
    input->pos_points = aux_points;
    input->mom_points = aux_points + input->pos_dims;
    input->pos_delta = aux_delta;
    input->mom_delta = aux_delta + input->pos_dims;


    for (int i = 0; i < input->pos_dims; ++i){
        input->pos_total *= input->pos_points[i];
    }
    for (int i = 0; i < input->mom_dims; ++i){
        input->mom_total *= input->mom_points[i];
    }

    REAL factor = input->pos_total*input->mom_total;
    input->grid_factor = malloc((input->pos_dims + input->mom_dims) * sizeof(int));
    input->space_factor = malloc(input->pos_dims * sizeof(int));
    for (int i = 0; i < input->pos_dims+input->mom_dims; ++i){
        factor /= input->pos_points[i];
        input->grid_factor[i] = factor;
    }
    for (int i = 0; i < input->pos_dims; ++i){
        input->space_factor[i] = input->grid_factor[i]/input->mom_total;
    }

    for (int i = 0; i < input->mom_dims; ++i){
        if(!strcmp(input->mom_bound_names[i],"dirichelet")){
            input->mom_bounds[i]=&dirichelet_bound;
        }
        if(!strcmp(input->mom_bound_names[i],"periodic")){
            input->mom_bounds[i]=&periodic_bound;
        }
    }
    for (int i = 0; i < input->pos_dims; ++i){
        if(!strcmp(input->pos_bound_names[i],"dirichelet")){
            input->pos_bounds[i]=&dirichelet_bound;
        }
        if(!strcmp(input->pos_bound_names[i],"periodic")){
            input->pos_bounds[i]=&periodic_bound;
        }
    }

    int field_param_count=0;
    for(int i = 0; i < input->n_fields; ++i){
        
        //Kernels
        if(!strcmp(input->kernel_names[i],"2Dpoisson")){
            input->kernels[i]=&poisson_2D;
        }
        if(!strcmp(input->kernel_names[i],"3Dpoisson")){
            input->kernels[i]=&poisson_3D;
        }

        //Forces
        if(!strcmp(input->force_names[i],"gradient")){
            input->forces[i]=&gradient_force;
        }
        if(!strcmp(input->force_names[i],"electricX")){
            input->forces[i]=&electric0_force;
        }
        if(!strcmp(input->force_names[i],"electricY")){
            input->forces[i]=&electric1_force;
        }
        if(!strcmp(input->force_names[i],"magnetic")){
            input->forces[i]=&magnetic_force;
        }

        //Field initial conditions
        if(!strcmp(input->field_init_names[i],"gaussian")){
            input->field_inits[i]=&gaussian;
            for (int j = 0; j < 3; ++j){
                input->field_init_params[i][j]=field_parameters[field_param_count];
                field_param_count++;
            }
        }
        if(!strcmp(input->field_init_names[i],"uniform")){
            input->field_inits[i]=&uniform;
            for (int j = 0; j < 3; ++j){
                input->field_init_params[i][j]=field_parameters[field_param_count];
                field_param_count++;
            }
        }
    }

    int pos_param_count=0;
    int mom_param_count=0;
    for(int i = 0; i < input->n_species; ++i){
        // Dispersion Relation
        if(!strcmp(input->dispersion_names[i],"linear")){
            input->dispersions[i]=&linear_dispersion;
        }
        if(!strcmp(input->dispersion_names[i],"quadratic")){
            input->dispersions[i]=&quadratic_dispersion;
        }

        // Position Initial Conditions
        if(!strcmp(input->pos_init_names[i],"gaussian")){
            input->pos_inits[i]=&gaussian;
            for (int j = 0; j < 2+input->pos_dims; ++j){
                input->pos_init_params[i][j]=pos_parameters[pos_param_count];
                pos_param_count++;
            }
        }
        if(!strcmp(input->pos_init_names[i],"double gaussian")){
            input->pos_inits[i]=&double_gaussian;
            for (int j = 0; j < 4+2*input->pos_dims; ++j){
                input->pos_init_params[i][j]=pos_parameters[pos_param_count];
                pos_param_count++;
            }
        }
        if(!strcmp(input->pos_init_names[i],"uniform")){
            input->pos_inits[i]=&uniform;
            for (int j = 0; j < 1+2*input->pos_dims; ++j){
                input->pos_init_params[i][j]=pos_parameters[pos_param_count];
                pos_param_count++;
            }
        }
        if(!strcmp(input->pos_init_names[i],"noise")){
            input->pos_inits[i]=&noise;
            for (int j = 0; j < 2; ++j){
                input->pos_init_params[i][j]=pos_parameters[pos_param_count];
                pos_param_count++;
            }
        }

        // Momentum Initial Condition
        if(!strcmp(input->mom_init_names[i],"gaussian")){
            input->mom_inits[i]=&gaussian;
            for (int j = 0; j < 2+input->mom_dims; ++j){
                input->mom_init_params[i][j]=mom_parameters[mom_param_count];
                mom_param_count++;
            }
        }
        if(!strcmp(input->mom_init_names[i],"double gaussian")){
            input->mom_inits[i]=&double_gaussian;
            for (int j = 0; j < 4+2*input->mom_dims; ++j){
                input->mom_init_params[i][j]=mom_parameters[mom_param_count];
                mom_param_count++;
            }
        }
        if(!strcmp(input->mom_init_names[i],"uniform")){
            input->mom_inits[i]=&uniform;
            for (int j = 0; j < 1+2*input->mom_dims; ++j){
                input->mom_init_params[i][j]=mom_parameters[mom_param_count];
                mom_param_count++;
            }
        }
        if(!strcmp(input->mom_init_names[i],"noise")){
            input->mom_inits[i]=&noise;
            for (int j = 0; j < 0; ++j){
                input->mom_init_params[i][j]=mom_parameters[mom_param_count];
                mom_param_count++;
            }
        }
        if(!strcmp(input->mom_init_names[i],"linear fermi")){
            input->mom_inits[i]=&linear_fermi;
            for (int j = 0; j < 2+input->mom_dims; ++j){
                input->mom_init_params[i][j]=mom_parameters[mom_param_count];
                mom_param_count++;
            }
        }
        if(!strcmp(input->mom_init_names[i],"quadratic fermi")){
            input->mom_inits[i]=&quad_fermi;
            for (int j = 0; j < 2+input->mom_dims; ++j){
                input->mom_init_params[i][j]=mom_parameters[mom_param_count];
                mom_param_count++;
            }
        }
        if(!strcmp(input->mom_init_names[i],"flattened quadratic fermi")){
            input->mom_inits[i]=&flattened_quad_fermi;
            for (int j = 0; j < 2+input->mom_dims; ++j){
                input->mom_init_params[i][j]=mom_parameters[mom_param_count];
                mom_param_count++;
            }
        }
        if(!strcmp(input->mom_init_names[i],"linear bose")){
            input->mom_inits[i]=&linear_bose;
            for (int j = 0; j < 2+input->mom_dims; ++j){
                input->mom_init_params[i][j]=mom_parameters[mom_param_count];
                mom_param_count++;
            }
        }
        if(!strcmp(input->mom_init_names[i],"quadratic bose")){
            input->mom_inits[i]=&quad_bose;
            for (int j = 0; j < 2+input->mom_dims; ++j){
                input->mom_init_params[i][j]=mom_parameters[mom_param_count];
                mom_param_count++;
            }
        }
        if(!strcmp(input->mom_init_names[i],"linear maxwell")){
            input->mom_inits[i]=&linear_maxwell;
            for (int j = 0; j < 2+input->mom_dims; ++j){
                input->mom_init_params[i][j]=mom_parameters[mom_param_count];
                mom_param_count++;
            }
        }
        if(!strcmp(input->mom_init_names[i],"quadratic maxwell")){
            input->mom_inits[i]=&quad_maxwell;
            for (int j = 0; j < 2+input->mom_dims; ++j){
                input->mom_init_params[i][j]=mom_parameters[mom_param_count];
                mom_param_count++;
            }
        }
        if(!strcmp(input->mom_init_names[i],"circle")){
            input->mom_inits[i]=&circle;
            for (int j = 0; j < 2+input->mom_dims; ++j){
                input->mom_init_params[i][j]=mom_parameters[mom_param_count];
                mom_param_count++;
            }
        }
        if(!strcmp(input->mom_init_names[i],"sphere")){
            input->mom_inits[i]=&sphere;
            for (int j = 0; j < 2+input->mom_dims; ++j){
                input->mom_init_params[i][j]=mom_parameters[mom_param_count];
                mom_param_count++;
            }
        }
    }
}

input_t* read_input(const char* name, int rank, int size){
    
    input_t* input = malloc(sizeof(input_t));

    input->rank = rank;
    input->size = size;
    input->procs = NULL;
    input->parallel_pos = NULL;
    input->parallel_factor = NULL;

    input->filename = malloc(200);
    input->pos_dims = 0;
    input->mom_dims = 0;
    input->pos_points = NULL;
    input->mom_points = NULL;
    input->pos_total = 1;
    input->mom_total = 1;
    input->pos_min = NULL;
    input->pos_max = NULL;
    input->mom_min = NULL;
    input->mom_max = NULL;
    input->n_species = 0;
    input->n_fields = 0;
    input->deltaT = 0;
    input->n_timesteps = 0;
    input->dispersion_names = NULL;
    input->dispersions = NULL;
    input->pusher = calloc(STR_SIZE, sizeof(char));
    input->operator = calloc(STR_SIZE, sizeof(char));
    input->kernel_names = NULL;
    input->kernels = NULL;
    input->force_names = NULL;
    input->forces = NULL;
    input->pos_init_names = NULL;
    input->mom_init_names = NULL;
    input->field_init_names = NULL;
    input->pos_bound_names = NULL;
    input->mom_bound_names = NULL;
    input->pos_inits = NULL;
    input->mom_inits = NULL;
    input->field_inits = NULL;
    input->pos_bounds = NULL;
    input->mom_bounds = NULL;
    input->pos_init_params = NULL;
    input->mom_init_params = NULL;
    input->field_init_params = NULL;
    REAL* pos_parameters = NULL;
    REAL* mom_parameters = NULL;
    REAL* field_parameters = NULL;
    input->padding = 1;
    input->charges = NULL;

    char line[500];
	FILE* file = fopen(name, "r");

	while(fgets(line, 500, file)!=NULL){

        char num_part[500];
        char txt_part[500];

        strcpy(num_part, "");
        strcpy(txt_part, "");

        int count = 0;
        int str_test=0;

		for(int i = 0; i<500; ++i) {

            if ((int)(line[i])==59) break; // Character (;) ends line

            if ((int)(line[i])==34){ // Character (") toggles string mode
                if(str_test) str_test=0;
                else str_test=1; 
            }
                

			if((isdigit(line[i]) || (int)(line[i]) == 46 || str_test || (int)(line[i]) == 45)&&line[i]!='"'){
				strncat(num_part,&line[i],1);
            }
            else{
			    if (((int)(line[i]) >= 65 && (int)(line[i]) <= 90) || ((int)(line[i]) >= 97 && (int)(line[i]) <= 122) || (int)(line[i]) == 95){
				    strncat(txt_part,&line[i],1);
                }
            }

            if ((int)(line[i])==44){
                if(!strcmp(txt_part,"POS_POINTS")){
                    (input->pos_points)[count] = atof(num_part);   // Character (,) pushes back new value of min to array
                    strcpy(num_part, "");
                    count++;
                }
                if(!strcmp(txt_part,"MOM_POINTS")){
                    (input->mom_points)[count] = atof(num_part);   // Character (,) pushes back new value of min to array
                    strcpy(num_part, "");
                    count++;
                }
                if(!strcmp(txt_part,"POS_MIN")){
                    (input->pos_min)[count] = atof(num_part);   // Character (,) pushes back new value of min to array
                    strcpy(num_part, "");
                    count++;
                }
                if(!strcmp(txt_part,"MOM_MIN")){
                    (input->mom_min)[count] = atof(num_part);   // Character (,) pushes back new value of min to array
                    strcpy(num_part, "");
                    count++;
                }
                if(!strcmp(txt_part,"POS_MAX")){
                    (input->pos_max)[count] = atof(num_part);     // Character (,) pushes back new value of max to array
                    strcpy(num_part, "");
                    count++;
                }   
                if(!strcmp(txt_part,"PROCS")){
                    (input->procs)[count] = atoi(num_part);     // Character (,) pushes back new value of max to array
                    strcpy(num_part, "");
                    count++;
                }   
                if(!strcmp(txt_part,"MOM_MAX")){
                    (input->mom_max)[count] = atof(num_part);   // Character (,) pushes back new value of min to array
                    strcpy(num_part, "");
                    count++;
                }
                if(!strcmp(txt_part,"POS_INIT")){
                    strcpy((input->pos_init_names)[count], num_part);     // Character (,) pushes back new value of max to array
                    strcpy(num_part, "");
                    count++;
                }   
                if(!strcmp(txt_part,"MOM_INIT")){
                    strcpy((input->mom_init_names)[count], num_part);     // Character (,) pushes back new value of max to array
                    strcpy(num_part, "");
                    count++;
                } 
                if(!strcmp(txt_part,"FIELD_INIT")){
                    strcpy((input->field_init_names)[count], num_part);     // Character (,) pushes back new value of max to array
                    strcpy(num_part, "");
                    count++;
                }   
                if(!strcmp(txt_part,"POS_PARAMS")){
                    pos_parameters[count] = atof(num_part);     // Character (,) pushes back new value of max to array
                    strcpy(num_part, "");
                    count++;
                }    
                if(!strcmp(txt_part,"MOM_PARAMS")){
                    mom_parameters[count] = atof(num_part);     // Character (,) pushes back new value of max to array
                    strcpy(num_part, "");
                    count++;
                }   
                if(!strcmp(txt_part,"FIELD_PARAMS")){
                    field_parameters[count] = atof(num_part);     // Character (,) pushes back new value of max to array
                    strcpy(num_part, "");
                    count++;
                }    
                if(!strcmp(txt_part,"POS_BOUND")){
                    strcpy((input->pos_bound_names)[count], num_part);     // Character (,) pushes back new value of max to array
                    strcpy(num_part, "");
                    count++;
                }   
                if(!strcmp(txt_part,"MOM_BOUND")){
                    strcpy((input->mom_bound_names)[count], num_part);     // Character (,) pushes back new value of max to array
                    strcpy(num_part, "");
                    count++;
                }   
                if(!strcmp(txt_part,"KERNEL")){
                    strcpy((input->kernel_names)[count], num_part);     // Character (,) pushes back new value of max to array
                    strcpy(num_part, "");
                    count++;
                }   
                if(!strcmp(txt_part,"FORCE")){
                    strcpy((input->force_names)[count], num_part);     // Character (,) pushes back new value of max to array
                    strcpy(num_part, "");
                    count++;
                }   
                if(!strcmp(txt_part,"DISP")){
                    strcpy((input->dispersion_names)[count], num_part);     // Character (,) pushes back new value of max to array
                    strcpy(num_part, "");
                    count++;
                }
                if(!strcmp(txt_part,"CHARGES")){
                    (input->charges[0])[count] = atof(num_part);   // Character (,) pushes back new value of min to array
                    strcpy(num_part, "");
                    count++;
                }    
            }
            
		}
		if(!strcmp(txt_part,"POS_DIMS")){

            input->pos_dims = atoi(num_part);

            input->pos_min=malloc((input->pos_dims)*sizeof(REAL));
            input->pos_max=malloc((input->pos_dims)*sizeof(REAL));

            input->pos_points=malloc((input->pos_dims)*sizeof(int));

            input->pos_bound_names=malloc((input->pos_dims)*sizeof(char*));
            char* aux_bound_names = calloc((input->pos_dims) * STR_SIZE, sizeof(char));
            
            for (int i = 0; i < input->pos_dims; ++i){
                input->pos_bound_names[i] = aux_bound_names + i*STR_SIZE;
            }

            input->pos_bounds = malloc((input->pos_dims)*sizeof(bound_t));

            input->procs = malloc(input->pos_dims * sizeof(int));
            input->parallel_pos = malloc(input->pos_dims * sizeof(int));
            input->parallel_factor = malloc(input->pos_dims * sizeof(int));
        }

        if(!strcmp(txt_part,"MOM_DIMS")){
            
            input->mom_dims = atoi(num_part);

            input->mom_min=malloc((input->mom_dims)*sizeof(REAL));
            input->mom_max=malloc((input->mom_dims)*sizeof(REAL));

            input->mom_points=malloc((input->mom_dims)*sizeof(int));

            input->mom_bound_names=malloc((input->mom_dims)*sizeof(char*));
            char* aux_bound_names = calloc((input->mom_dims) * STR_SIZE, sizeof(char));

            for (int i = 0; i < input->mom_dims; ++i){
                input->mom_bound_names[i] = aux_bound_names + i*STR_SIZE;
            }

            input->mom_bounds = malloc((input->mom_dims)*sizeof(bound_t));
        }

        if(!strcmp(txt_part,"N_SPECIES")){
            
            input->n_species = atoi(num_part);

            input->dispersion_names = malloc(input->n_species * sizeof(char*));
            for (int i = 0; i < input->n_species; ++i){
                input->dispersion_names[i] = calloc(STR_SIZE, sizeof(char));
            } 

            input->dispersions = malloc((input->n_species)*sizeof(disp_t));

            input->charges = malloc(input->n_species * sizeof(REAL*));

            input->pos_init_names=malloc((input->n_species)*sizeof(char*));
            input->mom_init_names=malloc((input->n_species)*sizeof(char*));
            char* aux_pos_init_names = calloc((input->n_species) * STR_SIZE, sizeof(char));
            char* aux_mom_init_names = calloc((input->n_species) * STR_SIZE, sizeof(char));


            pos_parameters = malloc(input->n_species * NPARAMS * sizeof(REAL));
            mom_parameters = malloc(input->n_species * NPARAMS * sizeof(REAL));
            input->pos_init_params = malloc((input->n_species)*sizeof(REAL*));
            input->mom_init_params = malloc((input->n_species)*sizeof(REAL*));
            REAL* aux_pos_init_params = calloc(input->n_species * NPARAMS, sizeof(REAL));
            REAL* aux_mom_init_params = calloc(input->n_species * NPARAMS, sizeof(REAL));


            for (int i = 0; i < input->n_species; ++i){
                input->pos_init_names[i] = aux_pos_init_names + i*STR_SIZE;
                input->mom_init_names[i] = aux_mom_init_names + i*STR_SIZE;
                input->pos_init_params[i] = aux_pos_init_params + i*NPARAMS;
                input->mom_init_params[i] = aux_mom_init_params + i*NPARAMS;
            }


            input->pos_inits = malloc((input->n_species)*sizeof(init_t));
            input->mom_inits = malloc((input->n_species)*sizeof(init_t));

        }

        if(!strcmp(txt_part,"N_FIELDS")){
            
            input->n_fields = atoi(num_part);

            input->kernel_names = malloc(input->n_fields * sizeof(char*));
            input->force_names = malloc(input->n_fields * sizeof(char*));

            for (int i = 0; i < input->n_fields; ++i){
                input->kernel_names[i] = calloc(STR_SIZE, sizeof(char));
                input->force_names[i] = calloc(STR_SIZE, sizeof(char));
            } 

            input->kernels = malloc((input->n_fields)*sizeof(kernel_t));
            input->forces = malloc((input->n_fields)*sizeof(force_t));

            input->charges[0] = malloc(input->n_species * input->n_fields * sizeof(REAL));
            for (int i = 1; i < input->n_species; ++i){
                input->charges[i] = input->charges[0] + i*input->n_fields;
            }

            input->field_init_names=malloc((input->n_fields)*sizeof(char*));
            char* aux_field_init_names = calloc((input->n_fields) * STR_SIZE, sizeof(char));

            field_parameters = malloc(input->n_fields * NPARAMS * sizeof(REAL));
            input->field_init_params = malloc((input->n_fields)*sizeof(REAL*));
            REAL* aux_field_init_params = calloc(input->n_fields * NPARAMS, sizeof(REAL));


            for (int i = 0; i < input->n_fields; ++i){
                input->field_init_names[i] = aux_field_init_names + i*STR_SIZE;
                input->field_init_params[i] = aux_field_init_params + i*NPARAMS;
            }


            input->field_inits = malloc((input->n_fields)*sizeof(init_t));

        }
        
        if(!strcmp(txt_part,"POS_POINTS")) input->pos_points[count] = atoi(num_part);
        if(!strcmp(txt_part,"MOM_POINTS")) input->mom_points[count] = atoi(num_part);
        if(!strcmp(txt_part,"N_TIMESTEPS")) input->n_timesteps = atoi(num_part);
        if(!strcmp(txt_part,"POS_MIN")) input->pos_min[count] = atof(num_part);
        if(!strcmp(txt_part,"POS_MAX")) input->pos_max[count] = atof(num_part);
        if(!strcmp(txt_part,"PROCS")) input->procs[count] = atoi(num_part);
        if(!strcmp(txt_part,"MOM_MIN")) input->mom_min[count] = atof(num_part);
        if(!strcmp(txt_part,"MOM_MAX")) input->mom_max[count] = atof(num_part);
        if(!strcmp(txt_part,"POS_INIT")) strcpy((input->pos_init_names)[count], num_part);
        if(!strcmp(txt_part,"MOM_INIT")) strcpy((input->mom_init_names)[count], num_part);
        if(!strcmp(txt_part,"FIELD_INIT")) strcpy((input->field_init_names)[count], num_part);
        if(!strcmp(txt_part,"POS_BOUND")) strcpy((input->pos_bound_names)[count], num_part);
        if(!strcmp(txt_part,"MOM_BOUND")) strcpy((input->mom_bound_names)[count], num_part);
        if(!strcmp(txt_part,"DISP")) strcpy(input->dispersion_names[count],num_part);
        if(!strcmp(txt_part,"KERNEL")) strcpy(input->kernel_names[count],num_part);
        if(!strcmp(txt_part,"FORCE")) strcpy(input->force_names[count],num_part);
        if(!strcmp(txt_part,"PUSHER")) strcpy(input->pusher,num_part);
        if(!strcmp(txt_part,"DELTAT")) input->deltaT = atof(num_part);
        if(!strcmp(txt_part,"POS_PARAMS")) pos_parameters[count] = atof(num_part);
        if(!strcmp(txt_part,"MOM_PARAMS")) mom_parameters[count] = atof(num_part);
        if(!strcmp(txt_part,"FIELD_PARAMS")) field_parameters[count] = atof(num_part);
        if(!strcmp(txt_part,"FILENAME")) strcpy(input->filename, num_part);
        if(!strcmp(txt_part,"OPER")) strcpy(input->operator, num_part);
        if(!strcmp(txt_part,"POS_DIAG_FREQ")) input->pos_diag_freq = atoi(num_part);
        if(!strcmp(txt_part,"MOM_DIAG_FREQ")) input->mom_diag_freq = atoi(num_part);
        if(!strcmp(txt_part,"PADDING")) input->padding = atoi(num_part);
        if(!strcmp(txt_part,"CHARGES")) input->charges[0][count] = atof(num_part);
	}

    initialize_input(input, pos_parameters, mom_parameters, field_parameters);

    free(pos_parameters);
    free(mom_parameters);
    free(field_parameters);

    fclose(file);

    return input;
}

void print_input(input_t* input){
    printf("Output filename = %s\n", input->filename);
    printf("Diagnostic Frequency = %d\n", input->pos_diag_freq);
    printf("Number of timesteps = %d\n", input->n_timesteps);
    printf("Length of timestep = %f\n", input->deltaT);
    printf("Number of configuration space dimentions = %d\n", input->pos_dims);
    printf("Number of momentum space dimentions = %d\n", input->mom_dims);
    for(int j = 0; j < input->pos_dims; ++j){
        printf("Number of configuration space points on axis %d = %d\n", j, input->pos_points[j]);
    }
    for(int j = 0; j < input->mom_dims; ++j){
        printf("Number of momentum space points on axis %d = %d\n", j, input->mom_points[j]);
    }
    for(int j = 0; j < input->pos_dims; ++j){
        printf("Initial distribution of position axis %d: %s\n", j, input->pos_init_names[j]);
    }
    for(int j = 0; j < input->mom_dims; ++j){
        printf("Initial distribution of momentum axis %d: %s\n", j, input->mom_init_names[j]);
    }
    for (int i = 0; i < input->pos_dims; ++i){
        printf("Limits of configuration space axis %d: [%f ; %f]\n", i, input->pos_min[i], input->pos_max[i]);
    }
    for (int i = 0; i < input->mom_dims; ++i){
        printf("Limits of momentum space axis %d: [%f ; %f]\n", i, input->mom_min[i], input->mom_max[i]);
    }
    for(int j = 0; j < input->pos_dims; ++j){
        printf("Boundary Condition of position axis %d: %s\n", j, input->pos_bound_names[j]);
    }
    for(int j = 0; j < input->mom_dims; ++j){
        printf("Boundary Condition of momentum axis %d: %s\n", j, input->mom_bound_names[j]);
    }
    printf("Number of Species = %d\n", input->n_species);
    printf("Number of Fields = %d\n", input->n_fields);
    for (int i = 0; i < input->n_fields; ++i){
        printf("Kernel for field %d: %s\n", i, input->kernel_names[i]);
    }
    for (int i = 0; i < input->n_species; ++i){
        printf("Dispersion Relation for species %d: %s\n", i, input->dispersion_names[i]);
    }
    printf("Time-advancing algorithm = %s\n", input->pusher);
    printf("Space-operator algorithm = %s\n", input->operator);
}

void free_input(input_t* input){
    
    free(input->filename);
    free(input->pos_points);
    free(input->pos_delta);
    free(input->lambda);
    free(input->pos_min);
    free(input->pos_max);
    free(input->mom_min);
    free(input->mom_max);
    
    for (int i = 0; i < input->n_species; ++i){
        free(input->dispersion_names[i]);
    }
    free(input->charges[0]);
    free(input->charges);
    free(input->dispersion_names);
    free(input->dispersions);
    free(input->pusher);
    free(input->operator);
    for (int i = 0; i < input->n_fields; ++i){
        free(input->kernel_names[i]);
        free(input->force_names[i]);
    }
    free(input->kernel_names);
    free(input->force_names);
    free(input->kernels);
    free(input->forces);

    free(input->pos_inits);
    free(input->mom_inits);
    free(input->field_inits);

    free(input->pos_init_params[0]);
    free(input->mom_init_params[0]);
    free(input->field_init_params[0]);
    free(input->pos_init_params);
    free(input->mom_init_params);
    free(input->field_init_params);

    free(input->pos_bound_names[0]);
    free(input->pos_init_names[0]);
    free(input->field_init_names[0]);
    free(input->pos_bound_names);
    free(input->pos_init_names);
    free(input->field_init_names);

    free(input->mom_bound_names[0]);
    free(input->mom_init_names[0]);
    free(input->mom_bound_names);
    free(input->mom_init_names);
    free(input->pos_bounds);
    free(input->mom_bounds);

    free(input->grid_factor);
    free(input->space_factor);
    free(input->dk);
    free(input->procs);
    free(input->parallel_factor);
    free(input->parallel_pos);

    free(input);
}