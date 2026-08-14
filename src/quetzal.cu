// Quetzal Kinetic Simulation

#include "quetzal.h"

void quetzal(const char* filename){

    // Initialize Parallelization
    int rank,size;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    input_t* input;

    for(int i = 0; i < size; ++i){
        if(rank==i) input = read_input(filename, rank, size);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    open_file(input);
    
    if(!rank){
        print_input(input);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    REAL** results = solve(input);

    free(results[0]);
    free(results);

    free_input(input);

    MPI_Finalize();
    
}
