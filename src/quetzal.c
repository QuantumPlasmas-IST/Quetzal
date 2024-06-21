#include "quetzal.h"

void quetzal(const char* filename){

    input_t* input;

    input = read_input(filename);
    print_input(input);
    open_file(input);

    REAL** results = solve(input);

    free(results[0]);
    free(results);

    free_input(input);
    
}
