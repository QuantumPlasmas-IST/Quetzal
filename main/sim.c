#include "stdlib.h"
#include "stdio.h"
#include "input.h"
#include "pusher.h"
#include "aux.h"
#include "output.h"
#include "quetzal.h"
#include "hdf5.h"
#include "utilities.h"


int main(int argc, char** argv){

    if(argc>1){
        quetzal(argv[1]);
    }
    else{
        printf("No input file given!\n");
        return 1;
    }

    return 0;
}
