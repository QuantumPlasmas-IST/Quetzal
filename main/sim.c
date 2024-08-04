#include "stdlib.h"
#include "stdio.h"
#include "input.h"
#include "pusher.h"
#include "aux.h"
#include "output.h"
#include "quetzal.h"
#include "hdf5.h"
#include "utilities.h"
#include "time.h"

int main(int argc, char** argv){

    double tstart = (double) time(NULL);

    if(argc>1){
        quetzal(argv[1]);
    }
    else{
        printf("No input file given!\n");
        return 1;
    }

    printf("Total Running Time: %fs\n", (double) time(NULL)-tstart);

    return 0;
}
