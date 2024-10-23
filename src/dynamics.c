#include "dynamics.h"

COMPLEX maxwell_2D_Ex(int fld, int qs, COMPLEX** fields, REAL* wavevector){

    return I * wavevector[1] * fields[fld+2][qs];

}

COMPLEX maxwell_2D_Ey(int fld, int qs, COMPLEX** fields, REAL* wavevector){

    return -I * wavevector[0] * fields[fld+1][qs];
    
}

COMPLEX maxwell_2D_Bz(int fld, int qs, COMPLEX** fields, REAL* wavevector){

    return I * (wavevector[1] * fields[fld-2][qs] - wavevector[0] * fields[fld-1][qs]);
}

COMPLEX null_dynamics(int fld, int qs, COMPLEX** fields, REAL* wavevector){

    return 0;
}