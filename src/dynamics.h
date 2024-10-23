#ifndef __DYNAMICS__
#define __DYNAMICS__

#include "aux.h"
#include "compile.h"

COMPLEX null_dynamics(int fld, int qs, COMPLEX** fields, REAL* wavevector);

COMPLEX maxwell_2D_Ex(int fld, int qs, COMPLEX** fields, REAL* wavevector);

COMPLEX maxwell_2D_Ey(int fld, int qs, COMPLEX** fields, REAL* wavevector);

COMPLEX maxwell_2D_Bz(int fld, int qs, COMPLEX** fields, REAL* wavevector);

#endif