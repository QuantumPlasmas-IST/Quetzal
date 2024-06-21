#include "algos.h"

void LUdecomp(REAL** A, REAL** L, REAL** U, int N){

    for (int i = 0; i<N; ++i){
        for (int j=0; j<N; ++j){
            L[i][j]=0;
            U[i][j]=0;
            if (i==j) L[i][j]=1;
        }
    }

    for (int i = 0; i<N;++i){
        for(int j = i; j<N; ++j){
            U[i][j]=A[i][j];
            for(int k = 0; k<i; ++k){
                U[i][j]=U[i][j]-L[i][k]*U[k][j];
            }
        }
        for(int j = i+1; j<N; ++j){
            L[j][i]=A[j][i];
            for(int k = 0; k<i; ++k){
                L[j][i]=L[j][i]-L[j][k]*U[k][i];
            }
            L[j][i]=L[j][i]/U[i][i];
        }
    }
}

void LUinverse(REAL** A, REAL** L, REAL** U, int N){
    
    REAL* Bcol=malloc(N*sizeof(REAL));
    REAL* Icol=malloc(N*sizeof(REAL));
    
    for(int i=0; i<N; ++i){ //loop for each column of the inverse matrix (Bcol)
        for(int j=0; j<N; ++j){ //loop to fill Icol
            Icol[j]=0;
            Bcol[j]=0;
            if(j==i) Icol[j]=1;
        }
        
        LUsolve(L, U, Icol, Bcol, N);

        // Put Bcol in A
        for(int j=0; j<N; ++j){
            A[j][i]=Bcol[j];
        }
    }

    free(Bcol);
    free(Icol);
}


void LUsolve(REAL** L, REAL** U, REAL* b, REAL* x, int N){

    REAL* z = malloc(N*sizeof(REAL));
    
    //Solve L.z=b
    for(int i=0; i<N; ++i){
        REAL aux=0;
        for(int j=0; j<i; ++j){
            aux+=L[i][j]*z[j];
        }
        z[i]=b[i]-aux;
    }
    //Solve U.x=z
    for(int i=N-1; i>=0; --i){
        REAL aux2=0;
        for(int j=i+1; j<N; ++j){
            aux2+=U[i][j]*x[j];
        }
        x[i]=(z[i]-aux2)/U[i][i];
    }

    free(z);
}

// Performs the discrete fourier transform of function f(t) to g(w)
void fastFourierAuxiliary(int N, double complex* ft, double complex* gw, int s){

    if(N==1){
        gw[0] = ft[0];
        return;
    }

    fastFourierAuxiliary(N/2, ft, gw, 2*s);
    fastFourierAuxiliary(N/2, ft+s, gw+N/2, 2*s);

    for (int k = 0; k < N/2; ++k){
        double complex p = gw[k];
        double complex q = gw[N/2+k] * cexp(- 2 * M_PI * I / N * k);
        gw[k] = p + q;
        gw[k+N/2] = p - q;
    }
}

void fastInverseFourierAuxiliary(int N, double complex* ft, double complex* gw, int s){

    if(N==1){
        gw[0] = ft[0];
        return;
    }

    fastInverseFourierAuxiliary(N/2, ft, gw, 2*s);
    fastInverseFourierAuxiliary(N/2, ft+s, gw+N/2, 2*s);

    for (int k = 0; k < N/2; ++k){
        double complex p = gw[k];
        double complex q = gw[N/2+k] * cexp(2 * M_PI * I / N * k);
        gw[k] = p + q;
        gw[k+N/2] = p - q;
    }
}

void fastFourier(int N, double complex* ft, double complex* gw, int s){
    fastFourierAuxiliary(N, ft, gw, s);
}

void fastInverseFourier(int N, double complex* ft, double complex* gw, int s){
    fastInverseFourierAuxiliary(N, ft, gw, s);
    
    for(int i = 0; i < N; ++i){
        gw[i] /= N;
    }
}

void NdFourier(int n_dims, int* N, double complex* ft, double complex* gw){

    if(n_dims == 1){
        fastFourier(N[0], ft, gw, 1);
        return;
    }

    int skip = 1;
    for (int i = 1; i < n_dims; ++i){
        skip *= N[i];
    }

    for (int j = 0; j < N[0]; ++j){
        NdFourier(n_dims-1, N+1, ft+j*skip, gw+j*skip);
    }

    for (int i = 0; i < skip; ++i){
        fastFourier(N[0], gw + i, ft + i*N[0], skip);
    }

    for (int i = 0; i < N[0]; ++i){
        for (int j = 0; j < skip; ++j){
            gw[i*skip + j] = ft[j*N[0] + i];
        }
    }
}

void NdInverseFourier(int n_dims, int* N, double complex* ft, double complex* gw){

    if(n_dims == 1){
        fastInverseFourier(N[0], ft, gw, 1);
        return;
    }

    int skip = 1;
    for (int i = 1; i < n_dims; ++i){
        skip *= N[i];
    }

    for (int j = 0; j < N[0]; ++j){
        NdInverseFourier(n_dims-1, N+1, ft+j*skip, gw+j*skip);
    }

    for (int i = 0; i < skip; ++i){
        fastInverseFourier(N[0], gw + i, ft + i*N[0], skip);
    }

    for (int i = 0; i < N[0]; ++i){
        for (int j = 0; j < skip; ++j){
            gw[i*skip + j] = ft[j*N[0] + i];
        }
    }
}

void fourierTransform(int n_dims, int* N, double complex* ft, double complex* gw, double* dt, double* t0){
    
    NdFourier(n_dims, N, ft, gw);

    int totalN = 1;
    for (int i = 0; i < n_dims; i++){
        totalN *= N[i];
    }
    int auxN = totalN;

    for (int i = 0; i < n_dims; i++){
        auxN/=N[i];
        for (int j = 0; j < totalN; ++j){
            gw[j] *= dt[i];
            gw[j] *= cexp(- 2*M_PI * I * t0[i]/(N[i]*dt[i]) * ((j/auxN)%N[i]));
        }
    } 
}

void inverseFourierTransform(int n_dims, int* N, double complex* ft, double complex* gw, double* dt, double* t0){
    
    int totalN = 1;
    for (int i = 0; i < n_dims; i++){
        totalN *= N[i];
    }
    int auxN = totalN;

    for (int i = 0; i < n_dims; i++){
        auxN/=N[i];
        for (int j = 0; j < totalN; ++j){
            gw[j] *= cexp(2*M_PI*I * t0[i]/(N[i]*dt[i]) * ((j/auxN)%N[i]));
        }
    }
    auxN=totalN;

    NdInverseFourier(n_dims, N, ft, gw);

    for (int i = 0; i < n_dims; i++){
        auxN/=N[i];
        for (int j = 0; j < totalN; ++j){
            gw[j] /= dt[i];
        }
    } 
}

void hadamardProduct(int n_dims, int* N, double complex* f1, double complex* f2, double complex* g0){

    int totalN = 1;
    for (int i = 0; i < n_dims; i++){
        totalN *= N[i];
    }

    for (int j = 0; j < totalN; ++j){
        g0[j] = f1[j] * f2[j];
    }

}

void kroneckerProduct(int n_dims, int* N1, int* N2, double complex* f1, double complex* f2, double complex* g0){

    if(!n_dims){
        g0[0] = f1[0]*f2[0];
        return;
    }

    int totalN1 = 1;
    int totalN2 = 1;
    for (int i = 1; i < n_dims; i++){
        totalN1 *= N1[i];
        totalN2 *= N2[i];
    }

    for (int j = 0; j < N1[0]; ++j){
        for (int k = 0; k < N2[0]; ++k){
            kroneckerProduct(n_dims-1, N1+1, N2+1, f1+j*totalN1, f2+k*totalN2, g0+j*totalN1*totalN2*N2[0]+k*totalN1*totalN2);
        }
    }

}

void FFTConvolution(int n_dims, int* N, double complex* f1, double complex* f2, double complex* g0){

    NdFourier(n_dims, N, f1, g0);
    NdFourier(n_dims, N, f2, f1);

    hadamardProduct(n_dims, N, g0, f1, f2);

    NdInverseFourier(n_dims, N, f2, g0);
}