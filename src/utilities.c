#include "utilities.h"

int factorial(int n){
    if(!n) return 1;
    return n*factorial(n-1);
}

int binomial_coeff(int n, int m){
    return factorial(n)/(factorial(m)*factorial(n-m));
}

int ipow(int base, int exp){
    return ((int)(pow(base,exp)));
}

int* list_combinations(int n, int m, int* list, int k){
    
    if(!m) return list;

    list[0] = 0;

    //printf("n = %d\n", n);
    //printf("m = %d\n", m);
    
    for (int i = 0; i < n; ++i){
        //printf("k = %d\n", k);
        //printf("comb = %d\n", binomial_coeff(n-1-list[0], m-1));
        if (k >= binomial_coeff(n-1-list[0], m-1)){
            k -= binomial_coeff(n-1-list[0], m-1);
            list[0]++;
        }
        else break;
    }

    //printf("\n");

    if(m>1) list_combinations(n-list[0]-1, m-1, list+1, k);
    
    for(int i = 1; i < m; ++i){
        list[i] += (list[0]+1);
    }

    return list;

}

int* list_combinations_all(int n, int m, int* list, int k){

    list_combinations(n, m, list, k);

    int* aux_list = malloc(n * sizeof(int));

    for (int i = 0; i < n; ++i){
        aux_list[i] = i;
    }
    for (int i = 0; i < m; ++i){
        aux_list[list[i]] = -1;
    }
    int counter = 0;
    for(int i = m; i < n; ++i){
        while(aux_list[counter] == -1){
            counter++;
        }
        list[i] = aux_list[counter];
        counter++;
    }

    free(aux_list);

    return list;
}

int list_index(int n_dims, int* Ns, int* is){

    int res = 0;

    for (int j = 0; j < n_dims; ++j){
        res *= Ns[j];
        res += is[j];
    }

    return res;
}

void axis_index(int n_dims, int* Ns, int* is, int index){

    int total = 1;

    for (int j = 0; j < n_dims; ++j){
        total *= Ns[j];
    }

    for (int j = 0; j < n_dims; ++j){
        total /= Ns[j];
        is[j] = (index / total) % Ns[j];
    }

}

REAL minmod3(REAL x1, REAL x2, REAL x3){
    if(x1*x2 < 0) return 0;
    if(x1*x3 < 0) return 0;
    if(fabs(x1) < fabs(x2)){
        if(fabs(x1) < fabs(x3)) return x1;
        return x3;
    }
    if(fabs(x2) < fabs(x3)) return x2;
    return x3;
}

REAL maxmod3(REAL x1, REAL x2, REAL x3){
    if(x1*x2 < 0) return 0;
    if(x1*x3 < 0) return 0;
    if(fabs(x1) > fabs(x2)){
        if(fabs(x1) > fabs(x3)) return x1;
        return x3;
    }
    if(fabs(x2) > fabs(x3)) return x2;
    return x3;
}

REAL minmod2(REAL x1, REAL x2){
    if(x1*x2<0) return 0;
    if(fabs(x1) < fabs(x2)) return x1;
    return x2;
}

REAL maxmod2(REAL x1, REAL x2){
    if(x1*x2<0) return 0;
    if(fabs(x1) < fabs(x2)) return x2;
    return x1;
}

REAL min2(REAL x1, REAL x2){
    if(x1 < x2) return x1;
    return x2;
}

REAL max2(REAL x1, REAL x2){
    if(x1 < x2) return x2;
    return x1;
}

REAL min3(REAL x1, REAL x2, REAL x3){
    if(x1 < x2){
        if(x1 < x3) return x1;
        return x3;
    }
    if(x2 < x3) return x2;
    return x3;
}

REAL max3(REAL x1, REAL x2, REAL x3){
    if(x1 > x2){
        if(x1 > x3) return x1;
        return x3;
    }
    if(x2 > x3) return x2;
    return x3;
}

REAL minmod_slopeLimiter(REAL u0, REAL u1, REAL u2, REAL Dx){
    return minmod2((u2-u1)/Dx, (u1-u0)/Dx);
}

REAL superbee_slopeLimiter(REAL u0, REAL u1, REAL u2, REAL Dx){
    return maxmod2(minmod2(2 * (u1 - u0)/Dx, (u2 - u1)/Dx), minmod2((u1 - u0)/Dx, 2 * (u2 - u1)/Dx));
}

REAL MC_slopeLimiter(REAL u0, REAL u1, REAL u2, REAL Dx){
    return minmod3(2 * (u2 - u1)/Dx, (u2 - u1)/(2 * Dx), 2 * (u1 - u0)/Dx);
}

REAL vanLeer_slopeLimiter(REAL u0, REAL u1, REAL u2, REAL Dx){
    REAL r = (u2-u1)/(u1-u0);
    return (r+fabs(r))/(1+fabs(r));
}

REAL minmod_fluxLimiter(REAL u0, REAL u1, REAL u2){
    return minmod2(1, (u1-u0)/(u2-u1));
}

REAL superbee_fluxLimiter(REAL u0, REAL u1, REAL u2){
    REAL theta = (u1-u0)/(u2-u1);
    return max3(0, min2(1, 2*theta), min2(2, theta));
}

REAL MC_fluxLimiter(REAL u0, REAL u1, REAL u2){
    REAL theta = (u1-u0)/(u2-u1);
    return max2(0, min3((1+theta)/2, 2, 2*theta));
}

REAL vanLeer_fluxLimiter(REAL u0, REAL u1, REAL u2){
    REAL theta = (u1-u0)/(u2-u1);
    return (theta+fabs(theta))/(1+fabs(theta));
}