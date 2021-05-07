#include <stdio.h>
#include <math.h> 
#include <stdlib.h>
#include "cminpack.h"
#define real __cminpack_real__

real *invDerPol;
int NN;
real *vanDerPol;

void preDerPol(int N, float h){
    
    int i, j, k;   
    
    int *pivotArray;
    int errorHandler;
    real p;

    NN = N;
    vanDerPol = malloc(N*N * sizeof(real));
    
    for (i=0; i < NN; i++) {
        for (j=0; j < NN; j++) {
           *(vanDerPol+NN*i+j) = pow(i*h, j); 
        }
        
    }
      
    invDerPol = malloc(N*N * sizeof(real));
    
    for (i=0; i < NN; i++) {
        for (j=0; j < NN; j++) {
            
           if (i == j){
               *(invDerPol+NN*i+j) = 1; 
           } else {
               *(invDerPol+NN*i+j) = 0;
            }
           
        }
        
    }

    // Gauss, descente, pas d'échange de pivot car VanderMonde
    for (j = 0; j < NN; j++){
        for (i = j+1; i < NN; i++){
            p = *(vanDerPol+NN*i+j)/(*(vanDerPol+NN*j+j));
            for (k = j; k < NN; k++){
                *(vanDerPol+NN*i+k) -= p*(*(vanDerPol+NN*j+k));
            }
            
            for (k = 0; k <= j; k++){
                *(invDerPol+NN*i+k) -= p*(*(invDerPol+NN*j+k));
            }
        }
    }
    
    // Gauss, normalisation
    for (i = 0; i < NN; i++){
        p = *(vanDerPol+NN*i+i);
        
        for (j = 0; j < NN; j++){
            *(vanDerPol+NN*i+j) /= p;
            *(invDerPol+NN*i+j) /= p;
            
        }
    }
    
    // Gauss, remontée, les coeffs diag sont à 1
    for (j = NN-1; j >= 0; j--){
        for (i = 0; i < j; i++){
            p = *(vanDerPol+NN*i+j);

            for (k = j; k < NN; k++){
                *(vanDerPol+NN*i+k) -= p*(*(vanDerPol+NN*j+k));
            }
            
            for (k = 0; k < NN; k++){
                *(invDerPol+NN*i+k) -= p*(*(invDerPol+NN*j+k));
            }
        }
    }
    
}

// Recherche le point le plus proche par dichotomie
int bissect(double *t, int N, double t0){
    int i, j, k;
    
    if (t0 < *t) {
        return(-1);
    }
    
    if (t0 > *(t+N-1)){
        return(-1);
    }
    
    i = 0;
    j = N-1;
    while ( i < j-1 ){
        k = (i+j) / 2;
        

        if (*(t+k) <= t0){
            i = k;
        } else {
            j = k;
        }
    }
     
    return(i);
}

// Interpole X(t) en t0 avec un degré ...
void interp(int dim, int N, double *t, double* X, double tItrp, double *XItrp){
    int i, j, k, m;
    const int k0 = 1;
    double h, kItrp;
    
    h = t[1]-t[0];
    
    // Recherche l'indice de t0 dans t
    k = bissect(t, N, tItrp);
    // On ramène t[k] dans l'intervalle 
    kItrp = (tItrp-t[k])/h+k0;
    
    /*if (k==-1) {
        throw("b cannot be zero!");
    }*/
    
    double *polItrp;
    
    polItrp = malloc(NN * sizeof(double));
    
    // Pour chaque dimension
    
    for (m=0; m < dim; m++){
        
        // Calcul le polynôme interpolant
        for (i=0; i < NN; i++){
            polItrp[i] = 0;
            
            for (j = 0; j< NN; j++){
                polItrp[i] += *(invDerPol+NN*i+j)*(*(X+dim*(j+k-k0)+m));
            }
        }
        
        // Puis évalue par le schéma de Horner
        XItrp[m] = 0;
        
        for (i=0; i < NN; i++){
            XItrp[m] = kItrp*XItrp[m]+polItrp[NN-1-i];
        }
        
    }
}

void interpDer(int dim, int N, double *t, double* X, double tItrp, double *XItrp, double *XpItrp){
    int i, j, k, m;
    const int k0 = 1;
    double h, kItrp;
    double leftDer, rightDer;
    
    h = t[1]-t[0];
    
    // Recherche l'indice de t0 dans t
    k = bissect(t, N, tItrp);
    // On ramène t[k] dans l'intervalle 
    kItrp = (tItrp-t[k])/h+k0;
    
    /*if (k==-1) {
        throw("b cannot be zero!");
    }*/
    
    double *polItrp;
    
    polItrp = malloc(NN * sizeof(double));
    
    // Pour chaque dimension
    
    for (m=0; m < dim; m++){
        
        // Calcul le polynôme interpolant
        for (i=0; i < NN; i++){
            polItrp[i] = 0;
            
            for (j = 0; j< NN; j++){
                polItrp[i] += *(invDerPol+NN*i+j)*(*(X+dim*(j+k-k0)+m));
            }
        }
        
        // Puis évalue par le schéma de Horner
        XItrp[m] = 0;
        XpItrp[m] = 0;

        
        for (i=NN-1; i >= 0; i--){
            XItrp[m] = kItrp*XItrp[m]+polItrp[i];
            
            if (i > 0 ){
                XpItrp[m] = (kItrp*XpItrp[m]+(i)*polItrp[i])/h;
            }
 
        } 
        
        // Si non dérivable, alors calcul à droite ou à gauche...
        leftDer = (*(X+dim*(k+1)+m)-*(X+dim*k+m))/h;
        rightDer = (*(X+dim*k+m)-*(X+dim*(k-1)+m))/h;
        
        if ( fabs(leftDer-rightDer) > fmin(fabs(leftDer), fabs(rightDer))) {
        
          // printf("Discontinuité de la dérivée à %f \n", tItrp);
            if (fabs(leftDer)< fabs(rightDer)) {
                XpItrp[0] = leftDer;
            } else {
                XpItrp[0] = rightDer;
            }
        }
        
    }
    
//   printf("%f, %f\n", *(X+dim*(j+k-k0)+0), XItrp[0]);

}


struct re_val {

    double *t;
    double *x;  

};

struct re_val rk4Neutral(int dim, int N0, double *t0, double *X0, float T, void (*f)(double, double *, double*,  double *, double *, double *), double (*tau)(double, double *, double *), int N, double *params){
    int i, j;
    double *t, *X, k1[dim], k2[dim], k3[dim], k4[dim], y[dim], h;
    double delayed[dim], delayedDer[dim];
    h = t0[1]-t0[0];
    
    t = malloc((N+N0) * sizeof(double));
    X = malloc((N+N0) * dim * sizeof(double));
    
    // Création du vecteur temps
    for (i = 0; i < N0; i++){
        t[i] = t0[i];
    }
    
    for (i = N0; i < N+N0; i++){
        t[i] = t[i-1]+h;
    }

    // Condition initiale
    for (j=0; j < dim*N0; j++){
        *(X+j) = *(X0+j);
    }

    // Boucle principale
    for (i = N0-1; i <= N+N0-2; i++){
        
        // Calcul de k1
          //retard        
        interpDer(dim, N+N0, t, X, tau(t[i], X+dim*i, params), delayed, delayedDer);
        f(*(t+i), X+dim*i, delayed, delayedDer, k1, params);
        
        
        // Calcul de k2
        for (j = 0; j < dim; j++){
            y[j] = *(X+j+dim*i)+h/2*k1[j];
        }
        interpDer(dim, N+N0, t, X, tau(t[i]+h/2, y, params), delayed, delayedDer);
        f(*(t+i)+h/2, y, delayed, delayedDer, k2, params);
        
        
        // Calcul de k3
        for (j = 0; j < dim; j++){
            y[j] = *(X+j+dim*i)+h/2*k2[j];
        }
        interpDer(dim, N+N0, t, X, tau(t[i]+h/2, y, params), delayed, delayedDer);
        f(*(t+i)+h/2, y, delayed, delayedDer, k3, params);
        
        
        // Calcul de k4
        for (j = 0; j < dim; j++){
            y[j] = *(X+j+dim*i)+h*k3[j];
        }
        interpDer(dim, N+N0, t, X, tau(t[i]+h/2, y, params), delayed, delayedDer);
        f(*(t+i)+h, y, delayed, delayedDer, k4, params);
        
        // Calcul du résultat
        for (j = 0; j < dim; j++){
            *(X+j+dim*(i+1)) = *(X+j+dim*i)+h/6*(k1[j]+2*k2[j]+2*k3[j]+k4[j]);
        }
    }
        
    struct re_val r;
    r.t = t;
    r.x = X;
    
    return(r);
}

int fcn2(void *p, int n, const real *x, real *fvec, int iflag){
    printf("test function\n");
    return(0);
};

struct re_val impTrNeutral(int dim, int N0, double *t0, double *X0, float T, void (*f)(double, double *, double*,  double *, double *, double *), double (*tau)(double, double *, double *), int N, double *params){
    int i, j;
    double *t, *X, k1[dim], k2[dim], k3[dim], k4[dim], y[dim], h;
    double delayed[dim], delayedDer[dim];
    
    
    int maxfev, ml, mu, mode, nprint, info, nfev, ldfjac, lr;
    real xtol, epsfcn, factor, fnorm;
    real x[dim], fvec[dim], diag[dim], fjac[dim*dim], r[dim*(dim+1)/2], qtf[dim], wa1[dim], wa2[dim], wa3[dim], wa4[dim];

    ldfjac = dim;
    lr = dim*(dim+1)/2;

    xtol = sqrt(__cminpack_func__(dpmpar)(1));
    
    maxfev = 2000;
    ml = dim;
    mu = dim;
    epsfcn = 0.;
    mode = 2;
    
    for (j=0; j< dim; j++) {
        diag[j] = 1.;
    } 
  

    factor = 1.e2;
    nprint = 0;
  
  
    
    h = t0[1]-t0[0];
    
    t = malloc((N+N0) * sizeof(double));
    X = malloc((N+N0) * dim * sizeof(double));
    
    // Création du vecteur temps
    for (i = 0; i < N0; i++){
        t[i] = t0[i];
    }
    
    for (i = N0; i < N+N0; i++){
        t[i] = t[i-1]+h;
    }
    
    // Condition initiale
    for (j=0; j < dim*N0; j++){
        X[j] = X0[j];
    }
   

    // Boucle principale
    for (i = N0-1; i <= N+N0-2; i++){ //N+N0-dim; i++){
                
          //retard        
        interpDer(dim, N+N0, t, X, tau(t[i], X+dim*i, params), delayed, delayedDer);
        f(*(t+i), X+dim*i, delayed, delayedDer, k1, params);
        
        // Starting point
        for (j = 0; j < dim; j++){
            x[j] = *(X+dim*i+j);
        }
        
        int fcn(void *p, int n, const real *x, real *fvec, int iflag){
            
            int k;
            double xtemp[dim];
            
            for (k=0; k < dim; k++){
                xtemp[k] = x[k];
            }
            f(*(t+i), xtemp, delayed, delayedDer, k2, params);                

            for (k = 0; k < dim; k++){
                fvec[k] = xtemp[k]-(*(X+dim*i+k)+h/2*(k1[k]+k2[k]));
                
            }
            
            return(0);
        }
        
        
        
        

        info = __cminpack_func__(hybrd)(fcn, 0, dim, x, fvec, xtol, maxfev, ml, mu, epsfcn,
                                        diag, mode, factor, nprint, &nfev,
                                        fjac, ldfjac, r, lr, qtf, wa1, wa2, wa3, wa4);
        
        // Calcul du résultat
        for (j = 0; j < dim; j++){
            X[dim*(i+1)+j] = x[j]; 
        }
        
    }
    
    
    struct re_val res;
    res.t = t;
    res.x = X;
    
    return(res);
}



struct re_val eulerImpNeutral(int dim, int N0, double *t0, double *X0, float T, void (*f)(double, double *, double*,  double *, double *, double *), double (*tau)(double, double *, double *), int N, double *params){
    int i, j;
    double *t, *X, k1[dim], k2[dim], k3[dim], k4[dim], y[dim], h;
    double delayed[dim], delayedDer[dim];
    
    
    int n, maxfev, ml, mu, mode, nprint, info, nfev, ldfjac, lr;
    real xtol, epsfcn, factor, fnorm;
    real x[dim], fvec[dim], diag[dim], fjac[dim*dim], r[dim*(dim+1)/2], qtf[dim], wa1[dim], wa2[dim], wa3[dim], wa4[dim];

    printf("test\n");
    n = dim;

    ldfjac = dim;
    lr = dim*(dim+1)/2;
    

    xtol = sqrt(__cminpack_func__(dpmpar)(1));
    
    maxfev = 2000;
    ml = dim;
    mu = dim;

    epsfcn = 0.;
    mode = 2;
  
    for (j=0; j< dim; j++) {
            diag[j] = 1.;
        } 
  

    factor = 1.e2;
    nprint = 0;
    
  
    
    h = t0[1]-t0[0];
    
    t = malloc((N+N0) * sizeof(double));
    X = malloc((N+N0) * dim * sizeof(double));

    
    // Création du vecteur temps
    for (i = 0; i < N0; i++){
        t[i] = t0[i];
    }
    
    for (i = N0; i < N+N0; i++){
        t[i] = t[i-1]+h;
    }
    
        
    
    // Condition initiale
    for (j=0; j < dim*N0; j++){
        X[j] = X0[j];
    }


    // Boucle principale
    for (i = N0-1; i <= N+N0-dim-2; i++){ 
        
        // Calcul de k1
          //retard        
        interpDer(dim, N+N0, t, X, tau(t[i], X+dim*i, params), delayed, delayedDer);
        f(*(t+i), X+dim*i, delayed, delayedDer, k1, params);

        int fcn(void *p, int n, const real *x, real *fvec, int iflag){
            
            int k;
            double xtemp[dim];
            
            for (k=0; k < dim; k++){
                xtemp[k] = x[k];
            }
            f(*(t+i), xtemp, delayed, delayedDer, k1, params);                
            
            for (k = 0; k < dim; k++){
                fvec[k] = xtemp[k]-(*(X+dim*i+k)+h*k1[k]);
                
            }
            
            return(0);
        }
        
        for (j = 0; j < dim; j++){
            x[j] = *(X+dim*i+j);
        }
        
        info = __cminpack_func__(hybrd)(fcn, 0, n, x, fvec, xtol, maxfev, ml, mu, epsfcn,
                                        diag, mode, factor, nprint, &nfev,
                                        fjac, ldfjac, r, lr, qtf, wa1, wa2, wa3, wa4);
        
        // Calcul du résultat
        for (j = 0; j < dim; j++){
            *(X+j+dim*(i+1)) = *(X+j+dim*i)+h/2*(k1[j]); 
        }
    }
        
    struct re_val res;
    res.t = t;
    res.x = X;
    
    return(res);
}

struct re_val rk4(int dim, int N0, double *t0, double *X0, float T, void (*f)(double, double *, double*,  double *, double *), double (*tau)(double, double *, double *), int N, double *params){
    int i, j;
    double *t, *X, k1[dim], k2[dim], k3[dim], k4[dim], y[dim], h;
    double delayed[dim];
    h = t0[1]-t0[0]; 
    
    
    t = malloc((N+N0) * sizeof(double));
    X = malloc((N+N0) * dim * sizeof(double));
    
    // Création du vecteur temps
    for (i = 0; i < N0; i++){
        t[i] = t0[i];
    }
    
    for (i = N0; i < N+N0; i++){
        t[i] = t[i-1]+h;
    }
    
    // Condition initiale
    for (j=0; j < dim*N0; j++){
        *(X+j) = *(X0+j);
    }
    
    // Boucle principale
    for (i = N0-1; i <= N+N0-2; i++){
        
        // Calcul de k1
          //retard        
        interp(dim, N+N0, t, X, tau(t[i], X+dim*i, params), delayed);
        f(*(t+i), X+dim*i, delayed, k1, params);
        
        
        // Calcul de k2
        for (j = 0; j < dim; j++){
            y[j] = *(X+j+dim*i)+h/2*k1[j];
        }
        interp(dim, N+N0, t, X, tau(t[i]+h/2, y, params), delayed);
        f(*(t+i)+h/2, y, delayed, k2, params);
        
        
        // Calcul de k3
        for (j = 0; j < dim; j++){
            y[j] = *(X+j+dim*i)+h/2*k2[j];
        }
        interp(dim, N+N0, t, X, tau(t[i]+h/2, y, params), delayed);
        f(*(t+i)+h/2, y, delayed, k3, params);
        
        
        // Calcul de k4
        for (j = 0; j < dim; j++){
            y[j] = *(X+j+dim*i)+h*k3[j];
        }
        interp(dim, N+N0, t, X, tau(t[i]+h/2, y, params), delayed);
        f(*(t+i)+h, y, delayed, k4, params);
        
        // Calcul du résultat
        for (j = 0; j < dim; j++){
            *(X+j+dim*(i+1)) = *(X+j+dim*i)+h/6*(k1[j]+2*k2[j]+2*k3[j]+k4[j]);
        }
    }
        
    struct re_val r;
    r.t = t;
    r.x = X;
    
    return(r);
}
