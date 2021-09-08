#include <math.h> 
#include <stdio.h>
#include <complex.h> 

#define N 5
double omegal = 697.05, Ql = 7.0, mu = 9.0909, H = 0.0005, a = 5.6e-2, b=9.5e-5, w = 0.012;
double complex S[N] = {-12.86 +238.63*I,-17.35  +697.05*I,-22.21+1061.27*I,
 -26.05+1437.42*I,-28.75+1828.96*I}; // ok
double complex C[N] = {585.62+1.327*I ,635.81 +3.877*I,691.38 +5.902*I,
 730.67 +7.995*I,927.80  +10.172*I};

#define rho 1.1851
#define R 0.012 // Rayon à l'embouchure
#define c0 346.19

double Zc = rho*c0/(M_PI*R*R); //Impedance caracteristique du résonateur

// fonction de heavisyde
double theta(double x){
    if (x > 0){
        return(x);
    } else {
        return(0);
    }
}

// Le vecteur d'etat est sous la forme complexe : h, h', real(p1), image(p1), ...

void fun(double t, double *x, double *der, double *params){
    
    double complex pn[N], s;
    double u, pb, h, p = 0;
    int i;
        
    pb = 8000;
    //omegal = 2*M_PI*params[2];
    h = x[0];
    /// Plus précis pour le H
    //H = a+2*M_PI*b/omegal;
        
    for (i=0; i < N; i++) { // ok
        pn[i] = x[2*i+2]+x[2*i+3]*I;
        p += 2*x[2*i+2];
    }
        
    // Expression exacte du débit

    if (pb > p) { // ok
        u = w*theta(h)*sqrt(2*fabs(pb-p)/rho);
    } else {
        u = -w*theta(h)*sqrt(2*fabs(pb-p)/rho);
    }

    /*! OU régularisation

        !q = zeta*h*THETA(h)*(pb-p)*((pb-p)**2+eta)**(-0.25)*/
    //zeta = params[2];
    
    der[0] = x[1];
    der[1] = -omegal/Ql*x[1]-omegal*omegal*(x[0]-H) + (pb-p)/mu; // ok
        
    for (i=0; i < N; i++){
        
        s = S[i]*pn[i]+Zc*C[i]*u;
        der[2*i+2] = creal(s); 
        der[2*i+3] = cimag(s);
    }
}

double tau(double t, double *x, double *params){
    
    return(0);
}
