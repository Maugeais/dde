#include <math.h> 
#include <stdio.h>

const int m = 2;

const double a[2] = { 10, 30 };
const double omega[2] = { 2764, 5510 };
const double Q[2] = { 100, 100 };
const double alpha = 340;



void fun(double t, double *x, double *xdelay, double *der, double *params){
    
    double ydelay = 0, y = 0, vdelay = 0;
    int k;
    
    // x = [v1, y1, ...]
    
    for (k = 0; k < m; k++){
        vdelay += xdelay[2*k];
        ydelay += xdelay[2*k+1];
        y += x[2*k+1];
    } 
        
    for (k=0; k < m; k++){
        der[2*k] = x[2*k+1];
        der[2*k+1] = alpha*a[k]/omega[0]*(1-tanh(vdelay)*tanh(vdelay))*ydelay-(omega[k]*omega[k])/(omega[0]*omega[0])*x[2*k]-omega[k]/(omega[0]*Q[0])*x[2*k+1];
    }    
}

double tau(double t, double *x, double *params){
    
    double epsilon = 100, bruit = 0;
    if (t < epsilon) {
        bruit = cos(20*M_PI*t/epsilon);
    }
    return(t-(0.1+(params[1]/2-fabs(t-params[1]/2))*params[0]/params[1]+bruit));
    
}
