#include <math.h> 
#include <stdio.h>



void fun(double t, double *x, double *xdelay, double *xdelayDer, double *der, double *params){
    
    der[0] = -0.5*x[0]-xdelayDer[0];
    
}

double tau(double t, double *x, double *params){
    
    return((t-0.5)*x[0]*x[0]);
    
}
