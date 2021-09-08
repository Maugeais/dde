#include <math.h> 
#include <stdio.h>

void fun(double t, double *x,double *der, double *params){
    
    der[0] = (params[0]-params[1]*x[1])*x[0];
    der[1] = (params[2]*x[0]-1)*x[1];
    
}

double tau(double t, double *x, double *params){
    
    /*double epsilon = 100, bruit = 0;
    if (t < epsilon) {
        bruit = cos(20*M_PI*t/epsilon);
    }*/
    return(0);
    
}
