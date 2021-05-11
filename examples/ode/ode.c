#include <math.h> 
#include <stdio.h>


void fun(double t, double *x, double *xdelay, double *xdelayDer, double *der, double *params){
    
    der[0] = -0.5*x[0];
    
    //printf("belln = %f, %f\n", x[0], der[0]);
}

double tau(double t, double *x, double *params){
    
    /*double epsilon = 100, bruit = 0;
    if (t < epsilon) {
        bruit = cos(20*M_PI*t/epsilon);
    }*/
    return(0);
    
}
