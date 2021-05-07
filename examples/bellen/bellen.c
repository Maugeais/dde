#include <math.h> 
#include <stdio.h>

const int m = 2;

const double a[2] = { 10, 30 };
const double omega[2] = { 2764, 5510 };
const double Q[2] = { 100, 100 };
const double alpha = 340;



void fun(double t, double *x, double *xdelay, double *xdelayDer, double *der, double *params){
    
    der[0] = -0.5*x[0]-xdelayDer[0];
    
    //printf("belln = %f, %f\n", x[0], der[0]);
}

double tau(double t, double *x, double *params){
    
    /*double epsilon = 100, bruit = 0;
    if (t < epsilon) {
        bruit = cos(20*M_PI*t/epsilon);
    }*/
    return((t-0.5)*x[0]*x[0]);
    
}
