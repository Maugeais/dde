#include <math.h> 
#include <stdio.h>

int dim = 2;


void fun(double t, double *x, double *xdelay, double *y){
    
    int i;
        
//    for (i=0; i < N; i++){
        //y[0] = x[0]*(0.1-0.2*x[1]);
        //y[1] = x[1]*(0.9*x[0]-0.4);
        
        y[0] = -M_PI/2*xdelay[0]/3;
        y[1] = 3*M_PI/2*xdelay[1]/3;
//    }
        
}

double tau(double t, double *x){
    
    return(3);
        
}
