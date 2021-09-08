#include <math.h> 
#include <stdio.h>

double omega1= 1440, omega2 = 2903, F1= 1322, F2= 2386, Q1= 36.6, Q2= 41.2;

double rampeGamma(double t, double *params) 
{
    double val;
    val = t*params[1]+params[0];
    return(val);
}


// Le vecteur d'etat est sousl a forme p1, p1', p2, p2', ...

void fun(double t, double *x, double *der, double *params){
    
    double gamma, zeta, up;
    
    gamma = rampeGamma(t, params);
    zeta = params[2];
    
    if (gamma-x[0]-x[2] < 1) {
        up = -(x[1]+x[3])*zeta*(1-3*(gamma-x[0]-x[2]))/(2*sqrt(gamma-x[0]-x[2]));
    } else {
        up = 0;
    }

    der[0] = x[1];  
    der[1] = -omega1/Q1*x[1]-omega1*omega1*x[0]+F1*up;

    der[2] = x[3];
    der[3] = -omega2/Q2*x[3]-omega2*omega2*x[2]+F2*up;
    //printf("belln = %f, %f, %f, %f\n", x[0], x[1], x[2], x[3]);
}

double tau(double t, double *x, double *params){
    
    return(0);
}
