#include <math.h> 
#include <stdio.h>


#define K 6

// Données page 187, flûte à bec alto Zen On Bressan,
const double a[K] = { 11.22, 22.36, 16.39, 12.64, 10.55, 10.32 };
const double omega[K] = { 1.6, 2511, 5113, 7570, 9719, 11910 };
const double Q[K] = { 3.31, 44.9, 59.64, 67.2, 73.57, 79.99 };


const double alpha = 340;
const double W = 0.00425; // Table 2.1 page 30
const double h = 0.001; // Table 2.1 page 30
const double rho = 1.2; // Table 2.1 page 30
const double Y0 = 0.0001; //Table 2.1 page 30
const double alphaI = 0.4/h; // equation 1.4 page 14
const double  b = 2*h/5; //equation 1.10 page 16
const double alphaVc = 0.6; // equation 1.12 page 16

double computeUj(double t, double *params){
    double Uj;
    if (params[0] > 0) {
        if (t < (omega[1]*params[2])/2) {
            Uj = sqrt((params[0]+t/(omega[1]*params[2])*(params[1]-params[0]))*2/rho); // equation 1.2 page 12
        } else {
            Uj = sqrt((params[0]+(omega[1]*params[2]-t)/(omega[1]*params[2])*(params[1]-params[0]))*2/rho);
        }
        
    } else {
        Uj = sqrt(params[1]*2/rho); // equation 1.2 page 12
    }
    return(Uj);
}


void fun(double t, double *x, double *xdelay, double *xdelayDer, double *der, double *params){
    
    double vAc, dVac, vAcmTau, dVAcmTau, ddVAcmTau, eta, dEta, ddEta, deltaP, dDeltaP;
    double Uj, tanH;
    int k;
    double deltaD;
    
    deltaD = 4/M_PI*sqrt(2*h*W); // equation 1.6 page 14  
    
    Uj = computeUj(t, params);

    vAc = 0;
    dVac = 0;
    
    vAcmTau = 0;
    dVAcmTau = 0;
    ddVAcmTau = 0; 
    
    for (k =0; k < K; k++){
        
        vAc += x[2*k];
        dVac += x[1+2*k];
        
        vAcmTau += xdelay[2*k];
        dVAcmTau += xdelay[1+2*k]; 
        ddVAcmTau += xdelayDer[1+2*k];
    }
    
    eta = h/Uj*exp(alphaI*W)*vAcmTau;
    dEta =  h/Uj*exp(alphaI*W)*dVAcmTau;
    ddEta = h/Uj*exp(alphaI*W)*ddVAcmTau;

    tanH = tanh((eta-Y0)/b);
        
    deltaP =  rho*deltaD*b*Uj/W*dEta/b*(1-tanH*tanH);
    dDeltaP = rho*deltaD*b*Uj/W*(ddEta/b-2*tanH*dEta*dEta/(b*b))*(1-tanH*tanH);

    if (vAc >= 0) {
        
         deltaP += -rho/2*(vAc/alphaVc)*(vAc/alphaVc)/omega[1];
         
         dDeltaP += -rho*dVac*vAc/(alphaVc*alphaVc)/omega[1];
         
    } else {
        deltaP -= -rho/2*(vAc/alphaVc)*(vAc/alphaVc)/omega[1];
        dDeltaP -= -rho*dVac*vAc/(alphaVc*alphaVc)/omega[1];
    }
        
    // x = [v1, y1, ...]
    
    der[0] = (a[0]*deltaP - Q[0]/omega[1]*x[0])/omega[0];
    der[1] = (a[0]*dDeltaP - Q[0]/omega[1]*x[1])/omega[0];

    

    
    for (k=1; k < K; k++){
        der[2*k] = x[2*k+1];
        der[2*k+1] = -omega[k]/omega[1]/Q[k]*x[2*k+1]-omega[k]/omega[1]*omega[k]/omega[1]*x[2*k] + a[k]*dDeltaP;
    }   
    
}

double tau(double t, double *x, double *params){
    
    double epsilon = 100, bruit = 0, tau, Uj;
    if (t < epsilon) {
        bruit = 0*cos(20*M_PI*t/epsilon);
    }
    
    Uj = computeUj(t, params);
    
        
    tau = omega[1]*W/(0.4*Uj); // equation 1.5 page 14 

    return(t-tau); //+(params[1]/2-fabs(t-params[1]/2))*params[0]/params[1]+bruit));
    
}
