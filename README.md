# dde, beta version 0.1
Delay differential equation solvers, special emphasis on neutral dde of the form

x' = f(t, x, x(t-tau), x'(t-tau))

with tau a delay function depending on t and x: tau(t, x)

* Dependencies : 

    needs cminpack-1.3.8, downloadable at https://github.com/devernay/cminpack (make and make install)

    Add to the pythonpath 

    export PYTHONPATH=$PYTHONPATH:/home/maugeais/Documents/Acoustique/Vents/Flute/dde/
    
    Inside dde/src, type make


* Model file :

    The C must contain two functions : 
        - void fun(double t, double *x, double *xdelay, double *xdelayDer, double *der, double *params)
        - double tau(double t, double *x, double *params)
        
        t = time (double)
        x = position (double *)
        xdelay = delayed position x(t-tau), (double *)
        xdelayDer = delayed speed x'(t-tau), (double *)
        der = speed x' (output)
        params = parameters given by the user 
        
* Python file :

    import dde
        
    t, X = dde.rk4Delay(t0, X0, T, cFileName, params, alg)
    
    Input
    
        t0 vector of length N, X0 of shape NxK, it stores the initial conditions
        T is the end time for computation. All the methods use a constant step, given by t0[1]-t0[0]
        cFileName, name of the C file containing fun and tau
        params are parameters defined by the user that are passed to fun and tau 
        alg is the choice of algorythm
            * rk4Neutral, explicit Runge Kutta 4 scheme, adapted for delay neutral equations
            * impTrNeutral, implicite trapzoid scheme, adapted for deley neutral equations
            * eulerImpNeutral, simple implicit euler scheme for delay neutral equations
            * rk4, explicit Runge Kutta 4 scheme, for delay equations (non neutral !!!!)
        
    

    Output
    
        t, X : results of the computation. 
        t is of length n = len(t0) + T/(t0[1]-t0[0])
        X is of shape nxK
