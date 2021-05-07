# dde
Delay differential equation solvers, special emphasis on neutral dde of the form

x' = f(t, x, x(t-tau), x'(t-tau))

with tau a delay function depending on t and x: tau(t, x)

* Dependencies : 

    needs cminpack-1.3.8, downloadable at https://github.com/devernay/cminpack

    Compilation and copie of libcminpackld.a in /usr/local/lib

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
    
        t0 vector of length N, X0 un vecteur de type NxK, ils représentent les conditions initiales
        cFileName est le nom du fichier C contenant fun et tau
        les paramètres sont des paramètres utilisateurs qui seront passés à fun et tau 
        alg est l'algorythme de résolution  utilisé, choisis parmi
            * rk4Neutral, schéma de type Runge Kutta 4 (explicite)
            * impTrNeutral, schéma de type trapèze implicite
            * eulerImpNeutral, schéma de type Euler implicite
            * rk4, de type Runge Kutta 4 pour les systèmes non neutres
        
    Toutes les méthodes sont à pas constant (déterminer par t0[1]-t0[0]
