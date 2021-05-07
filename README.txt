Programme de résolution d'équation à retard neutre, de la forme 
x' = f(t, x, x(t-tau), x'(t-tau))

tau étant une fonction retard, tau = tau(t, x)

* Installation : 

    Il faut télécharger cminpack-1.3.8 par exemple sur https://github.com/devernay/cminpack

    Compilation + copie de libcminpackld.a dans /usr/local/lib

    Ajout dans le pythonpath

    export PYTHONPATH=$PYTHONPATH:/home/maugeais/Documents/Acoustique/Vents/Flute/dde/


* Fichier de données :

    Le fichier en C doit contenir deux fonctions : 
        - void fun(double t, double *x, double *xdelay, double *xdelayDer, double *der, double *params)
        - double tau(double t, double *x, double *params)
        
        t le temps (double)
        x un vecteur (double)
        xdelay vecteur de même taille que x (= x(t-tau))
        xdelayDer vecteur de même taille que x (=x'(t-tau))
        der vecteur de même taille que x (sortie x')
        params sont des paramêtres transmis par l'utilisateur
        
* Fichier python :

    import dde
        
    t, X = dde.rk4Delay(t0, X0, T, cFileName, params, alg) 
    
    - t0 est un vecteur de longueur N, X0 un vecteur de type NxK, ils représentent les conditions initiales
    - cFileName est le nom du fichier C contenant fun et tau
    - les paramètres sont des paramètres utilisateurs qui seront passés à fun et tau 
    - alg est l'algorythme de résolution  utilisé, choisis parmi
        * rk4Neutral, schéma de type Runge Kutta 4 (explicite)
        * impTrNeutral, schéma de type trapèze implicite
        * eulerImpNeutral, schéma de type Euler implicite
        * rk4, de type Runge Kutta 4 pour les systèmes non neutres
        
    Toutes les méthodes sont à pas constant (déterminer par t0[1]-t0[0]
