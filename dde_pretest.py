import numpy as np
import matplotlib.pyplot as plt

from bisect import bisect_left

interp0, interp1 = 0, 4
vdmInv = []

def interp(T, X, t) :
    
    global vdmInv, P
    
    if len(vdmInv) == 0 :
        
        """ On se ramène à l'intervalle [-interp0*h, interp1*h] par translation"""
        h = T[1]-T[0]
        T0 = np.arange(-interp0*h, interp1*h, h)
        
        vdm = np.array([T0**i for i in range(len(T0))])
        vdmInv = np.linalg.inv(vdm)   
        
    print(vdm)
    # Identifie la position de t dans l'intervalle T
    
    pos = bisect_left(T, t)
        
    if pos < interp0 :
        raise ValueError("retard trop long")
        
    if pos > len(T)-2 :
        raise ValueError("retard négatif")
    
    # Construit le polynôme
    
    P = X[:, pos-interp0:pos+interp1].dot(vdmInv)
                
    Y =np.zeros_like(X[:, 0])

    for j in range(P.shape[0]) :
        for i in range(P.shape[1]) :
            Y[j] += P[j, i]*(t-T[pos])**i
        
#    print(X[:, pos], Y)    
                
    return(Y)

def rk4(f, t0, tN, N, X0) :
    
    T = np.linspace(t0, tN, N)
    h = (tN-t0)/N
    
    X = np.zeros((len(T), len(X0)))
    
    X[0] = X0
    
    for i, t in enumerate(T[1:]) :
    
        k1 = f(t, X[i])
        k2 = f(t+h/2, X[i]+h/2*k1)
        k3 = f(t+h/2, X[i]+h/2*k2)
        k4 = f(t+h, X[i]+h*k3)
 
        X[i+1] = X[i]+h/6*(k1+2*k2+2*k3+k4)    
        
    return(T, X)
    
def solve(f, tau, t0, tN, h, X0) :
    """ Résoud l'équation diff X'=f(t, X, X(t-tau(X))) 
    entre t0 et tN par pas constant h avec initialisation X0"""
    
    T = np.arange(t0-h*X0.shape[1], tN, h)
       
    X = np.zeros((X0.shape[0], len(T)))
        
    X[:, :X0.shape[1]] = X0
    
    for i in range(X0.shape[1]-1, len(T)-1) :
        
        t = T[i]
        
        k1 = f(t, X[:, i], interp(T, X, t-tau(X[:, i])))
        k2 = f(t+h/2, X[:, i]+h/2*k1, interp(T, X, t+h/2-tau(X[:, i]+h/2*k1)))
        k3 = f(t+h/2, X[:, i]+h/2*k2, interp(T, X, t+h/2-tau(X[:, i]+h/2*k2)))
        k4 = f(t+h, X[:, i]+h*k3, interp(T, X, t+h-tau(X[:, i]+h*k3)))
 
        X[:, i+1] = X[:, i]+h/6*(k1+2*k2+2*k3+k4)    
        
    return(T, X)

def fun(t, x) :
    return(np.array([x[1], -10*x[0]]))
    
def examplerk4() :    
    
    t, X = rk4(fun, 0, 1, 100, np.array([1, 0]))
    plt.plot(t, X[:, 0])
    
    plt.plot(t, np.cos(np.sqrt(10)*t))
    

    
def exampledde1() :    
    
    tau0 = 3
    def taudde(x) :
        return(tau0)
        
    def fundde(t, x, xtau) :
        return(-np.pi/2*xtau/tau0)

    
    h = 1e-2
    
    tmin = np.arange(-taudde(0)-2*h, 0, h)
    X0 = np.zeros((1, len(tmin)))
    X0[0, :] = np.cos(np.pi/2*tmin/tau0)
    
    plt.plot(tmin, X0[0, :])
    
    t, X = solve(fundde, taudde, 0., 30., h, X0)
    plt.plot(t[X0.shape[1]:], X[0, X0.shape[1]:])
    
    
def exampledde2() :    
        
    tau0 = 3
    def taudde(x) :
        return(-x[0]+tau0)
        
    def fundde(t, x, xtau) :
        return(np.array([-np.pi/2*xtau[0]/tau0, 3*np.pi/2*xtau[1]/tau0]))
    
    
    h = 1e-3
    
    tmin = np.arange(-taudde([0])-(interp0+1)*h, 0, h)
    X0 = np.zeros((2, len(tmin)))
    X0[0, :] = np.cos(np.pi/2*tmin/tau0)
    X0[1, :] = np.cos(3*np.pi/2*tmin/tau0)
    
    
    plt.plot(tmin, X0[0, :])
    plt.plot(tmin, X0[1, :])
    
    
    t, X = solve(fundde, taudde, 0., 10., h, X0)
    plt.plot(t[X0.shape[1]:], X[0, X0.shape[1]:])    
    plt.plot(t[X0.shape[1]:], X[1, X0.shape[1]:])    
#    :
#    plt.plot(t, np.cos(np.sqrt(10)*t))

exampledde2()
#exampledde2()

#
#T = np.linspace(0, 1, 7)
#X = np.random.random(len(T))
#
#
#t = np.linspace(0, 1, 100)
#
#plt.plot(T, X)
#
#interpTest(T, X, t)
