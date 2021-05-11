import dde
import numpy as np
import matplotlib.pyplot as plt

def bellen(TMax, h) :    
    
    
    t0 = np.arange(-1., 0, h)
    X0 = np.zeros((len(t0), 1))  
    
    X0[:, 0] = 1-t0
    
        
    t, X = dde.rk4Delay(t0, X0, TMax, 'ode.c', [], alg = 'rk4Neutral') 
        
    plt.plot(t, X, label="rk4")
    
    t, X = dde.rk4Delay(t0, X0, TMax, 'ode.c', [], alg = 'eulerNeutral') 
        
    plt.plot(t, X, label='euler')
    
    t, X = dde.rk4Delay(t0, X0, TMax, 'ode.c', [], alg = 'eulerImpNeutral') 
        
    plt.plot(t, X, label='eulerImp')
    
    t, X = dde.rk4Delay(t0, X0, TMax, 'ode.c', [], alg = 'impTrNeutral') 
        
    plt.plot(t, X, label='impTrNeutral')
    
    plt.legend()
    
bellen(5., 1e-1) 
#toyModel(10, 1e-2)   

plt.show()    
