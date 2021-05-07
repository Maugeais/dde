import dde
import numpy as np
import matplotlib.pyplot as plt

def bellen(TMax, h) :    
    
    
    t0 = np.arange(-1., 0, h)
    X0 = np.zeros((len(t0), 1))  
    
    X0[:, 0] = 1-t0
    
        
    t, X = dde.rk4Delay(t0, X0, TMax, 'bellen.c', [], alg = 'rk4Neutral') 
        
    plt.plot(t, X)
    
bellen(5., 1e-2) 
#toyModel(10, 1e-2)   

plt.show()    
