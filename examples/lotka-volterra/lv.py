import dde
import numpy as np
import matplotlib.pyplot as plt

def clarinet(t0, TMax, h) :    
    
    
    X0 = np.array([0.4, 0.2])
        
    t, X = dde.rk4(t0, X0, h, TMax, 'lv.c', [2/3, 4/3, 1]) 
            
    plt.plot(X[:, 0], X[:, 1], label="rk4")

clarinet(0., 10., 1e-1) 
#toyModel(10, 1e-2)   

plt.show()    
