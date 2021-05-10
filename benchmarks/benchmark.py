import dde
import yin
import numpy as np
import matplotlib.pyplot as plt

from time import process_time

""" Définition de stateArtModel """

K = 6
omega1 = 2511

def stateArtModel(TMax, h, params, alg = 'rk4Neutral', tref = [], pref = []) :  
    
    h /= 44100
    
    t0 = np.arange(-1, h, h)
    X0 = np.zeros((len(t0), 2*K+2))  
        
    for i in range(2*K+2) :
        
        if i % 2 == 0 :
            X0[:, i] = 1e-10*np.cos(2*np.pi*t0*omega1/3.5)
        else :
            X0[:, i] = -1e-10*np.sin(2*np.pi*t0*omega1/3.5)
       
    t0 *= omega1
    T = TMax*omega1
    
    start = process_time() 
            
    t, X = dde.rk4Delay(t0, X0, T, '../examples/stateArt/stateArt.c', params, alg = alg) 
    
    stop = process_time()
    
    p = np.sum(X[:, ::2], axis = 1)
    
    if (len(tref) > 0) :
        
        
        interp  = np.abs(p-np.interp(t, tref, pref))/max(abs(pref))     
        
        return([omega1*h, max(interp), stop-start])
    
    else :
        return(t, p)
    
""" Test de bellen """    
def bellenModel(TMax, h, params, alg = 'rk4Neutral', tref = [], pref = []) :
    
    h /= 10

    t0 = np.arange(-1., 0, h)
    X0 = np.zeros((len(t0), 1))  
    
    X0[:, 0] = 1-t0
      
    
    start = process_time() 
            
    t, p = dde.rk4Delay(t0, X0, TMax, '../examples/bellen/bellen.c', params, alg = alg) 
    
    stop = process_time()
    
    p = p.reshape((len(p)))
        
    if (len(tref) > 0) :
        
        interp  = np.abs(p-np.interp(t, tref, pref))/max(abs(pref))

        return([h, max(interp), stop-start])
    
    else :
        return(t, p)

def testModel(func = stateArtModel, H = []) :
        
    output = open('benchmark_'+func.__name__+'_.txt', 'w')
    
    tref, pref = func(1, H[-1], params=[0, 500, 10], alg = 'rk4Neutral') 
    
    for alg in ['rk4Neutral', 'impTrNeutral', 'eulerImpNeutral'] :
    
        output.write(alg+'\n')
        
        r1 = []
        
        for h in H[:-1] :
            
            r1.append(func(1, h, params=[0, 500, 10], alg = alg, tref = tref, pref = pref))
            
        output.write(str(r1)+'\n')
        
    
    output.close()
    


# testModel(stateArtModel, H=1/(2**np.arange(4)))
testModel(bellenModel, H=1/(2**np.arange(16)))