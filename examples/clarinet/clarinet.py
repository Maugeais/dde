import dde
import numpy as np
import matplotlib.pyplot as plt
import yin

def clarinet(t0, TMax, h) :    
    
    
    X0 = np.array([0.01, 0, 0, 0])
        
    t, X = dde.rk4(t0, X0, h, TMax, 'clarinet.c', [1.3, -0.2, 0.1]) 
            
    #plt.plot(t, X[:, 0]+X[:, 2], label="rk4")
    
    # Calcul les fréquences
    w_len=512*16
    w_step=256//2
    
    f0_min=200
    f0_max=2500
    harmo_thresh=0.1
   
    sr = int(44100)
    p = X[:, 0] + X[:, 2]
    plt.plot(t, p)
        
    pitches, harmonic_rates, argmins, times = yin.compute_yin(p, sr, None, w_len, w_step, f0_min, f0_max, harmo_thresh)

    pitches = np.array(pitches)
    times = np.array(times)
    
    Pm = 0.35+times*0.2
    
    plt.figure("Yin")
    ax = plt.gca()
    
    color=next(ax._get_lines.prop_cycler)['color']
    plt.plot(Pm, np.array(pitches), '-')
    
    plt.legend()
    

clarinet(0., 5, 1/44100)
#toyModel(10, 1e-2)   

plt.show()    
