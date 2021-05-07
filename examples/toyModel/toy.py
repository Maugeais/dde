import dde
import yin
import numpy as np
import matplotlib.pyplot as plt

def toyModel(TMax, h) :    
    
    omega1 = 2764.
    
    t0 = np.arange(-5, 0, h)
    X0 = np.zeros((len(t0), 4))  
    
    X0[:, 0] = 1e-1+0*np.cos(100*np.pi*t0)
    X0[:, 1] = 0*np.cos(3*np.pi/2*t0)
    
    T = TMax*omega1
        
    t, X = dde.rk4Delay(t0, X0, T, 'toy.c', [30, T]) 
    
    p = X[:, 0]+X[:, 2]
    
    #plt.plot(t, p)

    # Calcul les fréquences
    w_len=1024*4
    w_step=256*4
    f0_min=0.05
    f0_max=3
    harmo_thresh=0.85
    
    sr = int(1/(t[1]-t[0]))
    pitches, harmonic_rates, argmins, times = yin.compute_yin(p[len(t0):], sr, None, w_len, w_step, f0_min, f0_max, harmo_thresh)
    pitches = np.array(pitches)
    times = np.array(times)
    
    # Changements de signes
    """ Calcul bourrin
    I = np.where(p[1:]*p[:-1] < 0)[0]
    F = t[I]
    
    treal = 0.1+(T/2-abs(t[I[2:]]-T/2))*30/T
    Freal = omega1/(F[2:]-F[:-2])
    plt.plot(treal[:len(treal)//2], Freal[:len(treal)//2])
    plt.plot(treal[len(treal)//2:], Freal[len(Freal)//2:])"""
    
    tau = 0.1+(T/2-abs(times-T/2))*30/T
    #print(times[-1], t[-1])
    plt.plot(tau[:len(tau)//2], omega1*np.array(pitches[:len(tau)//2]))
    #plt.plot(tau[len(tau)//2:], omega1*np.array(pitches[len(tau)//2:]))

    
    from scipy.io.wavfile import write
    scaled = np.int16(p/np.max(np.abs(p)) * 32767)
    write('test.wav', int(sr*omega1), scaled)
    
    
toyModel(10, 5e-3) 
#toyModel(10, 1e-2)   

plt.show()    
