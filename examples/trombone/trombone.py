import dde
import numpy as np
import matplotlib.pyplot as plt
import yin

def trombone(t0, TMax, h) :    
    
    N = 5
    X0 = 1e-5*np.ones(5*N+2)
    X0[0] = 0*7.7863e-04
    
    a = 8000
    b = 0
        
    t, X = dde.rk4(t0, X0, h, TMax, 'trombone.c', [a, b, 30]) 
    
    p = 2*np.sum(X[:, 2::2], 1)
    plt.plot(t, p, label="pressure")

    #return(t, p)        
    #plt.plot(t, X[:, 0], label="h")

    # Calcul les fréquences
    w_len=1024*4
    w_step=256
    
    f0_min=40
    f0_max=2500
    harmo_thresh=0.9
   
    sr = int(44100)
 #   p = X[:, 0] + X[:, 2]
 #   plt.plot(t, p)
 #       
    pitches, harmonic_rates, argmins, times = yin.compute_yin(p, sr, None, w_len, w_step, f0_min, f0_max, harmo_thresh)

    pitches = np.array(pitches)
    times = np.array(times)
    
#    Pm = 0.35+times*0.2
#    
    plt.figure("Yin")
    ax = plt.gca()
    
    color=next(ax._get_lines.prop_cycler)['color']
    plt.plot(times, np.array(pitches), '-')
    harm = np.argmin(harmonic_rates)
    plt.figure("test")
    #plt.plot(p[harmonic_rates[harm]:harmonic_rates[harm+1]], '-')
    i = np.where(t >= times[harm])[0][0] #, t <= times[harm]))
    q = p[i: int(i+sr//pitches[harm])]
    j = 1+np.where(np.logical_and(q[1:]*q[:-1] < 0, q[1:] >= 0))[0][0]
    plt.plot(p[i+j: int(i+j+sr//pitches[harm])])

    

    
        
from scipy import signal   

def fltr(p, fs = 44100.0) :

    cutoff = 1000.0    # Desired cutoff frequency, Hz

    trans_width = 40  # Width of transition from pass band to stop band, Hz

    numtaps = 125    # Size of the FIR filter.

    taps = signal.remez(numtaps, [0, cutoff - trans_width, cutoff, 0.5*fs], [0, 1], Hz=fs)

    w, h = signal.freqz(taps, [1], worN=2000)
    
    #plt.plot(0.5*fs*w/np.pi, 20*np.log10(np.abs(h)))
    #plt.show()
    
    print(taps)
    
    pp = signal.lfilter(taps, 1, p)
    
    print(p.shape, pp.shape)
    
    return(np.real(pp))

    #plot_response(fs, w, h, "Band-pass Filter")
    
trombone(0., 1, 1/44100)

#plt.plot(t, abs(np.fft.fft(fltr(np.copy(p)))))
#fltr(p)
#toyModel(10, 1e-2)   

plt.show()    
