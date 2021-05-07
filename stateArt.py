import dde
import yin
import numpy as np
import matplotlib.pyplot as plt

K = 6
omega1 = 2511

def stateArtModel(TMax, h, params, alg = 'rk4Neutral') :    
    
    t0 = np.arange(-1, h, h)
    X0 = np.zeros((len(t0), 2*K+2))  
    
    label = 'h= ' + str(h) + ', ' + alg #str(h) + ', ' + str(params) + 
    
    for i in range(2*K+2) :
        
        if i % 2 == 0 :
            X0[:, i] = 1e-10*np.cos(2*np.pi*t0*omega1/3.5)
        else :
            X0[:, i] = -1e-10*np.sin(2*np.pi*t0*omega1/3.5)
       
    t0 *= omega1
    T = TMax*omega1
        
    t, X = dde.rk4Delay(t0, X0, T, 'stateArt.c', params, alg = alg) 
    
    p = np.sum(X[:, ::2], axis = 1)

    plt.figure("Wave")
    plt.plot(t/omega1, p, label=label)
    plt.legend()
    
    # Calcul les fréquences
    w_len=1024*4
    w_step=256//2
    
    f0_min=200
    f0_max=2500
    harmo_thresh=0.6
   
    sr = int(omega1/(t[1]-t[0]))
        
    pitches, harmonic_rates, argmins, times = yin.compute_yin(p[len(t0):], sr, None, w_len, w_step, f0_min, f0_max, harmo_thresh)
    pitches = np.array(pitches)
    times = np.array(times)
    
    Pm = params[0]+(TMax/2-abs(times-TMax/2))/TMax*(params[1]-params[0])
    
    plt.figure("Yin")
    ax = plt.gca()
    
    color=next(ax._get_lines.prop_cycler)['color']
    plt.plot(Pm[:len(Pm)//2], np.array(pitches[:len(Pm)//2]), '-', c = color, label = label)
    plt.plot(Pm[len(Pm)//2:], np.array(pitches[len(Pm)//2:]), ':', c = color)
    
    plt.legend()

    
#    from scipy.io.wavfile import write
#    scaled = np.int16(p/np.max(np.abs(p)) * 32767)
#    write('test'+'.'.join(str(e) for e in params)+'.wav', int(sr), scaled)
    
    
    
stateArtModel(20, 1/(1*44100), params=[30, 500, 20], alg = 'rk4Neutral')
stateArtModel(20, 1/(2*44100), params=[30, 500, 20], alg = 'rk4Neutral') 
stateArtModel(20, 1/(8*44100), params=[30, 500, 20], alg = 'rk4Neutral') 
stateArtModel(20, 1/(16*44100), params=[30, 500, 20], alg = 'rk4Neutral') 
#stateArtModel(20, 1/(32*44100), params=[30, 500, 20], alg = 'rk4Neutral') 


#stateArtModel(15, 1/(20*44100), params=[70, 400, 15], alg = 'rk4') 
#stateArtModel(30, 1/(44100), params=[50, 600, 30], alg = 'impTr') 
#stateArtModel(10, 1/(4*44100), params=[50, 600, 10], alg = 'rk4') 
#stateArtModel(10, 1/(4*44100), params=[50, 600, 10], alg = 'impTr') 
#stateArtModel(10, 1/(4*44100), params=[600, 10], alg = 'impTr') 

plt.show()
#toyModel(10, 1e-2)   

plt.show()    
