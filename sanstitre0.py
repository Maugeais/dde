import numpy as np
import matplotlib.pyplot as plt
import yin


t = np.arange(0, 5, 1/44100)

sig = np.cos((220+100*t)*2*np.pi*t)
sr = 44100

w_len=1024
w_step=256
f0_min=200
f0_max=1000
harmo_thresh=0.85

pitches, harmonic_rates, argmins, times = yin.compute_yin(sig, sr, None, w_len, w_step, f0_min, f0_max, harmo_thresh)


plt.plot(times, pitches)