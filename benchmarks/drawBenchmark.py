#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 15:28:36 2021

@author: maugeais
"""

import matplotlib.pyplot as plt
import ast

f = open('benchmark_stateArtModel_.txt', 'r')
# f = open('benchmark_bellenModel_.txt', 'r')

content = f.read().replace('inf', '10')

content = content.split('\n')

color = {'rk4Neutral' : 'r',
          'eulerNeutral' : 'b',
          'eulerImpNeutral' : 'g', 
          'impTrNeutral' : 'm'
          }

linestyle = {'1' : '-', '2' : ':', '3' : '--'}

for i, alg in enumerate(content[:-1:2]) :
    
    
    data = ast.literal_eval(content[2*i+1])
    
    H = [float(d[0]) for d in data]
    eps = [float(d[1]) for d in data]
    T = [float(d[2]) for d in data]
    
    label = content[2*i] 
    
    c = label.split(' ')
    
    plt.figure('error')
    plt.loglog(H, eps, label = label, c = color[c[0]], linestyle=linestyle[c[1]])
    plt.xlabel('Step (h)')
    plt.ylabel('Relative error (max)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    
    plt.figure('computation time')
    plt.loglog(H, T, label = label, c = color[c[0]], linestyle=linestyle[c[1]])
    plt.xlabel('Step (h)')
    plt.ylabel('Computation time (s)')    
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
