#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 15:28:36 2021

@author: maugeais
"""

import matplotlib.pyplot as plt
import ast

f = open('benchmark_bellenModel_.txt', 'r')

content = f.read().split('\n')

for i, alg in enumerate(content[:-1:2]) :
    
    
    # print(content[2*i+1])
    data = ast.literal_eval(content[2*i+1])
    
    H = [float(d[0]) for d in data]
    eps = [float(d[1]) for d in data]
    T = [float(d[2]) for d in data]
    
    label = content[2*i] 
    
    plt.figure('error')
    plt.loglog(H, eps, label = label)
    plt.legend()
    
    plt.figure('computation time')
    plt.loglog(H, T, label = label)
    plt.legend()
    