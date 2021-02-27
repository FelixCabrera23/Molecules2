#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 22:33:02 2021

@author: walberto

Funciones de graficas para el proyecto Moleculas2
"""

import matplotlib.pyplot as plt
import numpy as np

from Modulo_funciones import*

# Grafica de las energias:
            
def Graf_En():
    """
    Esto plotea la energia promedio de las particulas
    """
    # Variables necesarias
    t = []
    Ek = []
    V = []
    L = []
    
    tau = Get_data('t')
    
    # Leemos los datos del archivo
    with open('energias.dat') as energias:
        for line in energias:
            dat = line.split()
            t.append(float(dat[0]))
            Ek.append(float(dat[1]))
            V.append(float(dat[2]))
            L.append(float(dat[3]))
    
    # Generamos la grafica
    fig1 = plt.figure()
    ax1 = fig1.add_subplot()
    
    l1 = ax1.plot(t,Ek,label='T')
    l2 = ax1.plot(t,V,label='V')
    l3 = ax1.plot(t,L,label='H')
    
    ax1.legend(loc='upper right')
    ax1.set(xlabel=r'$t \times %f [s]$' % tau,
            ylabel=r'Energía [joules]')
    #plt.savefig('Energias.pdf')
    plt.show()   
