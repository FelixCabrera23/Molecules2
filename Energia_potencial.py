#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 01:34:26 2021

@author: walberto

Estudio de la energia potencial del sistema
para el caso de dos particulas

"""
import matplotlib.pyplot as plt
import numpy as np

from Modulo_funciones import*


U = []
d = []
F = []

y = [0,0]
set_start_cond(y,'y')

i = 0.5
while i < 1.5:
    

    x = [-i,i]
    set_start_cond(x,'x')
    
    Ui = Un(2)
    U.append(Ui)
    
    Fn = F_LJ(0,1)[0]
    F.append(Fn)
    
    d.append(i*2)
    
    i = i +0.01
    
# Generamos la grafica
fig3 = plt.figure()
ax3 = fig3.add_subplot()

l4 = ax3.plot(d,U,label='V')

ax3.legend()
ax3.set(xlabel=r'Separación $ [m]$',
        ylabel=r'Energía potencial [joules]')
#plt.savefig('Energias.pdf')
plt.show()   



# Grafica de Fuerza vz Distancia
plt.clf()

fig4 = plt.figure()
ax4 = fig4.add_subplot()

l5 = ax4.plot(d,F,label='F_LJ')

ax4.legend()
ax4.set(xlabel=r'Separación $ [m]$',
        ylabel=r'Fuerza [N]')
#plt.savefig('Energias.pdf')
plt.show()   