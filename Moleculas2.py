# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 11:58:39 2020

@author: User

Moleculas 2

Simualacion de particulas coloidales bajo el potencial de Lennard-Jones
utilizando el metodo de Leap Frog
"""

import matplotlib.pyplot as plt
import numpy as np
import random
import sys
from time import time

#Importanto modulo de funciones propias

from Modulo_funciones import*

"Constantes importantes"
Set_config(
    120.0, # tamaño de la caja
    1.0,  # parametro epsilon
    1.0, # parametro sigma
    30, #  Cantidad de particulas
    1.0, # Masa de las particulas
    0.01, # tiempo por paso
    0.5, # Radio de las particulas
    1999, # Semilla
    1.0 # Rango de las velocidades
    )

Ek = []  # Energia del sistema
V = []  # Esta es la energia potencial del sistema
"Las particulas estan definidas por su posicion y velocidad"
x = []
y = []
vx = []
vy = []


# Archivo donde vamos a ir escribiendo los resultados

posiciones = open('posisiones.dat','w')
velocidades = open('velocidades.dat','w')
energias = open('energias.dat','w')

####### Main


# Generamos el sistema inicial
# Start(Get_data('n'))

# Si queremos sobreescribir el sistema generado automaticamente:
# Esto pone coordenadas especificas
x = [-10,10]
set_start_cond(x,'x')
y = [0,0]
set_start_cond(y,'y')
vx = [0,0]
set_start_cond(vx,'vx')
vy = [0,0]
set_start_cond(vy,'vy')

# Conseguimos datos del programa principal
t = Get_data('t')
n = Get_data('n')

pasos = 400000 # Cantidad de pasos a simular
datos = 10000 # Cantidad de pasos por grafica a generar

Sum = []
tiempos = []
tf = 0

for k in range(pasos):
    # Esta parte hace graficas en tiempo real
    
    if (np.mod(k,datos)==0):
        
        x = Get_data('x')
        y = Get_data('y')
        vx = Get_data('vx')
        vy = Get_data('vy')
        # Genera la grafica
        Graf()
        
        # Guarda los datos de las energias
        E = En(n)
        U = Un(n)
        L = E+U
        Ek.append(E)
        V.append(U)
        Sum.append(L)
        tiempos.append(tf)
        
        tf = tf + (datos*t)
        
        # Escribiendo datos en archivo
        # Escribimos las posiciones y velocidades individuales
        for i in range(n):
            posiciones.write('%f    %f    ' %(x[i],y[i]))
            velocidades.write('%f    %f    ' %(vx[i],vy[i]))
        
        posiciones.write('\n') #Esto guarda un espacio entre cada movimiento
        velocidades.write('\n')
        
        energias.write('%f    %f    %f    %f \n' % (tf,E,U,L))
        
    # Esta parte realiza el movimiento    
    Mov()
    

        
posiciones.close()
velocidades.close()        
energias.close()

 
    
    
        