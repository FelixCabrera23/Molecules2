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

"Constantes importantes"

L = 50 # tamaño de la caja
e = 1  # parametro epsilon
s = 1 # parametro sigma
n = 100 #  Cantidad de particulas
m = 1 # Masa de las particulas
t = 0.0001 # tiempo por paso
R = 0.5 # Radio de las particulas
Ek = []  # Energia del sistema
V = []  # Esta es la energia potencial del sistema

"Las particulas estan definidas por su posicion y velocidad"
x = []
y = []
vx = []
vy = []

# Constantes necesarias para el programa

l = L/2 # para poner el centro de la caja en el origen
random.seed(1999) # Semilla para los procesos aleatorios
v0 = 1 # Rango de las velocidades

# Archivo donde vamos a ir escribiendo los resultados

posiciones = open('posisiones.dat','w')
velocidades = open('velocidades.dat','w')

# Definimos la funcion de fuerza
def F_LJ(pA,pB):
    """
    Esta función calcula la fuerza entre un par de particulas
    
    Devuelve la fuerza en componentes x y
    """
    x0 = x[pA]-x[pB]
    y0 = y[pA]-y[pB]
    r = np.sqrt(x0**2+y0**2)
    
    # Esta parte nos asegura que las particulas son solidas
    # Consideramos que la fuerza debido al potencial L-J es mayor que la
    # transferencia de momentum linal debido a las colisiones
    if r<(2*R):
        r = 2*R
    
    F = -(24*e*s**6*(r**6-2*s**6))/r**13
    
    Fx = F*(x0/r)
    Fy = F*(y0/r)
    
    return(Fx,Fy)
    
def Reb(p):
    """
    Esta función chequea si la particula rebota contra las paredes
    """
    #Para x
    x0 = abs(x[p])
    xs =x[p]/x0 # signo de la coordenada
    
    if x0 > l-R:
        dx = x0 - l+R
        xn = x0 - 2*dx
        
        x[p] = xn*xs # Esto cambia la coordenada reflejada
        
        vx[p] = vx[p]*(-1) # Esto refleja la velocidad
        
    #Para y
    y0 = abs(y[p])
    ys = y[p]/y0
    
    if y0 > l-R:
        dy = y0 - l+R
        yn = y0 -2*dy
        
        y[p] = yn*ys
        
        vy[p] = vy[p]*(-1)
        
    return()

def Start(n):
    """
    Esta función empiesa el sistema con posiciónes y velocidades 
    aletorias que dependen de un parametro v0 y r0 

    """
    for i in range(n):
        sig1 = (-1)**(random.randint(0,1))
        sig2 = (-1)**(random.randint(0,1))
        ang = random.random()*2*np.pi
        xi = sig1*random.random()*(l-R)
        yi = sig2*random.random()*(l-R)
        
        vi = v0*random.random()
        
        vxi = vi*np.cos(ang)
        vyi = vi*np.sin(ang)
        
        x.append(xi)
        y.append(yi)
        
        vx.append(vxi)
        vy.append(vyi)
        
    return()

def Graf():
    
    """
    Esta función grafica el sistema en tiempo real
    
    """
    
    plt.figure()
    plt.axis([-l,l,-l,l])
    ax = plt.gca()
    
    # plt.scatter(x, y)
    for i in range(len(x)):
        ax.add_patch( plt.Circle((x[i],y[i]),0.5))
    plt.show
    
    
def En(n):
    """
    Esta función calcula la energia cinetica del sistema

    """
    E = 0
    for i in range(n):
        E = E + (0.5*m*(vx[i]**2+vy[i]**2))
    # Ek.append(E)
    return(E)

def E_LJ(pA,pB):
    """
    Esta función  determina la energia entre dos particulas

    """
    x0 = x[pA]-x[pB]
    y0 = y[pA]-y[pB]
    r = np.sqrt(x0**2+y0**2)
    
    V = 4*e*((s/r)**12 - (s/r)**6)
    
    return(V)
    
def Un(n):
    U = 0
    
    for i in range(n):
        Ui = 0
        for j in range(n):
            if i == j:
                continue
            else:
                Ui = Ui + E_LJ(i,j)
                
        U = U + Ui

    Uf = U/2
    return(Uf)    

def Mov():
    """
   Esta función mueve un paso al sistema

    """
        
    for i in range(n):
        #####################################################################
        # Antes de mover guardamos la posicion de la particula en un archivo
        posiciones.write('%f      %f      ' % (x[i],y[i]))
        velocidades.write('%f      %f      ' %(vx[i],vy[i]))
        #####################################################################
        
        # Calculamos la sumatoria de fuerzas
        Fxn,Fyn = 0,0
        
        for j in range(n):
            if i==j:
                continue

            F = F_LJ(i,j)
            Fxn = Fxn + F[0]
            Fyn = Fyn + F[1]
            
         # Movemos la particula i
            
        x[i] = x[i] + vx[i]*t + (Fxn/m)*0.5*t**2
        y[i] = y[i] + vy[i]*t + (Fyn/m)*0.5*t**2
         
        vx[i] = vx[i] + (Fxn/m)*t
        vy[i] = vy[i] + (Fyn/m)*t
         
        # Checamos por colisiones
        Reb(i)
         
    return()

####### Main

n = 60
Start(n)

# x = [-1,1]
# y = [0,0]
# vx = [0,0]
# vy = [0,0]

pasos = 500000 # Cantidad de pasos a simular
datos = 1000 # Cantidad de pasos por grafica a generar

Sum = []
tiempos = []
tf = 0

for k in range(pasos):
    # Esta parte hace graficas en tiempo real
    
    if (np.mod(k,datos)==0):
        # Graf()
        E = En(n)
        U = Un(n)
        L = E+U
        Ek.append(E)
        V.append(U)
        Sum.append(L)
        tiempos.append(tf)
        tf = tf + (datos*t)
    # Esta parte realiza el movimiento    
    Mov()
    
    posiciones.write('\n') #Esto guarda un espacio entre cada movimiento
    velocidades.write('\n')
        
posiciones.close()
velocidades.close()        


# Esto plotea la energia promedio de las particulas
plt.figure()    
plt.plot(tiempos,Ek)
plt.plot(tiempos,V)
plt.plot(tiempos,Sum)
plt.show()    
    
    
        