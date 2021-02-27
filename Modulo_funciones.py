#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 17:36:03 2021

@author: walberto


Modulo de funciónes de Moleculas 2

"""

import matplotlib.pyplot as plt
import numpy as np
import random
import sys
from time import time

L = float # tamaño de la caja
e = float  # parametro epsilon
s = float # parametro sigma
n = int #  Cantidad de particulas
m = float # Masa de las particulas
t = float # tiempo por paso
R = float # Radio de las particulas
l = float # Lado medio de la caja
seed = int #semilla
v0 = float # velocidades

# ########## valores por defecto

L = 50.0 # tamaño de la caja
e = 1.0  # parametro epsilon
s = 1.0 # parametro sigma
n = 10.0 #  Cantidad de particulas
m = 1.0 # Masa de las particulas
t = 0.0001 # tiempo por paso
R = 0.5 # Radio de las particulas
v0 = 1.0 # Rango de las velocidades
seed = 1999 # Semilla


Ek = []  # Energia del sistema
V = []  # Esta es la energia potencial del sistema

"Las particulas estan definidas por su posicion y velocidad"
x = []
y = []
vx = []
vy = []


def Set_config(L_in,e_in,s_in,n_in,m_in,t_in,R_in,seed_in,v0_in):
    """
    Esta función adquiere los datos del problema del programa principal.

    """
    global L, e, s, n, m, t, R, seed, v0, l
    L = L_in
    e = e_in
    s = s_in
    n = n_in
    m = m_in
    t = t_in
    R = R_in
    seed = seed_in
    v0 = v0_in
    
    l = L/2.0 # para poner el centro de la caja en el origen
    random.seed(1999) # Semilla para los procesos aleatorios
    v0 = 1 # Rango de las velocidades

def Get_data(a):
    
    global L, e, s, n, m, t, R, seed, v0, l, x, y, vx, vy
    
    if (a == 'L'):
        return L
    if (a == 'e'):
        return e
    if (a == 't'):
        return t
    if (a == 'n'):
        return n
    if (a == 'x'):
        return x
    if (a == 'y'):
        return y
    if (a == 'vx'):
        return vx
    if (a == 'vy'):
        return vy
    
def set_start_cond(coord, a):
    
    global x,y,vx,vy,n
    
    if (a == 'x'):
        x = coord
    if (a == 'y'):
        y = coord
    if (a == 'vx'):
        vx = coord
    if (a == 'vy'):
        vy = coord    
    n = len(coord)
    
    


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
    if (x[p] != 0):        
        x0 = abs(x[p])
        xs =x[p]/x0 # signo de la coordenada
    
        if x0 > l-R:
            dx = x0 - l+R
            xn = x0 - 2*dx
            
            x[p] = xn*xs # Esto cambia la coordenada reflejada
            
            vx[p] = vx[p]*(-1) # Esto refleja la velocidad
        
    #Para y
    if (y[p] != 0):
        
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
    Potencial de Lenard-Jones
    Esta función  determina la energia entre dos particulas
    debida al potencial de Lenard-Jones

    """
    x0 = x[pA]-x[pB]
    y0 = y[pA]-y[pB]
    r = np.sqrt(x0**2+y0**2)
    
    V = 4*e*((s/r)**12 - (s/r)**6)
    
    return(V)
    
def Un(n):
    """
    Esta función calcula la energia potencial de todo el sistema

    """
    
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