#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 17:30:35 2021

@author: walberto

Modulo de funciónes 
Necesario para el proyecto de n-cuerpos chunk 


"""
import matplotlib.pyplot as plt
import numpy as np
import random
import sys
from time import time


##### VARIABLES GLOBALES Y PREDEFINIDAS ######

L = 50.0 # tamaño de la caja
e = 1.0  # parametro epsilon
s = 1.0 # parametro sigma
n = 10 # Cantidad de particulas
t = 0.0001 # tiempo por paso
v0 = 1.0 # Rango de las velocidades
seed = 1999 # Semilla



# Empezamos definiendo la clase particula,
# Este ente sera el que vivira y se vera afectado por los campos

class Particula(object):
    """
    Particulas que viven en los chunks, estas se ven atraidas y repelidas
    entre ellas segun la fuerza definida por el usuario.
    Poseen atributos como carga, masa, posición y momentum
    """
    
    def __init__(self, xy=[0,0],vxy=[0,0],r=0.5,m=1,chnk=0, n=0):
        """
        Propiedades intrinsecas de la particula

        """
        self.xy = xy        # Posición carteciana
        self.vxy = vxy      # Velocidad vector
        self.r = r          # Radio
        self.m = m          # Masa
        self.chnk = chnk    # Chunk
        self.n = n          # Numero
        
    
        
    def __repr__(self):
        return('Particula({0.xy!r},{0.vxy!r},{0.r!r},{0.m!r})'.format(self))
        
        
        
        
        
        
        
        
        
        
        
        