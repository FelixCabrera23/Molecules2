#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 17:18:42 2021

@author: walberto

Este programa resuelve el problema de N cuerpos por metodo de leap frog

METODO DE MALLAS: Se generará una malla que divida el espacio en chunks discretos donde las particuas vivan.
La fuerza solo se calculara dentro de una cierta cantidad de chuncs de distancia de cada particula. 

Se intentara incorporar un metodo orientado a objetos para que a futuro permita la paralelización del problema
"""

# Lybrays

import matplotlib.pyplot as plt
import numpy as np
import random
from time import time
import sys











