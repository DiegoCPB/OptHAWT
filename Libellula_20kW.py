#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 13:06:31 2018

@author: diego
"""
import numpy as np
import otimizacao as ot
import graficos

w = 46.23*2*np.pi/60 #rad/s
v = 5
eff = 0.81
pot = 6.2e3/eff

obj = ot.Otimizador(pot,v,w,1,2,['circular','S808','S805A','S806A','S806A'])
pa = obj.pa_opt

graficos.geometria(pa)
graficos.planta(pa)

#Ajuste no passo das pás
#pa.betas = pa.betas+2
#
#graficos.potência(pa,v_bound = [3,20],eff = eff)
#graficos.desempenho(pa)
