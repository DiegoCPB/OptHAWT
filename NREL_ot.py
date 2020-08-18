#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 13:06:31 2018

@author: diego
"""
import otimizacao as ot
import graficos

lamb = 7.69 # Razão de velocidades fixada
r = 61.63 #Raio do rotor NREL
v = 8 #Velocidade de máxima eficiência
pot = 1.796e6 #Potência requerida
c_raiz = 3.542 #Corda na raiz da pá

w = lamb*v/r #rad/s

obj = ot.Otimizador(pot,v,w,3,3,
                    ['circular','DU 93W405 AD','DU 93W350 AD',
                     'DU 93W300 LM','DU 93W250 LM','DU 93W210 LM',
                     'NACA 64618','NACA 64618'], c_circular = c_raiz)
pa = obj.pa_opt

graficos.geometria(pa)
graficos.planta(pa)
graficos.desempenho(pa,lambda_bound=[5,10])

#import macro_CATIA as mc
#mc.exportarPerfis(pa,"NREL_5MW_ot")
