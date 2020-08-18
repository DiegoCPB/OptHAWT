#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 13:06:31 2018

@author: diego
"""
import otimizacao as ot
import graficos

lamb = 7.5 # Razão de velocidades fixada
r = 60 #Raio do rotor Siemens
v = 5 #Velocidade de máxima eficiência
eff = 0.77 #Eficiência da turbina
pot = 310e3/eff #Potência requerida
c_raiz = 4.2 #Corda na raiz da pá

w = lamb*v/60 #rad/s

obj = ot.Otimizador(pot,v,w,2.5,3,
                    ['circular','NACA 63630','FW W3 211','FW W3 211',
                     'FW W1 152','FW W1 152'], c_circular = c_raiz)
pa = obj.pa_opt

graficos.potência(pa,v_bound = [3,25], maxP=3.6e6,fixa=False, eff = eff )
graficos.geometria(pa)
graficos.planta(pa)
graficos.desempenho(pa)