# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 19:19:22 2016

@author: Diego Chou
"""

import numpy as np
from time import clock
import os
from bisect import bisect_right, bisect_left

def tempoDeExecucao(func):
    """
    Calcula o tempo de execução da função especificada.
    """
    def nested(*args, **kwargs):
        def mostrador(string, tempo):        
            h = tempo//3600
            m = (tempo%3600)//60
            s = (tempo%3600)%60
            string_s = ' s'
            string_m = ' min '
            string_h = ' h '
            if h == 0:
                h = string_h = ''
            if m == 0:
                m = string_m = ''
            print("\n---> Tempo de %s : %s%s%s%s%f%s " %(string,
                                                         h,string_h,
                                                         m,string_m,
                                                         s,string_s))
        START_EXEC_TIME = clock()
        result = func(*args, **kwargs) #Função avaliada
        EXEC_TIME = clock()-START_EXEC_TIME
        mostrador('execucao',EXEC_TIME)
        
        return result
    return nested

def executarNaPasta(string):
   """
   Muda o diretório de execução para o endereço especifidado pela string
   """
   def nested(func):
       def wrapper(*args, **kwargs):
            _cwd = os.getcwd()
            try:
                os.chdir(string)
            except OSError:
                # O makedirs cria diretórios recursivamente
                # mkdir('dir1/dir2')    -> ERRO
                # makedirs('dir1/dir2') -> OK
                os.makedirs(string)
                os.chdir(string)
            try:    
                result = func(*args, **kwargs)
            except:
                os.chdir(_cwd)
                raise
            os.chdir(_cwd)
            return result 
       return wrapper
   return nested
   
def interplinear(x1,y1,x2,y2,x):   
    if x1<=x<=x2 or x2<=x<=x1:
        try:
            y = ((x-x1)*y2+(x2-x)*y1)/(x2-x1)
        except ZeroDivisionError:
            y = y1
    else:
        raise ValueError("O ponto %f está fora do intervalo: [%f,%f]" %(x,min(x1,x2),max(x1,x2)))
    return y 

def intervalo(lista, valor):
    lista = sorted(lista)
    if valor < lista[0] or valor > lista[-1]:
        raise ValueError('Input %f fora do intervalo de valores: %s' %(valor, [lista[0],lista[-1]]))
        
    def achar_maior_antes():
        """
        Acha o maior valor da lista menor ou igual a x
        """
        i = bisect_right(lista,valor)
        if lista[i-1] == lista[-1]:
            return i-2
        else:
            return i-1
        
    def achar_menor_depois():
        """
        Acha o menor valor da lista maior ou igual a x
        """
        i = bisect_left(lista,valor)
        if lista[i] == lista[0]:
            return i+1
        elif (valor in lista) and (valor != lista[-1]):
            return i+1
        else:            
            return i
    
    index_antes = achar_maior_antes()
    index_depois = achar_menor_depois()
            
    return index_antes, index_depois

def nan2zero(val):
    if np.isnan(val):
        return 0
    else:
        return val
    
def reset_dados_perfis():
    mypath = 'Perfis/'
    files = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
    
    for f in files:
        if '__' in f or 'XFoil' in f:
            os.remove(os.path.join(mypath, f))