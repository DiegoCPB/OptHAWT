# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 18:47:07 2019

@author: Diego Chou
"""
import perfil

import numpy as np
import openpyxl
import shutil
from scipy.signal import savgol_filter

def rot(x, y, alpha):
    alpha = np.deg2rad(alpha)
    xr = x*np.cos(alpha) - y*np.sin(alpha)
    yr = x*np.sin(alpha) + y*np.cos(alpha)
    return xr,yr

def alinhamento(pa):
    lista = pa._lista_nomes_perfis
    alignX = []
    alignY = []
    
    for string in lista:
        if '__' in string:
            p1 = perfil.PerfilInterpolado(string)
        else:
            p1 = perfil.Perfil(string)
        xc,yc = p1.centroide()
        alignX.append(xc)
        alignY.append(yc)
    
    return np.array(alignX), np.array(alignY)
        
def secoes(pa):
    output = []
    
    xc,yc = alinhamento(pa)
    
    wl = len(xc)//4
    if wl%2 == 0: wl+=1
    align_x = savgol_filter(xc,wl,3)
    align_y = savgol_filter(yc,wl,3)
    
    cordas = pa.cordas
    betas = pa.betas
    perfis = pa.coords
    x = pa.x
    
    for index,item in perfis.items():
        yp,zp = item
        xp = np.ones(len(yp))*x[index]*1000
        yp,zp = yp*cordas[index],zp*cordas[index]
        yp,zp = yp-align_x[index]*cordas[index],zp-align_y[index]*cordas[index]
        yp,zp = rot(yp,zp,betas[index]-90)
        yp,zp = yp*1000,zp*1000
        output.append([xp,yp,zp])
    
    return output

def exportarPerfis(pa, name):
    myfile = 'Catia/'+name+'.xlsm' 
    shutil.copyfile('Catia/modelo.xlsm', myfile)
    
    book = openpyxl.load_workbook(myfile,read_only=False, keep_vba=True)
    sheet = book['Feuil1']
    
    sheet.append(("StartLoft",))
    
    for aerofolio in secoes(pa):
        sheet.append(("StartCurve",))
        for i in range(len(aerofolio[0])):
            sheet.append((aerofolio[0][i],aerofolio[1][i],aerofolio[2][i]))
        sheet.append(("EndCurve",))
            
    sheet.append(("EndLoft",))
    sheet.append(("End",))
#    
    book.save(myfile)
    
    