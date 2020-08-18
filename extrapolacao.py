# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 23:40:31 2018

@author: Diego Chou
"""

import numpy as np

def viterna(limitesAlfa,funcaoCl,funcaoCd,AR):
    """
    Extrapolação das polares de sustentação e arrasto para o perfil
    segundo o metodo de viterna.
    
    Retorna, a partir de valores de Cl e Cd limitados, funções válidas
    para todos os ângulos.
    
    INPUTS:
    ------
    - limitesAlfa : Limites dos dados tabelados para o perfil [alfa_min,alfa_max]
    - funcaoCl :    Polar de sustentação do perfil
    - funcaoCd :    Polar de arrsto do perfil
    - AR :          Razão de aspecto (b**2/S) da pá em questão, pode ser aproximado para 11.
    
    OUTPUT:
    - extrapCl : Função extrapolada dos dados tabelados
    - extrapCd : Função extrapolada dos dados tabelados
    
    """ 
    stall_minus, stall_plus = np.deg2rad(limitesAlfa)
    
    
    def coeffs(alfa_stall):
        Cdmax = 1.11+0.018*17#AR #Fórmula válida para AR < 50
        A1 = Cdmax/2
        B1 = Cdmax
        
        Cl_stall = funcaoCl(np.rad2deg(alfa_stall))
        Cd_stall = funcaoCd(np.rad2deg(alfa_stall))
        
        A2 = (Cl_stall-Cdmax*np.sin(alfa_stall)*np.cos(alfa_stall))*np.sin(alfa_stall)/(np.cos(alfa_stall)**2)
        B2 = (Cd_stall-Cdmax*np.sin(alfa_stall)**2)/np.cos(alfa_stall)
        
        def Cl_viterna(alfa):
            alfa = np.deg2rad(alfa)
            return A1*np.sin(2*alfa)+A2*np.cos(alfa)**2/np.sin(alfa)
        def Cd_viterna(alfa):
            alfa = np.deg2rad(alfa)
            return B1*np.sin(alfa)**2+B2*np.cos(alfa)
        return Cl_viterna, Cd_viterna

    Cl_plus,Cd_plus = coeffs(stall_plus)
    Cl_minus,Cd_minus = coeffs(stall_minus)
    
    def extrapCl(alfa):
        """
        Alfa em graus
        """
        if alfa > limitesAlfa[0] and alfa < limitesAlfa[1]:
            return funcaoCl(alfa)
        elif alfa >= 360+limitesAlfa[0]:
            return extrapCl(alfa-360)
        elif alfa >= limitesAlfa[1]:
            return Cl_plus(alfa)
        elif alfa <= limitesAlfa[0]:
            return Cl_minus(alfa)
    
    def extrapCd(alfa):
        """
        Alfa em graus
        """
        if alfa > limitesAlfa[0] and alfa < limitesAlfa[1]:
            return funcaoCd(alfa)
        elif alfa >= 360+limitesAlfa[0]:
            return extrapCd(alfa-360)
        elif alfa >= limitesAlfa[1]:
            return Cd_plus(alfa)
        elif alfa <= limitesAlfa[0]:
            return Cd_minus(alfa)

    return extrapCl, extrapCd
    
if __name__ == '__main__':
    import perfil
    p = perfil.Perfil("S807",Re=7e6)
    x1,y1 = p.coords_padrao(50)
    fCl = p._funcaoCl
    fCd = p._funcaoCd
    lAlfa = p.limitesAng
    AR = 20
    
    eCl,eCd = viterna(lAlfa,fCl,fCd,AR)
    
    ang = np.linspace(-100,100,500)
    lCl = []
    lCd = []
    for i in ang:
        lCl.append(eCl(i))
        lCd.append(eCd(i))
    
    import matplotlib.pyplot as plt
    
#    plt.figure()
#    plt.grid(True)
#    plt.axis('equal')
#    plt.xlim([-0.1,1.1])
#    plt.plot(x1,y1,'k')#,label = '%s spline ' %(name))
#    plt.show()
    

    fig = plt.figure()
    fig.suptitle(p.name + ' - $Re=%d$' %(p.Re))
    ax1 = fig.add_subplot(2, 1, 1)
    ax1.plot(ang,lCl,'k:',label="VITERNA (AR = %.1f)" %AR)
    ax1.plot(np.linspace(lAlfa[0],lAlfa[1]),fCl(np.linspace(lAlfa[0],lAlfa[1])), label="XFoil")
    ax1.grid(True)
    ax1.set_ylabel(r'$C_l$')
    ax1.legend()
    
    ax2 = fig.add_subplot(2, 1, 2)
    ax2.plot(ang,lCd,'k:')
    ax2.plot(np.linspace(lAlfa[0],lAlfa[1]),fCd(np.linspace(lAlfa[0],lAlfa[1])))
    ax2.grid(True)
    ax2.set_ylabel(r'$C_d$')
    ax2.set_xlabel(r'$\alpha$ ($graus$)')
    plt.show()
    
    
    