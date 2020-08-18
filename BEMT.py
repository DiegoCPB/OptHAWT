# -*- coding: utf-8 -*-
"""
Created on Mon May 28 23:47:15 2018

@author: Diego Chou Pazo Blanco
"""

import numpy as np
import perfil as pr
import extrapolacao as ex
import apoio

class Elem_otimo():
    """
    Segundo a Teoria do Momento em um Elemento de Pá
    
    Condição ótima, solução analítica
    
    Para encontrar a condição ótima, as variáveis são inicializadas
    segundo os valores ótimos encontrados para a aproximação de Cl/Cd >> 1
    para o perfil.
    
    INPUTS:
    -------
    R0 :          Raio mínimo da pá (m)
    R :           Raio máximo da pá (m)
    N :           Numéro de pás
    v_inf :       Velocidade do vento (m/s)
    w :           Velocidade de rotação da pá (rad/s)
    r :           Posição do elemento na pá (0 < r < R; m)
    c :           Comprimento da corda (m)
    rho :         Densidade do ar (kg/m3)
    nome_perfil : String contendo o nome do perfil
    
    OUTPUTS:
    --------
    beta :   Ângulo do elemento em relação ao plano de rotação (graus)
    sigma :  Coeficiente de solidez
    aoa :    Ângulo de ataque do elemento
    a1 :     Fator a de indução translacional
    a2 :     Fator a' de indução rotacional
    cn :     Coeficiente aerodinâmico normal ao plano de rotação
    ct :     Coeficiente aerodinâmico tangencial ao plano de rotação
    phi :    Ângulo induzido no elemento
    dT :     Empuxo máximo
    dQ :     Torque máximo
    dP :     Potência máxima
    Fp :     Coeficiente de perda de Prandtl
    coords : Coordenadas do perfil
    """
   
    def __init__(self,R0,R,N,v_inf,w,r,
                 rho=1.225,nome_perfil="NACA 0011",Re=-1):
        # Argumentos de entrada
        
        # Essa correção esvita problemas de descontinuidade nas pontas da pá
        if R == r: 
            self.r = r-R*0.0001
        elif R0 == r: 
            self.r = r+R*0.0001
        else:
            self.r = r
            
        self.v1 = v_inf
        self.w = w
        self.R0 = R0
        self.R = R
        self.N = N
        self.rho = rho
        
        if "__" in nome_perfil:
            self.perfil = pr.PerfilInterpolado(nome_perfil,Re)
        else:
            self.perfil = pr.Perfil(nome_perfil,Re)
            
        self.coords = self.perfil.coords_padrao(50)
        
        # Loop principal de calculo
       
        self._main_calc()
        
        if self.a1 < 0.4:
            self.dT = 4*self.Fp*np.pi*self.r*self.rho*self.v1**2*(1-self.a1)*self.a1
        else:
            dCT = (50/9-4*self.Fp)*self.a1**2+(4*self.Fp-40/9)*self.a1+8/9
            self.dT = dCT*self.rho*self.v1*2*np.pi*self.r 
        self.dQ = 4*self.Fp*np.pi*self.r**3*self.rho*self.v1*self.w*(1-self.a1)*self.a2
        self.dP = self.dQ*self.w
        
    def _coeffs(self,phi, AR = -1):
        meuPerfil = self.perfil 
        aoa = self.aoa
        
        #Extrapolação de coeficientes pelo método de viterna
        fcd = meuPerfil._funcaoCd
        fcl = meuPerfil._funcaoCl
        limitesAlfa = meuPerfil.limitesAng
        if AR < 0:
            AR = 11 # AR é aproximado para 11
        fcz, fcx = ex.viterna(limitesAlfa,fcl,fcd,AR)
        
        cx = fcx(aoa)
        cz = fcz(aoa)
        
        cn = cz*np.cos(phi)+cx*np.sin(phi)
        ct = cz*np.sin(phi)-cx*np.cos(phi)
            
        return cn, ct
        
    def _Fp(self,phi):
        k = lambda r1,r2 : self.N/2.0*(r1-r2)/(self.r*np.sin(phi))
        F = lambda k : 2/np.pi*np.arccos(np.exp(-k))
        
        Ft = F(k(self.R,self.r))
        Fr = F(k(self.r,self.R0))
        
#        if np.isnan(Ft): Ft = 0
#        if np.isnan(Fr): Fr = 0
        return Ft*Fr
        
    def _main_calc(self):
        # Variáveis calculadas
        lambda_r = self.r*self.w/self.v1
        
        # Inicialização de variáveis na configuração ótima
        phi = 2./3*np.arctan(1/lambda_r)
        self.phi = np.rad2deg(phi)
        
        # Angulo de eficiencia máxima do perfil
        self.aoa = self.perfil.AlfaCl_Cdmax
        
         # Angulo de incidência
        self.beta = self.phi-self.aoa
        Cl = self.perfil.Cl(self.aoa)
        
        self.c = (1-np.cos(phi))*8*np.pi*self.r/(self.N*Cl)
        self.sigma = self.N*self.c/(2*np.pi*self.r)
        
        self.a1 = 1/(1+4*np.sin(phi)**2/(self.sigma*Cl*np.cos(phi)))
        self.a2 = (1-3*self.a1)/(4*self.a1-1) # a' - rotacional
        
        self.cn, self.ct = self._coeffs(phi)
        
        self.Fp = self._Fp(phi)
        
        
class Elem_generico(Elem_otimo):
    """
    Segundo a Teoria do Momento em um Elemento de Pá
    
    Condição qualquer, solução numérica
    
    INPUTS:
    -------
    R0 :          Raio mínimo da pá (m)
    R :           Raio máximo da pá (m)
    N :           Numéro de pás
    v_inf :       Velocidade do vento (m/s)
    w :           Velocidade de rotação da pá (rad/s)
    r :           Posição do elemento na pá (0 < r < R; m)
    c :           Comprimento da corda (m)
    beta :        Ângulo do elemento em relação ao plano de rotação (graus)
    rho :         Densidade do ar (kg/m3)
    nome_perfil : String contendo o nome do perfil
    
    OUTPUTS:
    --------
    sigma :  Coeficiente de solidez
    aoa :    Ângulo de ataque do elemento
    a1 :     Fator a de indução translacional
    a2 :     Fator a' de indução rotacional
    cn :     Coeficiente aerodinâmico normal ao plano de rotação
    ct :     Coeficiente aerodinâmico tangencial ao plano de rotação
    phi :    Ângulo induzido no elemento
    dT :     Empuxo
    dQ :     Torque
    dP :     Potência
    Cp :     Coeficiente de potência
    Fp :     Coeficiente de perda de Prandtl
    coords : Coordenadas do perfil
    """

    def __init__(self,R0,R,N,v_inf,w,r,c,beta,
                 rho=1.225,nome_perfil="NACA 0011",Re=-1):
        # Argumentos de entrada
        self.c = c
        self.beta = beta
    
        # Inicialização de variáveis
        self.sigma = N*c/(2*np.pi*r)
        self.Fp = 0
        self.dT = 0
        self.dQ = 0
        self.dP = 0
        self.phi = 0
        self.aoa = 0
        self.a1 = 0.25 # a - axial
        self.a2 = 0 # a' - rotacional
        self.cn = 0
        self.ct = 0
        
        #Inicialização dos parâmetros da classe de base
        super(Elem_generico,self).__init__(R0,R,N,v_inf,w,r,rho,nome_perfil,Re)
        self.c = c
        self.beta = beta
        
    
    def _main_calc(self, k_r = 0.25, pres=1e-4, max_count = 200):
        #Contador
        count = 0

        while True:
            count += 1
            if count > max_count:
#                print("\nBEMT: Número de iterações máximo alcançado.")
                self.a1 = np.nan
                self.a2 = np.nan
                break
        
            phi = np.arctan2(self.v1*(1-self.a1),self.r*self.w*(1+self.a2))

            self.phi = np.rad2deg(phi)
            self.aoa = self.phi-self.beta
            
            self.cn, self.ct = self._coeffs(phi, AR = (self.R-self.R0)/self.c)
            
            Fp = self._Fp(phi)
            
            CT = (self.sigma*(1-self.a1)**2*self.cn)/(np.sin(phi)**2)
            if CT > 0.96*Fp:
                a1 = (18*Fp-20-3*np.sqrt(CT*(50-36*Fp)+12*Fp*(3*Fp-4)))/(36*Fp-50)
            else:
                a1 = (4*Fp*np.sin(phi)**2/(self.sigma*self.cn) + 1)**-1
            
            a2 = (4*Fp*np.sin(phi)*np.cos(phi)/(self.sigma*self.ct) - 1)**-1
                
            k1 = abs(a1-self.a1)
            k2 = abs(a2-self.a2)
            erro = max(k1,k2)
            if erro < pres:
#                print("\nBEMT: convergência apos %d iteracoes." %(count))
                break   
                
            self.a1 = k_r*a1+(1-k_r)*self.a1
            self.a2 = k_r*a2+(1-k_r)*self.a2 
        
        self.Fp = Fp
        
        
if __name__ == "__main__":
    R = 0.5
    R0 = 0.05*R
    N = 3
    elem_perfil = "DU 93W300 LM"
    ar_rho = 1.225
    
    lamb = 20
    v_design = 7
    w_design = lamb*v_design/R
    
    r = R
            
    epo = Elem_otimo(R0,R,N,
                     v_design,w_design,r,
                     rho=ar_rho,nome_perfil=elem_perfil,Re=1e6)

    epg = Elem_generico(R0,R,N,
                        v_design,w_design,r,epo.c,epo.beta,
                        rho=ar_rho,nome_perfil=elem_perfil,Re= 1e6)
    print(epo.a1,epo.a2)
    print(epg.a1,epg.a2)
    print(epo.dP,epg.dP)
