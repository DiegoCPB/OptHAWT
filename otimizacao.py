# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 17:14:18 2018

@author: Diego Chou
"""

import numpy as np
from pso import pso
import Pa_BEMT

class Otimizador():
    """
    Classe que engloba todas as variáveis necessárias à otimização aerodinâmica
    das pás da turbina eólica horizontal.
    
    INPUT:
    ------
    pot :          Potência eólica requerida de projeto (W)
    v :            Velocidade do vento (m/s)
    w :            VElocidade de rotação (rad/s)
    lista_perfis : Lista dos perfis empregados na pá
    r_hub :        Raio do hub central
    N :            Número de pás
    N_elems :      Número de seções por pá
    rho :          Massa específica do ar
    c_circular:    Diametro da seção circular, se existente (m)
    """
    
    def __init__(self, pot, v, w, r_hub, N,lista_perfis, 
                 N_elems=20, rho = 1.225, c_circular = 0, tipo = 'tabelado'):
        self.tipo = tipo
        self.pot = pot
        self.v = v
        self.w = w
        
        self.lista_perfis = lista_perfis
        self.r_hub = r_hub
        
        self.rho = rho
        self.N_elems = N_elems
        self.N = N
        
        self.c_circular = c_circular
        
        self._var = self._otimo_teorico()
        
        self._ub,self._lb = self._limites()
        self._var_opt, self.pa_opt = self._otimo_pso()
    
    def _calc_blade(self, R, pos_perfis, tipo, final = False):
        if final == True:
            N_elems = 25
        else:
            N_elems = self.N_elems
        
        #Cálculo da pá
        blade = Pa_BEMT.Pa_otima(self.r_hub,R,self.N,N_elems,self.lista_perfis,pos_perfis,
                                 smooth=True, rho = self.rho, c_circular = self.c_circular)
        blade.calcular(self.v,self.w)
        
        blade = Pa_BEMT.Pa_generica(self.N,self.lista_perfis,pos_perfis,
                                    blade.x,blade.cordas,blade.betas, rho = self.rho,
                                    tipo = tipo)
        blade.calcular(self.v,self.w)
        
        return blade
        
    def _otimo_teorico(self, pres = 1e-3):
        """
        Método que define o raio da pá iterativamente.
        """
        print("\nCalculando pá ótima teórica para P = %.1f W..." %(self.pot))
        
        #Posições relativas dos aerofólios sobre a pá
        pos_perfis = np.linspace(0,1,len(self.lista_perfis))
        
        #Limites inferior e superior de R
        Cp_betz = 16/27
        R = np.sqrt(self.pot/(Cp_betz*0.5*self.rho*np.pi*self.v**3))
        
        count = 0
        while True:
            count += 1
            blade = self._calc_blade(R,pos_perfis,tipo='tabelado')
            Cp = blade.Cp
            Ri = np.sqrt(self.pot/(0.5*self.rho*self.v**3*np.pi*Cp)) 
            if np.isnan(Ri):
                raise ValueError("Reveja o valor do raio do cubo")
            elif abs(Ri-R) < pres:
                break
            else:
                R = Ri
       
        print("\tNúmero de iterações:\t %d" %(count))
        print("\tRaio de pá encontrado:\t %.3f m" %(R))
        print("\tPotência estimada:\t %.3f W" %(blade.P))
        
        var = np.concatenate([[R],pos_perfis[1:-1]])
        return var
        
    def _limites(self, delta = 0.05):
        R_max = self._var[0]
        pos_perfis = self._var[1:]
        
        #Limites para R em m
        ub_R = [R_max]
        lb_R = [(1-delta)*R_max]
        
        #Limites para pos_perfis. o primeiro perfil é fixado `0 e o ultimo `1
        ub_pos = [1-delta]*len(pos_perfis)
        lb_pos = [delta]*len(pos_perfis)
        
        return np.concatenate([ub_R,ub_pos]), np.concatenate([lb_R,lb_pos])
    
    def _restricao1(self,var):
        """
        Restrição imposta para que a ordem dos perfis não mude.
        """
        pos_perfis = list(var[1:])

        if pos_perfis == sorted(pos_perfis):
            return 1
        else:
            return -1
        
    def _restricao2(self,var,delta=0.05):
        """
        Restrição imposta para que a ordem dos perfis não mude.
        """
        pos_perfis = np.sort(var[1:])
        dif_pos = pos_perfis[1:]-pos_perfis[:-1]

        if np.all(dif_pos > delta):
            return 1
        else:
            return -1

    def _pontuacao(self,var,tipo = 'tabelado',final = False):
        """
        Função que avalia uma pá genérica.
        Recebe como entrada o argumento x, que possui a seguinte estrutura:
        
        var : [R, pos_perfis]
                       
        pos_perfis :   posições dos perfis, excluidos os dois das pontas, 
                       que variam entre 0 e 1. 
        """    
        R = var[0]
        
        pos_perfis = [0]+sorted(list(var[1:]))+[1]
                
        # O algoritmo minimiza essa função
        if final:
            return self._calc_blade(R,pos_perfis,tipo,final)
        else:
            blade = self._calc_blade(R,pos_perfis,tipo)
            k_pot = 1+abs(blade.P-self.pot)/self.pot
            pontuacao = 1000*(16/27-blade.Cp)*k_pot*blade.sigma
            
            return pontuacao
    
    def _otimo_pso(self,pres=1e-3):
        """
        Calculo principal
        """
        if len(self._var) == 1:
            var_opt = self._var
        else:
            var_opt = pso(self._pontuacao,self._lb,self._ub,
                          seed = self._var,ieqcons=[self._restricao1,self._restricao2],
                          omega=0.2, phip=0.4, phig=0.4,
                          swarmsize = 15, maxiter = 50, minstep=1e-3, minfunc=1e-4, debug = True)[0]
        
        print("\nRealizando cálculos finais...")
        if self.tipo == 'xfoil':
            print("\nCalculando polares no XFOIL") 
            pa_opt = self._pontuacao(var_opt,tipo = self.tipo,final = True)
            print("\tPotência alcancada XFoil: {:} W \n\tCp alcancado XFoil: {:}".format(pa_opt.P,pa_opt.Cp))
        pa_opt = self._pontuacao(var_opt,tipo = 'tabelado',final = True)
        print("\tPotência alcancada: {:} W \n\tCp alcancado: {:}".format(pa_opt.P,pa_opt.Cp))
        
        return var_opt, pa_opt
        
if __name__ == "__main__":
    import graficos
    import apoio
    #               pot, v,  c_lambda, r_hub, N,
    ot = Otimizador(30e3, 10, 7.6,  0.2,   3, 
                    ['DU 93W405 AD','DU 93W350 AD',
                     'DU 93W300 LM','DU 93W250 LM','DU 93W210 LM',
                     'NACA 64618','NACA 64618'])

    
#    graficos.geometria(ot.pa_opt)
    graficos.planta(ot.pa_opt)
#    graficos.potência(ot.pa_opt, fixa = False, maxP = 5e6)
    graficos.desempenho(ot.pa_opt,lambda_bound=[1,18])
    apoio.reset_dados_perfis()
    
        