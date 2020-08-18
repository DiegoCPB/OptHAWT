# -*- coding: utf-8 -*-
"""
Created on Sat Sep 15 16:11:51 2018

@author: Diego Chou
"""

import numpy as np
import BEMT
import apoio

class Pa_otima():
    """
    Objeto que constroi a pá ótima, calculada com a classe Elem_otimo,
    segunto o método BEMT.
    
    Paramêtros da classe:
    --------------------
    R0 :           Raio mínimo da pá (m)
    R :            Raio máximo da pá (m)
    N_pas :        Numéro de pás
    N_elems :      Número de elementos por pá
    lista_perfis : Lista com os nomes dos perfis aplicados à pá. ([string])
    pos_perfis :   Lista das posições dos perfis em porcentagem da envergadura (0 < pos_perfis[i] < 1)
    rho :          Massa específica do ar (kg/m^3)
    mu :           Viscosidade dinâmica do ar (m^2/s)
    tipo :         'tabelado' ou 'xfoil' (string)
    smooth :       Suavização das cordas e betas (Boolean)
    c_circular:    Diametro da seção circular, se existente (m)
    
    OBS: Obrigatoriamente len(lista_perfis) = len(pos_perfis)
    
    INPUTS da função calcular:
    -------------------------
    v_inf : Velocidade do vento (m/s)
    w :     Velocidade de rotação da pá (rad/s)
    
    OUTPUTS:
    --------
    area :   Area planiforme da pá (m^2)
    sigma :   Solidez da pá
    T :      Empuxo
    Q :      Torque
    P :      Potência
    Cp :     Coeficiente de potência
    x :      posições das seções na pá
    cordas : Lista das cordas ao longo da pá, segundo as posições do parametro x
    betas :  Lista de angulos ao longo da pá, segundo as posições do parametro x
    coords : Lista de coordenadas do perfil da seção, en função de x
    """
    
    def __init__(self,R0,R,N_pas,N_elems,lista_perfis,pos_perfis,
                 rho=1.225, mu = 18.2e-6, tipo='tabelado', smooth=True, c_circular = 0):
        self.R0 = R0
        self.R = R
        self.N = N_pas
        self.N_elems = N_elems
        self.rho = rho
        self.mu = mu
        self.tipo = tipo
        self.smooth = smooth
        self._lista_nomes_perfis = []
        
        if c_circular != 0:
            self.c_circular = c_circular
        else:
            self.c_circular = 0.05*R
                
        self.coords = {}
        
        if tipo not in ['tabelado','xfoil']:
            raise ValueError("O parametro tipo deve ter um dos dois valores: 'tabelado' ou 'xfoil'.")
        
        try:
            len(lista_perfis)
            len(pos_perfis)
        except:
            
            raise TypeError("Os paremetros lista_perfis e pos_perfis devem ser de um tipo iterável.")
            
        if len(lista_perfis) != len(pos_perfis):
            raise ValueError("Número de perfis incompatível com o número de posições.")
            
        if pos_perfis[0] != 0:
            raise ValueError("Não há perfil definido na raíz da pá. Reveja o parâmetro pos_perfis.")
        
        self.pos_perfis = pos_perfis
        self.lista_perfis = lista_perfis
        
    def _discretizacao(self):
        ang = np.linspace(0,np.pi, self.N_elems)
        self.x = self.R0 + 0.5*(self.R-self.R0)*(1-np.cos(ang))
        self.cordas = np.array([])
        self.betas = np.array([])
        
    def _string_perfil(self,pos):
        pos_perfis = self.pos_perfis
        lista_perfis = self.lista_perfis
    
        if len(lista_perfis) == 1:
            return lista_perfis[0]
        
        elif len(lista_perfis) > 1:
            #Transformação de porcentagem da envergadura para valores absolutos
            pos_perfis = self.R0 + (self.R-self.R0)*np.array(pos_perfis) 
            
            i_antes, i_depois = apoio.intervalo(pos_perfis,pos)
            pos_antes = pos_perfis[i_antes]
            pos_depois = pos_perfis[i_depois]
            perfil_antes = lista_perfis[i_antes]
            perfil_depois = lista_perfis[i_depois]
            
            porcentagem = (1-(pos-pos_antes)/(pos_depois-pos_antes))*100
            
            if int(np.round(porcentagem)) == 100:
                return perfil_antes
            elif int(np.round(porcentagem)) == 0:
                return perfil_depois
            else:
                return "%s__%s__%.1f" %(perfil_antes,perfil_depois,porcentagem)
        
        else:
            raise ValueError("Nome do perfil aerodinamico não informado.")
    
    def _integral(self,val):
        return 0.5*sum((val[1:]+val[:-1])*(self.x[1:]-self.x[:-1]))    
    
    def calcular(self,v_inf,w):
        self.v_inf = v_inf
        self.w = w
        c_var = False
        
        if 'circular' in self.lista_perfis:
            c_var = True
            if self.lista_perfis[0] != 'circular':
                raise ValueError("A seção circular deve sempre ser a primeira.")
            
            # Para o cálculo da geometria da pá, o perfil circular é substituido pelo seguinte.
            # Ele será restituído depois.
            self.lista_perfis[0] = self.lista_perfis[1] 
        
        #Realiza a discretização da pá e cria os parâmetros self.x, self.cordas, self.betas 
        self._discretizacao()
        
        T = np.array([])
        Q = np.array([])
        P = np.array([])
        v_rel = np.array([])
        
        for i in range(len(self.x)):
            string_perfil = self._string_perfil(self.x[i])
            self._lista_nomes_perfis.append(string_perfil)
            
            # Configuração ótima do elemento pela fórmula analítica
            e_analitico = BEMT.Elem_otimo(self.R0,self.R,self.N,self.v_inf,self.w,self.x[i],
                                          self.rho,nome_perfil=string_perfil)
            
            self.cordas = np.append(self.cordas,e_analitico.c)
            self.betas = np.append(self.betas,e_analitico.beta)
            self.coords[i] = e_analitico.coords
    
            if self.tipo == 'xfoil':
                v_rel = np.append(v_rel,np.sqrt((self.v_inf*(1-e_analitico.a1))**2+(self.x[i]*self.w*(1+e_analitico.a2))**2))
            
            elif self.tipo == 'tabelado':
                T = np.append(T,e_analitico.dT)
                Q = np.append(Q,e_analitico.dQ)
                P = np.append(P,e_analitico.dP)
        
        if self.tipo == 'xfoil':
            v_rel = np.mean(v_rel)
            S_asa = sum(0.5*(self.cordas[1:]+self.cordas[:-1])*(self.x[1:]-self.x[:-1]))
            c_aero = 1/S_asa*(sum((0.5*(self.cordas[1:]+self.cordas[:-1]))**2*(self.x[1:]-self.x[:-1])))
            Re_medio = c_aero*v_rel*self.rho/self.mu 
            
            self.cordas = np.array([])
            self.betas = np.array([])
            
            for i in range(len(self.x)):
                string_perfil = self._string_perfil(self.x[i])
                e_analitico = BEMT.Elem_otimo(self.R0,self.R,self.N,self.v_inf,self.w,self.x[i],
                                              self.rho,nome_perfil=string_perfil,Re = Re_medio)
                
                self.cordas = np.append(self.cordas,e_analitico.c)
                self.betas = np.append(self.betas,e_analitico.beta)
                # self.coords não precisa ser recalculado
                
                T = np.append(T,e_analitico.dT)
                Q = np.append(Q,e_analitico.dQ)
                P = np.append(P,e_analitico.dP)
        
        if c_var:
            # Restituição do perfil circular e ajuste da geometria da pá
            self.lista_perfis[0] = 'circular'
            corda_c = self.c_circular
            
            #Parametros do perfil seguinte
            pos_p = self.R0+self.pos_perfis[1]*(self.R-self.R0)
            index_p = np.searchsorted(self.x,pos_p)
            corda_p = self.cordas[index_p]
            
            for i in range(index_p):
                div = (pos_p+self.R0)/2
                if self.x[i] <= div:
                    self.cordas[i] = corda_c
                else:
                    self.cordas[i] = apoio.interplinear(div,corda_c,
                                                        pos_p,corda_p, self.x[i])
                
        # As cordas e betas são alisadas por um filtro de savgol de 2º grau
        if self.smooth:
            from scipy.signal import savgol_filter
            wl = self.N_elems//4
            if wl%2 == 0:
                wl+=1
            for i in range(2):
                self.cordas = savgol_filter(self.cordas, wl, 2, mode='mirror')
                self.betas = savgol_filter(self.betas, wl, 2, mode='mirror')
          
        self.area = self._integral(self.cordas)
        self.sigma = self.N*self.area/(np.pi*self.R**2)
        self.T = self._integral(T)
        self.Q = self._integral(Q)
        self.P = self._integral(P)
        self.Cp = self.P/(0.5*self.rho*np.pi*self.R**2*self.v_inf**3)
        
        
    
class Pa_generica(Pa_otima):
    """
    Objeto que constroi uma pá qualquer, calculada com a classe Elem_generico,
    segunto o método BEMT.
    
    INPUTS:
    -------
    R0 :           Raio mínimo da pá (m)
    R :            Raio máximo da pá (m)
    N :            Numéro de pás
    lista_perfis : Lista com os nomes dos perfis aplicados à pá.
    pos_perfis :   Lista das posições dos perfis em porcentagem da envergadura (0 < pos_perfis[i] < 1)
    x :            posições das seções na pá
    cordas :       Lista das cordas ao longo da pá, segundo as posições do parametro x
    betas :        Lista de angulos ao longo da pá, segundo as posições do parametro x
    
    OBS: Obrigatóriamente len(lista_perfis) = len(pos_perfis)
    
    INPUTS da função calcular:
    -------------------------
    v_inf : Velocidade do vento (m/s)
    w :     Velocidade de rotação da pá (rad/s)
    
    OUTPUTS:
    --------
    area :    Area planiforme da pá (m^2)
    sigma :   Solidez da pá
    T :       Empuxo
    Q :       Torque
    P :       Potência
    Cp :      Coeficiente de potência
    coords :  Lista de coordenadas do perfil da seção, en função de x
    N_elems : Número de elementos por pá
    """
    
    def __init__(self,N_pas,lista_perfis,pos_perfis,x,cordas,betas,
                 rho=1.225,mu = 18.2e-6,tipo='tabelado'):        
        try:
            len(x)
            len(cordas)
            len(betas)
        except:
            raise TypeError("Os paremetros x, cordas e betas devem ser de um tipo iterável.")
            
        if len(x) != len(cordas):
            raise ValueError("len(x) != len(cordas)")
            
        if len(x) != len(betas):
            raise ValueError("len(x) != len(betas)")
    
        self.x = x
        self.cordas = cordas
        self.betas = betas
        
        self.area = self._integral(cordas)
        self.sigma = N_pas*self.area/(np.pi*x[-1]**2)
        
        #Inicialização dos parâmetros da classe de base
        super(Pa_generica,self).__init__(x[0],x[-1],N_pas,len(x),lista_perfis,pos_perfis,
                                         rho,mu,tipo,False)
   
    def calcular(self,v_inf,w):
        self.v_inf = v_inf
        self.w = w
        
        T = np.array([])
        Q = np.array([])
        P = np.array([])
        v_rel = np.array([])
        
        for i in range(len(self.x)):
            string_perfil = self._string_perfil(self.x[i])
            self._lista_nomes_perfis.append(string_perfil)
            # Configuração ótima do elemento pela fórmula analítica
            e_numerico = BEMT.Elem_generico(self.R0,self.R,self.N,self.v_inf,self.w,self.x[i],
                                            self.cordas[i],self.betas[i],
                                            self.rho,nome_perfil=string_perfil)
            self.coords[i] = e_numerico.coords
    
            if self.tipo == 'xfoil':
                local_v = np.sqrt((self.v_inf*(1-e_numerico.a1))**2+(self.x[i]*self.w*(1+e_numerico.a2))**2)
                if local_v > self.v_inf:
                    v_rel = np.append(v_rel,local_v)
            
            elif self.tipo == 'tabelado':
                T = np.append(T,apoio.nan2zero(e_numerico.dT))
                Q = np.append(Q,apoio.nan2zero(e_numerico.dQ))
                P = np.append(P,apoio.nan2zero(e_numerico.dP))
        
        if self.tipo == 'xfoil':
            v_rel = np.mean(v_rel[~np.isnan(v_rel)]) #Os eventuais NaN são filtrados
            S_asa = sum(0.5*(self.cordas[1:]+self.cordas[:-1])*(self.x[1:]-self.x[:-1]))
            c_aero = 1/S_asa*(sum((0.5*(self.cordas[1:]+self.cordas[:-1]))**2*(self.x[1:]-self.x[:-1])))
            Re_medio = c_aero*v_rel*self.rho/self.mu 
        
            
            for i in range(len(self.x)):
                string_perfil = self._string_perfil(self.x[i])
                
                e_numerico = BEMT.Elem_generico(self.R0,self.R,self.N,self.v_inf,self.w,self.x[i],
                                                 self.cordas[i],self.betas[i],
                                                 self.rho,nome_perfil=string_perfil,Re = Re_medio)
                
                T = np.append(T,apoio.nan2zero(e_numerico.dT))
                Q = np.append(Q,apoio.nan2zero(e_numerico.dQ))
                P = np.append(P,apoio.nan2zero(e_numerico.dP))
                
        self.T = self._integral(T)
        self.Q = self._integral(Q)
        self.P = self._integral(P)
        self.Cp = self.P/(0.5*self.rho*np.pi*(self.R**2)*self.v_inf**3)
    
if __name__ == '__main__':   
    import graficos

    R = 51.33016479
    R0 = 2.5
    c_raiz = 4.2
    N = 3
    ar_rho = 1.225
    
    lamb = 6
    v = 5
    w =lamb*v/R
    
    N_elems = 25
    lista_perfis = ['circular','NACA 63630',
                    'FW W3 211','FW W1 152']
    pos_perfis = [0,0.16654043,0.21665741,0.29001587,0.73603106,0.95,1]
    
    a = Pa_otima(R0,R,N,N_elems,lista_perfis,pos_perfis,rho=ar_rho,smooth=True,
                 c_circular = c_raiz)
    a.calcular(v,w)
    
    b = Pa_generica(N,lista_perfis,pos_perfis,a.x,a.cordas,a.betas+2)
    b.calcular(v,w)

#    graficos.geometria(b)
    graficos.planta(b)
#    graficos.desempenho(b)
#    graficos.potência(b)
#    graficos.potência(b, v_bound=[1,20])
#    apoio.reset_dados_perfis()