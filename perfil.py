# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 22:34:49 2018

@author: Diego Chou
"""
import os
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline, splprep, splev
import apoio
import xfoil_module as xfm

class Perfil(): 
    """
    Classe que define as propriedades de um perfil aerodinâmico.
    
    INPUTS:
    -------
    name : Nome do arquivo de coordenadas do perfil sem extensão
    Re :   Número de Reynolds. Se negativo, o programa leva em conta a polar tabelada.
           Se positivo, a polar aerodinâmica sera calculada no referido Reynolds.
    
    OUTPUTS:
    --------
    limitesAng :    Limites superior e inferior dos angulos de ataque da polar.
    Clmax :         Coeficiente de sustentação máximo        
    AlfaClmax :     Angulo de Clmax
    Cdmin :         Coeficiente de arrasto mínimo
    AlfaCdmin :     Angulo de Cdmin
    Cl_Cdmax :      Cl/Cd máximo
    AlfaCl_Cdmax :  Angulo de Cl_Cdmax 
    dCl_dAlfa :     Variação do Cl em função do angulo de ataque em graus
    dCm_dAlfa :     Variação do Cm em função do angulo de ataque em graus
    ac :            Centro aerodinâmico do perfil em % da corda
    coords_padrao : Coordenadas do perfil 
    centroide :     
        
    FUNÇÕES:
    --------
    Cl(alfa) : Cl em função do ângulo de ataque em graus 
    Cd(alfa) : Cd em função do ângulo de ataque em graus
    Cm(alfa) : Cm em função do ângulo de ataque em graus
    
    """
    
    def __init__(self, name, Re=-1):
        self.name = name
        self.Re = Re
        
        if '__' not in name:
            #Verificação dos limites de ângulo de ataque tabelados
            if name == 'circular':
                self.limitesAng = [-180,180]
            else:
                self.limitesAng = self._dominioAng()
            
            # Calculo dos splines dos coeficientes aerodinâmicos
            self._funcaoCl, self._funcaoCd, self._funcaoCm = self._funcaoClCdCm()
            
            # Valores fixos calculados
            if name == 'circular':
                self.Clmax, self.AlfaClmax = 0,0
                self.Clmin, self.AlfaClmin = 0,0
                self.Cdmin, self.AlfaCdmin = self._funcaoCd(0),0
                self.Cl_Cdmax, self.AlfaCl_Cdmax = 0,0 
                self.dCl_dAlfa, self.dCm_dAlfa = 0,0
            else:
                self.Clmax, self.AlfaClmax = self._readClmax()
                self.Cdmin, self.AlfaCdmin = self._readCdmin()
                self.Cl_Cdmax, self.AlfaCl_Cdmax = self._readCl_Cdmax() 
                self.dCl_dAlfa, self.dCm_dAlfa = self._ajustelinear()
            
            # Centro aerodinâmico do perfil. Válido somente se o Cm for medido 
            # à 25% da corda.
            if name == 'circular':
                self.ac = 0.5
            else:
                self.ac = 0.25-self.dCm_dAlfa/self.dCl_dAlfa 
    
    @apoio.executarNaPasta('..')
    def _calcXfoilPolar(self, n =-8):
        """
        Função que calcula a polar do perfil pelo xfoil
        """
        minAlfa,maxAlfa = self.limitesAng
        print("Entre os ângulos de ataque %.1f e %.1f graus" %(minAlfa,maxAlfa))
        lista_alfas = np.arange(minAlfa,maxAlfa+0.5,0.5)
        
        while True:
            nRe = np.round(self.Re,n)#Arredondando o numero de Reynolds
            if nRe == 0: n+=1
            else: break
        
        xfm.save_polar("Perfis/"+self.name+".dat", alfas=lista_alfas, Reynolds=nRe)            
            
    @apoio.executarNaPasta('Perfis')
    def _readPoints(self):
        """
        Função que retorna os pontos do arquivo do perfil.
        O arquivo deve estar no formato aceito pelo XFLR5
        """
        name = self.name
        filename = '%s.dat' %(name)
        
        try:
            f = open(filename, 'r')
        except IOError as er:
            print(er)
            print('')
            raise
            
        flines = f.readlines()[1:]
        listaX = []
        listaY = []
        
        for i in range(len(flines)):
#            print flines[i]
            words = flines[i].split()
#            print words[1]
            try:
                x = float(words[0])
                y = float(words[1])
                listaX.append(x)
                listaY.append(y)
            except:
                pass    
            
        #normalização do perfil
        corda = max(listaX)
        listaX = np.array(listaX)/corda
        listaY = np.array(listaY)/corda
        
        return list(listaX),list(listaY)
    
    @apoio.executarNaPasta('Perfis')
    def _openFile(self,var = False,n=-8):
        name = self.name
        
        #Se o número de Reynolds for negativo, os valores tabelados são usados
        #Se positivo, os valores são recalculados no XFoil
        if self.Re < 0 or var == True:
            filename = '%s_polar.txt' %(name)
        else:
            # Calcula e salva a polar do perfil usando o XFoil 
            while True:
                nRe = np.round(self.Re,n)#Arredondando o numero de Reynolds
                if nRe == 0: n+=1
                else: break
            filename = '%s_polarXFoil_Re%.3f.txt' %(name, nRe*1e-6)
            if not os.path.isfile(filename):
                print("\nCalculando polar do perfil %s (Re = %d) pelo Xfoil..." %(self.name, nRe))
                self._calcXfoilPolar()
                print("OK")
            
        try:
            f = open(filename, 'r')
        except IOError:
            raise
            
        flines = f.readlines()
        return flines
    
    
        
    def _readAllClCdCm(self,var = False):
        flines = self._openFile(var)
        lista = []
        
        for i in range(11,len(flines)):
            words = flines[i].split()
            try:
                Cl = float(words[1])
                Cd = float(words[2])
                Cm = float(words[4])
                angulo = float(words[0])
                lista.append([angulo, Cl, Cd, Cm])
            except:
                pass
       
        return lista
    
    def _readClCdCm(self, alfa):       
        flines = self._openFile()

        for i in range(11,len(flines)):
            words = flines[i].split()
            try:
                Cl = float(words[1])
                Cd = float(words[2])
                Cm = float(words[4])
                angulo = float(words[0])
                if angulo == alfa:
                    lista = [Cl,Cd,Cm]
                    break
            except:
                pass
            
        return lista
                    
    def _readClmax(self):
        flines = self._openFile()
        CLmax = 0
        angulo_clmax = 0
        for i in range(11,len(flines)):
            words = flines[i].split()

            try:
                CL = float(words[1])
                angulo = float(words[0])
                if(CL>CLmax):
                    CLmax = CL
                    angulo_clmax = angulo 
            except:
                pass
                    
        return CLmax, angulo_clmax
    
    def _readCdmin(self):
        flines = self._openFile()
        CDmin = 100
        angulo_cdmin = 0
        for i in range(11,len(flines)):
            words = flines[i].split()
            
            try:
                CD = float(words[2])
                angulo = float(words[0])
                if(CD<CDmin):
                    CDmin = CD
                    angulo_cdmin = angulo 
            except:
                pass
                    
        return CDmin, angulo_cdmin
    
    def _readCl_Cdmax(self):
        flines = self._openFile()
        
        Cl_Cdmax = 0
        angulo_Cl_Cdmax = 0
        for i in range(11,len(flines)):
            words = flines[i].split()
            
            try:
                Cd = float(words[2])
                Cl = float(words[1])
                angulo = float(words[0])
                if(Cl/Cd>Cl_Cdmax):
                    Cl_Cdmax = Cl/Cd
                    angulo_Cl_Cdmax = angulo 
            except:
                pass
                    
        return Cl_Cdmax, angulo_Cl_Cdmax
    
    def _spline(self):
        """
        Gera 2 splines a partir dos pontos do perfil. 
        Um para o extradorso e outro para o intradorso.
        """ 
        
        x,y = self._readPoints()
        
        if y[1] < y[-2]:
            x = x[::-1]
            y = y[::-1]
            
        ba = min(x) # abscissa do bordo de ataque
        index = x.index(ba) # posição do ponto na lista
        
        if abs(y[0]-y[-1]) >= 1e-5:
            x_extra = x[:index+1]
            y_extra = y[:index+1]
            x_intra = x[index:]
            y_intra = y[index:]
        else:
            delta_y = abs(y[1]-y[-2])
            x_extra = x[:index+1]
            y_extra = np.append(y[0]+delta_y/2,y[1:index+1])
            x_intra = x[index:]
            y_intra = np.append(y[index:-1],y[-1]-delta_y/2)
        
        s_extra = splprep([x_extra,y_extra], s=0)[0]
        s_intra = splprep([x_intra,y_intra], s=0)[0]
        
        return s_extra, s_intra
    
    def _funcaoClCdCm(self):
        """
        Função que retorna os coeficientes do perfil para qualquer
        angulo dentro do intervalo de pontos,
        """
        if self.name == 'circular':
            sCl = lambda x: 0.0
            sCd = lambda x: 0.35
            sCm = lambda x: 0.0
        else:
            lista = np.array(self._readAllClCdCm())
            alfa = lista[:,0]
            Cl = lista[:,1]
            Cd = lista[:,2]
            Cm = lista[:,3]
            
            sCl = InterpolatedUnivariateSpline(alfa,Cl,k=3)
            sCd = InterpolatedUnivariateSpline(alfa,Cd,k=3)
            sCm = InterpolatedUnivariateSpline(alfa,Cm,k=3)
            
        return sCl, sCd, sCm
    
    def _dominioAng(self):
        try:
            lista = np.array(self._readAllClCdCm(True))
            alfa = lista[:,0]
        except:
            alfa = np.arange(-21,21)
        
        return min(alfa), max(alfa)
    
    def _verifLimites(self,alfa):
        if alfa < self.limitesAng[0] or alfa > self.limitesAng[1]:
            return False
        else:
            return True
    
    def _ajustelinear(self):
        alpha_CLmax = self.AlfaClmax
        sCL = self._funcaoCl
        sCM = self._funcaoCm
        
        alphas = np.linspace(0,alpha_CLmax-3.0)
        CL = sCL(alphas)        
        CM= sCM(alphas)
        
        a0, a1 = np.polyfit(alphas, CL, deg=1)
        m0, m1 = np.polyfit(alphas, CM, deg=1)
        return a0, m0   
    
    def _tri(self): 
        """
        Essa função triangula o perfil para que a área possa ser calculada
        pelo produto vetorial dos vetores dos triângulos gerados.
        """
        x,y = self.coords_padrao()
        tri = np.array([[0.0, 0.0, 0.0]])
        N = len(x)
        Ni = N//2 
        
        if Ni == N/2.0: # Caso com número par de pontos
            val = 1
        else:
            val = 0
       
        # Triangulacao do 1º perfil
        # Vale a regra da mão direita
        for i in range(0,Ni-1):
            tri = np.append(tri,[[i, i+1, N-i-2]], axis=0)
        
        if y[0] != y[-1]:
            tri = np.append(tri,[[0, N-2, N-1]], axis=0)
            
        for i in range(1,Ni-val):
            tri = np.append(tri,[[i, N-i-2, N-i-1]], axis=0)
        
        return tri[1:]      
    
    def centroide(self):
        """
        Calcula a posição do centróide do perfil em % da corda.
        """
        x,y = self.coords_padrao()
        x = np.array(x)
        y = np.array(y)
        x_centroide = 0.0 
        y_centroide = 0.0
        area_total = 0.0

        tri = self._tri()
        
        for i in range(len(tri)):
            # Calculando o vetor normal
            x0 = x[int(tri[i][0])]
            x1 = x[int(tri[i][1])]
            x2 = x[int(tri[i][2])]
            
            y0 = y[int(tri[i][0])]
            y1 = y[int(tri[i][1])]
            y2 = y[int(tri[i][2])]   
            
            xc = (x0+x1+x2)/3
            yc = (y0+y1+y2)/3
            AB = np.array([x1-x0, y1-y0,0.0])
            AC = np.array([x2-x0, y2-y0,0.0])
            p_vetorial = np.cross(AB, AC)
            area = np.linalg.norm(p_vetorial)
            x_centroide += area*xc
            y_centroide += area*yc
            area_total += area
        
        return x_centroide/area_total, y_centroide/area_total

    def coords_padrao(self,n=100):
        """
        A partir dos pontos lidos no arquivo .dat do perfil, gera um conjunto
        de coordenadas padrao, com pontos em posicoes especificas dos intra e
        extradorsos.
        
        Input:
        -----
        n : numero de pontos sobre o contorno do perfil
        """
        beta = np.linspace(0,np.pi,n//2)
        s1,s2 = self._spline()
        
        u1 = 0.5*(1-np.cos(beta))
        u2 = 1-u1
        
        x_extra,y_extra = splev(u1,s1)
        x_intra,y_intra = splev(u2,s2)
        
        x = np.append(x_extra[:-1],x_intra[::-1])
        y = np.append(y_extra[:-1],y_intra[::-1])
        
        return x,y
        
    def Cl(self,alfa):
        if self._verifLimites(alfa):
            return self._funcaoCl(alfa)
        else:
            print("\nATENCAO: AoA %.1f graus fora do envelope Cl de %s" %(alfa,self.name))
            return np.nan
    
    def Cd(self,alfa):
        if self._verifLimites(alfa):
            return self._funcaoCd(alfa)
        else:
            print("\nATENCAO: AoA %.1f graus fora do envelope Cd de %s" %(alfa,self.name))
            return np.nan
    
    def Cm(self,alfa):
        if self._verifLimites(alfa):
            return self._funcaoCm(alfa)
        else:
            print("\nATENCAO: AoA %.1f graus fora do envelope Cm de %s" %(alfa,self.name))
            return np.nan

        
class PerfilInterpolado(Perfil): 
    """
    Classe que define as propriedades de um perfil aerodinâmico.
    
    INPUTS:
    -------
    name : Nome do arquivo de coordenadas do perfil sem extensão, no formato
           "perfil1__perfil2__porcentagemInterpolacao"
    Re :   Número de Reynolds. Se negativo, o programa leva em conta a polar tabelada.
           Se positivo, a polar aerodinâmica sera calculada no referido Reynolds.
    
    OUTPUTS:
    --------
    name1 :         Nome do primeiro perfil
    name2 :         Nome do segundo perfil
    pct :           Porcentagem de interpolacao   
    limitesAng :    Limites superior e inferior dos angulos de ataque da polar.
    Clmax :         Coeficiente de sustentação máximo        
    AlfaClmax :     Angulo de Clmax
    Cdmin :         Coeficiente de arrasto mínimo
    AlfaCdmin :     Angulo de Cdmin
    Cl_Cdmax :      Cl/Cd máximo
    AlfaCl_Cdmax :  Angulo de Cl_Cdmax 
    dCl_dAlfa :     Variação do Cl em função do angulo de ataque em graus
    dCm_dAlfa :     Variação do Cm em função do angulo de ataque em graus
    ac :            Centro aerodinâmico do perfil em % da corda
    coords_padrao : Coordenadas do perfil interpolado
        
    FUNÇÕES:
    --------
    Cl(alfa) : Cl em função do ângulo de ataque em graus 
    Cd(alfa) : Cd em função do ângulo de ataque em graus
    Cm(alfa) : Cm em função do ângulo de ataque em graus
    
    """
    
    def __init__(self, name, Re=-1):
        if not "__" in name:
            raise NameError("A string de entrada nao respeita o formato necessário.")
        self.name = name
        self.name1, self.name2, self.pct = name.split("__")
        self.Re = Re
    
        self.pct = float(self.pct)/100
        
        if self.pct > 1:
            raise ValueError("Fator de interpolação de %.2f (> 1) não permitido." %(self.pct))
        
        self._perfil1 = Perfil(self.name1,self.Re)
        self._perfil2 = Perfil(self.name2,self.Re)
        
        if (Re < 0) or ('circular' in self.name):    
            if Re > 0:
                self._saveCoordsToDat()
            
            #Verificação dos limites de ângulo de ataque tabelados
            self.limitesAng = self._dominioAng()
            
            # Calculo dos splines dos coeficientes aerodinâmicos
            self._funcaoCl, self._funcaoCd, self._funcaoCm = self._funcaoClCdCm()
            
            #Lista de angulos de ataque em que os coeficientes sao definidos
            self._alfas = np.arange(self.limitesAng[0],self.limitesAng[1]+0.5,0.5)
            
            # Valores fixos calculados
            self.Clmax, self.AlfaClmax = self._getClmax()
            self.Cdmin, self.AlfaCdmin = self._getCdmin()
            self.Cl_Cdmax, self.AlfaCl_Cdmax = self._getCl_Cdmax() 
            self.dCl_dAlfa, self.dCm_dAlfa = self._ajustelinear()
        else:
            if self.name1 == self.name2:
                self.name = self.name1
            else:
                self._saveCoordsToDat()
            
            self._perfil = Perfil(self.name,self.Re)
            
            #Verificação dos limites de ângulo de ataque tabelados
            self.limitesAng = self._dominioAng()
            self._perfil.limitesAng = self.limitesAng 
        
            # Calculo dos splines dos coeficientes aerodinâmicos
            self._perfil._funcaoCl, self._perfil._funcaoCd, self._perfil._funcaoCm = self._perfil._funcaoClCdCm()
            self._funcaoCl = self._perfil._funcaoCl
            self._funcaoCd = self._perfil._funcaoCd 
            self._funcaoCm = self._perfil._funcaoCm
            
            # Valores fixos calculados
            self._perfil.Clmax, self._perfil.AlfaClmax = self._perfil._readClmax()
            self.Clmax = self._perfil.Clmax
            self.AlfaClmax = self._perfil.AlfaClmax
            
            self._perfil.Cdmin, self._perfil.AlfaCdmin = self._perfil._readCdmin()
            self.Cdmin = self._perfil.Cdmin,
            self.AlfaCdmin = self._perfil.AlfaCdmin
            
            self._perfil.Cl_Cdmax, self._perfil.AlfaCl_Cdmax = self._perfil._readCl_Cdmax() 
            self.Cl_Cdmax = self._perfil.Cl_Cdmax
            self.AlfaCl_Cdmax = self._perfil.AlfaCl_Cdmax
            
            self._perfil.dCl_dAlfa, self._perfil.dCm_dAlfa = self._perfil._ajustelinear()
            self.dCl_dAlfa = self._perfil.dCl_dAlfa
            self.dCm_dAlfa = self._perfil.dCm_dAlfa
        
        # Centro aerodinâmico do perfil. Válido somente se o Cm for medido 
        # à 25% da corda.
        self.ac = 0.25-self.dCm_dAlfa/self.dCl_dAlfa
    
    def _funcaoClCdCm(self):
        fCl1, fCd1, fCm1 = self._perfil1._funcaoClCdCm()
        fCl2, fCd2, fCm2 = self._perfil2._funcaoClCdCm()
        
        fCl = lambda alfa: self.pct*fCl1(alfa)+(1-self.pct)*fCl2(alfa)
        fCd = lambda alfa: self.pct*fCd1(alfa)+(1-self.pct)*fCd2(alfa)
        fCm = lambda alfa: self.pct*fCm1(alfa)+(1-self.pct)*fCm2(alfa)
        
        return fCl, fCd, fCm
        
    def _dominioAng(self):
        min1,max1 = self._perfil1.limitesAng
        min2,max2 = self._perfil2.limitesAng
        return max([min1,min2]), min([max1,max2])
    
    def _getClmax(self):
        Clmax = 0
        angulo_clmax = 0
        
        for alfa in self._alfas: 
            Cl = self._funcaoCl(alfa)
            if(Cl>Clmax):
                Clmax = Cl
                angulo_clmax = alfa
                
        return Clmax, angulo_clmax
    
    def _getCdmin(self):
        Cdmin = 100
        angulo_cdmin = 0
        
        for alfa in self._alfas: 
            Cd = self._funcaoCd(alfa)
            if(Cd<Cdmin):
                Cdmin = Cd
                angulo_cdmin = alfa
                
        return Cdmin, angulo_cdmin
        
    def _getCl_Cdmax(self):
        Cl_Cdmax = 0
        angulo_Cl_Cdmax = 0
        
        for alfa in self._alfas:
            Cl = self._funcaoCl(alfa)
            Cd = self._funcaoCd(alfa)
            if(Cl/Cd>Cl_Cdmax):
                Cl_Cdmax = Cl/Cd
                angulo_Cl_Cdmax = alfa
                    
        return Cl_Cdmax, angulo_Cl_Cdmax
    
    def coords_padrao(self,n=50):
        list_x1,list_y1 = self._perfil1.coords_padrao(n)
        list_x2, list_y2 = self._perfil2.coords_padrao(n) #O parametro list_x por padrão é o mesmo

        list_y = self.pct*list_y1+(1-self.pct)*list_y2
        list_x = self.pct*list_x1+(1-self.pct)*list_x2
       
        return list_x, list_y
    
    def _saveCoordsToDat(self,n = 100):
        nome = self.name        
        x,y = self.coords_padrao(n)
        
        filename = 'Perfis/'+nome+'.dat'
        if not os.path.isfile(filename):
            print("\nEscrevendo coordenadas do perfil %s..." %(self.name))
            myFile = open(filename,'w')
            myFile.write(nome+'\n')
            for i in range(len(x)):
                myFile.write("%f %f\n" %(x[i],y[i]))
            myFile.close()
                
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    
    name = 'NACA 64618'
    p1 = Perfil(name,Re=3e6)
    
    x1,y1 = p1.coords_padrao(100)
    
    w, h = plt.figaspect((max(y1)-min(y1))/(max(x1)-min(x1)))
    plt.figure(figsize=(w,h))
    plt.grid(True)
    plt.plot(x1,y1,'k-x')
    plt.xlim([-0.01,1.01])
    plt.tight_layout()
    plt.show()
    
#    apoio.reset_dados_perfis()
    
        
        
        
        
        
        
        
        