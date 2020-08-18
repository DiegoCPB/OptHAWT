# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 18:56:49 2018

@author: Diego Chou
"""
import numpy as np
import perfil
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

def geometria(pa):
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    
    xc = alinhamento(pa)[0]
    wl = len(xc)//4
    if wl%2 == 0: wl+=1
    align = savgol_filter(xc,wl,2)
    
    cordas = pa.cordas
    betas = pa.betas
    betas_rad = np.deg2rad(betas)
    perfis = pa.coords
    
    x = pa.x
    yba = -align*cordas*np.sin(betas_rad)
    zba = align*cordas*np.cos(betas_rad) 
    ybf = (1-align)*cordas*np.sin(betas_rad) 
    zbf = -(1-align)*cordas*np.cos(betas_rad)
    
    fig = plt.figure(figsize=plt.figaspect(0.4)*2)
#    fig.suptitle('Aspectos geométricos')
    
    ax1 = fig.add_subplot(2, 1, 1)
    ax1.plot(x, cordas, 'b')
    ax1.grid(True)
    ax1.set_ylabel(r'Corda ($m$)')
    
    ax2 = fig.add_subplot(2, 1, 2)
    ax2.plot(x,betas,'r')
    ax2.grid(True)
    ax2.set_ylabel(r'$\beta$ (graus)')
    ax2.set_xlabel(r'Raio ($m$)')
    
    fig2 = plt.figure(figsize=[11.69,8.27]) #folha A4
    ax3 = fig2.add_subplot(1, 1, 1, projection='3d',aspect='equal')
    ax3.azim = -135
    
    for i in range(pa.N):
        ang = i*360/pa.N
        xr1,zbar = rot(x,zba,ang)
        xr2,zbfr = rot(x,zbf,ang)
    
        ax3.plot(xr1,yba,zbar,'k')
        ax3.plot(xr2,ybf,zbfr,'k')
        ax3.plot(x,np.zeros(len(x)),'k--')
        
        for index,item in perfis.items():
            yp,zp = item
            xp = np.ones(len(yp))*x[index]
            yp,zp = yp*cordas[index],zp*cordas[index]
            yp,zp = yp-align[index]*cordas[index],zp
            yp,zp = rot(yp,zp,betas[index]-90)
            xp,zp = rot(np.ones(len(yp))*x[index],zp,ang)
            ax3.plot(xp,yp,zp,'k')
    
    ax3.set_xlim(ax3.get_xlim()-np.mean(ax3.get_xlim()))   
    ax3.set_ylim(ax3.get_xlim())
    ax3.set_zlim(ax3.get_xlim())    
    ax3.invert_xaxis()
    
    plt.show()
    
def planta(pa, secoes = 6):
    import matplotlib.pyplot as plt
    from scipy.interpolate import InterpolatedUnivariateSpline
    import perfil
    
    centroide = alinhamento(pa)[0]
    wl = len(centroide)//4
    if wl%2 == 0: wl+=1
    align = savgol_filter(centroide,wl,2)
    
    cordas = pa.cordas
    x = pa.x
    s_cordas = InterpolatedUnivariateSpline(x,cordas,k=1)
    perfis = pa.coords
    
    yba = align*cordas
    s_yba = InterpolatedUnivariateSpline(x,yba,k=1)    
    ybf = (align-1)*cordas
    
    
    plt.figure(figsize=[11.69,8.27]) #folha A4
#    plt.title("PLANTA")
    plt.axis("equal")
    plt.grid(True)
    
    for index,item in perfis.items():
        xp,yp = item
        xp,yp = xp*cordas[index],yp*cordas[index]
        xp,yp = rot(xp,yp,-90)
        xp += x[index]
        yp += yba[index]
        plt.plot(xp,yp,'#D3D3D3')
            
    plt.plot(x,yba,'k')
    plt.plot(x,ybf,'k')
    plt.show()
            
    pos_text = np.linspace(pa.R0,pa.R,len(pa.lista_perfis))
    for i in range(len(pa.lista_perfis)):
        x_abs = pa.R0+ (pa.R-pa.R0)*pa.pos_perfis[i]
        yba_abs = s_yba(x_abs)
        corda = s_cordas(x_abs)
        
        p1 = perfil.Perfil(pa.lista_perfis[i])
        xp,yp = p1.coords_padrao()
        xp,yp = xp*corda,yp*corda
        xp,yp = rot(xp,yp,-90)
        xp += x_abs 
        yp += s_yba(x_abs)
        plt.plot(xp,yp,'b')
        plt.annotate(p1.name, xy=(x_abs, yba_abs), 
                     xytext=(pos_text[i], 1.2*max(pa.cordas)), ha='center',
                     bbox=dict(boxstyle="round4", fc="w"),
                     arrowprops=dict(arrowstyle="->"))
        
    plt.plot(pa.x,yba-centroide*cordas, 'gx-', label = 'Centróide')
    
    plt.legend()
    plt.show()
    
def desempenho(pa, lambda_bound=[2,20]):
    import Pa_BEMT
    import matplotlib.pyplot as plt
    
    try:
        v = pa.v_inf
    except:
        ValueError("A pá ainda nao foi calculada.")
    
    R = pa.R
    N = pa.N
    lista_perfis = pa.lista_perfis
    pos_perfis = pa.pos_perfis
    ar_rho = pa.rho
    
    lamb = np.arange(lambda_bound[0],lambda_bound[1]+0.5,0.5)
    b1 = Pa_BEMT.Pa_generica(N,lista_perfis,pos_perfis,pa.x,pa.cordas,pa.betas,rho=ar_rho)
    b2 = Pa_BEMT.Pa_generica(N,lista_perfis,pos_perfis,pa.x,pa.cordas,pa.betas,rho=ar_rho,tipo='xfoil')
    
    Cp1 = [0]
    Cp2 = [0]
    
    for blade in [b1,b2]:
        if blade == b1:
            Cp = Cp1
        else:
            Cp = Cp2
        for l in lamb:
            blade.calcular(v,l*v/R)
#            if blade.Cp < 0 and Cp[-1] > blade.Cp: break
            Cp.append(blade.Cp)
            
    plt.figure()
#    plt.title("Coeficiente de Potência x Razão de velocidades")
    plt.grid(True)
    plt.plot(lamb[:len(Cp1)-1],Cp1[1:],label='Polares experimentais')
    plt.plot(lamb[:len(Cp2)-1],Cp2[1:],'--',label='Polares XFoil')
    plt.plot([5,6,7,8,9,10],[0.36160,0.44780,0.47774,0.48044,0.46870,0.44368], 'kv:', label='NREL 5MW')
    plt.xlabel(r'TSR ($\lambda$)')
    plt.ylabel(r'$C_p$')
    plt.legend()
    plt.show()

def potência(pa, v_bound=[2,20], maxP = 0, fixa = True, eff = 1):
    import Pa_BEMT
    import matplotlib.pyplot as plt
    
    try:
        w_ref = pa.w
        v_ref = pa.v_inf
    except:
        ValueError("A pá ainda nao foi calculada.")
        
    R = pa.R
    N = pa.N
    lista_perfis = pa.lista_perfis
    pos_perfis = pa.pos_perfis
    ar_rho = pa.rho
    lamb = w_ref*R/v_ref
    
    vel = np.arange(v_bound[0],v_bound[1])
    
    b1 = Pa_BEMT.Pa_generica(N,lista_perfis,pos_perfis,pa.x,pa.cordas,pa.betas,rho=ar_rho)
    P1 = []
    
    b2 = Pa_BEMT.Pa_generica(N,lista_perfis,pos_perfis,pa.x,pa.cordas,pa.betas,rho=ar_rho,tipo='xfoil')
    P2 = []
    
    for blade in [b1,b2]:
        if blade == b1:
            P = P1
        else:
            P = P2
        for v in vel:
            if fixa == True:
                w = w_ref
            else:
                w = lamb*v/R
            blade.calcular(v,w)
            if blade.P > maxP/eff and maxP != 0:
                P.append(maxP/eff)
            else:
                P.append(blade.P)
     
    if fixa:
        str_title = r"(rotação FIXA, $\omega = %.1f rpm$, $\varepsilon = %.1f$)" %(w_ref*60/(2*np.pi),eff)
    else:
        str_title = r"(rotação VARIÁVEL, $\lambda = %.1f$, $\varepsilon = %.1f$)" %(lamb,eff)
    
#    Pv = np.array([0.004955334066718109, 0.0006431496772627554, -0.0004393696404250136, -0.0015218889581092299, -0.002604408275793446, 0.02538005805245902, 0.09857983538548964, 0.17177961271852027, 0.6002425479462907, 1.3678203157099453, 2.071262044559898, 2.8608363265251953, 3.6252950668663573, 4.447730917573455, 5.2181107105461955, 6.0224450574938935, 6.8287110687825034, 7.646981446076154, 8.456370244422438, 9.232185127891494, 10.032579291972056, 10.801725668994948, 11.576793098649414, 12.345939475672306, 13.11672837659497, 13.907670413126743, 14.750583920926715, 15.555997428726686, 16.374141473048013, 17.18315736457739, 17.992173256106767, 18.562467693408788, 19.019997614282488, 19.43231222415141, 19.60563161870932, 19.766032352980154, 19.855380455672037, 19.893053917215596, 19.927497713687387, 19.677730983843386, 19.28262932576972, 18.877838672480742, 18.469818354119997, 18.07148703097456, 17.73451934419276, 17.484752614348757, 17.244674879720073, 16.994908149876068, 16.725763429601447, 16.4049440681785, 15.922641453167033, 15.414501517581407, 15.119519476732624, 14.96664269904173, 14.80407692613553, 14.599525507296313, 14.378825763098245, 14.187193004546112, 14.031086561783448, 13.882457554109394])*10**3
#    vv = np.array([0.10835398746172942, 0.48940078728677694, 0.870447587111824, 1.2514943869368715, 1.632541186761919, 2.013587986586966, 2.3946347864120137, 2.775681586237061, 3.1567283860621087, 3.5377751858871562, 3.88418136754629, 4.126665694707683, 4.317189094620208, 4.525032803615688, 4.7155562035282115, 4.923399912523692, 5.165884239685085, 5.425688875929436, 5.685493512173786, 5.927977839335181, 6.153141857413617, 6.360985566409099, 6.568829275404578, 6.776672984400058, 7.001837002478496, 7.2789619478058025, 7.5734072022160674, 7.867852456626331, 8.196938329202508, 8.543344510861642, 8.889750692520776, 9.253477183262866, 9.634523983087915, 10.015570782912961, 10.396617582738008, 10.777664382563056, 11.158711182388103, 11.539757982213152, 11.920804782038198, 12.301851581863247, 12.682898381688293, 13.06394518151334, 13.444991981338388, 13.826038781163435, 14.207085580988482, 14.588132380813532, 14.969179180638578, 15.350225980463625, 15.731272780288672, 16.11231958011372, 16.49336637993877, 16.874413179763813, 17.255459979588863, 17.636506779413907, 18.017553579238957, 18.398600379064007, 18.77964717888905, 19.1606939787141, 19.541740778539143, 19.90546726928124])

    
    v_siemens = np.array([0,3,5,7.5,10,12.5,15,20])

    p_siemens = np.array([0, 0, 309.96, 915.13, 
                          2095.94, 3365.31, 3600, 3600])*10**3
        
    plt.figure()
    plt.title("Curva de Potência " + str_title)
    plt.grid(True)
    plt.xlabel(r'v ($m/s$)')
    plt.ylabel(r'P ($W$)')
#    plt.plot(vel,np.array(P1)*eff,label='Polares experimentais')
    plt.plot(vel,np.array(P2)*eff,label='Polares XFoil')
    plt.plot(v_siemens,p_siemens,'kx',label='Siemens SWT-3.6-120')
#    plt.plot(vv,Pv,'k--',label='Libellula 20 kW')
    plt.legend()
    plt.show()