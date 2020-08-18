import numpy as np

import Pa_BEMT
import graficos

R0 = 2.8667
R = 61.633
B = 3
ar_rho = 1.225

lamb = 7.69
v = 8

x_nodes = np.array([R0,5.6,8.3333,11.75,15.85,19.95,24.05,
                    28.15,32.25,36.35,40.45,44.55,
                    48.65,52.75,56.1667,58.9,R])

c_nodes = np.array([3.542,3.854,4.167,4.557,4.652,4.458,4.249,
                    4.007,3.748,3.502,3.256,3.010,
                    2.764,2.518,2.313,2.086,1.419])

beta_nodes = np.array([13.308,13.308,13.308,13.308,11.480,10.162,9.011,
                       7.795,6.544,5.361,4.188,3.125,
                       2.319,1.526,0.863,0.370,0.106])

lista_perfis = ['circular',
                'circular',
                'circular',
                'DU 93W405 AD',
                'DU 93W350 AD',
                'DU 93W350 AD',
                'DU 93W300 LM',
                'DU 93W250 LM',
                'DU 93W250 LM',
                'DU 93W210 LM',
                'DU 93W210 LM',
                'NACA 64618',
                'NACA 64618',
                'NACA 64618',
                'NACA 64618',
                'NACA 64618',
                'NACA 64618']

pos_perfis = (x_nodes-R0)/(R-R0)

b = Pa_BEMT.Pa_generica(B,lista_perfis,pos_perfis,x_nodes,c_nodes,beta_nodes)
b.calcular(v,lamb*v/R)

#graficos.geometria(b)
#graficos.planta(b)
graficos.desempenho(b, lambda_bound=[5,10])

#import macro_CATIA as mc
#mc.exportarPerfis(b,"NREL_5MW")
