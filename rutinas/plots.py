#%%

#Abrimos librerias necesarias

import Oceano
import numpy as np

#%%

#Calculo escalas típicas del problema para las simulaciones

tau = 0.2         #Tension del viento [N/m^2]
L = 750000      # Longitud de la cuenca [m]
D = 2000        # Profundidad [m]
beta = 2e-11    # Coeficiente de Coriolis [1/(s*m)]
rho = 1025      # Densidad [kg/m^3]

A = 205         # Coeficiente de viscosidad lateral [m^2/s] (Lo pongo solo para correr esta función)
k = 1.16e-5     # Coeficiente de fricción de fondo [1/s] (Lo pongo solo para correr esta función)

# Los parámetros que cálculo acá no son los reales

U_sim1,Ro_sim1,Ef_sim1,Ev1_sim1,deltai_sim1,deltaf_sim1,deltav1_sim1,Wi_sim1,Wf_sim1,Wv1_sim1 = Oceano.parametros(tau,L,D,k,A,beta,rho)
U_sim2,Ro_sim2,Ef_sim2,Ev1_sim2,deltai_sim2,deltaf_sim2,deltav1_sim2,Wi_sim2,Wf_sim2,Wv1_sim2 = Oceano.parametros(tau,L,D,k,A,beta,rho)
U_sim3,Ro_sim2,Ef_sim3,Ev1_sim3,deltai_sim3,deltaf_sim3,deltav1_sim3,Wi_sim2,Wf_sim3,Wv1_sim3 = Oceano.parametros(tau,L,D,k,A,beta,rho)

#%%

#Abro los archivos correspondientes a las simulaciones
Lx = 750
Ly = 750
nx = 150
ny = 150

dir_salida_sim1 = '/home/daniu/Documentos/materias/Métodos Numéricos/monografia_final/modelo_QG/simulaciones/sim1/'
dir_salida_sim2 = '/home/daniu/Documentos/materias/Métodos Numéricos/monografia_final/modelo_QG/simulaciones/sim2/'
dir_salida_sim3 = '/home/daniu/Documentos/materias/Métodos Numéricos/monografia_final/modelo_QG/simulaciones/sim3/'

psi_temp_sim1,vort_temp_sim1,psiF_sim1,vortF_sim1,QG_diag_sim1,QG_curlw_sim1,X_sim1,Y_sim1,dx_sim1,dy_sim1 = Oceano.extractdata(dir_salida_sim1,Lx,Ly,nx,ny)
psi_temp_sim2,vort_temp_sim2,psiF_sim2,vortF_sim2,QG_diag_sim2,QG_curlw_sim2,X_sim2,Y_sim2,dx_sim2,dy_sim2 = Oceano.extractdata(dir_salida_sim2,Lx,Ly,nx,ny)
psi_temp_sim3,vort_temp_sim3,psiF_sim3,vortF_sim3,QG_diag_sim3,QG_curlw_sim2,X_sim3,Y_sim3,dx_sim3,dy_sim3 = Oceano.extractdata(dir_salida_sim3,Lx,Ly,nx,ny)

#%%
#Gráfico energía cinética en el centro del dominio

ecin = np.stack((QG_diag_sim1, QG_diag_sim2, QG_diag_sim3), axis=0)
Nom_sim = ['sim1','sim2','sim3']
dir_graf = '/home/daniu/Documentos/materias/Métodos Numéricos/monografia_final/modelo_QG/figuras/Ecin.png'

Oceano.graf_ecin(ecin,Nom_sim, dir_graf)

#%%

ecin = np.stack((QG_diag_sim1, QG_diag_sim2, QG_diag_sim3), axis=0)

TiempoEst = Oceano.Calc_TiempoEst(ecin)

#%%

#Gráfico función corriente de todas las simulaciones

psiF_todos = np.stack((psiF_sim1, psiF_sim2, psiF_sim3), axis=0)
X_todos = np.stack((X_sim1,X_sim2,X_sim3), axis=0)
Y_todos = np.stack((Y_sim1,Y_sim2,Y_sim3), axis=0)
U_todos = np.stack((U_sim1,U_sim2,U_sim3), axis=0)
L_todos = np.stack((L, L, L), axis=0)
clevs = [-30,-28,-26,-24,-22,-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,-1,-0.5,0]
dir_graf = '/home/daniu/Documentos/materias/Métodos Numéricos/monografia_final/modelo_QG/figuras/psiF.png'

Oceano.graf_psiF(psiF_todos, X_todos, Y_todos, U_todos, L_todos, clevs, Nom_sim, dir_graf)

#%%

#Gráfico vorticidad relativa de todas las simulaciones

vortF_todos = np.stack((vortF_sim1, vortF_sim2, vortF_sim3), axis=0)
clevs = [-50,-20,-10,-5,-2,-0.5,0,0.5,2,5,10,20,50]
dir_graf = '/home/daniu/Documentos/materias/Métodos Numéricos/monografia_final/modelo_QG/figuras/vortF.png'

Oceano.graf_vortF(vortF_todos, X_todos, Y_todos, U_todos, L_todos, clevs, Nom_sim, dir_graf)

#%%
#
# #Gráfico vorticidad total de todas las simulaciones
#
# vortF_todos = np.stack((psiF_sim1,psiF_sim2), axis=0)
# X_todos = np.stack((X_sim1,X_sim2), axis=0)
# Y_todos = np.stack((Y_sim1,Y_sim2), axis=0)
# U_todos = np.stack((U_sim1,U_sim2), axis=0)
# L_todos = np.stack((L,L), axis=0)
# clevs = np.arange(-8, -5.2,  0.1)
# lat = -30
# Nom_sim = ['sim1','sim2']
# dir_graf = '/home/daniu/Documentos/materias/Métodos Numéricos/monografia_final/modelo_QG/figuras/vortFf.png'
#
# Oceano.graf_vorttotalF(vortF_todos, X_todos, Y_todos, U_todos, L_todos, clevs, lat, Nom_sim, dir_graf)

#%%

#Calculo del parámetro de la corriente de borde oeste por efectos inerciales a partir del número de Roosby

Ro_sim1 = 0.005 # Número de Rossby
Ro_sim2 = 0.005
Ro_sim3 = 0.005

Ro = [Ro_sim1, Ro_sim2, Ro_sim3]

Ef_sim1 = 0.01 # Número de Ekman vertical
Ef_sim2 = 0.01
Ef_sim3 = 0.01

Ef = [Ef_sim1, Ef_sim2, Ef_sim3]

Ev1_sim1 = 0.001 # Número de Ekman horizontal
Ev1_sim2 = 0.001
Ev1_sim3 = 0.001

Ev1 = [Ev1_sim1, Ev1_sim2, Ev1_sim3]

# Coc_delta_i_f, Coc_delta_i_v1 = Oceano.Calc_Coc_delta(Ro, Ef, Ev1)

#%%
