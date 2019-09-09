#%%
# Abrimos librerias necesarias

import Oceano
import numpy as np

#%%
# Calculo escalas típicas del problema para las simulaciones

tau = 0.2       # Tension del viento [N/m^2]
L = 1500000     # Longitud de la cuenca [m]
D = 2000        # Profundidad [m]
beta = 2e-11    # Coeficiente de Coriolis [1/(s*m)]
rho = 1025      # Densidad [kg/m^3]

A = 205         # Coeficiente de viscosidad lateral [m^2/s] (Lo pongo solo para correr esta función)
k = 1.16e-5     # Coeficiente de fricción de fondo [1/s] (Lo pongo solo para correr esta función)

# Los parámetros que cálculo acá no son los reales

U_sim1,Ro_sim1,Ef_sim1,Ev1_sim1,deltai_sim1,deltaf_sim1,deltav1_sim1,Wi_sim1,Wf_sim1,Wv1_sim1 = Oceano.parametros(tau,L,D,k,A,beta,rho)
U_sim2,Ro_sim2,Ef_sim2,Ev1_sim2,deltai_sim2,deltaf_sim2,deltav1_sim2,Wi_sim2,Wf_sim2,Wv1_sim2 = Oceano.parametros(tau,L,D,k,A,beta,rho)
U_sim3,Ro_sim2,Ef_sim3,Ev1_sim3,deltai_sim3,deltaf_sim3,deltav1_sim3,Wi_sim2,Wf_sim3,Wv1_sim3 = Oceano.parametros(tau,L,D,k,A,beta,rho)
U_sim4,Ro_sim2,Ef_sim4,Ev1_sim4,deltai_sim4,deltaf_sim4,deltav1_sim4,Wi_sim2,Wf_sim4,Wv1_sim4 = Oceano.parametros(tau,L,D,k,A,beta,rho)
U_sim5,Ro_sim2,Ef_sim5,Ev1_sim5,deltai_sim5,deltaf_sim5,deltav1_sim5,Wi_sim2,Wf_sim5,Wv1_sim5 = Oceano.parametros(tau,L,D,k,A,beta,rho)
U_sim6,Ro_sim2,Ef_sim6,Ev1_sim6,deltai_sim6,deltaf_sim6,deltav1_sim6,Wi_sim2,Wf_sim6,Wv1_sim6 = Oceano.parametros(tau,L,D,k,A,beta,rho)
U_sim7,Ro_sim2,Ef_sim7,Ev1_sim7,deltai_sim7,deltaf_sim7,deltav1_sim7,Wi_sim2,Wf_sim7,Wv1_sim7 = Oceano.parametros(tau,L,D,k,A,beta,rho)

#%%
# Abro los archivos correspondientes a las simulaciones
Lx = 1500
Ly = 1500
nx = 150
ny = 150

dir_salida_sim1 = '/home/daniu/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/simulaciones/sim1/out_tmp/'
dir_salida_sim2 = '/home/daniu/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/simulaciones/sim2/out_tmp/'
dir_salida_sim3 = '/home/daniu/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/simulaciones/sim3/out_tmp/'
dir_salida_sim4 = '/home/daniu/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/simulaciones/sim4/out_tmp/'
dir_salida_sim5 = '/home/daniu/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/simulaciones/sim5/out_tmp/'
dir_salida_sim6 = '/home/daniu/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/simulaciones/sim6/out_tmp/'
dir_salida_sim7 = '/home/daniu/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/simulaciones/sim7/out_tmp/'

psi_temp_sim1,vort_temp_sim1,psiF_sim1,vortF_sim1,QG_diag_sim1,QG_curlw_sim1,X_sim1,Y_sim1,dx_sim1,dy_sim1 = Oceano.extractdata(dir_salida_sim1,Lx,Ly,nx,ny)
psi_temp_sim2,vort_temp_sim2,psiF_sim2,vortF_sim2,QG_diag_sim2,QG_curlw_sim2,X_sim2,Y_sim2,dx_sim2,dy_sim2 = Oceano.extractdata(dir_salida_sim2,Lx,Ly,nx,ny)
psi_temp_sim3,vort_temp_sim3,psiF_sim3,vortF_sim3,QG_diag_sim3,QG_curlw_sim3,X_sim3,Y_sim3,dx_sim3,dy_sim3 = Oceano.extractdata(dir_salida_sim3,Lx,Ly,nx,ny)
psi_temp_sim4,vort_temp_sim4,psiF_sim4,vortF_sim4,QG_diag_sim4,QG_curlw_sim4,X_sim4,Y_sim4,dx_sim4,dy_sim4 = Oceano.extractdata(dir_salida_sim4,Lx,Ly,nx,ny)
psi_temp_sim5,vort_temp_sim5,psiF_sim5,vortF_sim5,QG_diag_sim5,QG_curlw_sim5,X_sim5,Y_sim5,dx_sim5,dy_sim5 = Oceano.extractdata(dir_salida_sim5,Lx,Ly,nx,ny)
psi_temp_sim6,vort_temp_sim6,psiF_sim6,vortF_sim6,QG_diag_sim6,QG_curlw_sim6,X_sim6,Y_sim6,dx_sim6,dy_sim6 = Oceano.extractdata(dir_salida_sim6,Lx,Ly,nx,ny)
psi_temp_sim7,vort_temp_sim7,psiF_sim7,vortF_sim7,QG_diag_sim7,QG_curlw_sim7,X_sim7,Y_sim7,dx_sim7,dy_sim7 = Oceano.extractdata(dir_salida_sim7,Lx,Ly,nx,ny)

#%%
# Gráfico energía cinética en el centro del dominio

ecin = np.stack((QG_diag_sim1, QG_diag_sim2, QG_diag_sim3, QG_diag_sim4, QG_diag_sim5, QG_diag_sim6, QG_diag_sim7), axis=0)
Nom_sim = ['sim1','sim2','sim3','sim4','sim5','sim6','sim7']
dir_graf = '/home/daniu/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/figuras/Ecin.png'

Oceano.graf_ecin(ecin, Nom_sim, dir_graf)

#%%

ecin = np.stack((QG_diag_sim1, QG_diag_sim2, QG_diag_sim3, QG_diag_sim4, QG_diag_sim5, QG_diag_sim6, QG_diag_sim7), axis=0)

TiempoEst = Oceano.Calc_TiempoEst(ecin)

#%%

# Gráfico función corriente de todas las simulaciones
psiF_todos = np.stack((psiF_sim1, psiF_sim2, psiF_sim3, psiF_sim4, psiF_sim5, psiF_sim6, psiF_sim7), axis=0)
X_todos = np.stack((X_sim1,X_sim2,X_sim3,X_sim4,X_sim5,X_sim6,X_sim7), axis=0)
Y_todos = np.stack((Y_sim1,Y_sim2,Y_sim3,Y_sim4,Y_sim5,Y_sim6,Y_sim7), axis=0)
U_todos = np.stack((U_sim1,U_sim2,U_sim3,U_sim4,U_sim5,U_sim6,U_sim7), axis=0)
L_todos = np.stack((L, L, L, L, L, L, L), axis=0)
clevs = [-30,-28,-26,-24,-22,-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,-1,-0.5,0]
dir_graf = '/home/daniu/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/figuras/psiF.png'

Oceano.graf_psiF(psiF_todos, X_todos, Y_todos, U_todos, L_todos, clevs, Nom_sim, dir_graf)

#%%

# Gráfico vorticidad relativa de todas las simulaciones
vortF_todos = np.stack((vortF_sim1, vortF_sim2, vortF_sim3, vortF_sim4, vortF_sim5, vortF_sim6, vortF_sim7), axis=0)
clevs = [-50,-20,-10,-5,-2,-0.5,0,0.5,2,5,10,20,50]
dir_graf = '/home/daniu/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/figuras/vortF.png'

Oceano.graf_vortF(vortF_todos, X_todos, Y_todos, U_todos, L_todos, clevs, Nom_sim, dir_graf)

#%%
