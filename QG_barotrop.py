"""
Modelo cuasi-geostrófico de circulación superficial oceánica
forzada por el viento.
El código original está hecho en fortran, aquí lo reescribimos
para python 3

Dani risaro
Agosto 2019
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.insert(0, '/home/daniu/Documentos/materias/Métodos Numéricos/monografia_final/modelo_QG_python')
import func_externas
import params

# Note: the subscript a,b,c denotes time steps t+1,t,t-1 respectively.

directorio = '/home/daniu/Documentos/materias/Métodos Numéricos/monografia_final/modelo_QG_python/out_tmp/'

pi     = np.pi
zero   = 0.0
ele    = 2.*pi/(params.jm-1)
kx     = pi/(params.im-1)
dssqr  = params.ds**2
dssqri = 1./dssqr
dssqri4= dssqri/4.0
epsdt  = params.eps*params.dt
epsdti = 1./(1.+epsdt)
rober  = 1.e-3
nplts  = params.nlpt
c1     = epsdti
c2     = params.BFP*(params.dt/params.ds)*epsdti
c3     = 2.*params.dt*epsdti
c4     = epsdt*epsdti
c5     = 2.*params.dt*epsdti*params.Ro      #/3.0
c6     = 2.*params.dt*params.Ah*epsdti
c7     = 2.*params.dt*params.Bh*epsdti

bc1 = -(1.-params.gamma*params.ds*0.5)/(1.+params.gamma*params.ds*0.5)
bc2 =   1.-params.gamma*params.ds*0.5
bc3 = - 1./(1. + params.gamma*params.ds*0.5)

imm1 = params.im-1
imm2 = params.im-2
imm3 = params.im-3
imm4 = params.im-4
jmm1 = params.jm-1
jmm2 = params.jm-2
jmm3 = params.jm-3
jmm4 = params.jm-4

# number of time steps to be saved
xx = (params.nend-params.nst)/params.nlpt
nplots = int(xx) + 1
print('Time steps to be saved = ' + str(nplots))
ksf = 100
kv = 1000
kp = 0

# over-relaxation constants
# ncrit: # of steps allowed to do the relaxation in subroutine helm
# pcrit: criterium to stop the relaxation
# alfa: constant used for the sequential relaxation
# fxr:  constant used for the sequential relaxation

const1 = 1./(imm2-2)**2
const2 = 1./(jmm2-2.)**2
alfa = 2.-4.442883*np.sqrt(const1+const2)
fxr = dssqr/4.0

# initialize all arrays

nx = params.im
ny = params.jm
r = np.empty((nx, ny))
pa = np.empty((nx, ny))
pb = np.empty((nx, ny))
pc = np.empty((nx, ny))
psia = np.empty((nx, ny))
psib = np.empty((nx, ny))
psic = np.empty((nx, ny))
curlt = np.empty((nx, ny))
delpsi = np.empty((nx, ny))
delpsi4 = np.empty((nx, ny))

jpp = 0
jxp = 0
jpx = 0
jxx = 0

alpha = 1
beta = 0
gamma = 0
delta = 0

# wind stress curl
for i in range(params.im):
    for j in range(params.jm):
        if (params.GYR==1):
            curlt[i,j] = -params.HEM*np.sin(ele*0.5*(j-1))
        elif (params.GYR==2):
            curlt[i,j]= -np.sin(ele*0.5*(j-1))*np.sin(kx*(i-1))
        else:
            disp('Wind Gyre type not defined')
            break

np.save(directorio + 'QG_wind_stress', curlt)

#------------------------------------------------------------------------------
# Integracion temporal
#------------------------------------------------------------------------------

#for itime in range(params.nst,30):
for itime in range(100,200):

    print('step number = ' + str(itime))      # imprime cada 100 pasos

    for i in range(1,imm1):
        for j in range(1,jmm1):

            # horizontal mixing
            delpsi, delpsi4 = func_externas.horizontal_mixing(delpsi, delpsi4, psic, dssqri, i, j)

            # Arakawa's jacobian
            jpp, jxp, jpx, jxx = func_externas.arakawa_jacobian(jpp, jxp, jpx, jxx, pb, psib, dssqri4, i, j)

            # update vorticity
            psia = func_externas.vorticity(psia, psic, jpp, jxp, jpx, jxx, pb, curlt, delpsi, delpsi4, c1, c2, c3, c4, c5, c6, c7, alpha, beta, gamma, delta, i, j)

    # update stream function
    nrelax = zero

    while nrelax <= params.ncrit:
        # solve the Laplacian for the stream function by over-relaxation
        for i in range(2,imm2):
            for j in range(2,jmm2):
                r, pa = func_externas.laplacian_p(r, pa, psia, dssqri, alfa, fxr, i, j)

        # check convergence
        n1 = zero
        for i in range(2,imm2):
            for j in range (2,jmm2):
                rabs = np.sqrt(r[i,j]**2)
                if(rabs > params.pcrit):
                    n1 = n1+1

        if (n1!=0):
            nrelax = nrelax +1
        else:
            break

    if (nrelax > params.ncrit):
        print('Warning: The subroutine does not relax  ntime=' + str (itime))

    # set the boundary conditons on the stream function
    for j in range(1,jmm1):
        pa[0,j] = bc1*pa[2,j]
        pa[1,j] = zero
        pa[params.im-1,j] = bc1*pa[imm2,j]
        pa[imm1,j] = zero

    for i in range(1,imm1):
        pa[i,0] = bc1*pa[i,2]
        pa[i,1] = zero
        pa[i,params.jm-1] = bc1*pa[i,jmm2]
        pa[i,jmm1] = zero

    # diagonostic calculation of the vorticity on the walls
    for j in range(1,jmm1):
        psia[1,j] = (pa[2,j] + pa[0,j])*dssqri
        psia[imm1,j]= (pa[params.im-1,j] + pa[imm2,j])*dssqri

    for i in range(1,imm1):
        psia[i,1] = (pa[i,2] + pa[i,0])*dssqri
        psia[i,jmm1] = (pa[i,params.jm-1] + pa[i,jmm2])*dssqri

    # set the boundary conditons on the vorticity
    for j in range(1,jmm1):
        psia[0,j] = bc3*(bc2*psia[2,j]-4.*psia[1,j]+psia[1,j+1]+ psia[1,j-1])
        psia[params.im-1,j] = bc3*(bc2*psia[imm2,j]-4.*psia[imm1,j]+psia[imm1,j+1]+psia[imm1,j-1])

    for i in range(1,imm1):
        psia[i,0]  = bc3*(bc2*psia[i,2]-4.*psia[i,1]+psia[i+1,1]+psia[i-1,1])
        psia[i,params.jm-1] = bc3*(bc2*psia[i,jmm2]-4.*psia[i,jmm1]+psia[i+1,jmm1]+psia[i-1,jmm1])

    # time smoothing the stream function using a Robert's filter
    for i in range(params.im):
        for j in range(params.jm):
            pb[i,j] = pb[i,j]+rober*(pa[i,j]-2.*pb[i,j]+pc[i,j])
            psib[i,j] = psib[i,j]+rober*(psia[i,j]-2.*psib[i,j]+psic[i,j])

    # kinetic energy (entire formain)
    tke = 0.
    for i in range(1,imm1):
        for j in range(1,jmm1):
            gradpsi = ((pb[i+1,j]-pb[i-1,j])**2+(pb[i,j+1]-pb[i,j-1])**2)*dssqri
            tke = tke + 0.5*gradpsi

    for i in range(params.im):
        for j in range(params.jm):
            pc[i,j] = pb[i,j]
            pb[i,j] = pa[i,j]
            psic[i,j] = psib[i,j]
            psib[i,j] = psia[i,j]

#     if (itime >= nplts):
#         kp = kp + 1
#         ksf = ksf + kp
#         kv = kv + kp
#
#         num_file = kp
#
#         if (kp < 10):
#             idx = 3
#         elif (kp >= 10 and kp < 100):
#             idx = 2
#         else:
#             idx = 1
#
# #        np.save(directorio + '', curlt)
#         print('step number = ' + str(itime) + ' writing file = ' + str(num_file))
#
#         if (idx==3):
#             if (nplots<100):
#                 a = 5
#             else:
#                 a = 3
#
#         if (idx==2):
#             if (nplots<100):
#                 a = 5
#             else:
#                 a = 3
#
#         if (idx==1):
#             nombre_psi = directorio + 'psi' + str(int(idx/3))
#             nombre_vor = directorio + 'vor' + str(int(idx/3))
# #            np.save(nombre_psi, )
#
# #            open(ksf,file='./out_tmp/psi'//num_file(idx:3)//'.dat',form='formatted',status='unknown')
# #            open(kv,file='./out_tmp/vor'//num_file(idx:3)//'.dat',form='formatted',status='unknown')

#         nplts = nplts + params.nlpt
#
# #    else:
# #        print('step number = ' + str(itime))
