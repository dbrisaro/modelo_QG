#%%
"""
Funciones para trabajar las salidas del modelo QG_barotrop.f
Adaptadas de Anthony Schrapffer
Leandro Díaz
2018
"""

#%%

"""
Extraccion de los datos de salida del modelo

INPUTS
dir_salida: str, direccion de la salida
Lx: float, tamano de la cuenca (direccion X)
Ly: float, tamano de la cuenca (direccion Y)
nx: int, numero de punto de grilla (direccion X)
ny: int, numero de punto de grilla (direccion Y)

OUTPUTS
psi_temp: Campos de función de corriente de todos los tiempos
vort_temp: Campos de vorticidada de todos los tiempos
psiF: Campo de función de corriente del tiempo final
vortF: Campo de vorticidad del tiempo final
QG_diag: Información temporal de función corriente, vorticidad y energía cinética en el punto central del dominio
QG_curlw: Campo del rotor del esfuerzo del viento utilizado en la simulación
X: Vector con los puntos del eje X dimensionalizado
Y: Vector con los puntos del eje Y dimensionalizado
dx: Distancia entre puntos del eje X
dy: Distancia entre puntos del eje y

"""
def extractdata(dir_salida,Lx,Ly,nx,ny):

    #Cargamos las librerias necesarias
    import os
    import numpy as np

    archivos=os.listdir(dir_salida) # nombre de los archivos en el directorio actual + \output

    tiempos=0

    for name in archivos:
        if name[0:3]=='psi':
            tiempos=tiempos+1

    # Creacion de las matriz para recibir los datos
    psi_temp=np.empty(shape=[ny+2,nx+2,tiempos])
    psi_temp[:]=np.nan
    vort_temp=np.empty(shape=[ny+2,nx+2,tiempos])
    vort_temp[:]=np.nan

    # Extraccion
    for name in archivos:
        if name[0:3]=='psi':
            k1=name[3]+name[4]
            psi_temp[:,:,int(k1)-1]=np.loadtxt(dir_salida+name) #fromfile(name)
        if name[0:3]=='vor':
            k2=name[3]+name[4]
            vort_temp[:,:,int(k2)-1]=np.loadtxt(dir_salida+name) #fromfile(name)
        if name[0:7]=='QG_diag':
            QG_diag=np.loadtxt(dir_salida+name) #fromfile(name)
            # En el punto central del dominio
            # (Tiempo, Funcion Corriente,Vorticidad, EnCin)
        if name[0:7]=='QG_wind':
            QG_curlw=np.loadtxt(dir_salida+name) #fromfile(name)

    # Recorte de los datos
    a1=1; a2=np.size(vort_temp,0)-1;
    b1=1; b2=np.size(vort_temp,1)-1;

    psi_temp=psi_temp[a1:a2,b1:b2,:]
    psiF=psi_temp[:,:,np.size(psi_temp,2)-1]

    vort_temp=vort_temp[a1:a2,b1:b2,:]
    vortF=vort_temp[:,:,np.size(vort_temp,2)-1]

    X=np.linspace(0,Lx,num=nx)
    Y=np.linspace(0,Ly,num=ny)

    dx=Lx/(nx-1);
    dy=Ly/(ny-1);

    return psi_temp,vort_temp,psiF,vortF,QG_diag,QG_curlw,X,Y,dx,dy

#%%

"""
Calculo de los parametros para la simulación

INPUTS
tau: Tension del viento [N/m^2]
L: Longitud de la cuenca [m]
D: Profundidad [m]
k: Coeficiente de fricción de fondo [1/s]
A: Coeficiente de viscosidad lateral [m^2/s]
beta: Coeficiente de Coriolis [1/(s*m)]
rho: Densidad [kg/m^3]

OUTPUTS
U: Escala típica velocidad [m/s]
Ro: Número de Rossby
Ef: Número de Ekman vertical
Ev1: Número de Ekman horizontal
deltai: Parámetro de la corriente de borde oeste por efectos inerciales
Wi: Ancho de la corriente de borde oeste por efectos inerciales  [m]
deltaf: Parámetro de la corriente de borde oeste por fricción de fondo
Wf: Ancho de la corriente de borde oeste por fricción de fondo [m]
deltav1: Parámetro de la corriente de borde oeste por fricción lateral
Wv1: Ancho de la corriente de borde oeste por fricción lateral [m]

"""
def parametros(tau,L,D,k,A,beta,rho):

    #Cargamos las librerias necesarias
    import math

    U=(2*math.pi*tau)/(rho*D*beta*L) # Velocidad
    Ro=(2*math.pi*tau)/(rho*D*(math.pow(beta,2)*(math.pow(L,3)))) # Numero de Rossby
    Ef=k/(beta*L)  # Numero de Ekman vertical
    Ev1=A/(beta*(math.pow(L,3))) # Numero de Ekman Horizontal
    deltai=math.sqrt(Ro) # Efectos inerciales
    deltaf=Ef
    deltav1=math.pow(Ev1,1/3)
    Wi=deltai*L
    Wf=deltaf*L
    Wv1=deltav1*L

    return U,Ro,Ef,Ev1,deltai,deltaf,deltav1,Wi,Wf,Wv1

#%%

"""
Gráfico de energía cinética

INPUTS
ecin: array con QG_DIAG concatenado (en la primera dimensión) de las simulaciones
Nom_sim: Nombres de las simulaciones
dir_graf: str, direccion con nombre donde guardar el grafico

"""
def graf_ecin(ecin,Nom_sim,dir_graf):

    #Cargamos las librerias necesarias
    import numpy as np
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(5,3))
    plt.plot(np.squeeze(ecin[0,:,0]),np.squeeze(ecin[:,:,3]).T/100000,linewidth=2)
    plt.xlabel('Numero de iteraciones',fontsize=8)
    plt.ylabel('Energia cinetica $(*10^5)$',fontsize=8)
    plt.title('Energia cinetica en funcion del numero de iteraciones', fontsize=8.5)
    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)
    plt.legend(Nom_sim,fontsize=7,loc = 5)
    fig.subplots_adjust(bottom=0.15,top=0.85,left=0.1,right=0.97)
    plt.savefig(dir_graf,dpi=500)

#%%

"""
Gráfico de función corriente

INPUTS
psiF_todos: array con función corriente para el tiempo final concatenadas para las distintas simulaciones
X_todos:Coordenadas x concatenadas para las distintas simulaciones
Y_todos:Coordenadas y concatenadas para las distintas simulaciones
U_todos: Escala típica velocidad [m/s] concatenadas para las distintas simulaciones
L_todos: Longitud de la cuenca [m]
clevs: Niveles para graficar
Nom_sim: Nombres de las simulaciones
dir_graf: str, direccion con nombre donde guardar el grafico

"""
def graf_psiF(psiF_todos,X_todos,Y_todos,U_todos,L_todos,clevs,Nom_sim,dir_graf):

    #Cargamos las librerias necesarias
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import gridspec

    fig, ax = plt.subplots(figsize=(6*3*int(3*(np.size(psiF_todos,0)/3-(int(np.size(psiF_todos,0)-0.01)/3))),4*int((np.size(psiF_todos,0)+1)/3)))
    gs = gridspec.GridSpec(int((np.size(psiF_todos,0)+1)/3),3*int(3*(np.size(psiF_todos,0)/3-(int(np.size(psiF_todos,0)-0.01)/3))))

    for k in range(np.size(psiF_todos,0)):

        ax = plt.subplot(gs[k])
        im=ax.contourf(np.squeeze(X_todos[k]),np.squeeze(Y_todos[k]),np.squeeze(psiF_todos[k])*U_todos[k]*L_todos[k]/10000,clevs,cmap=plt.get_cmap('gist_earth_r'),extend='both')
        ax.contour(np.squeeze(X_todos[k]),np.squeeze(Y_todos[k]),np.squeeze(psiF_todos[k])*U_todos[k]*L_todos[k]/10000,clevs,colors='k',linewidths=0.4,linestyles='solid')
        cbar=plt.colorbar(im)
        cbar.ax.tick_params(labelsize=9,direction="in",length = 2,pad=2)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        plt.xlabel('Longitug (km)',fontsize=8)
        plt.ylabel('Longitug (km)',fontsize=8)
        plt.title(str(Nom_sim[k]),fontsize=12)
        fig.subplots_adjust(bottom=0.05,left=0.1,top=0.85,hspace=0.3, right =0.95)
        plt.suptitle('Función corriente $(*10^4 m^2/s)$', fontsize=14)

        plt.savefig(dir_graf,dpi=500)

#%%

"""
Gráfico de vorticidad

INPUTS
vortF_todos: array con función corriente para el tiempo final concatenadas para las distintas simulaciones
X_todos:Coordenadas x concatenadas para las distintas simulaciones
Y_todos:Coordenadas y concatenadas para las distintas simulaciones
U_todos: Escala típica velocidad [m/s] concatenadas para las distintas simulaciones
L_todos: Longitud de la cuenca [m]
clevs: Niveles para graficar
Nom_sim: Nombres de las simulaciones
dir_graf: str, direccion con nombre donde guardar el grafico

"""
def graf_vortF(vortF_todos,X_todos,Y_todos,U_todos,L_todos,clevs,Nom_sim,dir_graf):

    #Cargamos las librerias necesarias
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import gridspec

    fig, ax = plt.subplots(figsize=(6*3*int(3*(np.size(vortF_todos,0)/3-(int(np.size(vortF_todos,0)-0.01)/3))),4*int((np.size(vortF_todos,0)+1)/3)))
    gs = gridspec.GridSpec(int((np.size(vortF_todos,0)+1)/3),3*int(3*(np.size(vortF_todos,0)/3-(int(np.size(vortF_todos,0)-0.01)/3))))

    for k in range(np.size(vortF_todos,0)):

        ax = plt.subplot(gs[k])
        im=ax.contourf(np.squeeze(X_todos[k]),np.squeeze(Y_todos[k]),np.squeeze(vortF_todos[k])*U_todos[k]/L_todos[k]*10000000,clevs,cmap=plt.get_cmap('seismic'),extend='both')
        ax.contour(np.squeeze(X_todos[k]),np.squeeze(Y_todos[k]),np.squeeze(vortF_todos[k])*U_todos[k]/L_todos[k]*10000000,clevs,colors='k',linewidths=0.4,linestyles='solid')
        cbar=plt.colorbar(im)
        cbar.ax.tick_params(labelsize=9,direction="in",length = 2,pad=2)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        plt.xlabel('Longitug (km)',fontsize=8)
        plt.ylabel('Longitug (km)',fontsize=8)
        plt.title(str(Nom_sim[k]),fontsize=12)
        fig.subplots_adjust(bottom=0.05,left=0.1,top=0.85,hspace=0.3, right =0.95)
        plt.suptitle('Vorticidad Relativa $(*10^{-7} s^{-1})$', fontsize=14)

        plt.savefig(dir_graf,dpi=500)


#%%

"""
Gráfico de transporte meridional

INPUTS
psiF_todos: array con función corriente para el tiempo final concatenadas para las distintas simulaciones
X_todos:Coordenadas x concatenadas para las distintas simulaciones
Y_todos:Coordenadas y concatenadas para las distintas simulaciones
U_todos: Escala típica velocidad [m/s] concatenadas para las distintas simulaciones
L_todos: Longitud de la cuenca [m]
D_todos: Profundidad de la cuenca [m]
Lx: float, tamano de la cuenca (direccion X)
nx: int, numero de punto de grilla (direccion X)
clevs: Niveles para graficar
Nom_sim: Nombres de las simulaciones
dir_graf: str, direccion con nombre donde guardar el grafico

"""
def graf_TrasMer(psiF_todos,X_todos,Y_todos,U_todos,L_todos,D_todos,Lx,nx,clevs,Nom_sim,dir_graf):

    #Cargamos las librerias necesarias
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import gridspec

    #Calculamos transporte meridional promediado en la vertical (es la derivada zonal de la función corriente multiplicada por la profundidad)

    fig, ax = plt.subplots(figsize=(6*3*int(3*(np.size(psiF_todos,0)/3-(int(np.size(psiF_todos,0)-0.01)/3))),4*int((np.size(psiF_todos,0)+1)/3)))
    gs = gridspec.GridSpec(int((np.size(psiF_todos,0)+1)/3),3*int(3*(np.size(psiF_todos,0)/3-(int(np.size(psiF_todos,0)-0.01)/3))))

    for k in range(np.size(psiF_todos,0)):
        TrasMer=np.diff(np.squeeze(psiF_todos[k]),n=1,axis=1)*D_todos[k]
        # Reticula alternativa
        Xalt=np.linspace(0,Lx,num=nx-1)
        ax = plt.subplot(gs[k])
        im=ax.contourf(Xalt,np.squeeze(Y_todos[k]),TrasMer*U_todos[k]*L_todos[k]/1000000,clevs,cmap=plt.get_cmap('seismic'),extend='both')
        cbar=plt.colorbar(im)
        cbar.ax.tick_params(labelsize=9,direction="in",length = 2,pad=2)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        plt.xlabel('Longitug (km)',fontsize=8)
        plt.ylabel('Longitug (km)',fontsize=8)
        plt.title(str(Nom_sim[k]),fontsize=12)
        fig.subplots_adjust(bottom=0.05,left=0.1,top=0.85,hspace=0.3, right =0.95)
        plt.suptitle('Transporte meridional (Sv)', fontsize=14)

        plt.savefig(dir_graf,dpi=500)

#%%

"""
Gráfico de transporte meridional en la latitud central de la cuenca

INPUTS
psiF_todos: array con función corriente para el tiempo final concatenadas para las distintas simulaciones
X_todos:Coordenadas x concatenadas para las distintas simulaciones
U_todos: Escala típica velocidad [m/s] concatenadas para las distintas simulaciones
L_todos: Longitud de la cuenca [m]
D_todos: Profundidad de la cuenca [m]
Lx: float, tamano de la cuenca (direccion X)
nx: int, numero de punto de grilla (direccion X)
Nom_sim: Nombres de las simulaciones
dir_graf: str, direccion con nombre donde guardar el grafico

"""
def graf_TrasMer_LatCent(psiF_todos,X_todos,U_todos,L_todos,D_todos,Lx,nx,Nom_sim,dir_graf):

    #Cargamos las librerias necesarias
    import numpy as np
    import matplotlib.pyplot as plt

    TrasMer = np.empty(shape=[np.size(psiF_todos,0),np.size(psiF_todos,1),np.size(psiF_todos,2)-1])
    TrasMer[:]=np.nan
    for k in range(np.size(psiF_todos,0)):
        TrasMer[k]=np.diff(np.squeeze(psiF_todos[k]),n=1,axis=1)*D_todos[k]

    # Reticula alternativa
    Xalt=np.linspace(0,Lx,num=nx-1)
    fig = plt.figure(figsize=(5,3))
    plt.plot(Xalt,np.squeeze(U_todos[k]*L_todos[k]/1000000*TrasMer[:,int(np.size(TrasMer,1)/2),:]).T,linewidth=0.8)
    plt.xlabel('Distancia al borde oeste (Km)',fontsize=8)
    plt.ylabel('Transporte meridional (Sv)',fontsize=8)
    plt.axhline(y=0,linewidth=.5,linestyle='--',color='k')
    plt.title('Transporte meridional en la latitud central de la cuenca', fontsize=8.5)
    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)
    plt.legend(Nom_sim,fontsize=7,loc = 5)
    fig.subplots_adjust(bottom=0.15,top=0.85,left=0.1,right=0.97)
    plt.savefig(dir_graf,dpi=500)

#%%

"""
Gráfico de vorticidad relativa en la latitud central de la cuenca

INPUTS
vortF_todos: array con función corriente para el tiempo final concatenadas para las distintas simulaciones
X_todos:Coordenadas x concatenadas para las distintas simulaciones
U_todos: Escala típica velocidad [m/s] concatenadas para las distintas simulaciones
L_todos: Longitud de la cuenca [m]
Nom_sim: Nombres de las simulaciones
dir_graf: str, direccion con nombre donde guardar el grafico

"""
def graf_vortF_LatCent(vortF_todos,X_todos,U_todos,L_todos,Nom_sim,dir_graf):

    #Cargamos las librerias necesarias
    import numpy as np
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(5,3))
    plt.plot(X_todos[0],np.squeeze(vortF_todos[:,int(np.size(vortF_todos,1)/2),:]).T*U_todos[0]/L_todos[0]*10000000,linewidth=0.8)
    plt.xlabel('Distancia al borde oeste (Km)',fontsize=8)
    plt.ylabel('Vorticidad Relativa $(*10^{-7} s^{-1})$',fontsize=8)
    plt.axhline(y=0,linewidth=.5,linestyle='--',color='k')
    plt.title('Vorticidad Relativa en la latitud central de la cuenca', fontsize=8.5)
    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)
    plt.legend(Nom_sim,fontsize=7,loc = 5)
    fig.subplots_adjust(bottom=0.15,top=0.85,left=0.1,right=0.97)
    plt.savefig(dir_graf,dpi=500)

#%%

"""
Calculo de transporte total meridional de la Corriente de Borde Oeste en la latitud central

INPUTS
psiF_todos: array con función corriente para el tiempo final concatenadas para las distintas simulaciones
U_todos: Escala típica velocidad [m/s] concatenadas para las distintas simulaciones
L_todos: Longitud de la cuenca [m]
D_todos: Profundidad de la cuenca [m]

OUTPUTS
TrasMer_CBO_LatCent: transporte total meridional de la Corriente de Borde Oeste (Sv)
"""
def Calc_TrasMer_CBO_LatCent(psiF_todos,U_todos,L_todos,D_todos):

    #Cargamos las librerias necesarias
    import numpy as np

    TrasMer_CBO_LatCent=np.empty(shape=[np.size(psiF_todos,0)])
    TrasMer_CBO_LatCent[:]=np.nan

    for k in range(np.size(psiF_todos,0)):
        TrasMer=np.diff(np.squeeze(psiF_todos[k]),n=1,axis=1)*D_todos[k]
        TrasMer_LatCent=np.squeeze(U_todos[k]*L_todos[k]/1000000*TrasMer[int(np.size(TrasMer,0)/2),:])
        a=TrasMer_LatCent[1]/abs(TrasMer_LatCent[1])
        i=2
        m=0
        while i<len(TrasMer_LatCent):
            b=TrasMer_LatCent[i]/abs(TrasMer_LatCent[i])
            if b==a:
                i=i+1
            else:
                m=i
                break
        Limite_CBO_LatCent=m
        TrasMer_CBO_LatCent[k]=sum(TrasMer_LatCent[0:int(Limite_CBO_LatCent)])

    return TrasMer_CBO_LatCent


#%%

"""
Calculo de transporte total meridional en la latitud central

INPUTS
psiF_todos: array con función corriente para el tiempo final concatenadas para las distintas simulaciones
U_todos: Escala típica velocidad [m/s] concatenadas para las distintas simulaciones
L_todos: Longitud de la cuenca [m]
D_todos: Profundidad de la cuenca [m]

OUTPUTS
TrasMer_Total_LatCent: transporte total meridional de la Corriente de Borde Oeste (Sv)
"""
def Calc_TrasMer_Total_LatCent(psiF_todos,U_todos,L_todos,D_todos):

    #Cargamos las librerias necesarias
    import numpy as np

    TrasMer_Total_LatCent=np.empty(shape=[np.size(psiF_todos,0)])
    TrasMer_Total_LatCent[:]=np.nan

    for k in range(np.size(psiF_todos,0)):
        TrasMer=np.diff(np.squeeze(psiF_todos[k]),n=1,axis=1)*D_todos[k]
        TrasMer_LatCent=np.squeeze(U_todos[k]*L_todos[k]/1000000*TrasMer[int(np.size(TrasMer,0)/2),:])
        TrasMer_Total_LatCent[k]=sum(TrasMer_LatCent)
    return TrasMer_Total_LatCent

#%%

"""
Extensión de la Corriente de Borde Oeste en la latitud central

INPUTS
psiF_todos: array con función corriente para el tiempo final concatenadas para las distintas simulaciones
dx_todos: Distancia entre puntos del eje X

OUTPUTS
Ext_CBO_LatCent: Extensión de la Corriente de Borde Oeste (Km)
"""
def Calc_Ext_CBO_LatCent(psiF_todos,dx_todos):

    #Cargamos las librerias necesarias
    import numpy as np

    Ext_CBO_LatCent=np.empty(shape=[np.size(psiF_todos,0)])
    Ext_CBO_LatCent[:]=np.nan

    for k in range(np.size(psiF_todos,0)):
        TrasMer=np.diff(np.squeeze(psiF_todos[k]),n=1,axis=1)
        TrasMer_LatCent=np.squeeze(TrasMer[int(np.size(TrasMer,0)/2),:])
        a=TrasMer_LatCent[1]/abs(TrasMer_LatCent[1])
        i=2
        m=0
        while i<len(TrasMer_LatCent):
            b=TrasMer_LatCent[i]/abs(TrasMer_LatCent[i])
            if b==a:
                i=i+1
            else:
                m=i
                break
        Limite_CBO_LatCent=m
        Ext_CBO_LatCent[k]=dx_todos[k]*Limite_CBO_LatCent

    return Ext_CBO_LatCent

#%%

"""
Gráfico de los términos de la ecuación estacionaria del modelo de Stommel en la latitud central de la cuenca para una simulación

INPUTS
psiF: función corriente para el tiempo final
vortF: vorticidad relativa para el tiempo final
QG_curlw: rotor tensión del viento
Ef:  Número de Ekman vertical
X:Coordenadas x
Lx: float, tamano de la cuenca (direccion X)
nx: int, numero de punto de grilla (direccion X)
Nom_sim: Nombres de las simulación
dir_graf: str, direccion con nombre donde guardar el grafico

"""
def graf_terminos_ModStommel(psiF,vortF,QG_curlw,Ef,X,Lx,nx,Nom_sim,dir_graf):

    #Cargamos las librerias necesarias
    import numpy as np
    import matplotlib.pyplot as plt

    # Reticula alternativa
    Xalt=np.linspace(0,Lx,num=nx-1)
    Ter1=np.diff(psiF,n=1,axis=1)
    Ter1_LatCent=np.squeeze(Ter1[int(np.size(Ter1,0)/2),:])*Lx/(nx-1)

    fig = plt.figure(figsize=(5,3))
    plt.plot(Xalt,Ter1_LatCent,color='r',linewidth=0.8)
    plt.plot(X,-QG_curlw[int(np.size(Ter1,0)/2),1:201],color='b',linewidth=0.8)
    plt.plot(X,Ef*vortF[int(np.size(Ter1,0)/2),:],color='k',linewidth=0.8)
    plt.xlabel('Distancia al borde oeste (Km)',fontsize=8)
    plt.axhline(y=0,linewidth=.5,linestyle='--',color='k')
    plt.title('Términos ecuación Stommel adimensionalizados en la latitud central '+Nom_sim, fontsize=8.5)
    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)
    plt.legend(['Término 1','Término 2','Término 3'],fontsize=7,loc = 4)
    fig.subplots_adjust(bottom=0.15,top=0.85,left=0.1,right=0.97)
    plt.savefig(dir_graf,dpi=500)

#%%

"""
Número de iteraciones para alcanzar el estado estacionario

INPUTS
ecin: array con QG_DIAG concatenado (en la primera dimensión) de las simulaciones

OUTPUTS
TiempoEst: número de iteración en la que se alcanza el estado estacionario
"""

def Calc_TiempoEst(ecin):

    #Cargamos las librerias necesarias
    import numpy as np

    TiempoEst=np.empty(shape=[np.size(ecin,0)])
    TiempoEst[:]=np.nan

    for k in range(np.size(ecin,0)):
        EC=np.squeeze(ecin[k,:,3])
        i=len(EC)-1
        ECfinal=EC[len(EC)-1]
        deltaEC=((ECfinal-EC)/ECfinal)*100
        while i>=0:
            if deltaEC[i]<1:
                i=i-1
            else:
                m=i+2
                break
        TiempoEst[k]=m

    return TiempoEst
