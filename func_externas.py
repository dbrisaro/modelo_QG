"""
Funciones externas para el cómputo del modelo QG
Dani risaro
Julio 2019

"""

# def arakawa_jacobian(jpp, jxp, jpx, pb, psib, delta, ix, iy):
#     # Arakawa's jacobian
#
#     jpp = ((pb[ix+1,iy]-pb[ix-1,iy])*(psib[ix,iy+1]-psib[ix,iy-1])-
#         (pb[ix,iy+1]-pb[ix,iy-1])*(psib[ix+1,iy]-psib[ix-1,iy]))*delta
#
#     jxp = (psib[ix,iy+1]*(pb[ix+1,iy+1]-pb[ix-1,iy+1])-psib[ix,iy-1]
#         *(pb[ix+1,iy-1]-pb[ix-1,iy-1])-psib[ix+1,iy]*(pb[ix+1,iy+1]
#         -pb[ix+1,iy-1])+psib[ix-1,iy]*(pb[ix-1,iy+1]-pb[ix-1,iy-1]))*delta
#
#     jpx = (pb[ix+1,iy]*(psib[ix+1,iy+1]-psib[ix+1,iy-1])-pb[ix-1,iy]
#         *(psib[ix-1,iy+1]-psib[ix-1,iy-1])-pb[ix,iy+1]*(psib[ix+1,iy+1]
#         -psib[ix-1,iy+1])+pb[ix,iy-1]*(psib[ix+1,iy-1]-psib[ix-1,iy-1]))*delta
#
#     return jpp, jxp, jpx

def arakawa_jacobian(jpp, jxp, jpx, jxx, pb, psib, delta, ix, iy):
    # Arakawa's jacobian

    jpp = ((pb[ix+1,iy]-pb[ix-1,iy])*(psib[ix,iy+1]-psib[ix,iy-1])-
        (pb[ix,iy+1]-pb[ix,iy-1])*(psib[ix+1,iy]-psib[ix-1,iy]))*delta

    jxp = (psib[ix,iy+1]*(pb[ix+1,iy+1]-pb[ix-1,iy+1]) -
           psib[ix,iy-1]*(pb[ix+1,iy-1]-pb[ix-1,iy-1]) -
           psib[ix+1,iy]*(pb[ix+1,iy+1]-pb[ix+1,iy-1]) +
           psib[ix-1,iy]*(pb[ix-1,iy+1]-pb[ix-1,iy-1]))*delta

    jpx = (pb[ix+1,iy]*(psib[ix+1,iy+1]-psib[ix+1,iy-1])-
           pb[ix-1,iy]*(psib[ix-1,iy+1]-psib[ix-1,iy-1])-
           pb[ix,iy+1]*(psib[ix+1,iy+1]-psib[ix-1,iy+1])+
           pb[ix,iy-1]*(psib[ix+1,iy-1]-psib[ix-1,iy-1]))*delta

    jxx = ((pb[ix+1,iy+1]-pb[ix-1,iy-1])*(psib[ix-1,iy+1]-psib[ix+1,iy-1]) -
          (pb[ix-1,iy+1]-pb[ix+1,iy-1])*(psib[ix+1,iy+1]-psib[ix-1,iy-1]))*(1/2)*delta

    return jpp, jxp, jpx, jxx

def horizontal_mixing(delpsi, delpsi4, psic, delta, ix, iy):
    # horizontal mixing

    delpsi[ix,iy] = (psic[ix+1,iy]+psic[ix-1,iy]+psic[ix,iy+1]+psic[ix,iy-1]-4.*psic[ix,iy])*delta

    delpsi4[ix,iy] = (delpsi[ix+1,iy]+delpsi[ix-1,iy]+delpsi[ix,iy+1]+ \
                delpsi[ix,iy-1]-4.*delpsi[ix,iy])*delta

    return delpsi, delpsi4

def vorticity(psia, psic, jpp, jxp, jpx, jxx, pb, curlt, delpsi, delpsi4, c1, c2, c3, c4, c5, c6, c7, ix, iy):
    # update vorticity

    psia[ix,iy] = c1*psic[ix,iy] - c5*(jpp+jxp+jpx+jxx) - c2*(pb[ix+1,iy] -
            pb[ix-1,iy]) + c3*curlt[ix,iy] - c4*psic[ix,iy] + c6*delpsi[ix,iy] - c7*delpsi4[ix,iy]

    return psia

def laplacian_p(r, pa, psia, delta, alfa, fxr, ix, iy):
    # solve the Laplacian for the stream function by over-relaxation

    r[ix,iy] = ((pa[ix-1,iy]+pa[ix+1,iy]+pa[ix,iy-1]+pa[ix,iy+1]-4.*pa[ix,iy])*
                delta-psia[ix,iy])
    pa[ix,iy]= pa[ix,iy]+alfa*r[ix,iy]*fxr

    return r, pa

def cond_borde_x_p(pa, bc1, iy, nx):
    # set the boundary conditons on the stream function

    pa[0,j] = bc1*pa[2,j]
    pa[1,j] = zero
    pa[nx-1,j] = bc1*pa[imm2,j]
    pa[imm1,j] = zero

def cond_borde_y_p(pa, bc1, ix, ny):
    # set the boundary conditons on the stream function

    pa[i,0] = bc1*pa[i,2]
    pa[i,1] = zero
    pa[i,ny-1] = bc1*pa[i,jmm2]
    pa[i,jmm1] = zero

    return pa
