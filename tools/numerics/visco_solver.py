import numpy as np
import derivatives

# calculate viscous fluxes and update FX,FY,FZ
def calc_viscfluxes(model,ist,ien, visc,cond,o_dx,o_dy,o_dz, U, V, W, T, FX,FY,FZ):

    # CD
    if (model==0):
        
        # Compute FLUXES
        for i in range(ist[0], ien[0]+1):
            for j in range(ist[1], ien[1]+1):
                for k in range(ist[2], ien[2]+1):
                    # compute gradients using CD
                    dudx = derivatives.dfdx_cc(i,j,k,0,U)*o_dx   
                    dudy = derivatives.dfdx_cc(i,j,k,1,U)*o_dy   
                    dudz = derivatives.dfdx_cc(i,j,k,2,U)*o_dz
                    dvdx = derivatives.dfdx_cc(i,j,k,0,V)*o_dx   
                    dvdy = derivatives.dfdx_cc(i,j,k,1,V)*o_dy   
                    dvdz = derivatives.dfdx_cc(i,j,k,2,V)*o_dz
                    dwdx = derivatives.dfdx_cc(i,j,k,0,W)*o_dx   
                    dwdy = derivatives.dfdx_cc(i,j,k,1,W)*o_dy   
                    dwdz = derivatives.dfdx_cc(i,j,k,2,W)*o_dz
                    dTdx = derivatives.dfdx_cc(i,j,k,0,T)*o_dx   
                    dTdy = derivatives.dfdx_cc(i,j,k,1,T)*o_dy   
                    dTdz = derivatives.dfdx_cc(i,j,k,2,T)*o_dz
                    # divu
                    divu = dudx+ dvdy + dwdz
                    # shear stress
                    tau_xx = visc*(2*dudx - 2*divu/3) 
                    tau_xy = visc*(dudy + dvdx) 
                    tau_xz = visc*(dudz + dwdx)
                    tau_yy = visc*(2*dvdy - 2*divu/3) 
                    tau_yz = visc*(dvdz + dwdy)
                    tau_zz = visc*(2*dwdz - 2*divu/3) 

                    u = U[i,j,k]
                    v = V[i,j,k]
                    w = W[i,j,k]

                    # the minus sign is to be consistent as it is added to FLuxes  FT = F - Fvisc
                    # F1viscx =0
                    FX[i,j,k,1] -= tau_xx
                    FX[i,j,k,2] -= tau_xy
                    FX[i,j,k,3] -= tau_xz
                    FX[i,j,k,4] -= cond*dTdx + u*tau_xx + v*tau_xy + w*tau_xz
                    # F1viscy = 0
                    FY[i,j,k,1] -= tau_xy
                    FY[i,j,k,2] -= tau_yy
                    FY[i,j,k,3] -= tau_yz
                    FY[i,j,k,4] -= cond*dTdy + u*tau_xy + v*tau_yy + w*tau_yz
                    # F1viscz = 0
                    FZ[i,j,k,1] -= tau_xz
                    FZ[i,j,k,2] -= tau_yz
                    FZ[i,j,k,3] -= tau_zz
                    FZ[i,j,k,4] -= cond*dTdz + u*tau_xz + v*tau_yz + w*tau_zz
     # FV
    if (model==1):

        # Compute Viscous Fluxes  X
        for i in range(ist[0], ien[0]+2):
            for j in range(ist[1], ien[1]+1):
                for k in range(ist[2], ien[2]+1):
                    # compute gradients using FV based on face  dU[i+1/2]
                    # FLUX[i]  face[i] is between cell [i+1] and cell[i]
                    dudx = derivatives.dfdx_fv(i,j,k,0,U)*o_dx   
                    dvdx = derivatives.dfdx_fv(i,j,k,0,V)*o_dx   
                    dwdx = derivatives.dfdx_fv(i,j,k,0,W)*o_dx   
                    dTdx = derivatives.dfdx_fv(i,j,k,0,T)*o_dx   
                    
                    # tangent deriv
                    dudy = derivatives.dfcross_fv(i,j,k,0,1,U)*o_dy  
                    dudz = derivatives.dfcross_fv(i,j,k,0,2,U)*o_dz
                    dvdy = derivatives.dfcross_fv(i,j,k,0,1,V)*o_dy   
                    dvdz = derivatives.dfcross_fv(i,j,k,0,2,V)*o_dz
                    dwdy = derivatives.dfcross_fv(i,j,k,0,1,W)*o_dy   
                    dwdz = derivatives.dfcross_fv(i,j,k,0,2,W)*o_dz

                    # divu
                    divu = dudx+ dvdy + dwdz
                    # shear stress
                    tau_xx = visc*(2*dudx - 2*divu/3) 
                    tau_xy = visc*(dudy + dvdx) 
                    tau_xz = visc*(dudz + dwdx)

                    # Shall I interpolate velocity in face
                    # u = U[i,j,k]
                    # v = V[i,j,k]
                    # w = W[i,j,k]
                    u = derivatives.interp_fv_phi(i,j,k,0,U)
                    v = derivatives.interp_fv_phi(i,j,k,0,V)
                    w = derivatives.interp_fv_phi(i,j,k,0,W)

                    # the minus sign is to be consistent as it is added to FLuxes  FT = F - Fvisc
                    FX[i,j,k,1] -= tau_xx
                    FX[i,j,k,2] -= tau_xy
                    FX[i,j,k,3] -= tau_xz
                    FX[i,j,k,4] -= cond*dTdx + u*tau_xx + v*tau_xy + w*tau_xz

        # Compute Viscous Fluxes  Y
        idir = 1 
        for i in range(ist[0], ien[0]+1):
            for j in range(ist[1], ien[1]+2):
                for k in range(ist[2], ien[2]+1):
                    # compute gradients using FV based on face  dU[i+1/2]
                    # FLUX[i]  face[i] is between cell [i+1] and cell[i]
                    dudy = derivatives.dfdx_fv(i,j,k,idir,U)*o_dy   
                    dvdy = derivatives.dfdx_fv(i,j,k,idir,V)*o_dy   
                    dwdy = derivatives.dfdx_fv(i,j,k,idir,W)*o_dy   
                    dTdy = derivatives.dfdx_fv(i,j,k,idir,T)*o_dy   
                    
                    # tangent deriv
                    dudx = derivatives.dfcross_fv(i,j,k,idir,0,U)*o_dx  
                    dudz = derivatives.dfcross_fv(i,j,k,idir,2,U)*o_dz
                    dvdx = derivatives.dfcross_fv(i,j,k,idir,0,V)*o_dx   
                    dvdz = derivatives.dfcross_fv(i,j,k,idir,2,V)*o_dz
                    dwdx = derivatives.dfcross_fv(i,j,k,idir,0,W)*o_dx   
                    dwdz = derivatives.dfcross_fv(i,j,k,idir,2,W)*o_dz

                    # divu
                    divu = dudx+ dvdy + dwdz
                    # shear stress
                    tau_xy = visc*(dudy + dvdx) 
                    tau_yy = visc*(2*dvdy - 2*divu/3) 
                    tau_yz = visc*(dvdz + dwdy)

                    # u = U[i,j,k]
                    # v = V[i,j,k]
                    # w = W[i,j,k]
                    u = derivatives.interp_fv_phi(i,j,k,idir,U)
                    v = derivatives.interp_fv_phi(i,j,k,idir,V)
                    w = derivatives.interp_fv_phi(i,j,k,idir,W)

                    # the minus sign is to be consistent as it is added to FLuxes  FT = F - Fvisc
                    FY[i,j,k,1] -= tau_xy
                    FY[i,j,k,2] -= tau_yy
                    FY[i,j,k,3] -= tau_yz
                    FY[i,j,k,4] -= cond*dTdy + u*tau_xy + v*tau_yy + w*tau_yz

        idir = 2
        # Compute Viscous Fluxes  Z
        for i in range(ist[0], ien[0]+1):
            for j in range(ist[1], ien[1]+1):
                for k in range(ist[2], ien[2]+2):
                    # compute gradients using FV based on face  dU[i+1/2]
                    # FLUX[i]  face[i] is between cell [i+1] and cell[i]
                    dudz = derivatives.dfdx_fv(i,j,k,idir,U)*o_dz  
                    dvdz = derivatives.dfdx_fv(i,j,k,idir,V)*o_dz   
                    dwdz = derivatives.dfdx_fv(i,j,k,idir,W)*o_dz   
                    dTdz = derivatives.dfdx_fv(i,j,k,idir,T)*o_dz   
                    
                    # tangent deriv
                    dudx = derivatives.dfcross_fv(i,j,k,idir,0,U)*o_dx  
                    dudy = derivatives.dfcross_fv(i,j,k,idir,1,U)*o_dy
                    dvdx = derivatives.dfcross_fv(i,j,k,idir,0,V)*o_dx   
                    dvdy = derivatives.dfcross_fv(i,j,k,idir,1,V)*o_dy
                    dwdx = derivatives.dfcross_fv(i,j,k,idir,0,W)*o_dx   
                    dwdy = derivatives.dfcross_fv(i,j,k,idir,1,W)*o_dy

                    # divu
                    divu = dudx+ dvdy + dwdz
                    # shear stress
                    tau_xz = visc*(dudz + dwdx)
                    tau_yz = visc*(dvdz + dwdy)
                    tau_zz = visc*(2*dwdz - 2*divu/3) 

                    u = derivatives.interp_fv_phi(i,j,k,idir,U)
                    v = derivatives.interp_fv_phi(i,j,k,idir,V)
                    w = derivatives.interp_fv_phi(i,j,k,idir,W)

                    # F1viscz = 0
                    FZ[i,j,k,1] -= tau_xz
                    FZ[i,j,k,2] -= tau_yz
                    FZ[i,j,k,3] -= tau_zz
                    FZ[i,j,k,4] -= cond*dTdz + u*tau_xz + v*tau_yz + w*tau_zz