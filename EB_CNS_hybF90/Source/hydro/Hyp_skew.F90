module advection_module2
  use amrex_fort_module, only : rt=>amrex_real
  implicit none
  private
  public :: hyp_skew_3d

contains

  subroutine hyp_skew_3d(q, qd_lo, qd_hi, &
                     lo, hi, dx, flux1, flux2, flux3)

    use amrex_mempool_module, only : amrex_allocate, amrex_deallocate
    use cns_module, only : urho, umx, umy, umz, ueden, ueint, utemp, nvar, &
         nvarsolve,qrho,qu,qv,qw,qp,qc,qeint,qtemp,qvar, smallp, smallr

    use skew_module

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in   ) ::     q( qd_lo(1): qd_hi(1), qd_lo(2): qd_hi(2), qd_lo(3): qd_hi(3),QVAR)
    real(rt), intent(  out) :: flux1(lo(1):hi(1)+1,lo(2):hi(2)  ,lo(3):hi(3)  ,nvarsolve) !<---
    real(rt), intent(  out) :: flux2(lo(1):hi(1)  ,lo(2):hi(2)+1,lo(3):hi(3)  ,nvarsolve) !<--
    real(rt), intent(  out) :: flux3(lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3)+1,nvarsolve) !<--

    integer :: i, j, k,nv,is,ii,jj,kk
    integer :: qtlo(3), qthi(3), fd1_lo(3), fd2_lo(3), fd3_lo(3), fd1_hi(3), fd2_hi(3), fd3_hi(3)
    ! SNM: new variables
    real(rt) :: kin,sp,sr,RR
    real(rt),  allocatable :: fr(:)
    real(rt),  allocatable :: UU(:,:),UUP(:,:),VV(:)
    real(rt),  allocatable :: PP(:)
    !----  

    qtlo = lo - 1
    qthi = hi + 1

    fd1_lo = lo
    fd1_hi = hi; fd1_hi(1) = hi(1)+1
    fd2_lo = lo
    fd2_hi = hi; fd2_hi(2) = hi(2)+1
    fd3_lo = lo
    fd3_hi = hi; fd3_hi(3) = hi(3)+1

    ! allocate arrays
    allocate( UU(iskew1-1:iskew2,nvarsolve),UUP(iskew1-1:iskew2,nvarsolve),VV(iskew1-1:iskew2) )
    allocate( PP(iskew1-1:iskew2))
    allocate( sk(iskew1:iskew2,iskew1:iskew2))
    !
    allocate(fr(nvarsolve))
    
    ! detector here  in X 
    ! jump in x, create variable
    !
 

    !-------------------------------------------------------   X
    call bl_proffortfuncstart_int(1)
    ! compute flux in x 
    do  k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)+1
          ! prepare arrays
          do is = iskew1-1,iskew2
            ii  = i + is
            UU(is,URHO)  =  q(ii,j,k,QRHO)
            UU(is,UMX)   =  q(ii,j,k,QRHO)*q(i,j,k,QU)
            UU(is,UMY)   =  q(ii,j,k,QRHO)*q(i,j,k,QV)
            UU(is,UMZ)   =  q(ii,j,k,QRHO)*q(i,j,k,QW)
    kin = UU(ii,UMX)*q(i,j,k,QU) + UU(ii,UMY)*q(i,j,k,QV) + UU(ii,UMZ)*q(i,j,k,QW)
            UU(is,UEDEN) =  q(ii,j,k,QRHO)*q(i,j,k,QEINT) + 0.5d0*kin
            ! species go here NRY
            !UU(is,ispec)   =  q(ii,j,k,QRHO)*q(i,j,k,ispec)
            VV(is) = q(ii,j,k,QU)
            PP(is) = q(ii,j,k,QP)
            ! enthalpy term
            UUP(is,UEDEN)   =  UU(is,UEDEN)  + PP(is)
          end do         
                
          !  skew fluxes
          do nv = 1, nvarsolve  
            fr(nv) = flux_advskew(UUP(iskew1-1:iskew2,nv),VV(iskew1-1:iskew2))
          end do       
          ! pressure gradient (linear flux)
          fr(nvv)  = fr(nvv)  + flux_pres(PP(iskew1-1:iskew2))
          ! shock capturing
          ! ------- compute wave speed
          sp   =  max( abs(  q(i,j,k,QU) +   q(i,j,k,QC)), abs(  q(i,j,k,QU) - q(i,j,k,QC))) 
          sr   =  max( abs(q(i+1,j,k,QU) + q(i+1,j,k,QC)), abs(q(i+1,j,k,QU) - q(i+1,j,k,QC)) )
          RR   =  0.5d0*(sp + sr) 

! #if shock_capt   
!           psip = sensor_dir(i,j,k,1)    ! needs to change
!           psir = sensor_dir(i+1,j,k,1)    
!           e2 = c2skew*abs(RR)*max(psip,psir)
!           do nv = 1,nvarsolve        
!             fr(nv) = fr(nv) - flux_shock(UU(iskew1:iskew2,nv),e2) 
!           end do
! #endif
          !---------- high frequency damping
#if skew_damp 
          e4 = max(0.0_wp,c4skew*abs(RR)-e2)
          do nv = 1,nvarsolve  
            fr(nv) = fr(nv) + flux_damp(UU(iskew1:iskew2,nv),e4)
          end do 
#endif
          do nv = 1,nvarsolve
            flux1(i,j,k,nv) = fr(nv)
          end do
        end do
       end do
    end do
    call bl_proffortfuncstop_int(1)
    !--------------------------------------------- 
    ! detector here  in Y

    call bl_proffortfuncstart_int(2)
    do  k = lo(3), hi(3)
      do  j = lo(2), hi(2)+1
        do i = lo(1), hi(1)

          ! prepare arrays
          do is = iskew1-1,iskew2
            jj  = j + is
            UU(is,URHO)  =  q(i,jj,k,QRHO)
            UU(is,UMX)   =  q(ii,jj,k,QRHO)*q(i,j,k,QU)
            UU(is,UMY)   =  q(ii,jj,k,QRHO)*q(i,j,k,QV)
            UU(is,UMZ)   =  q(ii,jj,k,QRHO)*q(i,j,k,QW)
    kin = UU(is,UMX)*q(i,jj,k,QU) + UU(is,UMY)*q(i,jj,k,QV) + UU(is,UMZ)*q(i,jj,k,QW)
            UU(is,UEDEN) =  q(ii,jj,k,QRHO)*q(i,jj,k,QEINT) + 0.5d0*kin
            ! species go here NRY
            !UU(is,ispec)   =  q(ii,j,k,QRHO)*q(i,j,k,ispec)
            VV(is) = q(i,jj,k,QV)
            PP(is) = q(i,jj,k,QP)
            ! enthalpy term
            UUP(is,UEDEN)   =  UU(is,UEDEN)  + PP(is)
          end do      
          !  skew fluxes
          do nv = 1, nvarsolve  
            fr(nv) = flux_advskew(UUP(iskew1-1:iskew2,nv),VV(iskew1-1:iskew2))
          end do       
          ! pressure gradient (linear flux)
          fr(nvv)  = fr(nvv)  + flux_pres(PP(iskew1-1:iskew2))
          
        end do   
      end do
    end do
    call bl_proffortfuncstop_int(2)
    !--------------------------------------------- 

    !--------------------------------------------- 
    ! detector here  in z


    call bl_proffortfuncstart_int(3)
    do  k = lo(3), hi(3)+1
      do  j = lo(2), hi(2)
        do i = lo(1), hi(1)
        end do  
      end do
    end do
    call bl_proffortfuncstop_int(3)

    ! deallocate arrays

  end subroutine hyp_skew_3d

end module advection_module2
