module skew_module

  use amrex_fort_module, only : rt=>amrex_real
  
  implicit none
  
  !integer, public :: iskew1,iskew2
#if (flux_skew ==2)
  integer, parameter, public  :: iskew1 = 0,iskew2 = 1  
#elif (flux_skew ==6)
  integer, parameter, public  :: iskew1 = -2 ,iskew2 = 3
#else
  integer, parameter, public  :: iskew1 = -1 ,iskew2 = 2
  
#endif    

  real(rt), parameter :: c2skew = 0.2d0, c4skew = 0.01d0


  real(rt), private, allocatable :: sk(:,:)

  ! Some common fractions
  real(rt), parameter :: r12 = 1.0d0/2.0d0, r13 = 1.0d0/3.0d0, &
                          r14 = 1.0d0/4.0d0, r23 = 2.0d0/3.0d0, &
                          r34 = 3.0d0/4.0d0, r35 = 3.0d0/5.0d0, &
                          r43 = 4.0d0/3.0d0, r16 = 1.0d0/6.0d0, &
                          r56 = 5.0d0/6.0d0, r76 = 7.0d0/6.0d0, &
                          r32 = 3.0d0/2.0d0, r15 = 1.0d0/5.0d0, &
                          r52 = 5.0d0/2.0d0, r45 = 4.0d0/5.0d0
   real(rt), parameter :: r11_6 =  11.0d0/6.0d0, r13_12 = 13.0d0/12.0d0,  &
                          r13_10 = 13.0d0/10.0d0, r1_12 = 1.0d0/12.0d0, &
                          r1_24 = 1.0d0/24.0d0, r7_5760 = 7.0d0/5760d0, &
                          r1_16 = 1.0d0/16.0d0, r7_12 = 7.0d0/12.0d0


  public :: flux_advskew,flux_shock,flux_damp

  contains

  !***************************************************************************!
  ! p at i+1/2 (second order) 
  function flux_pres2(Pi,Pip1) result(flux)
    
    implicit none

    real(rt), intent(in) ::Pi,Pip1
    real(rt) :: flux

    flux = r12*( Pi + Pip1 ) 

  end function flux_pres2  
  !***************************************************************************!
  ! p at i+1/2 (fourth order) 
  function flux_pres4(Pim1,Pi,Pip1,Pip2) result(flux)
    
    implicit none

    real(rt), intent(in) ::Pi,Pip1,pip2,pim1
    real(rt) :: flux

    flux = r7_12*(Pi + Pip1) - r1_12*(Pip2+Pim1)

  end function flux_pres4  
  !***************************************************************************!
  ! p at i+1/2 
  function flux_pres(P) result(flux)
    
    implicit none

    real(rt), intent(in) :: P(iskew1:)
    real(rt) :: flux

#if flux_skew==2
    flux = r12*( P(0) + P(1) ) 
#else    
    flux = r7_12*(P(0) + P(1)) - r1_12*(P(2)+P(-1))
#endif

  end function flux_pres  
  !***************************************************************************!
  ! F(U,V) at i+1/2 (second order) divergence
  function flux_div2(Ui,Uip1,Vi,Vip1) result(flux)
    
    implicit none

    real(rt), intent(in) :: Ui,Uip1,Vi,Vip1
    real(rt) :: flux

    flux = r12*( Ui*Vi + Uip1*Vip1 ) 

  end function flux_div2  
  !***************************************************************************!
  ! F(U,V) at i+1/2 (second order) skew (stencil i  to i+1)
  function flux_skew2(Ui,Uip1,Vi,Vip1) result(flux)
    
    implicit none

    real(rt), intent(in) :: Ui,Uip1,Vi,Vip1
    real(rt) :: flux

    flux = r14*( Ui + Uip1 )  * ( Vi + Vip1)

  end function flux_skew2  
  !***************************************************************************! 
  ! F(U,V) at i+1/2 (fourth order) skew (stencil i-1  to i+2)
  function flux_skew4(Uim1,Ui,Uip1,Uip2,Vim1,Vi,Vip1,Vip2) result(flux)
    
    implicit none

    real(rt), intent(in) :: Uim1,Ui,Uip1,Uip2,Vim1,Vi,Vip1,Vip2
    real(rt) :: flux

    flux = r13*( Ui + Uip1 )  * ( Vi + Vip1) - & 
    &  r1_24* (Uim1*Vim1 + Uim1*Vip1 + Ui*Vi + Ui*Vip2 + Uip1*Vip1 + Uip1*Vim1 + Uip2*Vi +Uip2*Vip2 )
   
  end function flux_skew4 
  !***************************************************************************! 
  ! F(U,V) at i+1/2 (sixth order) skew (stencil i-2  to i+3)
  function flux_skew6(U,V) result(flux)
    
    implicit none
    real(rt), intent(in) :: U(iskew1:),V(iskew1:)
    real(rt) :: flux
  
    flux = sk(0,0)*U(0)*V(0) + sk(1,0)*U(1)*V(0) + & 
      &    sk(2,0)*U(2)*V(0)+sk(3,0)*U(3)*V(0) + & 
      &    sk(-2,-2)*U(-2)*V(-2) + sk(1,-2)*U(1)*V(-2) + & 
      &    sk( -1,-1)*U(-1)*V(-1) + sk(1,-1)*U(1)*V(-1) + sk(2,-1)*U(2)*V(-1)+ &
      &    sk(0,1)*U(0)*V(1)+sk(-2,1)*U(-2)*V(1)+sk(-1,1)*U(-1)*V(1) +  & 
      &    sk(1,1)*U(1)*V(1) + sk(0,2)*U(0)*V(2)+sk(-1,2)*U(-1)*V(2) +  & 
      &    sk(2,2)*U(2)*V(2)+ sk(0,3)*U(0)*V(3)+sk(3,3)*U(3)*V(3) 

  end function flux_skew6 
  !***************************************************************************! 
  ! skew F(U,V) at i+1/2 
  function flux_advskew(U,V) result(flux)
    
    implicit none
    real(rt), intent(in) :: U(iskew1:),V(iskew1:)
    real(rt) :: flux
  
#if flux_skew==2
    flux = r14*( U(0) + U(1) )  * ( V(0) + V(1))
#elif flux_skew==6    
    flux = sk(0,0)*U(0)*V(0) + sk(1,0)*U(1)*V(0) + & 
      &    sk(2,0)*U(2)*V(0)+sk(3,0)*U(3)*V(0) + & 
      &    sk(-2,-2)*U(-2)*V(-2) + sk(1,-2)*U(1)*V(-2) + & 
      &    sk( -1,-1)*U(-1)*V(-1) + sk(1,-1)*U(1)*V(-1) + sk(2,-1)*U(2)*V(-1)+ &
      &    sk(0,1)*U(0)*V(1)+sk(-2,1)*U(-2)*V(1)+sk(-1,1)*U(-1)*V(1) +  & 
      &    sk(1,1)*U(1)*V(1) + sk(0,2)*U(0)*V(2)+sk(-1,2)*U(-1)*V(2) +  & 
      &    sk(2,2)*U(2)*V(2)+ sk(0,3)*U(0)*V(3)+sk(3,3)*U(3)*V(3) 
#else
    flux = r13*( U(0) + U(1) )  * ( V(0) + V(1)) - & 
      &  r1_24* (U(-1)*V(-1) + U(-1)*V(1) + U(0)*V(0) + U(0)*V(2) +  & 
      &          U(1)*V(1) + U(1)*V(-1) + U(2)*V(0) +U(2)*V(2) )
#endif
  end function flux_advskew 

  !***************************************************************************! 
  ! art visc correction to the fluxes
  function flux_dis(Uim1,Ui,Uip1,Uip2,e2,e4) result(flux)
    
    implicit none
  
    real(rt), intent(in) :: Uim1,Ui,Uip1,Uip2,e2,e4
    real(rt) :: flux
  
    flux = e2*( Uip1 - Ui)  - e4*(Uip2 - 3.0d0*Uip1 + 3.0d0*Ui - Uim1)
      
  end function flux_dis 
  !***************************************************************************! 
  ! art visc correction to the fluxes
  function flux_shock(W,e2) result(flux)
    
    implicit none
  
    real(rt), intent(in) :: W(iskew1:)
    real(rt), intent(in) :: e2
    real(rt) :: flux
  
    flux = e2*( W(1) - W(0))     
      
  end function flux_shock
  !***************************************************************************! 
  ! high frequency damping correction to the fluxes 
  ! (to apply in smooth parts of the flow)
  function flux_damp(W,e4) result(flux)
    
    implicit none
    
    real(rt), intent(in) :: W(iskew1:)
    real(rt), intent(in) :: e4
    real(rt) :: flux
    
#if (flux_skew==2)
    flux = 0.0d0 
#elif (flux_skew==6)   
    flux   =  r1_12*e4*(-W(3) + 17.0d0*W(2)-46.0d0*W(1) + & 
    & 46.0d0*W(0) -17.0d0*W(-1) + W(-2) )  
#else
    flux   =  e4*( W(2) - 3.0d0*W(1) + 3.0d0*W(0)- W(-1) )
#endif
     
  end function flux_damp 
  !***************************************************************************!
  ! calculates directonal sensor 
  !  NEW *********** NEED TO REWRITE
  subroutine calculate_sensor_skew(qvar,sens,solido)

    implicit none
  
    real(rt), intent(in)    :: qvar(nxlo:nxup,nylo:nyup,nzlo:nzup,nvar)
    real(rt), intent(inout) :: sens(nxlo:nxup,nylo:nyup,nzlo:nzup,ndim)
    
    ! local vars
    integer :: i,j,k,ip,is,js,je,ie,ks,ke
    real(rt) :: Ml   !local mach number
    real(rt) :: Mref !reference mach number
    real(rt) :: vel(nxlo:nxup,nylo:nyup,nzlo:nzup)

    do k = ks,ke
    do j = js,je
    do i = is,ie
        sens(i,j,k,ip) = detect_ad(ip,i,j,k,qvar)
    end do
    end do
    end do
    end do  

    !scale sensor by mach number and sensor on vel also
#if sensor_scale
#if twodim 
    vel  = (qvar(:,:,:,nvux)**2 + qvar(:,:,:,nvuy)**2 )**r12
#elif threedim 
    vel  = (qvar(:,:,:,nvux)**2 + qvar(:,:,:,nvuy)**2 + qvar(:,:,:,nvuz)**2)**r12
#endif
    Mref = artdif_param%Mref
    
    do k = ks,ke
    do j = js,je
    do i = is,ie
      call q2ths(qvar(i,j,k,:),q_ths)
      Ml = vel(i,j,k)/ths_sound(q_ths) 
      sens(i,j,k,:) = sens(i,j,k,:)*Ml/Mref
    end do
    end do
    end do  
#endif
    
  end subroutine calculate_sensor_skew
  !***************************************************************************! 
  ! detector P, rho and Y
  ! psi = sum(psi_k^2) / sum(psi_k)
  function detect_ad(ip,i,j,k,qvar) result(dd)

    implicit none
  
    integer,  intent(in) :: ip,i,j,k,nv
    real(rt), intent(in) :: qvar(nxlo:nxup,nylo:nyup,nzlo:nzup,nvar)
    real(rt) :: auxd, sum1,sum2,dd
    sum1  = detect_discon(ip,i,j,k,qvar(:,:,:,nvpr)) + 1.0e-8d0
    sum2  = sum1*sum1
#if (skewsensor_pres == 0)  
    auxd  = detect_discon(ip,i,j,k,qvar(:,:,:,nvrh))
    sum1  = sum1  + auxd
    sum2  = sum2  + auxd*auxd
#endif 
    ! scalars 
#if complexchem
    do nv = nvsp1,nvsp2
      auxd = detect_discon(ip,i,j,k,qvar(:,:,:,nv))
      sum1 = sum1 + auxd
      sum2 = sum2 + auxd*auxd
    end do 
#endif

    dd  =  sum2/sum1
     
  end function detect_ad

  !***************************************************************************!
  ! initialise skew coefficients (called  from init_flow )
  !
  subroutine init_skew

    implicit none
 

#if flux_skew==6
    sk(:,:) = 0.0d0

    sk(0,0) = 37.0d0/120.0d0
    sk(1,0) = 45.0d0/120.0d0
    sk(2,0) = -9.0d0/120.0d0
    sk(3,0) = 1.0d0/120.0d0

    sk(-2,-2) = 1.0d0/120.0d0
    sk( 1,-2) = 1.0d0/120.0d0

    sk(-1,-1) = -8.0d0/120.0d0
    sk(+1,-1) = -9.0d0/120.0d0
    sk(+2,-1) = 1.0d0/120.0d0

    sk(+0,+1)  = 45.0d0/120.0d0
    sk(-2,+1) = 1.0d0/120.0d0
    sk(-1,+1) = -9.0d0/120.0d0
    sk(+1,+1) = 37.0d0/120.0d0

    sk(+0,+2) = -9.0d0/120.0d0
    sk(-1,+2) = 1.0d0/120.0d0
    sk(+2,+2) = -8.0d0/120.0d0
    
    sk(+0,+3) = 1.0d0/120.0d0
    sk(+3,+3) = 1.0d0/120.0d0
#endif

  end subroutine init_skew  

!----------------------------------------------------------------------------------------
end module skew_module  
