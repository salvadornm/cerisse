module cns_tagging_module
  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: cns_tag_denerror,cns_tag_presserror

contains

  subroutine cns_tag_denerror (lo, hi, tag, tlo, thi, rho, rlo, rhi, flag, flo, fhi, &
       dengrad, tagval, clearval) bind(c,name='cns_tag_denerror')
    use iso_c_binding, only : c_char
    use amrex_ebcellflag_module, only : is_regular_cell
    integer, dimension(3), intent(in) :: lo, hi, tlo, thi, rlo, rhi, flo, fhi
    character(kind=c_char), intent(inout) :: tag(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3))
    real(kind=amrex_real),  intent(in)    :: rho(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    integer,                intent(in)   :: flag(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    real(kind=amrex_real), intent(in) :: dengrad
    character(kind=c_char), intent(in) :: tagval, clearval

    integer :: i,j,k
    real(amrex_real) :: ax, ay, az

   !  write(*,*) " [cns_tag_denerror] lo ",lo
   !  write(*,*) "         ",hi
   !  write(*,*) "tlo=",tlo
   !  write(*,*) "thi=",thi
   !  write(*,*) "rlo=",rlo
   !  write(*,*) "rhi=",rhi
   !  write(*,*) "tag(1,1,1) = ",tag(1,1,1)
   !  write(*,*) "  tagval clearval  = ",tagval,clearval
   !  write(*,*) " dengrad (from input)=",dengrad
    
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (is_regular_cell(flag(i,j,k))) then
                ax = ABS(rho(i+1,j,k) - rho(i,j,k))
                ay = ABS(rho(i,j+1,k) - rho(i,j,k))
                az = ABS(rho(i,j,k+1) - rho(i,j,k))
                ax = MAX(ax,ABS(rho(i,j,k) - rho(i-1,j,k)))
                ay = MAX(ay,ABS(rho(i,j,k) - rho(i,j-1,k)))
                az = MAX(az,ABS(rho(i,j,k) - rho(i,j,k-1)))
                if ( MAX(ax,ay,az) .ge. dengrad ) then
                   tag(i,j,k) = tagval
                endif
             end if
          end do
       end do
    end do

  end subroutine cns_tag_denerror

  ! -------------------------------------------------------------
  ! create one for pressure 
  ! -------------------------------------------------------------
  subroutine cns_tag_presserror (lo, hi, tag, tlo, thi, pres, rlo, rhi, flag, flo, fhi, &
    pgrad, tagval, clearval) bind(c,name='cns_tag_presserror')
 use iso_c_binding, only : c_char
 use amrex_ebcellflag_module, only : is_regular_cell
 integer, dimension(3), intent(in) :: lo, hi, tlo, thi, rlo, rhi, flo, fhi
 character(kind=c_char), intent(inout) :: tag(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3))
 real(kind=amrex_real),  intent(in)    :: pres(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
 integer,                intent(in)   :: flag(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
 real(kind=amrex_real), intent(in) :: pgrad
 character(kind=c_char), intent(in) :: tagval, clearval

 integer :: i,j,k
 real(amrex_real) :: ax, ay, az

!  write(*,*) " [cns_tag_preserror] lo ",lo
!  write(*,*) "         ",hi
!  write(*,*) "tlo=",tlo
!  write(*,*) "thi=",thi
!  write(*,*) "rlo=",rlo
!  write(*,*) "rhi=",rhi
!  write(*,*) "tag(1,1,1) = ",tag(1,1,1)
!  write(*,*) "  tagval clearval  = ",tagval,clearval
!  write(*,*) " dengrad (from input)=",dengrad
 
 do       k = lo(3), hi(3)
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (is_regular_cell(flag(i,j,k))) then
             ax = ABS(pres(i+1,j,k) - pres(i,j,k))
             ay = ABS(pres(i,j+1,k) - pres(i,j,k))
             az = ABS(pres(i,j,k+1) - pres(i,j,k))
             ax = MAX(ax,ABS(pres(i,j,k) - pres(i-1,j,k)))
             ay = MAX(ay,ABS(pres(i,j,k) - pres(i,j-1,k)))
             az = MAX(az,ABS(pres(i,j,k) - pres(i,j,k-1)))
             if ( MAX(ax,ay,az) .ge. pgrad ) then
                tag(i,j,k) = tagval
             endif
          end if
       end do
    end do
 end do

end subroutine cns_tag_presserror

end module cns_tagging_module
