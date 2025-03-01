module rayleigh_friction

!---------------------------------------------------------------------------------
! Module to apply rayleigh friction in region of model top.
! We specify a decay rate profile that is largest at the model top and
! drops off vertically using a hyperbolic tangent profile.
! We compute the tendencies in u and v using an Euler backward scheme.
! We then apply the negative of the kinetic energy tendency to "s", the dry
! static energy.
!---------------------------Code history--------------------------------
! This is a new routine written by Art Mirin in collaboration with Phil Rasch.
! Initial coding for this version:  Art Mirin, May 2007.
! First CCPP conversion by Andrew Gettelman, Sep 2021
! Final CCPP completion by Kate T-C, Oct 2024
!---------------------------------------------------------------------------------

  use ccpp_kinds, only:  kind_phys

  implicit none

  save
  private                         ! Make default type private to the module
!
! PUBLIC: interfaces
!
  public rayleigh_friction_init
  public rayleigh_friction_run

! PRIVATE: module variables
integer         :: rayk0 = 2           ! vertical level at which rayleigh friction term is centered
real(kind_phys) :: raykrange = 0._kind_phys   ! range of rayleigh friction profile 
                                ! if 0, range is set to satisfy x=2 (see below)
real(kind_phys) :: raytau0 = 0._kind_phys     ! approximate value of decay time at model top (days)
                                ! if 0., no rayleigh friction is applied
! Local
real (kind_phys) :: krange         ! range of rayleigh friction profile 
real (kind_phys) :: tau0           ! approximate value of decay time at model top
real (kind_phys) :: otau0          ! inverse of tau0
real (kind_phys), allocatable :: otau(:)     ! inverse decay time versus vertical level

! We apply a profile of the form otau0 * [1 + tanh (x)] / 2 , where
! x = (k0 - k) / krange. The default is for x to equal 2 at k=1, meaning
! krange = (k0 - 1) / 2. The default is applied when raykrange is set to 0.
! If otau0 = 0, no term is applied.

contains

  !> \section arg_table_rayleigh_friction_init Argument Table
  !! \htmlinclude rayleigh_friction_init.html
  subroutine rayleigh_friction_init(pver, raytau0_nl, raykrange_nl, rayk0_nl, masterproc, iulog, errmsg, errflg)

    integer, intent(in)                  :: pver
    real (kind_phys), intent(in)         :: raytau0_nl
    real (kind_phys), intent(in)         :: raykrange_nl
    integer, intent(in)                  :: rayk0_nl
    logical, intent(in)                  :: masterproc
    integer, intent(in)                  :: iulog
    character(len=512), intent(out)      :: errmsg
    integer, intent(out)                 :: errflg

    ! local variables
    character(len=*), parameter :: subname = 'rayleigh_friction_init'
    real (kind_phys) x
    integer k, ierr

    errmsg = ''
    errflg = 0
    raytau0 = raytau0_nl
    raykrange = raykrange_nl
    rayk0 = rayk0_nl

    if ( raytau0 > 0._kind_phys) then  ! Rayleigh Friction enabled
       ! Allocate module data
       allocate(otau(pver), stat=ierr)
       if (ierr /= 0) then
          errflg = ierr
          errmsg = trim(subname)//': Allocate of otau failed'
          return
       end if

       !-----------------------------------------------------------------------
       ! Compute tau array
       !-----------------------------------------------------------------------

       krange = raykrange
       if (raykrange .eq. 0._kind_phys) krange = (rayk0 - 1) / 2._kind_phys

       tau0 = (86400._kind_phys) * raytau0   ! convert to seconds
       otau0 = 0._kind_phys
       if (tau0 .ne. 0._kind_phys) otau0 = 1._kind_phys/tau0

       do k = 1, pver
          x = (rayk0 - k) / krange
          otau(k) = otau0 * (1 + tanh(x)) / (2._kind_phys)
       enddo
    end if ! Rayleigh Friction enabled

    if ( masterproc ) then
       if (raytau0 > 0._kind_phys) then
          write(iulog,*) 'tuning parameters rayleigh_friction_init: rayk0',rayk0
          write(iulog,*) 'tuning parameters rayleigh_friction_init: raykrange',raykrange
          write(iulog,*) 'tuning parameters rayleigh_friction_init: raytau0',raytau0
       else
          write(iulog,*) 'Rayleigh Friction not enabled.'
       endif
    endif

  end subroutine rayleigh_friction_init

  !=========================================================================================
  !> \section arg_table_rayleigh_friction_run Argument Table
  !! \htmlinclude rayleigh_friction_run.html
  subroutine rayleigh_friction_run(ncol, pver, ztodt, u, v, dudt, dvdt, dsdt, errmsg, errflg)

    !------------------------------Arguments--------------------------------
    integer,             intent(in) :: ncol  
    integer,             intent(in) :: pver
    real(kind_phys),     intent(in) :: ztodt   !physics timestep
    real(kind_phys),     intent(in) :: u(:,:)     
    real(kind_phys),     intent(in) :: v(:,:)     
    real(kind_phys),    intent(out) :: dudt(:,:) !tendency_of_x_wind    
    real(kind_phys),    intent(out) :: dvdt(:,:) !tendency_of_y_wind    
    real(kind_phys),    intent(out) :: dsdt(:,:)  !heating_rate 
      
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg
   
    !---------------------------Local storage-------------------------------
    character(len=*),     parameter :: subname = 'rayleigh_friction_run'
    integer :: k                                   ! level
    real(kind_phys) :: rztodt                      ! 1./ztodt
    real(kind_phys) :: c1, c2, c3                  ! temporary variables
    !-----------------------------------------------------------------------
    
    ! initialize values
    errmsg = ''
    errflg = 0
    dudt(:,:)=0._kind_phys
    dvdt(:,:)=0._kind_phys
    dsdt(:,:) =0._kind_phys

    if (raytau0 .eq. 0._kind_phys) return

    rztodt = 1._kind_phys/ztodt

    ! u, v and s are modified by rayleigh friction

    do k = 1, pver
       c2 = 1._kind_phys / (1._kind_phys + otau(k)*ztodt)
       c1 = -otau(k) * c2
       c3 = 0.5_kind_phys * (1._kind_phys - c2*c2) * rztodt
       dudt(:ncol,k) = c1 * u(:ncol,k)
       dvdt(:ncol,k) = c1 * v(:ncol,k)
       dsdt(:ncol,k) = c3 * (u(:ncol,k)**2 + v(:ncol,k)**2)
    enddo
    
  end subroutine rayleigh_friction_run

end module rayleigh_friction

