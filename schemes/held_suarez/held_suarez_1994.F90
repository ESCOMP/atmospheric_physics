!< \section arg_table_held_suarez_1994
!! \htmlinclude held_suarez_1994.html
module held_suarez_1994
  !-----------------------------------------------------------------------
  !
  ! Purpose: Implement idealized Held-Suarez forcings
  !    Held, I. M., and M. J. Suarez, 1994: 'A proposal for the
  !    intercomparison of the dynamical cores of atmospheric general
  !    circulation models.'
  !    Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
  !
  !-----------------------------------------------------------------------

  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  public :: held_suarez_1994_init
  public :: held_suarez_1994_run

  !!
  !! Forcing parameters
  !!
  real(kind_phys), parameter :: efoldf  =  1._kind_phys  ! efolding time for wind dissipation
  real(kind_phys), parameter :: efolda  = 40._kind_phys  ! efolding time for T dissipation
  real(kind_phys), parameter :: efolds  =  4._kind_phys  ! efolding time for T dissipation
  real(kind_phys), parameter :: sigmab  =  0.7_kind_phys ! threshold sigma level
  real(kind_phys), parameter :: t00     = 200._kind_phys ! minimum reference temperature
  real(kind_phys), parameter :: kf      = 1._kind_phys/(86400._kind_phys*efoldf) ! 1./efolding_time for wind dissipation

  real(kind_phys), parameter :: onemsig = 1._kind_phys - sigmab ! 1. - sigma_reference

  real(kind_phys), parameter :: ka      = 1._kind_phys/(86400._kind_phys * efolda) ! 1./efolding_time for temperature diss.
  real(kind_phys), parameter :: ks      = 1._kind_phys/(86400._kind_phys * efolds)

  !!
  !! Model constants, reset in init call
  !!
  real(kind_phys)            :: pref    = 0.0_kind_phys  ! Surface pressure



!=======================================================================
contains
!=======================================================================

!> \section arg_table_held_suarez_1994_init Argument Table
!! \htmlinclude held_suarez_1994_init.html
  subroutine held_suarez_1994_init(pref_in, errmsg, errflg)

    !! Dummy arguments
    real(kind_phys),    intent(in)  :: pref_in

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ' '
    errflg = 0

    pref   = pref_in

  end subroutine held_suarez_1994_init

!> \section arg_table_held_suarez_1994_run Argument Table
!! \htmlinclude held_suarez_1994_run.html
  subroutine held_suarez_1994_run(pver, ncol, pref_mid_norm, clat, cappa, &
       cpair, pmid, uwnd, vwnd, temp, du, dv, ds, scheme_name, errmsg, errflg)

    !
    ! Input arguments
    !
    integer,  intent(in)  :: pver                    ! Num vertical levels
    integer,  intent(in)  :: ncol                    ! Num active columns
    real(kind_phys), intent(in)  :: pref_mid_norm(:) ! reference pressure normalized by surface pressure
    real(kind_phys), intent(in)  :: clat(:)          ! latitudes(radians) for columns
    real(kind_phys), intent(in)  :: cappa(:,:)       ! ratio of dry air gas constant to specific heat at constant pressure
    real(kind_phys), intent(in)  :: cpair(:,:)       ! specific heat of dry air at constant pressure
    real(kind_phys), intent(in)  :: pmid(:,:)        ! mid-point pressure
    real(kind_phys), intent(in)  :: uwnd(:,:)        ! Zonal wind (m/s)
    real(kind_phys), intent(in)  :: vwnd(:,:)        ! Meridional wind (m/s)
    real(kind_phys), intent(in)  :: temp(:,:)        ! Temperature (K)
    !
    ! Output arguments
    !
    real(kind_phys),   intent(out) :: du(:,:)   ! Zonal wind tend
    real(kind_phys),   intent(out) :: dv(:,:)   ! Meridional wind tend
    real(kind_phys),   intent(out) :: ds(:,:)   ! Heating rate
    character(len=64), intent(out) :: scheme_name
    character(len=512),intent(out):: errmsg
    integer,           intent(out):: errflg
    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i, k          ! Longitude, level indices

    real(kind_phys) :: kv            ! 1./efolding_time (normalized) for wind
    real(kind_phys) :: kt            ! 1./efolding_time for temperature diss.
    real(kind_phys) :: trefa         ! "radiative equilibrium" T
    real(kind_phys) :: trefc         ! used in calc of "radiative equilibrium" T
    real(kind_phys) :: cossq(ncol)   ! coslat**2
    real(kind_phys) :: cossqsq(ncol) ! coslat**4
    real(kind_phys) :: sinsq(ncol)   ! sinlat**2
    real(kind_phys) :: coslat(ncol)  ! cosine(latitude)
    !
    !-----------------------------------------------------------------------
    !

    errmsg = ' '
    errflg = 0
    scheme_name = "HELD_SUAREZ"

    do i = 1, ncol
      coslat (i) = cos(clat(i))
      sinsq  (i) = sin(clat(i))*sin(clat(i))
      cossq  (i) = coslat(i)*coslat(i)
      cossqsq(i) = cossq (i)*cossq (i)
    end do

    !
    !-----------------------------------------------------------------------
    !
    ! Held/Suarez IDEALIZED physics algorithm:
    !
    !   Held, I. M., and M. J. Suarez, 1994: A proposal for the
    !   intercomparison of the dynamical cores of atmospheric general
    !   circulation models.
    !   Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
    !
    !-----------------------------------------------------------------------
    !
    ! Compute idealized radiative heating rates (as dry static energy)
    !
    !
    do k = 1, pver
      if (pref_mid_norm(k) > sigmab) then
        do i = 1, ncol
          kt      = ka + (ks - ka)*cossqsq(i)*(pref_mid_norm(k) - sigmab)/onemsig
          trefc   = 315._kind_phys - (60._kind_phys * sinsq(i))
          trefa   = (trefc - 10._kind_phys*cossq(i)*log((pmid(i,k)/pref)))*(pmid(i,k)/pref)**cappa(i,k)
          trefa   = max(t00,trefa)
          ds(i,k) = (trefa - temp(i,k))*kt*cpair(i,k)
        end do
      else
        do i = 1, ncol
          trefc   = 315._kind_phys - 60._kind_phys*sinsq(i)
          trefa   = (trefc - 10._kind_phys*cossq(i)*log((pmid(i,k)/pref)))*(pmid(i,k)/pref)**cappa(i,k)
          trefa   = max(t00,trefa)
          ds(i,k) = (trefa - temp(i,k))*ka*cpair(i,k)
        end do
      end if
    end do
    !
    ! Add diffusion near the surface for the wind fields
    !
    du(:,:) = 0._kind_phys
    dv(:,:) = 0._kind_phys

    !
    do k = 1, pver
      if (pref_mid_norm(k) > sigmab) then
        kv  = kf*(pref_mid_norm(k) - sigmab)/onemsig
        do i = 1, ncol
          du(i,k) = -kv*uwnd(i,k)
          dv(i,k) = -kv*vwnd(i,k)
        end do
      end if
    end do

  end subroutine held_suarez_1994_run

end module held_suarez_1994
