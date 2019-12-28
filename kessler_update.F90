module kessler_update

! use ccpp_kinds, only: kind_phys => kind_phys
use ccpp_kinds, only: kind_phys
use geopotential, only: geopotential_t

implicit none
private
public :: kessler_update_timestep_init, kessler_update_run, kessler_update_finalize

contains

!> \section arg_table_kessler_update_timestep_init  Argument Table
!! \htmlinclude kessler_update_timestep_init.html
  subroutine kessler_update_timestep_init(temp, temp_prev, ttend_t, errmsg, errflg)

    real(kind_phys), intent(in)    :: temp(:,:)
    real(kind_phys), intent(out)   :: temp_prev(:,:)
    real(kind_phys), intent(out)   :: ttend_t(:,:)
    character(len=512),      intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

!   Initialize the previous temperature and its tendency to zero
    temp_prev(:,:)  = temp(:,:)
    ttend_t(:,:)    = 0._kind_phys

    errmsg = ' '
    errflg = 0

  end subroutine kessler_update_timestep_init

!> \section arg_table_kessler_update_run  Argument Table
!! \htmlinclude kessler_update_run.html
  subroutine kessler_update_run(nz, pcols, ncol, gravit, cpair, rair, zvir, phis, temp, &
                 lnpint, lnpmid, pint, pmid, pdel, rpdel, qc, theta, exner, dt,  &
                 zi, zm, temp_prev, ttend_t, st_energy, errmsg, errflg )

    integer, intent(in)     :: nz
    integer, intent(in)     :: pcols
    integer, intent(in)     :: ncol
    real(kind_phys), intent(in)    :: gravit
    real(kind_phys), intent(in)    :: cpair
    real(kind_phys), intent(in)    :: rair
    real(kind_phys), intent(in)    :: zvir
    real(kind_phys), intent(in)    :: phis(:)
    real(kind_phys), intent(in)    :: temp(:,:)   ! temperature
    real(kind_phys), intent(in)    :: lnpint(:,:)
    real(kind_phys), intent(in)    :: lnpmid(:,:)
    real(kind_phys), intent(in)    :: pint(:,:)
    real(kind_phys), intent(in)    :: pmid(:,:)
    real(kind_phys), intent(in)    :: pdel(:,:)
    real(kind_phys), intent(in)    :: rpdel(:,:)
    real(kind_phys), intent(in)    :: qc(:,:)
    real(kind_phys), intent(in)    :: exner(:,:)
    real(kind_phys), intent(in)    :: dt

    real(kind_phys), intent(inout) :: theta(:,:)
    real(kind_phys), intent(inout) :: zi(:,:)
    real(kind_phys), intent(inout) :: zm(:,:)
    real(kind_phys), intent(inout) :: temp_prev(:,:)
    real(kind_phys), intent(inout) :: ttend_t(:,:)
    real(kind_phys), intent(inout) :: st_energy(:,:)

    character(len=512),      intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    integer                 :: k, rk
    integer                 :: vert_surf, vert_toa
    real(kind_phys)                :: rairv(pcols,nz)
    real(kind_phys)                :: zvirv(pcols,nz)
    real(kind_phys)                :: ptend_s(pcols,nz)

    errmsg = ' '
    errflg = 0

    rairv(:,:) = rair
    zvirv(:,:) = zvir

    ! Back out tendencies from updated fields
    do k = 1, nz
      ptend_s(:ncol,k) = (theta(:ncol,k) * exner(:ncol,k) - temp_prev(:ncol,k)) * cpair / dt
      ttend_t(:ncol,k) = ttend_t(:ncol,k) + ptend_s(:ncol,k)/cpair
    end do
 
    ! Kessler is bottom to top
    vert_toa = 30
    vert_surf = 1

    call geopotential_t  ( nz, nz+1, .true., vert_surf, vert_toa, &
            lnpint,    lnpmid,    pint  , pmid  , pdel  , rpdel  , &
            temp     , qc, rairv, gravit  , zvirv              , &
            zi    , zm      , ncol         )

    do k = 1, nz
      st_energy(:ncol,k) = temp(:ncol,k)*cpair+gravit*zm(:ncol,k)+phis(:ncol)
    end do


!    surf_state%precl(:ncol) = surf_state%precl(:ncol) + precl(:ncol)  ! KEEPING THIS HERE AS A REMINDER

    ! Set the previous q values to the current q
    temp_prev(:,:)     = temp(:,:)

  end subroutine kessler_update_run

!> \section arg_table_kessler_update_finalize  Argument Table
!! \htmlinclude kessler_update_finalize.html
  subroutine kessler_update_finalize(errmsg, errflg)

    character(len=512),      intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    errmsg = ' '
    errflg = 0

  end subroutine kessler_update_finalize
end module kessler_update
