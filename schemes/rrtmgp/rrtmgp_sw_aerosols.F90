!> \file rrtmgp_sw_aerosols.F90
!!

!> This module sets the RRTMGP aerosol shortwave optical properties
module rrtmgp_sw_aerosols

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public rrtmgp_sw_aerosols_run

  ! Mapping from RRTMG shortwave bands to RRTMGP
  integer, parameter, dimension(14) :: rrtmg_to_rrtmgp_swbands = &
     [ 14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 ]

  real(kind_phys), parameter :: tiny = 1.0e-80_kind_phys

contains

!> \section arg_table_rrtmgp_sw_aerosols_run Argument Table
!! \htmlinclude rrtmgp_sw_aerosols_run.html
!!
   subroutine rrtmgp_sw_aerosols_run( &
     doswrad, ncol, pver, nswbands, nday, idxday, ktopcam, ktoprad, &
     aer_tau, aer_tau_w, aer_tau_w_g, &
     aer_sw, errmsg, errflg)

    use ccpp_optical_props, only: ty_optical_props_2str_ccpp

    ! Inputs
    logical,          intent(in)    :: doswrad                    !< do shortwave calculation
    integer,          intent(in)    :: ncol                       !< number of columns
    integer,          intent(in)    :: pver                       !< number of vertical layers
    integer,          intent(in)    :: nswbands                   !< number of SW bands
    integer,          intent(in)    :: nday                       !< daytime columns count
    integer,          intent(in)    :: idxday(:)                  !< daytime column indices
    integer,          intent(in)    :: ktopcam                    !< top CAM level for RRTMGP
    integer,          intent(in)    :: ktoprad                    !< top RRTMGP level

    real(kind_phys),  intent(in)    :: aer_tau(:, :, :)           !< SW extinction OD (ncol,pver,nswbands)
    real(kind_phys),  intent(in)    :: aer_tau_w(:, :, :)         !< SW ssa*tau (ncol,pver,nswbands)
    real(kind_phys),  intent(in)    :: aer_tau_w_g(:, :, :)       !< SW asy*ssa*tau (ncol,pver,nswbands)

    ! Outputs
    class(ty_optical_props_2str_ccpp), intent(inout) :: aer_sw    !< Aerosol optical properties object

    character(len=*), intent(out) :: errmsg                       !< CCPP error message
    integer,          intent(out) :: errflg                       !< CCPP error flag

    ! Local variables
    integer :: i
    real(kind_phys) :: tau_reordered(ncol, pver, nswbands)
    real(kind_phys) :: tau_w_reordered(ncol, pver, nswbands)
    real(kind_phys) :: tau_w_g_reordered(ncol, pver, nswbands)

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doswrad .or. nday == 0) return

    ! Reorder from RRTMG band order to RRTMGP band order
    tau_reordered(:ncol, :, :)     = aer_tau(:ncol, :, rrtmg_to_rrtmgp_swbands)
    tau_w_reordered(:ncol, :, :)   = aer_tau_w(:ncol, :, rrtmg_to_rrtmgp_swbands)
    tau_w_g_reordered(:ncol, :, :) = aer_tau_w_g(:ncol, :, rrtmg_to_rrtmgp_swbands)

    ! Initialize defaults: tau=0, ssa=1, g=0
    aer_sw%optical_props%tau = 0.0_kind_phys
    aer_sw%optical_props%ssa = 1.0_kind_phys
    aer_sw%optical_props%g   = 0.0_kind_phys

    ! Subset to daytime columns and map CAM vertical to RRTMGP vertical
    do i = 1, nday
      ! Set aerosol optical depth, clip to zero
      aer_sw%optical_props%tau(i, ktoprad:, :) = &
        max(tau_reordered(idxday(i), ktopcam:, :), 0._kind_phys)

      ! Set single scattering albedo: ssa = tau_w / tau
      where (tau_reordered(idxday(i), ktopcam:, :) > 0._kind_phys)
        aer_sw%optical_props%ssa(i, ktoprad:, :) = &
          tau_w_reordered(idxday(i), ktopcam:, :) &
          / tau_reordered(idxday(i), ktopcam:, :)
      elsewhere
        aer_sw%optical_props%ssa(i, ktoprad:, :) = 1._kind_phys
      end where

      ! Set asymmetry parameter: g = tau_w_g / tau_w
      where (tau_w_reordered(idxday(i), ktopcam:, :) > tiny)
        aer_sw%optical_props%g(i, ktoprad:, :) = &
          tau_w_g_reordered(idxday(i), ktopcam:, :) &
          / tau_w_reordered(idxday(i), ktopcam:, :)
      elsewhere
        aer_sw%optical_props%g(i, ktoprad:, :) = 0._kind_phys
      end where
    end do

    ! Impose limits on the components
    aer_sw%optical_props%ssa = min(max(aer_sw%optical_props%ssa, 0._kind_phys), 1._kind_phys)
    aer_sw%optical_props%g   = min(max(aer_sw%optical_props%g, -1._kind_phys), 1._kind_phys)

  end subroutine rrtmgp_sw_aerosols_run
end module rrtmgp_sw_aerosols
