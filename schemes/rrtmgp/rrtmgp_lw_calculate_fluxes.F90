module rrtmgp_lw_calculate_fluxes

   use ccpp_kinds, only:  kind_phys
   implicit none
   private

   public :: rrtmgp_lw_calculate_fluxes_run  ! main routine


CONTAINS

   !> \section arg_table_rrtmgp_lw_calculate_fluxes_run  Argument Table
   !! \htmlinclude rrtmgp_lw_calculate_fluxes_run.html
   subroutine rrtmgp_lw_calculate_fluxes_run(num_diag_subcycles, icall, ncol, pverp, nlay, ktopcam, ktoprad, &
      active_calls, flw, flwc, flns, flnt, flwds, fnl, fcnl, errmsg, errflg)

      use ccpp_fluxes,        only: ty_fluxes_broadband_ccpp
      use ccpp_fluxes_byband, only: ty_fluxes_byband_ccpp
      !------------------------------------------------
      !   Input / output parameters
      !------------------------------------------------
      integer,                        intent(in) :: num_diag_subcycles  ! Number of diagnostics subcycles
      integer,                        intent(in) :: icall               ! Current diagnostic subcycle
      integer,                        intent(in) :: pverp               ! Number of vertical layer interfaces
      integer,                        intent(in) :: ncol                ! Number of horizontal grid points
      integer,                        intent(in) :: nlay                ! Number of vertical layers in RRTMGP
      integer,                        intent(in) :: ktopcam             ! Index in CAM arrays of top level (layer or interface) at which RRTMGP is active
      integer,                        intent(in) :: ktoprad             ! Index in RRTMGP array corresponding to top layer or interface of CAM arrays
      logical,                        intent(in) :: active_calls(:)     ! Logical array of flags for whether a specified subcycle is active
      type(ty_fluxes_byband_ccpp),    intent(in) :: flw                 ! Longwave all-sky flux object
      type(ty_fluxes_broadband_ccpp), intent(in) :: flwc                ! Longwave clear-sky flux object
      ! Output variables
      real(kind_phys),               intent(out) :: fnl(:,:)            ! Longwave net radiative flux [W m-2]
      real(kind_phys),               intent(out) :: fcnl(:,:)           ! Longwave net radiative clear-sky flux [W m-2]
      real(kind_phys),               intent(out) :: flns(:)             ! Longwave net upward flux at surface [W m-2]
      real(kind_phys),               intent(out) :: flnt(:)             ! Longwave net outgoing flux at model top [W m-2]
      real(kind_phys),               intent(out) :: flwds(:)            ! Longwave downward radiative flux at surface [W m-2]

      
      ! CCPP error handling variables
      character(len=*),   intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables
      integer :: diag_index, idx

      errmsg = ''
      errflg = 0

      diag_index = num_diag_subcycles - icall + 1

      ! Don't do anything if this subcycle is inactive
      if (.not. active_calls(diag_index)) then
         return
      end if

      fnl = 0.0_kind_phys
      fcnl = 0.0_kind_phys

      ! RTE-RRTMGP convention for net is (down - up) **CAM assumes (up - down) !!
      fnl( :,ktopcam:) = -1._kind_phys * flw%fluxes%flux_net(    :, ktoprad:)
      fcnl(:,ktopcam:) = -1._kind_phys * flwc%fluxes%flux_net(   :, ktoprad:)

      flns(:ncol) = fnl(:ncol, pverp)
      flnt(:ncol) = fnl(:ncol, ktopcam)
      flwds(:ncol) = flw%fluxes%flux_dn(:, nlay+1)

   end subroutine rrtmgp_lw_calculate_fluxes_run

end module rrtmgp_lw_calculate_fluxes
