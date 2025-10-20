module rrtmgp_sw_calculate_fluxes

   use ccpp_kinds, only:  kind_phys
   implicit none
   private
   save

   public :: rrtmgp_sw_calculate_fluxes_run  ! main routine


CONTAINS

   !> \section arg_table_rrtmgp_sw_calculate_fluxes_run  Argument Table
   !! \htmlinclude rrtmgp_sw_calculate_fluxes_run.html
   subroutine rrtmgp_sw_calculate_fluxes_run(num_diag_subcycles, icall, ncol, pverp, nlay, nday, idxday, ktopcam, ktoprad, &
      active_calls, fsw, fswc, fns, fcns, fsns, fsnt, soll, sols, solld, solsd, errmsg, errflg)

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
      integer,                        intent(in) :: nday                ! Daytime points dimension
      integer,                        intent(in) :: ktopcam             ! Index in CAM arrays of top level (layer or interface) at which RRTMGP is active
      integer,                        intent(in) :: ktoprad             ! Index in RRTMGP array corresponding to top layer or interface of CAM arrays
      integer,                        intent(in) :: idxday(:)           ! Daytime points indices
      logical,                        intent(in) :: active_calls(:)     ! Logical array of flags for whether a specified subcycle is active
      type(ty_fluxes_byband_ccpp),    intent(in) :: fsw                 ! Shortwave all-sky flux object
      type(ty_fluxes_broadband_ccpp), intent(in) :: fswc                ! Shortwave clear-sky flux object
      ! Output variables
      real(kind_phys),               intent(out) :: fns(:,:)            ! Shortwave net radiative flux [W m-2]
      real(kind_phys),               intent(out) :: fcns(:,:)           ! Shortwave net radiative clear-sky flux [W m-2]
      real(kind_phys),               intent(out) :: fsns(:)             ! Shortwave net upward flux at surface [W m-2]
      real(kind_phys),               intent(out) :: fsnt(:)             ! Shortwave net outgoing flux at model top [W m-2]
      real(kind_phys),               intent(out) :: soll(:)             ! Direct solar radiative flux at surface >= 700nm [W m-2]
      real(kind_phys),               intent(out) :: sols(:)             ! Direct solar radiative flux at surface < 700nm [W m-2]
      real(kind_phys),               intent(out) :: solld(:)            ! Diffuse solar radiative flux at surface >= 700nm [W m-2]
      real(kind_phys),               intent(out) :: solsd(:)            ! Diffuse solar radiative flux at surface < 700nm [W m-2]

      
      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables
      integer :: diag_index, idx
      real(kind_phys), dimension(size(fsw%fluxes%bnd_flux_dn,1), &
                                 size(fsw%fluxes%bnd_flux_dn,2), &
                                 size(fsw%fluxes%bnd_flux_dn,3)) :: flux_dn_diffuse

      errmsg = ''
      errflg = 0

      diag_index = num_diag_subcycles - icall + 1

      ! Don't do anything if this subcycle is inactive
      if (.not. active_calls(diag_index)) then
         return
      end if

      ! Initialize to provide 0.0 values for night columns.
      fns = 0.0_kind_phys
      fcns = 0.0_kind_phys
      fsns = 0.0_kind_phys
      fsnt = 0.0_kind_phys
      soll = 0.0_kind_phys
      sols = 0.0_kind_phys
      solld = 0.0_kind_phys
      solsd = 0.0_kind_phys

      do idx = 1, nday
         fns(idxday(idx), ktopcam:) = fsw%fluxes%flux_net(idx, ktoprad:)
         fcns(idxday(idx), ktopcam:) = fswc%fluxes%flux_net(idx, ktoprad:)
      end do

      fsns(:) = fns(:, pverp)
      fsnt(:) = fns(:, ktopcam)

      ! Export surface fluxes
      ! sols(pcols)      Direct solar rad on surface (< 0.7)
      ! soll(pcols)      Direct solar rad on surface (>= 0.7)
      ! RRTMGP: Near-IR bands (1-10), 820-16000 cm-1, 0.625-12.195 microns
      ! Put half of band 10 in each of the UV/visible and near-IR values,
      ! since this band straddles 0.7 microns:
      ! UV/visible bands 10-13, 16000-50000 cm-1, 0.200-0.625 micron
      ! Calculate diffuse flux from total and direct
      flux_dn_diffuse = fsw%fluxes%bnd_flux_dn - fsw%fluxes%bnd_flux_dn_dir

      do idx = 1, nday
         soll(idxday(idx)) = sum(fsw%fluxes%bnd_flux_dn_dir(idx,nlay+1,1:9)) &
                             + 0.5_kind_phys * fsw%fluxes%bnd_flux_dn_dir(idx,nlay+1,10)

         sols(idxday(idx)) = 0.5_kind_phys * fsw%fluxes%bnd_flux_dn_dir(idx,nlay+1,10)   &
                             + sum(fsw%fluxes%bnd_flux_dn_dir(idx,nlay+1,11:14))

         solld(idxday(idx)) = sum(flux_dn_diffuse(idx,nlay+1,1:9))         &
                              + 0.5_kind_phys * flux_dn_diffuse(idx,nlay+1,10)

         solsd(idxday(idx)) = 0.5_kind_phys * flux_dn_diffuse(idx, nlay+1, 10)    &
                              + sum(flux_dn_diffuse(idx,nlay+1,11:14))
      end do
   end subroutine rrtmgp_sw_calculate_fluxes_run

end module rrtmgp_sw_calculate_fluxes
