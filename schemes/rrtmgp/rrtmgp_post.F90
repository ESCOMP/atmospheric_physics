module rrtmgp_post

  implicit none
  private

  public :: rrtmgp_post_run

contains
!> \section arg_table_rrtmgp_post_run Argument Table
!! \htmlinclude rrtmgp_post_run.html
!!
subroutine rrtmgp_post_run(qrs_prime, qrl_prime, fsns, pdel, atm_optics_sw, cloud_sw, aer_sw,  &
                  fsw, fswc, atm_optics_lw, sources_lw, cloud_lw, aer_lw, flw, flwc, qrs, qrl, &
                  netsw, errmsg, errflg)
   use ccpp_kinds,             only: kind_phys
   use ccpp_optical_props,     only: ty_optical_props_1scl_ccpp, ty_optical_props_2str_ccpp
   use ccpp_source_functions,  only: ty_source_func_lw_ccpp
   use ccpp_fluxes,            only: ty_fluxes_broadband_ccpp
   use ccpp_fluxes_byband,     only: ty_fluxes_byband_ccpp
   real(kind_phys), dimension(:,:),  intent(in)    :: pdel           ! Layer thickness [Pa]
   real(kind_phys), dimension(:),    intent(in)    :: fsns           ! Surface net shortwave flux [W m-2]
   real(kind_phys), dimension(:,:),  intent(in)    :: qrs_prime      ! Shortwave heating rate [J kg-1 s-1]
   real(kind_phys), dimension(:,:),  intent(in)    :: qrl_prime      ! Longwave heating rate [J kg-1 s-1]
   type(ty_optical_props_1scl_ccpp), intent(inout) :: atm_optics_lw  ! Atmosphere optical properties object (longwave)
   type(ty_optical_props_2str_ccpp), intent(inout) :: atm_optics_sw  ! Atmosphere optical properties object (shortwave)
   type(ty_optical_props_1scl_ccpp), intent(inout) :: aer_lw         ! Aerosol optical properties object (longwave)
   type(ty_optical_props_2str_ccpp), intent(inout) :: aer_sw         ! Aerosol optical properties object (shortwave)
   type(ty_optical_props_1scl_ccpp), intent(inout) :: cloud_lw       ! Cloud optical properties object (longwave)
   type(ty_optical_props_2str_ccpp), intent(inout) :: cloud_sw       ! Cloud optical properties object (shortwave)
   type(ty_fluxes_broadband_ccpp),   intent(inout) :: fswc           ! Shortwave clear-sky flux object
   type(ty_fluxes_broadband_ccpp),   intent(inout) :: flwc           ! Longwave clear-sky flux object
   type(ty_fluxes_byband_ccpp),      intent(inout) :: fsw            ! Shortwave all-sky flux object
   type(ty_fluxes_byband_ccpp),      intent(inout) :: flw            ! Longwave all-sky flux object
   type(ty_source_func_lw_ccpp),     intent(inout) :: sources_lw     ! Longwave sources object
   real(kind_phys), dimension(:,:),  intent(out)   :: qrs            ! Shortwave heating rate adjusted by air pressure thickness [J Pa kg-1 s-1]
   real(kind_phys), dimension(:,:),  intent(out)   :: qrl            ! Longwave heating rate adjusted by air pressure thickness [J Pa kg-1 s-1]
   real(kind_phys), dimension(:),    intent(out)   :: netsw          ! Net shortwave flux to be sent to coupler [W m-2]
   character(len=*),                 intent(out)   :: errmsg
   integer,                          intent(out)   :: errflg

   ! Set error varaibles
   errflg = 0
   errmsg = ''
   ! The radiative heating rates are maintained across multiple physics timesteps
   ! as Q*dp (for energy conservation).
   qrs(:,:) = qrs_prime(:,:) * pdel(:,:)
   qrl(:,:) = qrl_prime(:,:) * pdel(:,:)

   ! Set the netsw to be sent to the coupler
   netsw(:) = fsns(:)

   call free_optics_sw(atm_optics_sw, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   call free_optics_sw(cloud_sw, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   call free_optics_sw(aer_sw, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   call free_fluxes_byband(fsw)
   call free_fluxes_broadband(fswc)

   call sources_lw%sources%finalize()
   call free_optics_lw(atm_optics_lw, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   call free_optics_lw(cloud_lw, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   call free_optics_lw(aer_lw, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   call free_fluxes_byband(flw)
   call free_fluxes_broadband(flwc)

end subroutine rrtmgp_post_run

  !=========================================================================================

subroutine free_optics_sw(optics, errmsg, errflg)
   use ccpp_optical_props,      only: ty_optical_props_2str_ccpp

   type(ty_optical_props_2str_ccpp), intent(inout) :: optics
   character(len=*),                 intent(out)   :: errmsg
   integer,                          intent(out)   :: errflg

   errmsg = ''
   errflg = 0

   errmsg = optics%optical_props%finalize_2str()
   if (len_trim(errmsg) /= 0) then
      errflg = 1
      return
   end if

   call optics%optical_props%finalize()

end subroutine free_optics_sw

!=========================================================================================

subroutine free_optics_lw(optics, errmsg, errflg)
   use ccpp_optical_props,      only: ty_optical_props_1scl_ccpp

   type(ty_optical_props_1scl_ccpp), intent(inout) :: optics
   character(len=*),                 intent(out)   :: errmsg
   integer,                          intent(out)   :: errflg

   errmsg = optics%optical_props%finalize_1scl()
   if (len_trim(errmsg) /= 0) then
      errflg = 1
      return
   end if

   call optics%optical_props%finalize()

end subroutine free_optics_lw

!=========================================================================================

subroutine free_fluxes_broadband(fluxes)
   use ccpp_fluxes,            only: ty_fluxes_broadband_ccpp

   class(ty_fluxes_broadband_ccpp), intent(inout) :: fluxes

   if (associated(fluxes%fluxes%flux_up)) deallocate(fluxes%fluxes%flux_up)
   if (associated(fluxes%fluxes%flux_dn)) deallocate(fluxes%fluxes%flux_dn)
   if (associated(fluxes%fluxes%flux_net)) deallocate(fluxes%fluxes%flux_net)
   if (associated(fluxes%fluxes%flux_dn_dir)) deallocate(fluxes%fluxes%flux_dn_dir)

end subroutine free_fluxes_broadband

!=========================================================================================

subroutine free_fluxes_byband(fluxes)
   use ccpp_fluxes_byband,     only: ty_fluxes_byband_ccpp

   class(ty_fluxes_byband_ccpp), intent(inout) :: fluxes

   if (associated(fluxes%fluxes%flux_up)) deallocate(fluxes%fluxes%flux_up)
   if (associated(fluxes%fluxes%flux_dn)) deallocate(fluxes%fluxes%flux_dn)
   if (associated(fluxes%fluxes%flux_net)) deallocate(fluxes%fluxes%flux_net)
   if (associated(fluxes%fluxes%flux_dn_dir)) deallocate(fluxes%fluxes%flux_dn_dir)

   if (associated(fluxes%fluxes%bnd_flux_up)) deallocate(fluxes%fluxes%bnd_flux_up)
   if (associated(fluxes%fluxes%bnd_flux_dn)) deallocate(fluxes%fluxes%bnd_flux_dn)
   if (associated(fluxes%fluxes%bnd_flux_net)) deallocate(fluxes%fluxes%bnd_flux_net)
   if (associated(fluxes%fluxes%bnd_flux_dn_dir)) deallocate(fluxes%fluxes%bnd_flux_dn_dir)

end subroutine free_fluxes_byband

end module rrtmgp_post
