module rrtmgp_sw_mcica_subcol_gen

implicit none
private
save

public :: rrtmgp_sw_mcica_subcol_gen_run

!==================================================================================================
contains
!==================================================================================================

subroutine rrtmgp_sw_mcica_subcol_gen_run(dosw, kdist_sw, nswbands, nswgpts, nday, nlay, &
              pver, tiny, idxday, ktopcam, ktoprad, cldfprime, c_cld_tau,       &
              c_cld_tau_w, c_cld_tau_w_g, cloud_sw, pmid, errmsg, errflg)
   use ccpp_kinds,              only: kind_phys
   use ccpp_gas_concentrations, only: ty_gas_concs_ccpp
   use ccpp_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp_ccpp
   use ccpp_optical_props,      only: ty_optical_props_2str_ccpp

   ! Compute combined cloud optical properties.
   ! Create MCICA stochastic arrays for cloud SW optical properties.
   ! Initialize optical properties object (cloud_sw) and load with MCICA columns.

   ! arguments
   class(ty_gas_optics_rrtmgp_ccpp), intent(in)  :: kdist_sw             ! shortwave gas optics object
   integer,                          intent(in)  :: nswbands
   integer,                          intent(in)  :: nswgpts
   integer,                          intent(in)  :: nlay                 ! number of layers in radiation calculation (may include "extra layer")
   integer,                          intent(in)  :: nday                 ! number of daylight columns
   integer,                          intent(in)  :: pver
   integer,                          intent(in)  :: ktopcam
   integer,                          intent(in)  :: ktoprad
   integer,                          intent(in)  :: idxday(:)            ! indices of daylight columns in the chunk
   real(kind_phys),                  intent(in)  :: tiny
   real(kind_phys),                  intent(in)  :: c_cld_tau(:,:,:)
   real(kind_phys),                  intent(in)  :: c_cld_tau_w(:,:,:)
   real(kind_phys),                  intent(in)  :: c_cld_tau_w_g(:,:,:)
   real(kind_phys),                  intent(in)  :: cldfprime(:,:)       ! combined cloud fraction
   real(kind_phys),                  intent(in)  :: pmid(:,:)
   logical,                          intent(in)  :: dosw

   type(ty_optical_props_2str_ccpp), intent(inout) :: cloud_sw  ! SW cloud optical properties object
   character(len=512),               intent(out)   :: errmsg
   integer,                          intent(out)   :: errflg

   ! Local variables

   integer :: i, k, ncol
   integer :: igpt, nver
   integer :: istat
   integer, parameter :: changeseed = 1
   character(len=256) :: alloc_errmsg

   ! Arrays for converting from CAM chunks to RRTMGP inputs.
   real(kind_phys), allocatable :: cldf(:,:)
   real(kind_phys), allocatable :: tauc(:,:,:)
   real(kind_phys), allocatable :: ssac(:,:,:)
   real(kind_phys), allocatable :: asmc(:,:,:)
   real(kind_phys), allocatable :: taucmcl(:,:,:)
   real(kind_phys), allocatable :: ssacmcl(:,:,:)
   real(kind_phys), allocatable :: asmcmcl(:,:,:)
   real(kind_phys), allocatable :: day_cld_tau(:,:,:)
   real(kind_phys), allocatable :: day_cld_tau_w(:,:,:)
   real(kind_phys), allocatable :: day_cld_tau_w_g(:,:,:)
   !--------------------------------------------------------------------------------

   ! if no daylight columns the cloud_sw object isn't initialized
   if (nday > 0 .and. dosw) then

      ! number of CAM's layers in radiation calculation.  Does not include the "extra layer".
      nver = pver - ktopcam + 1

      allocate( &
         cldf(nday,nver),                      &
         day_cld_tau(nswbands,nday,nver),      &
         day_cld_tau_w(nswbands,nday,nver),    &
         day_cld_tau_w_g(nswbands,nday,nver),  &
         tauc(nswbands,nday,nver), taucmcl(nswgpts,nday,nver), &
         ssac(nswbands,nday,nver), ssacmcl(nswgpts,nday,nver), &
         asmc(nswbands,nday,nver), asmcmcl(nswgpts,nday,nver), stat=istat, errmsg=alloc_errmsg)
      if (istat /= 0) then
         errflg = 1
         write(errmsg,*) 'rrtmgp_sw_mcica_subcol_gen_run: failed to allocate variable(s) - message: ', alloc_errmsg
      end if

      ! Subset "chunk" data so just the daylight columns, and the number of CAM layers in the
      ! radiation calculation are used by MCICA to produce subcolumns.
      cldf            = cldfprime(       idxday(1:nday), ktopcam:)
      day_cld_tau     = c_cld_tau(    :, idxday(1:nday), ktopcam:)
      day_cld_tau_w   = c_cld_tau_w(  :, idxday(1:nday), ktopcam:)
      day_cld_tau_w_g = c_cld_tau_w_g(:, idxday(1:nday), ktopcam:)

      ! Compute the optical properties needed for the 2-stream calculations.  These calculations
      ! are the same as the RRTMG version.

      ! set cloud optical depth, clip @ zero
      tauc = merge(day_cld_tau, 0.0_kind_phys, day_cld_tau > 0.0_kind_phys)
      ! set value of asymmetry
      asmc = merge(day_cld_tau_w_g / max(day_cld_tau_w, tiny), 0.0_kind_phys, day_cld_tau_w > 0.0_kind_phys)
      ! set value of single scattering albedo
      ssac = merge(max(day_cld_tau_w, tiny) / max(tauc, tiny), 1.0_kind_phys , tauc > 0.0_kind_phys)
      ! set asymmetry to zero when tauc = 0
      asmc = merge(asmc, 0.0_kind_phys, tauc > 0.0_kind_phys)

      ! MCICA uses spectral data (on bands) to construct subcolumns (one per g-point)
      call mcica_subcol_sw( &
         kdist_sw%gas_props, nswbands, nswgpts, nday, nlay, &
         nver, changeseed, pmid, cldf, tauc,     &
         ssac, asmc, taucmcl, ssacmcl, asmcmcl)
   
      ! If there is an extra layer in the radiation then this initialization
      ! will provide the optical properties there.
      cloud_sw%optical_props%tau = 0.0_kind_phys
      cloud_sw%optical_props%ssa = 1.0_kind_phys
      cloud_sw%optical_props%g   = 0.0_kind_phys

      ! Set the properties on g-points.
      do igpt = 1,nswgpts
         cloud_sw%optical_props%g  (:,ktoprad:,igpt) = asmcmcl(igpt,:,:)
         cloud_sw%optical_props%ssa(:,ktoprad:,igpt) = ssacmcl(igpt,:,:)
         cloud_sw%optical_props%tau(:,ktoprad:,igpt) = taucmcl(igpt,:,:)
      end do

      ! validate checks that: tau > 0, ssa is in range [0,1], and g is in range [-1,1].
      errmsg = cloud_sw%optical_props%validate()
      if (len_trim(errmsg) > 0) then
         errflg = 1
         return
      end if

      ! delta scaling adjusts for forward scattering
      errmsg = cloud_sw%optical_props%delta_scale()
      if (len_trim(errmsg) > 0) then
         errflg = 1
         return
      end if

      ! All information is in cloud_sw, now deallocate local vars.
      deallocate( &
         cldf, tauc, ssac, asmc, &
         taucmcl, ssacmcl, asmcmcl,&
         day_cld_tau, day_cld_tau_w, day_cld_tau_w_g )

   end if

end subroutine rrtmgp_sw_mcica_subcol_gen_run

end module rrtmgp_sw_mcica_subcol_gen
