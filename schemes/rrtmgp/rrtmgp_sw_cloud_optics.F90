module rrtmgp_sw_cloud_optics
use ccpp_kinds, only: kind_phys

!--------------------------------------------------------------------------------
! Transform data for inputs from CAM's data structures to those used by
! RRTMGP.  Subset the number of model levels if CAM's top exceeds RRTMGP's
! valid domain.  Add an extra layer if CAM's top is below 1 Pa.
! The vertical indexing increases from top to bottom of atmosphere in both
! CAM and RRTMGP arrays.   
!--------------------------------------------------------------------------------

implicit none
private
save

public :: rrtmgp_sw_cloud_optics_run

! Mapping from RRTMG shortwave bands to RRTMGP.  Currently needed to continue using
! the SW optics datasets from RRTMG (even thought there is a slight mismatch in the
! band boundaries of the 2 bands that overlap with the LW bands).
integer, parameter, dimension(14) :: rrtmg_to_rrtmgp_swbands = &
   [ 14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 ]

!==================================================================================================
contains
!==================================================================================================

subroutine rrtmgp_sw_cloud_optics_run(dosw, ncol, pver, ktopcam, ktoprad, nlay, nswgpts, nday, idxday, fillvalue, &
   nswbands, iulog, pgam, lamc, nnite, idxnite, cld, cldfsnow, cldfgrau, cldfprime, cld_tau, grau_tau, &
   snow_tau, degrau, dei, des, iclwpth, iciwpth, icswpth, icgrauwpth, tiny_in, idx_sw_diag, do_graupel, &
   do_snow, kdist_sw, c_cld_tau, c_cld_tau_w, c_cld_tau_w_g, tot_cld_vistau, tot_icld_vistau,          &
   liq_icld_vistau, ice_icld_vistau, snow_icld_vistau, grau_icld_vistau, errmsg, errflg)
   use ccpp_gas_optics_rrtmgp,    only: ty_gas_optics_rrtmgp_ccpp
   use ccpp_optical_props,        only: ty_optical_props_2str_ccpp
   use rrtmgp_cloud_optics_setup, only: g_mu, g_lambda, nmu, nlambda, g_d_eff, n_g_d
   use rrtmgp_cloud_optics_setup, only: ext_sw_liq, asm_sw_liq, ssa_sw_liq
   use rrtmgp_cloud_optics_setup, only: ext_sw_ice, asm_sw_ice, ssa_sw_ice

   ! Compute combined cloud optical properties.

   ! arguments
   integer,  intent(in) :: nlay                      ! Number of layers in radiation calculation (may include "extra layer")
   integer,  intent(in) :: ncol                      ! Total number of columns
   integer,  intent(in) :: nday                      ! Number of daylight columns
   integer,  intent(in) :: idxday(:)                 ! Indices of daylight columns
   integer,  intent(in) :: nswgpts                   ! Number of shortwave g-points
   integer,  intent(in) :: pver                      ! Number of vertical layers
   integer,  intent(in) :: ktopcam                   ! Index in CAM arrays of top level (layer or interface) at which RRTMGP is active
   integer,  intent(in) :: ktoprad                   ! Index in RRTMGP array corresponding to top layer or interface of CAM arrays
   integer,  intent(in) :: nswbands                  ! Number of shortwve bands
   integer,  intent(in) :: nnite                     ! Number of night columns
   integer,  intent(in) :: idxnite(:)                ! Indices of night columns in the chunk
   integer,  intent(in) :: iulog                     ! Logging unit
   integer,  intent(in) :: idx_sw_diag               ! Index for band that contains 500-nm wave

   logical,  intent(in) :: do_snow                   ! Flag to include snow in radiation calculation
   logical,  intent(in) :: do_graupel                ! Flag to include graupel in radiation calculation
   logical,  intent(in) :: dosw                      ! Flag to do shortwave radiation this timestep

   real(kind_phys), intent(in) :: fillvalue          ! Fill value for night columns
   real(kind_phys), intent(in) :: tiny_in            ! Definition of tiny for RRTMGP

   real(kind_phys), intent(in) :: lamc(:,:)          ! Prognosed value of lambda for cloud [1]
   real(kind_phys), intent(in) :: pgam(:,:)          ! Prognosed value of mu for cloud [1]
   real(kind_phys), intent(in) :: dei(:,:)           ! Mean effective radius for ice cloud [micron]
   real(kind_phys), intent(in) :: des(:,:)           ! Mean effective radius for snow [micron]
   real(kind_phys), intent(in) :: degrau(:,:)        ! Mean effective radius for graupel [micron]
   real(kind_phys), intent(in) :: iclwpth(:,:)       ! In-cloud liquid water path [kg m-2]
   real(kind_phys), intent(in) :: iciwpth(:,:)       ! In-cloud ice water path [kg m-2]
   real(kind_phys), intent(in) :: icswpth(:,:)       ! In-cloud snow water path [kg m-2]
   real(kind_phys), intent(in) :: icgrauwpth(:,:)    ! In-cloud graupel water path [kg m-2]
   real(kind_phys), intent(in) :: cld(:,:)           ! Cloud fraction (liq+ice) [fraction]
   real(kind_phys), intent(in) :: cldfsnow(:,:)      ! Cloud fraction of just "snow clouds" [fraction]
   real(kind_phys), intent(in) :: cldfgrau(:,:)      ! Cloud fraction of just "graupel clouds" [fraction]
   real(kind_phys), intent(in) :: cldfprime(:,:)     ! Combined cloud fraction [fraction]

   class(ty_gas_optics_rrtmgp_ccpp), intent(in)  :: kdist_sw             ! shortwave gas optics object
   real(kind_phys),                  intent(out) :: cld_tau(:,:,:)       ! liquid + ice optical depth
   real(kind_phys),                  intent(out) :: snow_tau(:,:,:)      ! snow optical depth
   real(kind_phys),                  intent(out) :: grau_tau(:,:,:)      ! graupel optical depth
   real(kind_phys),                  intent(out) :: c_cld_tau(:,:,:)     ! combined cloud extinction optical depth
   real(kind_phys),                  intent(out) :: c_cld_tau_w  (:,:,:) ! combined cloud single scattering albedo * tau
   real(kind_phys),                  intent(out) :: c_cld_tau_w_g(:,:,:) ! combined cloud asymmetry parameter * w * tau

   ! Diagnostic outputs
   real(kind_phys), intent(out) :: tot_cld_vistau(:,:)   ! gbx total cloud optical depth
   real(kind_phys), intent(out) :: tot_icld_vistau(:,:)  ! in-cld total cloud optical depth
   real(kind_phys), intent(out) :: liq_icld_vistau(:,:)  ! in-cld liq cloud optical depth
   real(kind_phys), intent(out) :: ice_icld_vistau(:,:)  ! in-cld ice cloud optical depth
   real(kind_phys), intent(out) :: snow_icld_vistau(:,:) ! snow in-cloud visible sw optical depth
   real(kind_phys), intent(out) :: grau_icld_vistau(:,:) ! Graupel in-cloud visible sw optical depth

   ! Error variables
   character(len=512), intent(out) :: errmsg
   integer,            intent(out) :: errflg

   ! Local variables

   integer :: i, k
   integer :: igpt, nver
   integer :: istat
   integer, parameter :: changeseed = 1

   ! cloud radiative parameters are "in cloud" not "in cell"
   real(kind_phys) :: liq_tau    (nswbands,ncol,pver)  ! liquid extinction optical depth
   real(kind_phys) :: liq_tau_w  (nswbands,ncol,pver)  ! liquid single scattering albedo * tau
   real(kind_phys) :: liq_tau_w_g(nswbands,ncol,pver)  ! liquid asymmetry parameter * tau * w
   real(kind_phys) :: ice_tau    (nswbands,ncol,pver)  ! ice extinction optical depth
   real(kind_phys) :: ice_tau_w  (nswbands,ncol,pver)  ! ice single scattering albedo * tau
   real(kind_phys) :: ice_tau_w_g(nswbands,ncol,pver)  ! ice asymmetry parameter * tau * w
   real(kind_phys) :: snow_tau_w (nswbands,ncol,pver)  ! snow single scattering albedo * tau
   real(kind_phys) :: snow_tau_w_g(nswbands,ncol,pver) ! snow asymmetry parameter * tau * w
   real(kind_phys) :: cld_tau_w  (nswbands,ncol,pver)  ! cloud single scattering albedo * tau
   real(kind_phys) :: cld_tau_w_g(nswbands,ncol,pver)  ! cloud asymmetry parameter * w * tau
   real(kind_phys) :: grau_tau_w  (nswbands,ncol,pver) ! graupel single scattering albedo * tau
   real(kind_phys) :: grau_tau_w_g(nswbands,ncol,pver) ! graupel asymmetry parameter * tau * w

   ! RRTMGP does not use this property in its 2-stream calculations.
   real(kind_phys) :: sw_tau_w_f(nswbands,ncol,pver) ! Forward scattered fraction * tau * w.

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

   character(len=*), parameter :: sub = 'rrtmgp_set_cloud_sw'
   !--------------------------------------------------------------------------------

   if (.not. dosw) then
      return
   end if

   ! Combine the cloud optical properties.

   ! gammadist liquid optics
   call get_liquid_optics_sw(ncol, pver, nswbands, ext_sw_liq, asm_sw_liq, ssa_sw_liq, lamc, pgam, g_lambda, g_mu, iclwpth, liq_tau, liq_tau_w, liq_tau_w_g, sw_tau_w_f, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   ! Mitchell ice optics
   call interpolate_ice_optics_sw(ncol, pver, nswbands, tiny_in, ext_sw_ice, asm_sw_ice, ssa_sw_ice, iciwpth, dei, g_d_eff, ice_tau, ice_tau_w, ice_tau_w_g, sw_tau_w_f)

   cld_tau(:,:ncol,:)     =  liq_tau(:,:ncol,:)     + ice_tau(:,:ncol,:)
   cld_tau_w(:,:ncol,:)   =  liq_tau_w(:,:ncol,:)   + ice_tau_w(:,:ncol,:)
   cld_tau_w_g(:,:ncol,:) =  liq_tau_w_g(:,:ncol,:) + ice_tau_w_g(:,:ncol,:)

   ! add in snow
   if (do_snow) then
      call interpolate_ice_optics_sw(ncol, pver, nswbands, tiny_in, ext_sw_ice, asm_sw_ice, ssa_sw_ice, icswpth, des, g_d_eff, snow_tau, snow_tau_w, snow_tau_w_g, sw_tau_w_f)
      do i = 1, ncol
         do k = 1, pver
            if (cldfprime(i,k) > 0._kind_phys) then
               c_cld_tau(:,i,k)     = ( cldfsnow(i,k)*snow_tau(:,i,k) &
                                      + cld(i,k)*cld_tau(:,i,k) )/cldfprime(i,k)
               c_cld_tau_w(:,i,k)   = ( cldfsnow(i,k)*snow_tau_w(:,i,k)  &
                                      + cld(i,k)*cld_tau_w(:,i,k) )/cldfprime(i,k)
               c_cld_tau_w_g(:,i,k) = ( cldfsnow(i,k)*snow_tau_w_g(:,i,k) &
                                      + cld(i,k)*cld_tau_w_g(:,i,k) )/cldfprime(i,k)
            else
               c_cld_tau(:,i,k)     = 0._kind_phys
               c_cld_tau_w(:,i,k)   = 0._kind_phys
               c_cld_tau_w_g(:,i,k) = 0._kind_phys
            end if
         end do
      end do
   else
      c_cld_tau(:,:ncol,:)     = cld_tau(:,:ncol,:)
      c_cld_tau_w(:,:ncol,:)   = cld_tau_w(:,:ncol,:)
      c_cld_tau_w_g(:,:ncol,:) = cld_tau_w_g(:,:ncol,:)
   end if

   ! add in graupel
   if (do_graupel) then
      call get_grau_optics_sw(ncol, pver, nswbands, tiny_in, g_d_eff, ext_sw_ice, asm_sw_ice, ssa_sw_ice, iulog, icgrauwpth, degrau, idx_sw_diag, grau_tau, grau_tau_w, grau_tau_w_g, sw_tau_w_f)
      do i = 1, ncol
         do k = 1, pver
            if (cldfprime(i,k) > 0._kind_phys) then
               c_cld_tau(:,i,k)     = ( cldfgrau(i,k)*grau_tau(:,i,k) &
                                      + cld(i,k)*c_cld_tau(:,i,k) )/cldfprime(i,k)
               c_cld_tau_w(:,i,k)   = ( cldfgrau(i,k)*grau_tau_w(:,i,k)  &
                                      + cld(i,k)*c_cld_tau_w(:,i,k) )/cldfprime(i,k)
               c_cld_tau_w_g(:,i,k) = ( cldfgrau(i,k)*grau_tau_w_g(:,i,k) &
                                      + cld(i,k)*c_cld_tau_w_g(:,i,k) )/cldfprime(i,k)
            else
               c_cld_tau(:,i,k)     = 0._kind_phys
               c_cld_tau_w(:,i,k)   = 0._kind_phys
               c_cld_tau_w_g(:,i,k) = 0._kind_phys
            end if
         end do
      end do
   end if

   ! cloud optical properties need to be re-ordered from the RRTMG spectral bands
   ! (assumed in the optics datasets) to RRTMGP's
   ice_tau(:,:ncol,:)       = ice_tau(rrtmg_to_rrtmgp_swbands,:ncol,:)
   liq_tau(:,:ncol,:)       = liq_tau(rrtmg_to_rrtmgp_swbands,:ncol,:)
   c_cld_tau(:,:ncol,:)     = c_cld_tau(rrtmg_to_rrtmgp_swbands,:ncol,:)
   c_cld_tau_w(:,:ncol,:)   = c_cld_tau_w(rrtmg_to_rrtmgp_swbands,:ncol,:)
   c_cld_tau_w_g(:,:ncol,:) = c_cld_tau_w_g(rrtmg_to_rrtmgp_swbands,:ncol,:)
   if (do_snow) then
      snow_tau(:,:ncol,:)   = snow_tau(rrtmg_to_rrtmgp_swbands,:ncol,:)
   else
      snow_tau(:,:ncol,:)   = 0._kind_phys
   end if
   if (do_graupel) then
      grau_tau(:,:ncol,:)   = grau_tau(rrtmg_to_rrtmgp_swbands,:ncol,:)
   else
      grau_tau(:,:ncol,:)   = 0._kind_phys
   end if

   ! Set arrays for diagnostic output.
   ! cloud optical depth fields for the visible band
   tot_icld_vistau(:ncol,:) = c_cld_tau(idx_sw_diag,:ncol,:)
   liq_icld_vistau(:ncol,:) = liq_tau(idx_sw_diag,:ncol,:)
   ice_icld_vistau(:ncol,:) = ice_tau(idx_sw_diag,:ncol,:)
   if (do_snow) then
      snow_icld_vistau(:ncol,:) = snow_tau(idx_sw_diag,:ncol,:)
   else
      snow_icld_vistau(:ncol,:) = 0._kind_phys
   endif
   if (do_graupel) then
      grau_icld_vistau(:ncol,:) = grau_tau(idx_sw_diag,:ncol,:)
   else
      grau_icld_vistau(:ncol,:) = 0._kind_phys
   endif

   ! multiply by total cloud fraction to get gridbox value
   tot_cld_vistau(:ncol,:) = c_cld_tau(idx_sw_diag,:ncol,:)*cldfprime(:ncol,:)

   ! overwrite night columns with fillvalue
   do i = 1, nnite
      tot_cld_vistau(idxnite(i),:)   = fillvalue
      tot_icld_vistau(idxnite(i),:)  = fillvalue
      liq_icld_vistau(idxnite(i),:)  = fillvalue
      ice_icld_vistau(idxnite(i),:)  = fillvalue
      snow_icld_vistau(idxnite(i),:) = fillvalue
      grau_icld_vistau(idxnite(i),:) = fillvalue
   end do

end subroutine rrtmgp_sw_cloud_optics_run

!==============================================================================

subroutine get_grau_optics_sw(ncol, pver, nswbands, tiny_in, g_d_eff, ext_sw_ice, asm_sw_ice, ssa_sw_ice, &
                iulog, icgrauwpth, degrau, idx_sw_diag, tau, tau_w, tau_w_g, tau_w_f)

   integer, intent(in)  :: ncol
   integer, intent(in)  :: pver
   integer, intent(in)  :: nswbands
   integer, intent(in)  :: iulog
   integer, intent(in)  :: idx_sw_diag
   real(kind_phys), intent(in) :: tiny_in
   real(kind_phys), intent(in) :: ext_sw_ice(:,:)
   real(kind_phys), intent(in) :: asm_sw_ice(:,:)
   real(kind_phys), intent(in) :: ssa_sw_ice(:,:)
   real(kind_phys), intent(in) :: degrau(:,:)
   real(kind_phys), intent(in) :: g_d_eff(:)
   real(kind_phys), intent(in) :: icgrauwpth(:,:)
   
   real(kind_phys),intent(out) :: tau    (:,:,:) ! extinction optical depth
   real(kind_phys),intent(out) :: tau_w  (:,:,:) ! single scattering albedo * tau
   real(kind_phys),intent(out) :: tau_w_g(:,:,:) ! asymmetry parameter * tau * w
   real(kind_phys),intent(out) :: tau_w_f(:,:,:) ! forward scattered fraction * tau * w

   integer :: i,k

   ! This does the same thing as get_ice_optics_sw, except with a different
   ! water path and effective diameter.
   call interpolate_ice_optics_sw(ncol, pver, nswbands, tiny_in, ext_sw_ice, asm_sw_ice, ssa_sw_ice, icgrauwpth, degrau, g_d_eff, tau, tau_w, &
        tau_w_g, tau_w_f)
   do i = 1, ncol
      do k = 1, pver
         if (tau(idx_sw_diag,i,k).gt.100._kind_phys) then
            write(iulog,*) 'WARNING: SW Graupel Tau > 100  (i,k,icgrauwpth,degrau,tau):'
            write(iulog,*) i,k,icgrauwpth(i,k), degrau(i,k), tau(idx_sw_diag,i,k)
         end if
      enddo
   enddo

end subroutine get_grau_optics_sw

!==============================================================================

subroutine get_liquid_optics_sw(ncol, pver, nswbands, tiny_in, ext_sw_liq, asm_sw_liq, ssa_sw_liq, lamc, pgam, g_lambda, &
    g_mu, iclwpth, tau, tau_w, tau_w_g, tau_w_f, errmsg, errflg)
   integer, intent(in)  :: ncol
   integer, intent(in)  :: pver
   integer, intent(in)  :: nswbands
   real(kind_phys), intent(in)  :: tiny_in
   real(kind_phys), intent(in)  :: g_lambda(:,:)
   real(kind_phys), intent(in)  :: g_mu(:)
   real(kind_phys), intent(in)  :: ext_sw_liq(:,:,:)
   real(kind_phys), intent(in)  :: asm_sw_liq(:,:,:)
   real(kind_phys), intent(in)  :: ssa_sw_liq(:,:,:)
   real(kind_phys), intent(in)  :: iclwpth(:,:)
   real(kind_phys), intent(in)  :: lamc(:,:)
   real(kind_phys), intent(in)  :: pgam(:,:)

   real(kind_phys), intent(out) :: tau    (:,:,:) ! extinction optical depth
   real(kind_phys), intent(out) :: tau_w  (:,:,:) ! single scattering albedo * tau
   real(kind_phys), intent(out) :: tau_w_g(:,:,:) ! asymmetry parameter * tau * w
   real(kind_phys), intent(out) :: tau_w_f(:,:,:) ! forward scattered fraction * tau * w
   character(len=512), intent(out) :: errmsg
   integer,            intent(out) :: errflg

   real(kind_phys), dimension(ncol,pver) :: kext
   integer i,k,swband

   do k = 1,pver
      do i = 1,ncol
         if(lamc(i,k) > 0._kind_phys) then ! This seems to be clue from microphysics of no cloud
            call gam_liquid_sw(nswbands, tiny_in, g_lambda, g_mu, ext_sw_liq, asm_sw_liq, ssa_sw_liq, iclwpth(i,k), &
                 lamc(i,k), pgam(i,k), tau(1:nswbands,i,k), tau_w(1:nswbands,i,k), tau_w_g(1:nswbands,i,k),         &
                 tau_w_f(1:nswbands,i,k), errmsg, errflg)
         else
            tau(1:nswbands,i,k) = 0._kind_phys
            tau_w(1:nswbands,i,k) = 0._kind_phys
            tau_w_g(1:nswbands,i,k) = 0._kind_phys
            tau_w_f(1:nswbands,i,k) = 0._kind_phys
         endif
      enddo
   enddo

end subroutine get_liquid_optics_sw

!==============================================================================

subroutine interpolate_ice_optics_sw(ncol, pver, nswbands, tiny_in, ext_sw_ice, asm_sw_ice, ssa_sw_ice, &
     iciwpth, dei, g_d_eff, tau, tau_w, tau_w_g, tau_w_f)
  ! SIMA-specific interpolation routines
  use interpolate_data, only: interp_type, lininterp, lininterp_init, lininterp_finish, extrap_method_bndry

  integer, intent(in) :: ncol
  integer, intent(in) :: pver
  integer, intent(in) :: nswbands
  real(kind_phys), intent(in) :: tiny_in
  real(kind_phys), intent(in) :: iciwpth(:,:)
  real(kind_phys), intent(in) :: dei(:,:)
  real(kind_phys), intent(in) :: g_d_eff(:)
  real(kind_phys), intent(in) :: ext_sw_ice(:,:)
  real(kind_phys), intent(in) :: asm_sw_ice(:,:)
  real(kind_phys), intent(in) :: ssa_sw_ice(:,:)

  real(kind_phys),intent(out) :: tau    (:,:,:) ! extinction optical depth
  real(kind_phys),intent(out) :: tau_w  (:,:,:) ! single scattering albedo * tau
  real(kind_phys),intent(out) :: tau_w_g(:,:,:) ! asymmetry parameter * tau * w
  real(kind_phys),intent(out) :: tau_w_f(:,:,:) ! forward scattered fraction * tau * w

  type(interp_type) :: dei_wgts

  integer :: i, k, swband
  integer :: n_g_d
  real(kind_phys) :: ext(nswbands), ssa(nswbands), asm(nswbands)

  n_g_d = size(g_d_eff)

  do k = 1,pver
     do i = 1,ncol
        if( iciwpth(i,k) < tiny_in .or. dei(i,k) == 0._kind_phys) then
           ! if ice water path is too small, OD := 0
           tau    (:,i,k) = 0._kind_phys
           tau_w  (:,i,k) = 0._kind_phys
           tau_w_g(:,i,k) = 0._kind_phys
           tau_w_f(:,i,k) = 0._kind_phys
        else
           ! for each cell interpolate to find weights in g_d_eff grid.
           call lininterp_init(g_d_eff, n_g_d, dei(i:i,k), 1, &
                extrap_method_bndry, dei_wgts)
           ! interpolate into grid and extract radiative properties
           do swband = 1, nswbands
              call lininterp(ext_sw_ice(:,swband), n_g_d, &
                   ext(swband:swband), 1, dei_wgts)
              call lininterp(ssa_sw_ice(:,swband), n_g_d, &
                   ssa(swband:swband), 1, dei_wgts)
              call lininterp(asm_sw_ice(:,swband), n_g_d, &
                   asm(swband:swband), 1, dei_wgts)
           end do
           tau    (:,i,k) = iciwpth(i,k) * ext
           tau_w  (:,i,k) = tau(:,i,k) * ssa
           tau_w_g(:,i,k) = tau_w(:,i,k) * asm
           tau_w_f(:,i,k) = tau_w_g(:,i,k) * asm
           call lininterp_finish(dei_wgts)
        endif
     enddo
  enddo

end subroutine interpolate_ice_optics_sw

!==============================================================================

subroutine gam_liquid_sw(nswbands, tiny_in, g_lambda, g_mu, ext_sw_liq, asm_sw_liq, ssa_sw_liq, clwptn, lamc, pgam, tau, tau_w, tau_w_g, tau_w_f, errmsg, errflg)
  ! SIMA-specific interpolation routines
  use interpolate_data,         only: interp_type, lininterp, lininterp_finish
  use radiation_utils,          only: get_mu_lambda_weights_ccpp
  use rrtmgp_cloud_optics_setup, only: nmu, nlambda

  integer,         intent(in)  :: nswbands
  real(kind_phys), intent(in)  :: tiny_in
  real(kind_phys), intent(in)  :: ext_sw_liq(:,:,:)
  real(kind_phys), intent(in)  :: asm_sw_liq(:,:,:)
  real(kind_phys), intent(in)  :: ssa_sw_liq(:,:,:)
  real(kind_phys), intent(in)  :: g_mu(:)
  real(kind_phys), intent(in)  :: g_lambda(:,:)
  real(kind_phys), intent(in)  :: lamc
  real(kind_phys), intent(in)  :: pgam
  real(kind_phys), intent(in)  :: clwptn ! cloud water liquid path new (in cloud) [kg m-2]
  real(kind_phys), intent(out) :: tau(:), tau_w(:), tau_w_f(:), tau_w_g(:)

  character(len=512), intent(out) :: errmsg
  integer,            intent(out) :: errflg

  integer :: swband ! sw band index

  real(kind_phys) :: ext(nswbands), ssa(nswbands), asm(nswbands)

  type(interp_type) :: mu_wgts
  type(interp_type) :: lambda_wgts

  ! Set error variables
  errmsg = ''
  errflg = 0

  if (clwptn < tiny_in) then
    tau = 0._kind_phys
    tau_w = 0._kind_phys
    tau_w_g = 0._kind_phys
    tau_w_f = 0._kind_phys
    return
  endif

  call get_mu_lambda_weights_ccpp(nmu, nlambda, g_mu, g_lambda, lamc, pgam, &
                  mu_wgts, lambda_wgts, errmsg, errflg)
  if (errflg /= 0) then
     return
  end if

  do swband = 1, nswbands
     call lininterp(ext_sw_liq(:,:,swband), nmu, nlambda, &
          ext(swband:swband), 1, mu_wgts, lambda_wgts)
     call lininterp(ssa_sw_liq(:,:,swband), nmu, nlambda, &
          ssa(swband:swband), 1, mu_wgts, lambda_wgts)
     call lininterp(asm_sw_liq(:,:,swband), nmu, nlambda, &
          asm(swband:swband), 1, mu_wgts, lambda_wgts)
  enddo

  ! compute radiative properties
  tau = clwptn * ext
  tau_w = tau * ssa
  tau_w_g = tau_w * asm
  tau_w_f = tau_w_g * asm

  call lininterp_finish(mu_wgts)
  call lininterp_finish(lambda_wgts)

end subroutine gam_liquid_sw

!==============================================================================

end module rrtmgp_sw_cloud_optics
