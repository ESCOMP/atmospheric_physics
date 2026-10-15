!> \file rrtmgp_sw_cloud_optics.F90
!!
module rrtmgp_sw_cloud_optics

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

!> \section arg_table_rrtmgp_sw_cloud_optics_run Argument Table
!! \htmlinclude rrtmgp_sw_cloud_optics_run.html
!!
subroutine rrtmgp_sw_cloud_optics_run(dosw, ncol, pver, ktopcam, ktoprad,  nswgpts, nday, idxday, fillvalue, &
   nswbands, iulog, pgam, lamc, rel, rei, nnite, idxnite, cld, cldfsnow, cldfgrau, cldfprime, &
   degrau, dei, des, iclwpth, iciwpth, icswpth, icgrauwpth, tiny_in, idx_sw_diag, do_graupel, &
   do_snow, kdist_sw, cld_tau, grau_tau, snow_tau, c_cld_tau, c_cld_tau_w, c_cld_tau_w_g, tot_cld_vistau,    &
   tot_icld_vistau, liq_icld_vistau, ice_icld_vistau, snow_icld_vistau, grau_icld_vistau, errmsg, errflg)
   use ccpp_gas_optics_rrtmgp,    only: ty_gas_optics_rrtmgp_ccpp
   use ccpp_optical_props,        only: ty_optical_props_2str_ccpp
   use rrtmgp_cloud_optics_setup, only: g_mu, g_lambda, nmu, nlambda, g_d_eff, n_g_d
   use rrtmgp_cloud_optics_setup, only: ext_sw_liq, asm_sw_liq, ssa_sw_liq
   use rrtmgp_cloud_optics_setup, only: ext_sw_ice, asm_sw_ice, ssa_sw_ice
   use rrtmgp_cloud_optics_setup, only: liq_cld_optics, ice_cld_optics
   use ccpp_kinds,                only: kind_phys

   ! Compute combined cloud optical properties.

   ! arguments
   integer,  intent(in) :: ncol                      ! Total number of columns
   integer,  intent(in) :: nday                      ! Number of daylight columns
   integer,  intent(in) :: idxday(:)                 ! Indices of daylight columns
   integer,  intent(in) :: nswgpts                   ! Number of shortwave g-points
   integer,  intent(in) :: pver                      ! Number of vertical layers
   integer,  intent(in) :: ktopcam                   ! Index in host model arrays of top level (layer or interface) at which RRTMGP is active
   integer,  intent(in) :: ktoprad                   ! Index in RRTMGP array corresponding to top layer or interface of host model arrays
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
   real(kind_phys), intent(in) :: rel(:,:)           ! Effective radius of stratiform cloud liquid water droplet [um]
   real(kind_phys), intent(in) :: rei(:,:)           ! Effective radius of stratiform cloud ice crystal [um]
   real(kind_phys), intent(in) :: dei(:,:)           ! Mean effective radius for ice cloud [um]
   real(kind_phys), intent(in) :: des(:,:)           ! Mean effective radius for snow [um]
   real(kind_phys), intent(in) :: degrau(:,:)        ! Mean effective radius for graupel [um]
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
   character(len=*),   intent(out) :: errmsg
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

   select case (trim(liq_cld_optics))
   case ('slingo')
      ! Slingo (1989) liquid optics
      call slingo_liq_optics_sw(ncol, pver, nswbands, cld, rel, iclwpth, liq_tau, liq_tau_w, liq_tau_w_g, sw_tau_w_f, errmsg, errflg)
   case ('gammadist')
      ! gammadist liquid optics
      call get_liquid_optics_sw(ncol, pver, nswbands, tiny_in, ext_sw_liq, asm_sw_liq, ssa_sw_liq, lamc, pgam, g_lambda, g_mu, iclwpth, liq_tau, liq_tau_w, liq_tau_w_g, sw_tau_w_f, errmsg, errflg)
   case default
      write(errmsg,'(a,a)') sub, ': liq_cld_optics must be either slingo or gammadist'
      errflg = 1
   end select
   if (errflg /= 0) then
      return
   end if

   select case (trim(ice_cld_optics))
   case ('ebertcurry')
      ! Ebert and Curry (1992) ice optics
      call ec_ice_optics_sw(ncol, pver, nswbands, cld, rei, iciwpth, ice_tau, ice_tau_w, ice_tau_w_g, sw_tau_w_f, errmsg, errflg)
   case ('mitchell')
      ! Mitchell ice optics
      call interpolate_ice_optics_sw(ncol, pver, nswbands, tiny_in, ext_sw_ice, asm_sw_ice, ssa_sw_ice, iciwpth, dei, g_d_eff, ice_tau, ice_tau_w, ice_tau_w_g, sw_tau_w_f)
   case default
      write(errmsg,'(a,a)') sub, ': ice_cld_optics must be either ebertcurry or mitchell'
      errflg = 1
   end select
   if (errflg /= 0) then
      return
   end if

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
   use ccpp_kinds,                only: kind_phys

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
   use ccpp_kinds,                only: kind_phys

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
   character(len=*),   intent(out) :: errmsg
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
  use ccpp_kinds,       only: kind_phys
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
  use interpolate_data,          only: interp_type, lininterp, lininterp_finish
  use radiation_utils,           only: get_mu_lambda_weights_ccpp
  use rrtmgp_cloud_optics_setup, only: nmu, nlambda
  use ccpp_kinds,                only: kind_phys

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

  character(len=*),   intent(out) :: errmsg
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

!==============================================================================

subroutine slingo_liq_optics_sw(ncol, pver, nswbands, cldn, rel, iclwpth, liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f, errmsg, errflg)
   ! Slingo (1989) (https://doi.org/10.1175/1520-0469(1989)046<1419:AGPFTS>2.0.CO;2) shortwave liquid cloud optics,
   ! using the in-cloud liquid water path.
   use radiation_utils, only: get_sw_spectral_boundaries_ccpp
   use ccpp_kinds,      only: kind_phys

   integer, intent(in)  :: ncol
   integer, intent(in)  :: pver
   integer, intent(in)  :: nswbands
   real(kind_phys), intent(in) :: cldn(:,:)          ! cloud fraction [fraction]
   real(kind_phys), intent(in) :: rel(:,:)           ! liquid effective drop radius [um]
   real(kind_phys), intent(in) :: iclwpth(:,:)       ! in-cloud liquid water path [kg m-2]

   real(kind_phys), intent(out) :: liq_tau    (:,:,:) ! extinction optical depth
   real(kind_phys), intent(out) :: liq_tau_w  (:,:,:) ! single scattering albedo * tau
   real(kind_phys), intent(out) :: liq_tau_w_g(:,:,:) ! asymmetry parameter * tau * w
   real(kind_phys), intent(out) :: liq_tau_w_f(:,:,:) ! forward scattered fraction * tau * w
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   ! Minimum cloud amount (as a fraction of the grid-box area) to
   ! distinguish from clear sky
   real(kind_phys), parameter :: cldmin = 1.0e-80_kind_phys

   ! Decimal precision of cloud amount (0 -> preserve full resolution;
   ! 10^-n -> preserve n digits of cloud amount)
   real(kind_phys), parameter :: cldeps = 0.0_kind_phys

   real(kind_phys), dimension(nswbands) :: wavmin_gp  ! boundaries in RRTMGP band order
   real(kind_phys), dimension(nswbands) :: wavmax_gp  ! boundaries in RRTMGP band order
   real(kind_phys), dimension(nswbands) :: wavmin     ! boundaries in RRTMG band order
   real(kind_phys), dimension(nswbands) :: wavmax     ! boundaries in RRTMG band order

   ! A. Slingo's data for cloud particle radiative properties (from 'A GCM
   ! Parameterization for the Shortwave Properties of Water Clouds' JAS
   ! vol. 46 may 1989 pp 1419-1427)
   real(kind_phys) :: abarl(4) = &  ! A coefficient for extinction optical depth
      [ 2.817e-02_kind_phys, 2.682e-02_kind_phys,2.264e-02_kind_phys,1.281e-02_kind_phys]
   real(kind_phys) :: bbarl(4) = &  ! B coefficient for extinction optical depth
      [ 1.305_kind_phys    , 1.346_kind_phys    ,1.454_kind_phys    ,1.641_kind_phys    ]
   real(kind_phys) :: cbarl(4) = &  ! C coefficient for single scat albedo
      [-5.62e-08_kind_phys ,-6.94e-06_kind_phys ,4.64e-04_kind_phys ,0.201_kind_phys    ]
   real(kind_phys) :: dbarl(4) = &  ! D coefficient for single  scat albedo
      [ 1.63e-07_kind_phys , 2.35e-05_kind_phys ,1.24e-03_kind_phys ,7.56e-03_kind_phys ]
   real(kind_phys) :: ebarl(4) = &  ! E coefficient for asymmetry parameter
      [ 0.829_kind_phys    , 0.794_kind_phys    ,0.754_kind_phys    ,0.826_kind_phys    ]
   real(kind_phys) :: fbarl(4) = &  ! F coefficient for asymmetry parameter
      [ 2.482e-03_kind_phys, 4.226e-03_kind_phys,6.560e-03_kind_phys,4.353e-03_kind_phys]

   real(kind_phys) :: abarli        ! A coefficient for current spectral band
   real(kind_phys) :: bbarli        ! B coefficient for current spectral band
   real(kind_phys) :: cbarli        ! C coefficient for current spectral band
   real(kind_phys) :: dbarli        ! D coefficient for current spectral band
   real(kind_phys) :: ebarli        ! E coefficient for current spectral band
   real(kind_phys) :: fbarli        ! F coefficient for current spectral band

   ! Caution... A. Slingo recommends no less than 4.0 micro-meters nor
   ! greater than 20 micro-meters

   integer :: ns, i, k, indxsl
   real(kind_phys) :: tmp1l, tmp2l, tmp3l, g

   ! Set error variables
   errmsg = ''
   errflg = 0

   ! get_sw_spectral_boundaries_ccpp returns the boundaries in RRTMGP band
   ! order; the optics in this scheme are computed in RRTMG band order and
   ! reordered afterwards (see rrtmg_to_rrtmgp_swbands), so permute the
   ! boundaries into RRTMG band order here.
   call get_sw_spectral_boundaries_ccpp(wavmin_gp, wavmax_gp, 'microns', errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   wavmin(rrtmg_to_rrtmgp_swbands) = wavmin_gp
   wavmax(rrtmg_to_rrtmgp_swbands) = wavmax_gp

   do ns = 1, nswbands
      ! Set index for cloud particle properties based on the wavelength,
      ! according to A. Slingo (1989) equations 1-3:
      ! Use index 1 (0.25 to 0.69 micrometers) for visible
      ! Use index 2 (0.69 - 1.19 micrometers) for near-infrared
      ! Use index 3 (1.19 to 2.38 micrometers) for near-infrared
      ! Use index 4 (2.38 to 4.00 micrometers) for near-infrared
      if(wavmax(ns) <= 0.7_kind_phys) then
         indxsl = 1
      else if(wavmax(ns) <= 1.25_kind_phys) then
         indxsl = 2
      else if(wavmax(ns) <= 2.38_kind_phys) then
         indxsl = 3
      else if(wavmax(ns) > 2.38_kind_phys) then
         indxsl = 4
      end if

      ! Set cloud extinction optical depth, single scatter albedo,
      ! asymmetry parameter, and forward scattered fraction:
      abarli = abarl(indxsl)
      bbarli = bbarl(indxsl)
      cbarli = cbarl(indxsl)
      dbarli = dbarl(indxsl)
      ebarli = ebarl(indxsl)
      fbarli = fbarl(indxsl)

      do k=1,pver
         do i=1,ncol

            ! note that optical properties for liquid valid only
            ! in range of 4.2 > rel > 16 micron (Slingo 89)
            if (cldn(i,k) >= cldmin .and. cldn(i,k) >= cldeps) then
               tmp1l = abarli + bbarli/min(max(4.2_kind_phys,rel(i,k)),16._kind_phys)
               liq_tau(ns,i,k) = 1000._kind_phys*iclwpth(i,k)*tmp1l
            else
               liq_tau(ns,i,k) = 0.0_kind_phys
            end if

            tmp2l = 1._kind_phys - cbarli - dbarli*min(max(4.2_kind_phys,rel(i,k)),16._kind_phys)
            tmp3l = fbarli*min(max(4.2_kind_phys,rel(i,k)),16._kind_phys)
            ! Do not let single scatter albedo be 1.  Delta-eddington solution
            ! for non-conservative case has different analytic form from solution
            ! for conservative case, and raddedmx is written for non-conservative case.
            liq_tau_w(ns,i,k) = liq_tau(ns,i,k) * min(tmp2l,0.999999_kind_phys)
            g = ebarli + tmp3l
            liq_tau_w_g(ns,i,k) = liq_tau_w(ns,i,k) * g
            liq_tau_w_f(ns,i,k) = liq_tau_w(ns,i,k) * g * g

         end do ! End do i=1,ncol
      end do    ! End do k=1,pver
   end do ! nswbands

end subroutine slingo_liq_optics_sw

!==============================================================================

subroutine ec_ice_optics_sw(ncol, pver, nswbands, cldn, rei, iciwpth, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f, errmsg, errflg)
   ! Ebert and Curry (1992) (https://doi.org/10.1029/91JD02472) shortwave ice cloud optics,
   ! using the in-cloud ice water path.
   use radiation_utils, only: get_sw_spectral_boundaries_ccpp
   use ccpp_kinds,      only: kind_phys

   integer, intent(in)  :: ncol
   integer, intent(in)  :: pver
   integer, intent(in)  :: nswbands
   real(kind_phys), intent(in) :: cldn(:,:)          ! cloud fraction [fraction]
   real(kind_phys), intent(in) :: rei(:,:)           ! ice effective drop size [um]
   real(kind_phys), intent(in) :: iciwpth(:,:)       ! in-cloud ice water path [kg m-2]

   real(kind_phys), intent(out) :: ice_tau    (:,:,:) ! extinction optical depth
   real(kind_phys), intent(out) :: ice_tau_w  (:,:,:) ! single scattering albedo * tau
   real(kind_phys), intent(out) :: ice_tau_w_g(:,:,:) ! asymmetry parameter * tau * w
   real(kind_phys), intent(out) :: ice_tau_w_f(:,:,:) ! forward scattered fraction * tau * w
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   real(kind_phys), parameter :: scalefactor = 1._kind_phys

   ! Minimum cloud amount (as a fraction of the grid-box area) to
   ! distinguish from clear sky
   real(kind_phys), parameter :: cldmin = 1.0e-80_kind_phys

   ! Decimal precision of cloud amount (0 -> preserve full resolution;
   ! 10^-n -> preserve n digits of cloud amount)
   real(kind_phys), parameter :: cldeps = 0.0_kind_phys

   real(kind_phys), dimension(nswbands) :: wavmin_gp  ! boundaries in RRTMGP band order
   real(kind_phys), dimension(nswbands) :: wavmax_gp  ! boundaries in RRTMGP band order
   real(kind_phys), dimension(nswbands) :: wavmin     ! boundaries in RRTMG band order
   real(kind_phys), dimension(nswbands) :: wavmax     ! boundaries in RRTMG band order

   ! ice water coefficients (Ebert and Curry,1992, JGR, 97, 3831-3836)
   real(kind_phys) :: abari(4) = &     ! a coefficient for extinction optical depth
      [ 3.448e-03_kind_phys, 3.448e-03_kind_phys,3.448e-03_kind_phys,3.448e-03_kind_phys]
   real(kind_phys) :: bbari(4) = &     ! b coefficient for extinction optical depth
      [ 2.431_kind_phys    , 2.431_kind_phys    ,2.431_kind_phys    ,2.431_kind_phys    ]
   real(kind_phys) :: cbari(4) = &     ! c coefficient for single scat albedo
      [ 1.00e-05_kind_phys , 1.10e-04_kind_phys ,1.861e-02_kind_phys,.46658_kind_phys   ]
   real(kind_phys) :: dbari(4) = &     ! d coefficient for single scat albedo
      [ 0.0_kind_phys      , 1.405e-05_kind_phys,8.328e-04_kind_phys,2.05e-05_kind_phys ]
   real(kind_phys) :: ebari(4) = &     ! e coefficient for asymmetry parameter
      [ 0.7661_kind_phys   , 0.7730_kind_phys   ,0.794_kind_phys    ,0.9595_kind_phys   ]
   real(kind_phys) :: fbari(4) = &     ! f coefficient for asymmetry parameter
      [ 5.851e-04_kind_phys, 5.665e-04_kind_phys,7.267e-04_kind_phys,1.076e-04_kind_phys]

   real(kind_phys) :: abarii           ! A coefficient for current spectral band
   real(kind_phys) :: bbarii           ! B coefficient for current spectral band
   real(kind_phys) :: cbarii           ! C coefficient for current spectral band
   real(kind_phys) :: dbarii           ! D coefficient for current spectral band
   real(kind_phys) :: ebarii           ! E coefficient for current spectral band
   real(kind_phys) :: fbarii           ! F coefficient for current spectral band

   integer :: ns, i, k, indxsl
   real(kind_phys) :: tmp1i, tmp2i, tmp3i, g

   ! Set error variables
   errmsg = ''
   errflg = 0

   ! get_sw_spectral_boundaries_ccpp returns the boundaries in RRTMGP band
   ! order; the optics in this scheme are computed in RRTMG band order and
   ! reordered afterwards (see rrtmg_to_rrtmgp_swbands), so permute the
   ! boundaries into RRTMG band order here.
   call get_sw_spectral_boundaries_ccpp(wavmin_gp, wavmax_gp, 'microns', errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   wavmin(rrtmg_to_rrtmgp_swbands) = wavmin_gp
   wavmax(rrtmg_to_rrtmgp_swbands) = wavmax_gp

   do ns = 1, nswbands

      if(wavmax(ns) <= 0.7_kind_phys) then
         indxsl = 1
      else if(wavmax(ns) <= 1.25_kind_phys) then
         indxsl = 2
      else if(wavmax(ns) <= 2.38_kind_phys) then
         indxsl = 3
      else if(wavmax(ns) > 2.38_kind_phys) then
         indxsl = 4
      end if

      abarii = abari(indxsl)
      bbarii = bbari(indxsl)
      cbarii = cbari(indxsl)
      dbarii = dbari(indxsl)
      ebarii = ebari(indxsl)
      fbarii = fbari(indxsl)

      do k=1,pver
         do i=1,ncol

            ! note that optical properties for ice valid only
            ! in range of 13 > rei > 130 micron (Ebert and Curry 92)
            if (cldn(i,k) >= cldmin .and. cldn(i,k) >= cldeps) then
               tmp1i = abarii + bbarii/max(13._kind_phys,min(scalefactor*rei(i,k),130._kind_phys))
               ice_tau(ns,i,k) = 1000.0_kind_phys*iciwpth(i,k)*tmp1i
            else
               ice_tau(ns,i,k) = 0.0_kind_phys
            end if

            tmp2i = 1._kind_phys - cbarii - dbarii*min(max(13._kind_phys,scalefactor*rei(i,k)),130._kind_phys)
            tmp3i = fbarii*min(max(13._kind_phys,scalefactor*rei(i,k)),130._kind_phys)
            ! Do not let single scatter albedo be 1.  Delta-eddington solution
            ! for non-conservative case has different analytic form from solution
            ! for conservative case, and raddedmx is written for non-conservative case.
            ice_tau_w(ns,i,k) = ice_tau(ns,i,k) * min(tmp2i,0.999999_kind_phys)
            g = ebarii + tmp3i
            ice_tau_w_g(ns,i,k) = ice_tau_w(ns,i,k) * g
            ice_tau_w_f(ns,i,k) = ice_tau_w(ns,i,k) * g * g

         end do ! End do i=1,ncol
      end do    ! End do k=1,pver
   end do ! nswbands

end subroutine ec_ice_optics_sw

end module rrtmgp_sw_cloud_optics
