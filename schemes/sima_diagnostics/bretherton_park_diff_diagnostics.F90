! Diagnostics for Bretherton-Park (UW) PBL scheme
module bretherton_park_diff_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: bretherton_park_diff_diagnostics_init
   public :: bretherton_park_diff_diagnostics_run

contains

!> \section arg_table_bretherton_park_diff_diagnostics_init  Argument Table
!! \htmlinclude bretherton_park_diff_diagnostics_init.html
   subroutine bretherton_park_diff_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! History add field calls - surface and horizontal-only fields
      call history_add_field('WGUSTD', 'Wind gusts from turbulence', horiz_only, 'avg', 'm s-1')
      call history_add_field('UW_errorPBL', 'Error function of UW PBL', horiz_only, 'avg', 'm2 s-1')
      call history_add_field('UW_pblh', 'PBL height', horiz_only, 'avg', 'm')
      call history_add_field('UW_pblhp', 'PBL top pressure', horiz_only, 'avg', 'Pa')
      call history_add_field('UW_tpert', 'Convective temperature excess', horiz_only, 'avg', 'K')
      call history_add_field('UW_qpert', 'Convective moisture excess', horiz_only, 'avg', 'kg kg-1')
      call history_add_field('UW_wpert', 'Convective velocity excess', horiz_only, 'avg', 'm s-1')
      call history_add_field('UW_ustar', 'Surface friction velocity', horiz_only, 'avg', 'm s-1')
      call history_add_field('UW_tkes', 'Surface TKE', horiz_only, 'avg', 'm2 s-2')
      call history_add_field('UW_minpblh', 'Minimum PBL height', horiz_only, 'avg', 'm')
      call history_add_field('UW_ncvfin_o', 'Initial total number of CL regimes', horiz_only, 'avg', '1')
      call history_add_field('UW_ncvfin_mg', 'Number of CLs after merging', horiz_only, 'avg', '1')
      call history_add_field('UW_ncvfin_f', 'Final number of CLs with SRCL', horiz_only, 'avg', '1')

      ! Level-based (midpoint) fields
      call history_add_field('UW_n2', 'Buoyancy frequency', 'lev', 'avg', 's-2')
      call history_add_field('UW_s2', 'Shear frequency', 'lev', 'avg', 's-2')
      call history_add_field('UW_ri', 'Interface Richardson number', 'lev', 'avg', '1')
      call history_add_field('UW_sfuh', 'Upper-half saturation fraction', 'lev', 'avg', '1')
      call history_add_field('UW_sflh', 'Lower-half saturation fraction', 'lev', 'avg', '1')
      call history_add_field('UW_cldn', 'Cloud fraction', 'lev', 'avg', '1')
      call history_add_field('UW_qrl', 'LW cooling rate', 'lev', 'avg', 'K s-1')
      call history_add_field('UW_ql', 'Liquid water content', 'lev', 'avg', 'kg kg-1')

      ! Interface-based fields
      call history_add_field('BPROD', 'Buoyancy production', 'ilev', 'avg', 'm2 s-3')
      call history_add_field('SFI', 'Interface-layer saturation fraction', 'ilev', 'avg', '1')
      call history_add_field('SPROD', 'Shear production', 'ilev', 'avg', 'm2 s-3')
      call history_add_field('UW_sfi', 'Interface saturation fraction', 'ilev', 'avg', '1')
      call history_add_field('UW_chu', 'Buoyancy coefficient chu', 'ilev', 'avg', 'm s-2 kg J-1')
      call history_add_field('UW_chs', 'Buoyancy coefficient chs', 'ilev', 'avg', 'm s-2 kg J-1')
      call history_add_field('UW_cmu', 'Buoyancy coefficient cmu', 'ilev', 'avg', 'm s-2 kg-1 kg')
      call history_add_field('UW_cms', 'Buoyancy coefficient cms', 'ilev', 'avg', 'm s-2 kg-1 kg')
      call history_add_field('UW_tke', 'Turbulent kinetic energy', 'ilev', 'avg', 'm2 s-2')
      call history_add_field('UW_wcap', 'Normalized TKE', 'ilev', 'avg', 'm2 s-2')
      call history_add_field('UW_bprod', 'Buoyancy production', 'ilev', 'avg', 'm2 s-3')
      call history_add_field('UW_sprod', 'Shear production', 'ilev', 'avg', 'm2 s-3')
      call history_add_field('UW_kvh', 'Eddy diffusivity of heat', 'ilev', 'avg', 'm2 s-1')
      call history_add_field('UW_kvm', 'Eddy diffusivity of momentum', 'ilev', 'avg', 'm2 s-1')
      call history_add_field('UW_turbtype', 'Interface turbulence type identifier', 'ilev', 'avg', '1')
      call history_add_field('UW_gh', 'Normalized buoyancy production at all interfaces', 'ilev', 'avg', '1')
      call history_add_field('UW_sh', 'Galperin instability function for heat-moisture at all interfaces', 'ilev', 'avg', '1')
      call history_add_field('UW_sm', 'Galperin instability function for momentum at all interfaces', 'ilev', 'avg', '1')
      call history_add_field('UW_ria', 'Richardson number at all interfaces', 'ilev', 'avg', '1')
      call history_add_field('UW_leng', 'Turbulence length scale at all interfaces', 'ilev', 'avg', 'm')

      ! Convective layer (CL) regime fields - these have dimension ncvmax
      ! Note: These need special handling in metadata for the ncvmax dimension;
      ! for history output, they are output throughout the vertical (and padded to zero otherwise)
      !
      ! Note: this logic is 'overbuilt' in this diagnostics scheme, because for now
      ! ncvmax is hardcoded to pver in the UW PBL scheme.
      ! But it is better to overbuild this now rather than risking an out-of-bounds error
      ! down the line if ncvmax is set to less than pver.
      call history_add_field('UW_kbase_o', 'Initial CL base external interface index', 'lev', 'avg', '1')
      call history_add_field('UW_ktop_o', 'Initial CL top external interface index', 'lev', 'avg', '1')
      call history_add_field('UW_kbase_mg', 'CL base after merging', 'lev', 'avg', '1')
      call history_add_field('UW_ktop_mg', 'CL top after merging', 'lev', 'avg', '1')
      call history_add_field('UW_kbase_f', 'Final CL base with SRCL', 'lev', 'avg', '1')
      call history_add_field('UW_ktop_f', 'Final CL top with SRCL', 'lev', 'avg', '1')
      call history_add_field('UW_wet', 'Entrainment rate at CL top', 'lev', 'avg', 'm s-1')
      call history_add_field('UW_web', 'Entrainment rate at CL base', 'lev', 'avg', 'm s-1')
      call history_add_field('UW_jtbu', 'Buoyancy jump across CL top', 'lev', 'avg', 'm s-2')
      call history_add_field('UW_jbbu', 'Buoyancy jump across CL base', 'lev', 'avg', 'm s-2')
      call history_add_field('UW_evhc', 'Evaporative enhancement factor at CL top', 'lev', 'avg', '1')
      call history_add_field('UW_jt2slv', 'Liquid water virtual static energy jump for evhc', 'lev', 'avg', 'J kg-1')
      call history_add_field('UW_n2ht', 'Buoyancy frequency at just below CL top interface', 'lev', 'avg', 's-2')
      call history_add_field('UW_n2hb', 'Buoyancy frequency at just above CL base interface', 'lev', 'avg', 's-2')
      call history_add_field('UW_lwp', 'Liquid water path in CL top layer', 'lev', 'avg', 'kg m-2')
      call history_add_field('UW_optdepth', 'Optical depth of CL top layer', 'lev', 'avg', '1')
      call history_add_field('UW_radfrac', 'Fraction of radiative cooling confined in CL top', 'lev', 'avg', '1')
      call history_add_field('UW_radf', 'Buoyancy production at CL top by radiative cooling', 'lev', 'avg', 'm2 s-3')
      call history_add_field('UW_wstar', 'Convective velocity in each CL', 'lev', 'avg', 'm s-1')
      call history_add_field('UW_wstar3fact', 'Enhancement of wstar3 due to entrainment', 'lev', 'avg', '1')
      call history_add_field('UW_ebrk', 'Net mean TKE of CL including entrainment effect', 'lev', 'avg', 'm2 s-2')
      call history_add_field('UW_wbrk', 'Net mean normalized TKE of CL', 'lev', 'avg', 'm2 s-2')
      call history_add_field('UW_lbrk', 'Energetic internal thickness of CL', 'lev', 'avg', 'm')
      call history_add_field('UW_ricl', 'CL-averaged Richardson number', 'lev', 'avg', '1')
      call history_add_field('UW_ghcl', 'CL-averaged normalized buoyancy production', 'lev', 'avg', '1')
      call history_add_field('UW_shcl', 'CL-averaged Galperin instability function for heat-moisture', 'lev', 'avg', '1')
      call history_add_field('UW_smcl', 'CL-averaged Galperin instability function for momentum', 'lev', 'avg', '1')
      call history_add_field('UW_wsed', 'Sedimentation velocity at top of each CL', 'lev', 'avg', 'm s-1')

   end subroutine bretherton_park_diff_diagnostics_init

!> \section arg_table_bretherton_park_diff_diagnostics_run  Argument Table
!! \htmlinclude bretherton_park_diff_diagnostics_run.html
   subroutine bretherton_park_diff_diagnostics_run( &
      ncol, pver, pverp, ncvmax, &
      ! Surface/horizontal-only fields
      errorPBL, pblh, pblhp, &
      tpert, qpert, wpert, &
      ustar, tkes, minpblh, &
      ncvfin_o, ncvfin_mg, ncvfin_f, &
      ! Level fields
      n2, s2, ri, &
      sfuh, sflh, &
      cldn, qrl, qlfd, &
      ! Interface fields
      bprod, sprod, sfi, &
      chu, chs, cmu, cms, &
      tke, wcap, &
      kvh, kvm, &
      turbtype, &
      ghi, shi, smi, rii, lengi, &
      ! Convective layer fields (ncvmax dimension)
      kbase_o, ktop_o, &
      kbase_mg, ktop_mg, &
      kbase_f, ktop_f, &
      wet, web, &
      jtbu, jbbu, &
      evhc, jt2slv, &
      n2ht, n2hb, &
      lwp, opt_depth, radinvfrac, radf, &
      wstar, wstar3fact, &
      ebrk, wbrk, lbrk, &
      ricl, ghcl, shcl, smcl, &
      wsed, &
      errmsg, errflg)

      use cam_history, only: history_out_field

      integer,            intent(in)  :: ncol
      integer,            intent(in)  :: pver
      integer,            intent(in)  :: pverp
      integer,            intent(in)  :: ncvmax              ! # of convective layers [count]

      real(kind_phys),    intent(in)  :: errorPBL(:)         ! Error function of UW PBL [m2 s-1]
      real(kind_phys),    intent(in)  :: pblh(:)             ! PBL height [m]
      real(kind_phys),    intent(in)  :: pblhp(:)            ! PBL top pressure [Pa]
      real(kind_phys),    intent(in)  :: tpert(:)            ! Convective temperature excess [K]
      real(kind_phys),    intent(in)  :: qpert(:)            ! Convective moisture excess [kg kg-1]
      real(kind_phys),    intent(in)  :: wpert(:)            ! Convective velocity excess [m s-1]
      real(kind_phys),    intent(in)  :: ustar(:)            ! Surface friction velocity [m s-1]
      real(kind_phys),    intent(in)  :: tkes(:)             ! Surface TKE [m2 s-2]
      real(kind_phys),    intent(in)  :: minpblh(:)          ! Minimum PBL height [m]
      real(kind_phys),    intent(in)  :: ncvfin_o(:)         ! Initial total number of CL regimes [1]
      real(kind_phys),    intent(in)  :: ncvfin_mg(:)        ! Number of CLs after merging [1]
      real(kind_phys),    intent(in)  :: ncvfin_f(:)         ! Final number of CLs with SRCL [1]

      real(kind_phys),    intent(in)  :: n2(:,:)             ! Buoyancy frequency [s-2]
      real(kind_phys),    intent(in)  :: s2(:,:)             ! Shear frequency [s-2]
      real(kind_phys),    intent(in)  :: ri(:,:)             ! Richardson number [1]
      real(kind_phys),    intent(in)  :: sfuh(:,:)           ! Upper-half saturation fraction [1]
      real(kind_phys),    intent(in)  :: sflh(:,:)           ! Lower-half saturation fraction [1]
      real(kind_phys),    intent(in)  :: cldn(:,:)           ! Cloud fraction [1]
      real(kind_phys),    intent(in)  :: qrl(:,:)            ! LW cooling rate [K s-1]
      real(kind_phys),    intent(in)  :: qlfd(:,:)           ! Liquid water content [kg kg-1]

      real(kind_phys),    intent(in)  :: bprod(:,:)          ! Buoyancy production [m2 s-3]
      real(kind_phys),    intent(in)  :: sprod(:,:)          ! Shear production [m2 s-3]
      real(kind_phys),    intent(in)  :: sfi(:,:)            ! Interface saturation fraction [1]
      real(kind_phys),    intent(in)  :: chu(:,:)            ! Heat buoyancy coefficient for dry states [kg J-1]
      real(kind_phys),    intent(in)  :: chs(:,:)            ! Heat buoyancy coefficient for saturated states [kg J-1]
      real(kind_phys),    intent(in)  :: cmu(:,:)            ! Moisture buoyancy coefficient for dry states [kg-1 kg]
      real(kind_phys),    intent(in)  :: cms(:,:)            ! Moisture buoyancy coefficient for saturated states [kg-1 kg]
      real(kind_phys),    intent(in)  :: tke(:,:)            ! Turbulent kinetic energy [m2 s-2]
      real(kind_phys),    intent(in)  :: wcap(:,:)           ! Normalized TKE [m2 s-2]
      real(kind_phys),    intent(in)  :: kvh(:,:)            ! Eddy diffusivity of heat [m2 s-1]
      real(kind_phys),    intent(in)  :: kvm(:,:)            ! Eddy diffusivity of momentum [m2 s-1]
      integer,            intent(in)  :: turbtype(:,:)       ! Turbulence type identifier [1]
      real(kind_phys),    intent(in)  :: ghi(:,:)            ! Normalized buoyancy production at all interfaces [1]
      real(kind_phys),    intent(in)  :: shi(:,:)            ! Galperin instability function for heat-moisture [1]
      real(kind_phys),    intent(in)  :: smi(:,:)            ! Galperin instability function for momentum [1]
      real(kind_phys),    intent(in)  :: rii(:,:)            ! Richardson number at all interfaces [1]
      real(kind_phys),    intent(in)  :: lengi(:,:)          ! Turbulence length scale [m]

      ! Convective layer input parameters (ncvmax dimension)
      real(kind_phys),    intent(in)  :: kbase_o(:,:)        ! Initial CL base external interface index [1]
      real(kind_phys),    intent(in)  :: ktop_o(:,:)         ! Initial CL top external interface index [1]
      real(kind_phys),    intent(in)  :: kbase_mg(:,:)       ! CL base after merging [1]
      real(kind_phys),    intent(in)  :: ktop_mg(:,:)        ! CL top after merging [1]
      real(kind_phys),    intent(in)  :: kbase_f(:,:)        ! Final CL base with SRCL [1]
      real(kind_phys),    intent(in)  :: ktop_f(:,:)         ! Final CL top with SRCL [1]
      real(kind_phys),    intent(in)  :: wet(:,:)            ! Entrainment rate at CL top [m s-1]
      real(kind_phys),    intent(in)  :: web(:,:)            ! Entrainment rate at CL base [m s-1]
      real(kind_phys),    intent(in)  :: jtbu(:,:)           ! Buoyancy jump across CL top [m s-2]
      real(kind_phys),    intent(in)  :: jbbu(:,:)           ! Buoyancy jump across CL base [m s-2]
      real(kind_phys),    intent(in)  :: evhc(:,:)           ! Evaporative enhancement factor at CL top [1]
      real(kind_phys),    intent(in)  :: jt2slv(:,:)         ! Liquid water virtual static energy jump [J kg-1]
      real(kind_phys),    intent(in)  :: n2ht(:,:)           ! Buoyancy frequency at just below CL top [s-2]
      real(kind_phys),    intent(in)  :: n2hb(:,:)           ! Buoyancy frequency at just above CL base [s-2]
      real(kind_phys),    intent(in)  :: lwp(:,:)            ! Liquid water path in CL top layer [kg m-2]
      real(kind_phys),    intent(in)  :: opt_depth(:,:)      ! Optical depth of CL top layer [1]
      real(kind_phys),    intent(in)  :: radinvfrac(:,:)     ! Fraction of radiative cooling in CL top [1]
      real(kind_phys),    intent(in)  :: radf(:,:)           ! Buoyancy production at CL top by radiation [m2 s-3]
      real(kind_phys),    intent(in)  :: wstar(:,:)          ! Convective velocity in each CL [m s-1]
      real(kind_phys),    intent(in)  :: wstar3fact(:,:)     ! Enhancement of wstar3 due to entrainment [1]
      real(kind_phys),    intent(in)  :: ebrk(:,:)           ! Net mean TKE of CL [m2 s-2]
      real(kind_phys),    intent(in)  :: wbrk(:,:)           ! Net mean normalized TKE of CL [m2 s-2]
      real(kind_phys),    intent(in)  :: lbrk(:,:)           ! Energetic internal thickness of CL [m]
      real(kind_phys),    intent(in)  :: ricl(:,:)           ! CL-averaged Richardson number [1]
      real(kind_phys),    intent(in)  :: ghcl(:,:)           ! CL-averaged normalized buoyancy production [1]
      real(kind_phys),    intent(in)  :: shcl(:,:)           ! CL-averaged Galperin function for heat-moisture [1]
      real(kind_phys),    intent(in)  :: smcl(:,:)           ! CL-averaged Galperin function for momentum [1]
      real(kind_phys),    intent(in)  :: wsed(:,:)           ! Sedimentation velocity at top of each CL [m s-1]

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variable for turbtype conversion from integer to real.
      real(kind_phys) :: turbtype_real(ncol, pverp)

      ! Temporary array for ncvmax to pver conversion (reusable)
      ! Note that ncvmax is up to pver, so the maximum this temp array
      ! is up to pver (not pverp).
      real(kind_phys) :: tmp_lev(ncol, pver)

      errmsg = ''
      errflg = 0

      ! Convert turbtype from integer to real for output
      turbtype_real = real(turbtype, kind_phys)

      ! Surface/horizontal-only fields
      call history_out_field('WGUSTD', wpert)
      call history_out_field('UW_errorPBL', errorPBL)
      call history_out_field('UW_pblh', pblh)
      call history_out_field('UW_pblhp', pblhp)
      call history_out_field('UW_tpert', tpert)
      call history_out_field('UW_qpert', qpert)
      call history_out_field('UW_wpert', wpert)
      call history_out_field('UW_ustar', ustar)
      call history_out_field('UW_tkes', tkes)
      call history_out_field('UW_minpblh', minpblh)
      call history_out_field('UW_ncvfin_o', ncvfin_o)
      call history_out_field('UW_ncvfin_mg', ncvfin_mg)
      call history_out_field('UW_ncvfin_f', ncvfin_f)

      ! Level fields
      call history_out_field('UW_n2', n2)
      call history_out_field('UW_s2', s2)
      call history_out_field('UW_ri', ri)
      call history_out_field('UW_sfuh', sfuh)
      call history_out_field('UW_sflh', sflh)
      call history_out_field('UW_cldn', cldn)
      call history_out_field('UW_qrl', qrl)
      call history_out_field('UW_ql', qlfd)

      ! Interface fields
      call history_out_field('BPROD', bprod)
      call history_out_field('SFI', sfi)
      call history_out_field('SPROD', sprod)
      call history_out_field('UW_sfi', sfi)
      call history_out_field('UW_chu', chu)
      call history_out_field('UW_chs', chs)
      call history_out_field('UW_cmu', cmu)
      call history_out_field('UW_cms', cms)
      call history_out_field('UW_tke', tke)
      call history_out_field('UW_wcap', wcap)
      call history_out_field('UW_bprod', bprod)
      call history_out_field('UW_sprod', sprod)
      call history_out_field('UW_kvh', kvh)
      call history_out_field('UW_kvm', kvm)
      call history_out_field('UW_turbtype', turbtype_real)
      call history_out_field('UW_gh', ghi)
      call history_out_field('UW_sh', shi)
      call history_out_field('UW_sm', smi)
      call history_out_field('UW_ria', rii)
      call history_out_field('UW_leng', lengi)

      ! Convective layer fields (convert from ncvmax to pver dimension)
      ! Copy data up to ncvmax then pad with zeroes.
      !
      ! The following code is overbuilt for when ncvmax < pver
      if (ncvmax < pver) then
         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = kbase_o(:, :ncvmax)
         call history_out_field('UW_kbase_o', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = ktop_o(:, :ncvmax)
         call history_out_field('UW_ktop_o', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = kbase_mg(:, :ncvmax)
         call history_out_field('UW_kbase_mg', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = ktop_mg(:, :ncvmax)
         call history_out_field('UW_ktop_mg', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = kbase_f(:, :ncvmax)
         call history_out_field('UW_kbase_f', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = ktop_f(:, :ncvmax)
         call history_out_field('UW_ktop_f', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = wet(:, :ncvmax)
         call history_out_field('UW_wet', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = web(:, :ncvmax)
         call history_out_field('UW_web', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = jtbu(:, :ncvmax)
         call history_out_field('UW_jtbu', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = jbbu(:, :ncvmax)
         call history_out_field('UW_jbbu', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = evhc(:, :ncvmax)
         call history_out_field('UW_evhc', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = jt2slv(:, :ncvmax)
         call history_out_field('UW_jt2slv', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = n2ht(:, :ncvmax)
         call history_out_field('UW_n2ht', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = n2hb(:, :ncvmax)
         call history_out_field('UW_n2hb', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = lwp(:, :ncvmax)
         call history_out_field('UW_lwp', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = opt_depth(:, :ncvmax)
         call history_out_field('UW_optdepth', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = radinvfrac(:, :ncvmax)
         call history_out_field('UW_radfrac', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = radf(:, :ncvmax)
         call history_out_field('UW_radf', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = wstar(:, :ncvmax)
         call history_out_field('UW_wstar', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = wstar3fact(:, :ncvmax)
         call history_out_field('UW_wstar3fact', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = ebrk(:, :ncvmax)
         call history_out_field('UW_ebrk', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = wbrk(:, :ncvmax)
         call history_out_field('UW_wbrk', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = lbrk(:, :ncvmax)
         call history_out_field('UW_lbrk', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = ricl(:, :ncvmax)
         call history_out_field('UW_ricl', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = ghcl(:, :ncvmax)
         call history_out_field('UW_ghcl', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = shcl(:, :ncvmax)
         call history_out_field('UW_shcl', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = smcl(:, :ncvmax)
         call history_out_field('UW_smcl', tmp_lev)

         tmp_lev = 0.0_kind_phys
         tmp_lev(:, :ncvmax) = wsed(:, :ncvmax)
         call history_out_field('UW_wsed', tmp_lev)
      else
         ! If ncvmax is pver, we do not need additional copying.
         call history_out_field('UW_kbase_o', kbase_o)
         call history_out_field('UW_ktop_o', ktop_o)
         call history_out_field('UW_kbase_mg', kbase_mg)
         call history_out_field('UW_ktop_mg', ktop_mg)
         call history_out_field('UW_kbase_f', kbase_f)
         call history_out_field('UW_ktop_f', ktop_f)
         call history_out_field('UW_wet', wet)
         call history_out_field('UW_web', web)
         call history_out_field('UW_jtbu', jtbu)
         call history_out_field('UW_jbbu', jbbu)
         call history_out_field('UW_evhc', evhc)
         call history_out_field('UW_jt2slv', jt2slv)
         call history_out_field('UW_n2ht', n2ht)
         call history_out_field('UW_n2hb', n2hb)
         call history_out_field('UW_lwp', lwp)
         call history_out_field('UW_optdepth', opt_depth)
         call history_out_field('UW_radfrac', radinvfrac)
         call history_out_field('UW_radf', radf)
         call history_out_field('UW_wstar', wstar)
         call history_out_field('UW_wstar3fact', wstar3fact)
         call history_out_field('UW_ebrk', ebrk)
         call history_out_field('UW_wbrk', wbrk)
         call history_out_field('UW_lbrk', lbrk)
         call history_out_field('UW_ricl', ricl)
         call history_out_field('UW_ghcl', ghcl)
         call history_out_field('UW_shcl', shcl)
         call history_out_field('UW_smcl', smcl)
         call history_out_field('UW_wsed', wsed)
      end if

   end subroutine bretherton_park_diff_diagnostics_run

end module bretherton_park_diff_diagnostics
