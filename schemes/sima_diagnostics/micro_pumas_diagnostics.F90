! Diagnostics for the PUMAS cloud microphysics scheme
!
! This is the CAM-SIMA conversion of the history output that was done
! in CAM's micro_pumas_cam.F90 (micro_pumas_cam_init / micro_pumas_cam_tend).
!
! There are two diagnostic schemes in this module:
!   micro_pumas_pre_diagnostics    - must be placed BEFORE PUMAS (micro_pumas_ccpp) in the SDF;
!                                    outputs the in-cloud water paths before microphysics
!   micro_pumas_diagnostics        - must be placed AFTER micro_pumas_ccpp_optics_limiter and
!                                    BEFORE apply_heating_rate in the SDF (the air temperature is
!                                    expected to be the value before the microphysics heating is
!                                    applied, as in pumas_post_main)
!
! NOTE: The macrophysics liquid condensation rate (cmeliq) and the shallow convective
!       precipitation rate (prec_sh) are not yet calculated in CAM-SIMA for CAM7 (they come from
!       CLUBB), so they are currently read from the snapshot (registry variables).
!
! Differences from CAM:
!   - Subcolumns are not supported, so all fields are output on the grid.
!   - The 'trop_cld_lev' history vertical coordinate (the model levels from the top level for
!     cloud physics to the bottom, i.e. the PUMAS levels) is registered in
!     micro_pumas_diagnostics_init, as CAM did in micro_pumas_cam_register.
!   - The constituent mixing ratios are output by sima_state_diagnostics.
!   - The binned TAU / emulated warm rain diagnostics (amk*, ank*) are not output, as the CAM-SIMA
!     history does not yet support 3-D fields.
!
! FOR_REVIEWER_FEEDBACK: The following fields from CAM's micro_pumas_cam.F90 are NOT included in
!    this conversion.  Are any of them needed?
!     - TAU / emulated warm rain binned diagnostics (only when micro_mg_warm_rain is 'tau' or 'emulated'):
!         amk_c, ank_c, amk_r, ank_r, amk, ank, amk_out, ank_out
!           These have an extra size bin dimension ('bins_ncd' in CAM), and the CAM-SIMA history does
!           not yet support 3-D fields, so their registration and output are commented out
!           (marked !UNCOMMENT WHEN 3D IMPLEMENTED).  The other TAU / emulated fields are output.
!     - HPROGCLD: CAM registers this field in micro_pumas_cam.F90 but never outputs it when PUMAS
!         is used.  If wanted, it could be output from the microphysics heating rate (as is done in
!         rk_stratiform_diagnostics).
!     - The constituent surface fluxes and the "before physics"/"after physics" constituent fields
!         (sflxnam, apcnst, bpcnst): CAM registers these in micro_pumas_cam.F90, but outputs them
!         in cam_diagnostics.F90 (general physics diagnostics), not in the microphysics.  The
!         constituent mixing ratios themselves are output by sima_state_diagnostics.
!     - FICE_SCOL, MPDLIQ_SCOL, MPDICE_SCOL: subcolumns are not supported in CAM-SIMA.
module micro_pumas_diagnostics

   use ccpp_kinds,  only: kind_phys
   use pumas_kinds, only: pumas_r8 => kind_r8

   implicit none
   private
   save

   public :: micro_pumas_diagnostics_init         ! init routine
   public :: micro_pumas_pre_diagnostics_run      ! before microphysics routine
   public :: micro_pumas_diagnostics_run          ! main routine
   public :: micro_pumas_diagnostics_final        ! finalize routine

   ! Value of the tropopause level index where no tropopause was found
   ! (must match NOTFOUND in tropopause_find.F90)
   integer, parameter :: NOTFOUND = -1

   ! Parameters for the rainbow diagnostic
   real(kind_phys), parameter :: rb_rmin  = 1.e-6_kind_phys                                 ! Strat Rain threshold (mixing ratio)
   real(kind_phys), parameter :: rb_rcmin = 5._kind_phys/(86400._kind_phys*1000._kind_phys) ! Conv Rain Threshold (mm/d--> m/s)
   real(kind_phys), parameter :: rb_pmin  = 85000._kind_phys                                ! Minimum pressure for surface layer

   ! Accumulators for the precipitation efficiency calculation which persist across
   ! timesteps (these were the ACPRECL, ACGCME and ACNUM physics buffer fields in CAM)
   !
   ! FOR_REVIEWER_FEEDBACK: In CAM, ACPRECL, ACGCME and ACNUM were 'global' physics buffer
   !    fields, which CAM writes to and reads from the restart file, so PE, APRL and PEFRAC were
   !    restart exact.  These module variables are NOT written to the restart file, so a restart
   !    run will reset the accumulators to zero and PE, APRL and PEFRAC will differ from a
   !    continuous run until each column reaches its next "end of cloud" reset.
   !    Once restart capability is available in CAM-SIMA, these should be changed to:
   !     1. Add three host model (registry) variables which are written to and read from the
   !        restart file and initialized to zero for initial runs, e.g.:
   !          accumulated_stratiform_precipitation_rate_for_precipitation_efficiency (m s-1)
   !          accumulated_cloud_liquid_condensation_rate_for_precipitation_efficiency (m s-1)
   !          number_of_timesteps_accumulated_for_precipitation_efficiency (count, integer)
   !        (the standard names are suggestions and need to be added to the standard names
   !        dictionary)
   !     2. Pass them to micro_pumas_diagnostics_run as intent(inout) arguments with
   !        dimensions (horizontal_loop_extent) and use them in place of acprecl, acgcme and acnum.
   !     3. Remove these module variables, the allocation and zeroing in
   !        micro_pumas_diagnostics_init, the deallocation in micro_pumas_diagnostics_final
   !        (the final phase can then be removed), and the col_start/col_end arguments and the
   !        icol indexing in micro_pumas_diagnostics_run.
   !    Is this the right approach, and are there preferred standard names?
   real(kind_phys), allocatable :: acprecl(:) ! accumulated precip across timesteps (m s-1)
   real(kind_phys), allocatable :: acgcme(:)  ! accumulated condensation across timesteps (m s-1)
   integer,         allocatable :: acnum(:)   ! counter for number of timesteps accumulated

   ! Values of the 'trop_cld_lev' history coordinate.  This must persist for the whole run as the
   ! history coordinate points to it, so it is intentionally not deallocated in the final phase.
   real(kind_phys), allocatable, target :: trop_cld_levs(:)

CONTAINS

   !> \section arg_table_micro_pumas_diagnostics_init  Argument Table
   !! \htmlinclude micro_pumas_diagnostics_init.html
   subroutine micro_pumas_diagnostics_init(pcols, pver, top_lev, micro_mg_warm_rain, micro_mg_adjust_cpt, &
                                           errmsg, errflg)

      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only, add_vert_coord, hist_coords, get_hist_coord_index
!     use cam_history_support, only: add_hist_coord                          !UNCOMMENT WHEN 3D IMPLEMENTED
!     use pumas_stochastic_collect_tau, only: ncd                            !UNCOMMENT WHEN 3D IMPLEMENTED

      integer,            intent(in)  :: pcols              ! horizontal dimension
      integer,            intent(in)  :: pver
      integer,            intent(in)  :: top_lev            ! top model level for cloud physics (index)
      character(len=*),   intent(in)  :: micro_mg_warm_rain ! warm rain method (kk2000, sb2001, tau, emulated)
      logical,            intent(in)  :: micro_mg_adjust_cpt ! flag for the cold point tropopause temperature adjustment

      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables
      integer :: lev_idx
      integer :: ierr

      errmsg = ''
      errflg = 0

      ! Allocate and initialize the precipitation efficiency accumulators
      ! FOR_REVIEWER_FEEDBACK: Zeroing here is only correct for initial runs.  When restarts are
      !    implemented, remove this allocation and initialization (see the comment at the
      !    acprecl, acgcme and acnum declarations at the top of this module).
      allocate(acprecl(pcols), acgcme(pcols), acnum(pcols), stat=ierr)
      if (ierr /= 0) then
         errflg = 1
         errmsg = 'micro_pumas_diagnostics_init: error allocating precipitation efficiency accumulators'
         return
      end if
      acprecl(:) = 0._kind_phys
      acgcme(:)  = 0._kind_phys
      acnum(:)   = 0

      ! Add the 'trop_cld_lev' history vertical coordinate for the PUMAS levels (top_lev to pver).
      ! This coordinate is specific to the PUMAS output, so it is registered here (CAM registers it
      ! in micro_pumas_cam_register).
      lev_idx = get_hist_coord_index('lev')
      if (lev_idx <= 0) then
         errflg = 1
         errmsg = "micro_pumas_diagnostics_init: unable to find history coordinate 'lev'"
         return
      end if
      allocate(trop_cld_levs(pver-top_lev+1), stat=ierr)
      if (ierr /= 0) then
         errflg = 1
         errmsg = 'micro_pumas_diagnostics_init: error allocating trop_cld_levs'
         return
      end if
      trop_cld_levs(:) = hist_coords(lev_idx)%real_values(top_lev:pver)
      call add_vert_coord('trop_cld_lev', pver-top_lev+1,                            &
           'troposphere hybrid level at midpoints (1000*(A+B))', 'hPa', trop_cld_levs, &
           positive='down')

      ! NOTE -subcyc comments below mean that the original CAM code had "sampled_on_subcycle=.true." in the old outfld call
      !        Will use this for guidance during future development

      call history_add_field ('CME',      'Rate of cond-evap within the cloud',          'lev',         'avg','kg kg-1 s-1') !subcyc
      call history_add_field ('PRODPREC', 'Rate of conversion of condensate to precip',  'lev',         'avg','kg kg-1 s-1') !subcyc
      call history_add_field ('EVAPPREC', 'Rate of evaporation of falling precip',       'lev',         'avg','kg kg-1 s-1') !subcyc
      call history_add_field ('EVAPSNOW', 'Rate of evaporation of falling snow',         'trop_cld_lev','avg','kg kg-1 s-1') !subcyc
      call history_add_field ('FICE',     'Fractional ice content within cloud',         'lev',         'avg','fraction') !subcyc
      call history_add_field ('CLDFSNOW', 'Cloud fraction adjusted for snow',            'lev',         'avg','1') !subcyc
      call history_add_field ('ICWMRST',  'Prognostic in-stratus water mixing ratio',    'lev',         'avg','kg kg-1') !subcyc
      call history_add_field ('ICIMRST',  'Prognostic in-stratus ice mixing ratio',      'lev',         'avg','kg kg-1') !subcyc

      ! MG microphysics diagnostics
      call history_add_field ('QCSEVAP',  'Rate of evaporation of falling cloud water',  'trop_cld_lev','avg','kg kg-1 s-1') !subcyc
      call history_add_field ('QISEVAP',  'Rate of sublimation of falling cloud ice',    'trop_cld_lev','avg','kg kg-1 s-1') !subcyc
      call history_add_field ('QVRES',    'Rate of residual condensation term',          'trop_cld_lev','avg','kg kg-1 s-1') !subcyc
      call history_add_field ('CMEIOUT',  'Rate of deposition/sublimation of cloud ice', 'trop_cld_lev','avg','kg kg-1 s-1') !subcyc
      call history_add_field ('VTRMC',    'Mass-weighted cloud water fallspeed',         'trop_cld_lev','avg','m s-1') !subcyc
      call history_add_field ('VTRMI',    'Mass-weighted cloud ice fallspeed',           'trop_cld_lev','avg','m s-1') !subcyc

      call history_add_field ('QCSEDTEN', &
            'Cloud water mixing ratio tendency from sedimentation',                      'trop_cld_lev','avg','kg kg-1 s-1') !subcyc
      call history_add_field ('QISEDTEN', &
            'Cloud ice mixing ratio tendency from sedimentation',                        'trop_cld_lev','avg','kg kg-1 s-1') !subcyc

      call history_add_field ('PRAO',     'Accretion of cloud water by rain',            'lev', 'avg', 'kg kg-1 s-1') !subcyc
      call history_add_field ('PRCO',     'Autoconversion of cloud water',               'lev', 'avg', 'kg kg-1 s-1') !subcyc
      call history_add_field ('MNUCCCO',  'Immersion freezing of cloud water',           'lev', 'avg', 'kg kg-1 s-1') !subcyc
      call history_add_field ('MNUCCTO',  'Contact freezing of cloud water',             'lev', 'avg', 'kg kg-1 s-1') !subcyc

      call history_add_field ('MNUCCDO', &
            'Homogeneous and heterogeneous nucleation from vapor',                       'trop_cld_lev','avg','kg kg-1 s-1') !subcyc

      call history_add_field ('MNUCCDOhet', 'Heterogeneous nucleation from vapor',               'lev','avg','kg kg-1 s-1') !subcyc
      call history_add_field ('MSACWIO', 'Conversion of cloud water from rime-splintering',      'lev','avg','kg kg-1 s-1') !subcyc
      call history_add_field ('PSACWSO', 'Accretion of cloud water by snow',                     'lev','avg','kg kg-1 s-1') !subcyc
      call history_add_field ('BERGSO',  'Conversion of cloud water to snow from bergeron',      'lev','avg','kg kg-1 s-1') !subcyc
      call history_add_field ('BERGO',   'Conversion of cloud water to cloud ice from bergeron', 'lev','avg','kg kg-1 s-1') !subcyc

      call history_add_field ('MELTO',    'Melting of cloud ice',                         'lev',         'avg','kg kg-1 s-1')!subcyc
      call history_add_field ('MELTSTOT', 'Melting of snow',                              'trop_cld_lev','avg','kg kg-1 s-1')!subcyc
      call history_add_field ('MNUDEPO',  'Deposition Nucleation',                        'trop_cld_lev','avg','kg kg-1 s-1')!subcyc
      call history_add_field ('HOMOO',    'Homogeneous freezing of cloud water',          'lev',         'avg','kg kg-1 s-1')!subcyc
      call history_add_field ('QCRESO',   'Residual condensation term for cloud water',   'lev',         'avg','kg kg-1 s-1')!subcyc
      call history_add_field ('PRCIO',    'Autoconversion of cloud ice to snow',          'lev',         'avg','kg kg-1 s-1')!subcyc
      call history_add_field ('PRAIO',    'Accretion of cloud ice to snow',               'lev',         'avg','kg kg-1 s-1')!subcyc
      call history_add_field ('QIRESO',   'Residual deposition term for cloud ice',       'lev',         'avg','kg kg-1 s-1')!subcyc
      call history_add_field ('MNUCCRO',  'Heterogeneous freezing of rain to snow',       'trop_cld_lev','avg','kg kg-1 s-1')!subcyc
      call history_add_field ('MNUCCRIO', 'Heterogeneous freezing of rain to ice',        'trop_cld_lev','avg','kg kg-1 s-1')!subcyc
      call history_add_field ('PRACSO',   'Accretion of rain by snow',                    'trop_cld_lev','avg','kg kg-1 s-1')!subcyc
      call history_add_field ('VAPDEPSO', 'Vapor deposition onto snow',                   'trop_cld_lev','avg','kg kg-1 s-1')!subcyc
      call history_add_field ('MELTSDT',  'Latent heating rate due to melting of snow',   'trop_cld_lev','avg','W kg-1')!subcyc

      call history_add_field ('FRZRDT', &
            'Latent heating rate due to homogeneous freezing of rain',                    'trop_cld_lev','avg','W kg-1')!subcyc
      call history_add_field ('QRSEDTEN', &
            'Rain mixing ratio tendency from sedimentation',                              'trop_cld_lev','avg','kg kg-1 s-1')!subcyc
      call history_add_field ('QSSEDTEN', &
            'Snow mixing ratio tendency from sedimentation',                              'trop_cld_lev','avg','kg kg-1 s-1')!subcyc
      call history_add_field ('NNUCCCO', &
            'Number Tendency due to Immersion freezing of cloud water',                   'trop_cld_lev','avg','kg-1 s-1')!subcyc
      call history_add_field ('NNUCCTO', &
            'Number Tendency due to Contact freezing of cloud water',                     'trop_cld_lev','avg','kg-1 s-1')!subcyc
      call history_add_field ('NNUCCDO',  'Number Tendency due to Ice nucleation',        'trop_cld_lev','avg','kg-1 s-1')!subcyc
      call history_add_field ('NNUDEPO',  'Number Tendency due to Deposition Nucleation', 'trop_cld_lev','avg','kg-1 s-1')!subcyc
      call history_add_field ('NHOMO', &
            'Number Tendency due to Homogeneous freezing of cloud water',                 'trop_cld_lev','avg','kg-1 s-1')!subcyc
      call history_add_field ('NNUCCRO', &
            'Number Tendency due to heterogeneous freezing of rain to snow',              'trop_cld_lev','avg','kg-1 s-1')!subcyc
      call history_add_field ('NNUCCRIO', &
            'Number Tendency due to Heterogeneous freezing of rain to ice',               'trop_cld_lev','avg','kg-1 s-1')!subcyc
      call history_add_field ('NSACWIO', &
            'Number Tendency due to Ice Multiplication- Rime-splintering',                'trop_cld_lev','avg','kg-1 s-1')!subcyc
      call history_add_field ('NPRAO', &
            'Number Tendency due to Accretion of cloud water by rain',                    'trop_cld_lev','avg','kg-1 s-1')!subcyc
      call history_add_field ('NPSACWSO', &
            'Number Tendency due to Accretion of cloud water by snow',                    'trop_cld_lev','avg','kg-1 s-1')!subcyc
      call history_add_field ('NPRAIO', &
            'Number Tendency due to Accretion of cloud ice to snow',                      'trop_cld_lev','avg','kg-1 s-1')!subcyc
      call history_add_field ('NPRACSO', &
            'Number Tendency due to Accretion of rain by snow',                           'trop_cld_lev','avg','kg-1 s-1')!subcyc
      call history_add_field ('NPRCO', &
            'Number Tendency due to Autoconversion of cloud water [to rain]',             'trop_cld_lev','avg','kg-1 s-1')!subcyc
      call history_add_field ('NPRCIO', &
            'Number Tendency due to Autoconversion of cloud ice to snow',                 'trop_cld_lev','avg','kg-1 s-1')!subcyc
      call history_add_field ('NCSEDTEN', &
            'Number Tendency due to cloud liquid sedimentation',                          'trop_cld_lev','avg','kg-1 s-1')!subcyc

      call history_add_field ('NISEDTEN', 'Number Tendency due to cloud ice sedimentation', 'trop_cld_lev','avg','kg-1 s-1') !subcyc

      call history_add_field ('NRSEDTEN', 'Number Tendency due to rain sedimentation',      'trop_cld_lev','avg','kg-1 s-1') !subcyc
      call history_add_field ('NSSEDTEN', 'Number Tendency due to snow sedimentation',      'trop_cld_lev','avg','kg-1 s-1') !subcyc
      call history_add_field ('NMELTO',   'Number Tendency due to Melting of cloud ice',    'trop_cld_lev','avg','kg-1 s-1') !subcyc
      call history_add_field ('NMELTS',   'Number Tendency due to Melting of snow',         'trop_cld_lev','avg','kg-1 s-1') !subcyc

      if (trim(micro_mg_warm_rain) == 'kk2000') then
         call history_add_field ('qctend_KK2000', &
            'cloud liquid mass tendency due to autoconversion accretion from KK2000', 'trop_cld_lev', 'avg', 'kg kg-1 s-1') !subcyc
         call history_add_field ('nctend_KK2000', &
            'cloud number mass tendency due to autoconversion accretion from KK2000', 'trop_cld_lev', 'avg', 'kg-1 s-1') !subcyc
         call history_add_field ('qrtend_KK2000', &
            'rain mass tendency due to autoconversion accretion from KK2000',         'trop_cld_lev', 'avg', 'kg kg-1 s-1') !subcyc
         call history_add_field ('nrtend_KK2000', &
            'rain number tendency due to autoconversion accretion from KK2000',       'trop_cld_lev', 'avg', 'kg-1 s-1') !subcyc
      end if

      if (trim(micro_mg_warm_rain) == 'sb2001') then
         call history_add_field ('qctend_SB2001', &
            'cloud liquid mass tendency due to autoconversion  accretion from SB2001',  'trop_cld_lev','avg','kg kg-1 s-1') !subcyc
         call history_add_field ('nctend_SB2001', &
            'cloud liquid number tendency due to autoconversion accretion from SB2001', 'trop_cld_lev','avg','kg-1 s-1') !subcyc
         call history_add_field ('qrtend_SB2001', &
            'rain mass tendency due to autoconversion accretion from SB2001',           'trop_cld_lev','avg','kg kg-1 s-1') !subcyc
         call history_add_field ('nrtend_SB2001', &
            'rain number tendency due to autoconversion accretion from SB2001',         'trop_cld_lev','avg','kg-1 s-1') !subcyc
      end if

      ! TAU / emulated stochastic collection (machine learning) warm rain diagnostics
      ! (registered in CAM's stochastic_tau_cam_init / stochastic_emulated_init_cam)
      ! FOR_REVIEWER_FEEDBACK: The binned fields (amk*, ank*) have a size bin dimension, and CAM-SIMA
      !    history does not yet support 3-D fields, so they are commented out.  CAM also registers
      !    qctend_TAU_diag, nctend_TAU_diag, qrtend_TAU_diag and nrtend_TAU_diag, but never outputs them,
      !    so they are not included.  The units of nc_fixer and nr_fixer are copied from CAM (kg/kg), but
      !    they are number changes (kg-1?).  Should these be corrected?
      if (trim(micro_mg_warm_rain) == 'tau' .or. trim(micro_mg_warm_rain) == 'emulated') then
!        call add_hist_coord('bins_ncd', ncd, 'bins for TAU microphysics') !UNCOMMENT WHEN 3D IMPLEMENTED
!        call history_add_field ('amk_c',   'cloud liquid mass from bins',                 & !UNCOMMENT WHEN 3D IMPLEMENTED
!                                (/'trop_cld_lev', 'bins_ncd    '/), 'avg', 'kg') !UNCOMMENT WHEN 3D IMPLEMENTED
!        call history_add_field ('ank_c',   'cloud liquid number concentration from bins', & !UNCOMMENT WHEN 3D IMPLEMENTED
!                                (/'trop_cld_lev', 'bins_ncd    '/), 'avg', 'kg-1') !UNCOMMENT WHEN 3D IMPLEMENTED
!        call history_add_field ('amk_r',   'rain mass from bins',                         & !UNCOMMENT WHEN 3D IMPLEMENTED
!                                (/'trop_cld_lev', 'bins_ncd    '/), 'avg', 'kg') !UNCOMMENT WHEN 3D IMPLEMENTED
!        call history_add_field ('ank_r',   'rain number concentration from bins',         & !UNCOMMENT WHEN 3D IMPLEMENTED
!                                (/'trop_cld_lev', 'bins_ncd    '/), 'avg', 'kg-1') !UNCOMMENT WHEN 3D IMPLEMENTED
!        call history_add_field ('amk',     'all liquid mass from bins',                   & !UNCOMMENT WHEN 3D IMPLEMENTED
!                                (/'trop_cld_lev', 'bins_ncd    '/), 'avg', 'kg') !UNCOMMENT WHEN 3D IMPLEMENTED
!        call history_add_field ('ank',     'all liquid number concentration from bins',   & !UNCOMMENT WHEN 3D IMPLEMENTED
!                                (/'trop_cld_lev', 'bins_ncd    '/), 'avg', 'kg-1') !UNCOMMENT WHEN 3D IMPLEMENTED
!        call history_add_field ('amk_out', 'all liquid mass from bins',                   & !UNCOMMENT WHEN 3D IMPLEMENTED
!                                (/'trop_cld_lev', 'bins_ncd    '/), 'avg', 'kg') !UNCOMMENT WHEN 3D IMPLEMENTED
!        call history_add_field ('ank_out', 'all liquid number concentration from bins',   & !UNCOMMENT WHEN 3D IMPLEMENTED
!                                (/'trop_cld_lev', 'bins_ncd    '/), 'avg', 'kg-1') !UNCOMMENT WHEN 3D IMPLEMENTED

         call history_add_field ('scale_nc',   'scaling factor for nc',                'trop_cld_lev', 'avg', '1') !subcyc
         call history_add_field ('scale_nr',   'scaling factor for nr',                'trop_cld_lev', 'avg', '1') !subcyc
         call history_add_field ('scale_qc',   'scaling factor for qc',                'trop_cld_lev', 'avg', '1') !subcyc
         call history_add_field ('scale_qr',   'scaling factor for qr',                'trop_cld_lev', 'avg', '1') !subcyc

         call history_add_field ('QC_TAU_in',  'qc in TAU',                            'trop_cld_lev', 'avg', 'kg kg-1') !subcyc
         call history_add_field ('NC_TAU_in',  'nc in TAU',                            'trop_cld_lev', 'avg', 'kg-1') !subcyc
         call history_add_field ('QR_TAU_in',  'qr in TAU',                            'trop_cld_lev', 'avg', 'kg kg-1') !subcyc
         call history_add_field ('NR_TAU_in',  'nr in TAU',                            'trop_cld_lev', 'avg', 'kg-1') !subcyc
         call history_add_field ('QC_TAU_out', 'qc out TAU',                           'trop_cld_lev', 'avg', 'kg kg-1') !subcyc
         call history_add_field ('NC_TAU_out', 'nc out TAU',                           'trop_cld_lev', 'avg', 'kg-1') !subcyc
         call history_add_field ('QR_TAU_out', 'qr out TAU',                           'trop_cld_lev', 'avg', 'kg kg-1') !subcyc
         call history_add_field ('NR_TAU_out', 'nr out TAU',                           'trop_cld_lev', 'avg', 'kg-1') !subcyc

         call history_add_field ('qctend_TAU', 'qc tendency due to TAU bin code',      'trop_cld_lev', 'avg', 'kg kg-1 s-1') !subcyc
         call history_add_field ('nctend_TAU', 'nc tendency due to TAU bin code',      'trop_cld_lev', 'avg', 'kg-1 s-1') !subcyc
         call history_add_field ('qrtend_TAU', 'qr tendency due to TAU bin code',      'trop_cld_lev', 'avg', 'kg kg-1 s-1') !subcyc
         call history_add_field ('nrtend_TAU', 'nr tendency due to TAU bin code',      'trop_cld_lev', 'avg', 'kg-1 s-1') !subcyc

         call history_add_field ('gmnnn_lmnnn_TAU', &
            'sum of mass gain and loss from bin code',                                 'trop_cld_lev', 'avg', '1') !subcyc
         call history_add_field ('ML_fixer',   'frequency that ML fixer is activated', 'trop_cld_lev', 'avg', '1') !subcyc
         call history_add_field ('qc_fixer',   'delta qc due to ML fixer',             'trop_cld_lev', 'avg', 'kg kg-1') !subcyc
         call history_add_field ('nc_fixer',   'delta nc due to ML fixer',             'trop_cld_lev', 'avg', 'kg kg-1') !subcyc
         call history_add_field ('qr_fixer',   'delta qr due to ML fixer',             'trop_cld_lev', 'avg', 'kg kg-1') !subcyc
         call history_add_field ('nr_fixer',   'delta nr due to ML fixer',             'trop_cld_lev', 'avg', 'kg kg-1') !subcyc
      end if

      ! FOR_REVIEWER_FEEDBACK: The units for LAMC, LAMR and N0R are copied from CAM ('unitless'), but these
      !    are size distribution slopes and an intercept, which in SI units would be m-1 (LAMC, LAMR) and
      !    m-4 (N0R).  PUMAS does not document the units.  Should the units be corrected here, or should the
      !    CAM units be kept so that the output matches CAM?
      call history_add_field ('LAMC', 'Size distribution parameter lambda for liquid',    'trop_cld_lev', 'avg', '1') !subcyc
      call history_add_field ('LAMR', 'Size distribution parameter lambda for rain',      'trop_cld_lev', 'avg', '1') !subcyc
      call history_add_field ('PGAM', 'Size distribution parameter mu (pgam) for liquid', 'trop_cld_lev', 'avg', '1') !subcyc
      call history_add_field ('N0R',  'Size distribution parameter n0 for rain',          'trop_cld_lev', 'avg', '1') !subcyc

      ! Graupel/hail diagnostics (PUMAS in CAM-SIMA always carries graupel; these rates are
      ! zero if neither graupel nor hail is turned on)
      call history_add_field ('NMELTG',   'Number Tendency due to Melting of graupel',    'trop_cld_lev','avg','kg-1 s-1')!subcyc
      call history_add_field ('NGSEDTEN', 'Number Tendency due to graupel sedimentation', 'trop_cld_lev','avg','kg-1 s-1')!subcyc
      call history_add_field ('PSACRO', &
            'Collisions between rain & snow (Graupel collecting snow)',                   'lev',         'avg','kg kg-1 s-1')!subcyc
      call history_add_field ('PRACGO',   'Change in q collection rain by graupel',       'lev',         'avg','kg kg-1 s-1')!subcyc
      call history_add_field ('PSACWGO',  'Change in q collection droplets by graupel',   'lev',         'avg','kg kg-1 s-1')!subcyc
      call history_add_field ('PGSACWO', &
            'Q conversion to graupel due to collection droplets by snow',                 'lev',         'avg','kg kg-1 s-1')!subcyc
      call history_add_field ('PGRACSO', &
            'Q conversion to graupel due to collection rain by snow',                     'lev',         'avg','kg kg-1 s-1')!subcyc

      call history_add_field ('PRDGO',    'Deposition of graupel',                        'lev',         'avg','kg kg-1 s-1')!subcyc
      call history_add_field ('QMULTGO',  'Q change due to ice mult droplets/graupel',    'lev',         'avg','kg kg-1 s-1')!subcyc
      call history_add_field ('QMULTRGO', 'Q change due to ice mult rain/graupel',        'lev',         'avg','kg kg-1 s-1')!subcyc
      call history_add_field ('QGSEDTEN', &
            'Graupel/Hail mixing ratio tendency from sedimentation',                      'trop_cld_lev','avg','kg kg-1 s-1')!subcyc
      call history_add_field ('NPRACGO',  'Change N collection rain by graupel',          'lev',         'avg','kg-1 s-1')!subcyc
      call history_add_field ('NSCNGO', &
            'Change N conversion to graupel due to collection droplets by snow',          'lev',         'avg','kg-1 s-1')!subcyc
      call history_add_field ('NGRACSO', &
            'Change N conversion to graupel due to collection rain by snow',              'lev',         'avg','kg-1 s-1')!subcyc
      call history_add_field ('NMULTGO',  'Ice mult due to acc droplets by graupel',      'lev',         'avg','kg-1 s-1')!subcyc
      call history_add_field ('NMULTRGO', 'Ice mult due to acc rain by graupel',          'lev',         'avg','kg-1 s-1')!subcyc
      call history_add_field ('NPSACWGO', 'Change N collection droplets by graupel',      'lev',         'avg','kg-1 s-1')!subcyc
      call history_add_field ('CLDFGRAU', 'Cloud fraction adjusted for graupel',          'lev',         'avg','1')!subcyc
      call history_add_field ('MELTGTOT', 'Melting of graupel',                           'trop_cld_lev','avg','kg kg-1 s-1')!subcyc


      call history_add_field ('RBFRAC',  'Fraction of sky covered by a potential rainbow',     horiz_only,'avg','fraction') !subcyc
      call history_add_field ('RBFREQ',  'Potential rainbow frequency',                        horiz_only,'avg','1') !subcyc
      call history_add_field( 'rbSZA',   'solar zenith angle',                                 horiz_only,'inst','degrees') !subcyc

      ! History variables for CAM5 microphysics
      call history_add_field ('MPDT',    'Heating tendency - Morrison microphysics',           'lev', 'avg', 'W kg-1') !subcyc
      call history_add_field ('MPDQ',    'Q tendency - Morrison microphysics',                 'lev', 'avg', 'kg kg-1 s-1') !subcyc
      call history_add_field ('MPDLIQ',  'CLDLIQ tendency - Morrison microphysics',            'lev', 'avg', 'kg kg-1 s-1') !subcyc
      call history_add_field ('MPDICE',  'CLDICE tendency - Morrison microphysics',            'lev', 'avg', 'kg kg-1 s-1') !subcyc
      call history_add_field ('MPDNLIQ', 'NUMLIQ tendency - Morrison microphysics',            'lev', 'avg', 'kg-1 s-1') !subcyc
      call history_add_field ('MPDNICE', 'NUMICE tendency - Morrison microphysics',            'lev', 'avg', 'kg-1 s-1') !subcyc
      call history_add_field ('MPDW2V',  'Water <--> Vapor tendency - Morrison microphysics',  'lev', 'avg', 'kg kg-1 s-1') !subcyc
      call history_add_field ('MPDW2I',  'Water <--> Ice tendency - Morrison microphysics',    'lev', 'avg', 'kg kg-1 s-1') !subcyc
      call history_add_field ('MPDW2P',  'Water <--> Precip tendency - Morrison microphysics', 'lev', 'avg', 'kg kg-1 s-1') !subcyc
      call history_add_field ('MPDI2V',  'Ice <--> Vapor tendency - Morrison microphysics',    'lev', 'avg', 'kg kg-1 s-1') !subcyc
      call history_add_field ('MPDI2W',  'Ice <--> Water tendency - Morrison microphysics',    'lev', 'avg', 'kg kg-1 s-1') !subcyc
      call history_add_field ('MPDI2P',  'Ice <--> Precip tendency - Morrison microphysics',   'lev', 'avg', 'kg kg-1 s-1') !subcyc
      call history_add_field ('ICWNC',   'Prognostic in-cloud water number conc',              'lev', 'avg', 'm-3') !subcyc
      call history_add_field ('ICINC',   'Prognostic in-cloud ice number conc',                'lev', 'avg', 'm-3') !subcyc

      call history_add_field ('EFFLIQ_IND', 'Prognostic droplet effective radius (indirect effect)', 'lev', 'avg', 'um') !subcyc

      call history_add_field ('CDNUMC', 'Vertically-integrated droplet concentration', horiz_only, 'avg', 'm-2') !subcyc
      call history_add_field ('MPICLWPI', &
            'Vertically-integrated in-cloud Initial Liquid WP (Before Micro)',         horiz_only, 'avg', 'kg m-2') !subcyc
      call history_add_field ('MPICIWPI', &
            'Vertically-integrated in-cloud Initial Ice WP (Before Micro)',            horiz_only, 'avg', 'kg m-2') !subcyc


      ! This is only if the coldpoint temperatures are being adjusted.
      ! NOTE: Some fields related to these and output later are added in tropopause.F90.
      ! (bug fix: CAM had units of 'K' and 'Hunidity' for TROPF_RHADJ, which is in percent)
      if (micro_mg_adjust_cpt) then
        call history_add_field ('TROPF_TADJ', 'Temperatures after cold point adjustment',      'lev',        'avg', 'K') !subcyc
        call history_add_field ('TROPF_RHADJ','Relative Humidity after cold point adjustment', 'lev',        'avg', 'percent') !subcyc
        call history_add_field ('TROPF_CDT',  'Cold point temperature adjustment',              horiz_only,  'avg', 'K') !subcyc
        call history_add_field ('TROPF_CDZ',  'Distance of coldpoint from coldest model level', horiz_only,  'avg', 'm') !subcyc
      end if


      ! Averaging for cloud particle number and size
      call history_add_field ('AWNC',   'Average cloud water number conc',             'lev', 'avg', 'm-3') !subcyc
      call history_add_field ('AWNI',   'Average cloud ice number conc',               'lev', 'avg', 'm-3') !subcyc
      call history_add_field ('AREL',   'Average droplet effective radius',            'lev', 'avg', 'um') !subcyc
      call history_add_field ('AREI',   'Average ice effective radius',                'lev', 'avg', 'um') !subcyc
      ! Frequency arrays for above
      call history_add_field ('FREQL',  'Fractional occurrence of liquid',             'lev', 'avg', 'fraction') !subcyc
      call history_add_field ('FREQI',  'Fractional occurrence of ice',                'lev', 'avg', 'fraction') !subcyc


      ! Average cloud top particle size and number (liq, ice) and frequency
      call history_add_field ('ACTREL', 'Average Cloud Top droplet effective radius',  horiz_only, 'avg', 'um') !subcyc
      call history_add_field ('ACTREI', 'Average Cloud Top ice effective radius',      horiz_only, 'avg', 'um') !subcyc
      call history_add_field ('ACTNL',  'Average Cloud Top droplet number',            horiz_only, 'avg', 'm-3') !subcyc
      call history_add_field ('ACTNI',  'Average Cloud Top ice number',                horiz_only, 'avg', 'm-3') !subcyc

      call history_add_field ('FCTL',   'Fractional occurrence of cloud top liquid',   horiz_only, 'avg', 'fraction') !subcyc
      call history_add_field ('FCTI',   'Fractional occurrence of cloud top ice',      horiz_only, 'avg', 'fraction') !subcyc

      ! New frequency arrays for mixed phase and supercooled liquid (only and mixed) for (a) Cloud Top and (b) everywhere..
      call history_add_field ('FREQM',   'Fractional occurrence of mixed phase',                  'lev', 'avg', 'fraction') !subcyc
      call history_add_field ('FREQSL',  'Fractional occurrence of only supercooled liquid',      'lev', 'avg', 'fraction') !subcyc
      call history_add_field ('FREQSLM', 'Fractional occurrence of super cooled liquid with ice', 'lev', 'avg', 'fraction') !subcyc

      call history_add_field ('FCTM',        'Fractional occurrence of cloud top mixed phase',  horiz_only,'avg','fraction') !subcyc
      call history_add_field ('FCTSL', &
            'Fractional occurrence of cloud top only supercooled liquid',                       horiz_only,'avg','fraction') !subcyc
      call history_add_field ('FCTSLM', &
            'Fractional occurrence of cloud top super cooled liquid with ice',                  horiz_only,'avg','fraction') !subcyc

      call history_add_field ('LS_FLXPRC',   'ls stratiform gbm interface rain+snow flux',      'ilev', 'avg', 'kg m-2 s-1') !subcyc
      call history_add_field ('LS_FLXSNW',   'ls stratiform gbm interface snow flux',           'ilev', 'avg', 'kg m-2 s-1') !subcyc

      call history_add_field ('REL',         'MG REL stratiform cloud effective radius liquid', 'lev', 'avg', 'um') !subcyc
      call history_add_field ('REI',         'MG REI stratiform cloud effective radius ice',    'lev', 'avg', 'um') !subcyc
      call history_add_field ('LS_REFFRAIN', 'ls stratiform rain effective radius',             'lev', 'avg', 'um') !subcyc
      call history_add_field ('LS_REFFSNOW', 'ls stratiform snow effective radius',             'lev', 'avg', 'um') !subcyc
      call history_add_field ('CV_REFFLIQ',  'convective cloud liq effective radius',           'lev', 'avg', 'um') !subcyc
      call history_add_field ('CV_REFFICE',  'convective cloud ice effective radius',           'lev', 'avg', 'um') !subcyc
      call history_add_field ('MG_SADICE',   'MG surface area density ice',                     'lev', 'avg', 'cm2 cm-3') !subcyc
      call history_add_field ('MG_SADSNOW',  'MG surface area density snow',                    'lev', 'avg', 'cm2 cm-3') !subcyc

      ! diagnostic precip
      call history_add_field ('QRAIN',       'Diagnostic grid-mean rain mixing ratio',          'lev', 'avg', 'kg kg-1') !subcyc
      call history_add_field ('QSNOW',       'Diagnostic grid-mean snow mixing ratio',          'lev', 'avg', 'kg kg-1') !subcyc
      call history_add_field ('NRAIN',       'Diagnostic grid-mean rain number conc',           'lev', 'avg', 'm-3') !subcyc
      call history_add_field ('NSNOW',       'Diagnostic grid-mean snow number conc',           'lev', 'avg', 'm-3') !subcyc

      ! size of precip
      ! FOR_REVIEWER_FEEDBACK: PUMAS computes rercld in meters (calc_rercld in micro_pumas_v1.F90), and the
      !    'm' units here are correct.  However, micro_pumas_ccpp.meta (PUMAS) and
      !    micro_pumas_ccpp_dimensions_post.meta declare it in 'um'.  The rercld entry in
      !    micro_pumas_diagnostics.meta uses 'um' to match, so that the CCPP framework does not apply a
      !    unit conversion.  Should the units be fixed to 'm' in the PUMAS and dimensions_post metadata
      !    (and then here)?
      call history_add_field ('RERCLD',    'Diagnostic effective radius of Liquid Cloud and Rain',    'lev', 'avg', 'm') !subcyc
      ! FOR_REVIEWER_FEEDBACK: The DSNOW units are copied from CAM ('m'), but the values are in microns
      !    (CAM multiplied des by 1.e6 before output, and pumas_post_main does the same).  Should the units
      !    be corrected to 'um', or should the CAM units be kept so that the output matches CAM?
      call history_add_field ('DSNOW',     'Diagnostic grid-mean snow diameter',                      'lev', 'avg', 'm') !subcyc

      ! diagnostic radar reflectivity, cloud-averaged
      call history_add_field ('REFL',      '94 GHz radar reflectivity',                               'lev','avg','dBZ')!subcyc
      call history_add_field ('AREFL',     'Average 94 GHz radar reflectivity',                       'lev','avg','dBZ')!subcyc
      call history_add_field ('FREFL',     'Fractional occurrence of radar reflectivity',             'lev','avg','fraction')!subcyc

      call history_add_field ('CSRFL',     '94 GHz radar reflectivity (CloudSat thresholds)',         'lev', 'avg', 'dBZ') !subcyc
      call history_add_field ('ACSRFL',    'Average 94 GHz radar reflectivity (CloudSat thresholds)', 'lev', 'avg', 'dBZ') !subcyc

      call history_add_field ('FCSRFL', &
            'Fractional occurrence of radar reflectivity (CloudSat thresholds)',                      'lev','avg','fraction')!subcyc

      call history_add_field ('AREFLZ',    'Average 94 GHz radar reflectivity',                       'lev','avg','mm6 m-3') !subcyc

      ! 10cm (rain) radar reflectivity
      call history_add_field ('REFL10CM',  '10cm (Rain) radar reflectivity (Dbz)',                    'lev','avg','dBZ') !subcyc
      call history_add_field ('REFLZ10CM', '10cm (Rain) radar reflectivity (Z units)',                'lev','avg','mm6 m-3') !subcyc

      ! Aerosol information
      call history_add_field ('NCAL',      'Number Concentation Activated for Liquid',                'lev', 'avg', 'm-3') !subcyc
      call history_add_field ('NCAI',      'Number Concentation Activated for Ice',                   'lev', 'avg', 'm-3') !subcyc

      ! Average rain and snow mixing ratio (Q), number (N) and diameter (D), with frequency
      call history_add_field ('AQRAIN',    'Average rain mixing ratio',                               'lev','avg','kg kg-1')!subcyc
      call history_add_field ('AQSNOW',    'Average snow mixing ratio',                               'lev','avg','kg kg-1')!subcyc
      call history_add_field ('ANRAIN',    'Average rain number conc',                                'lev','avg','m-3')!subcyc
      call history_add_field ('ANSNOW',    'Average snow number conc',                                'lev','avg','m-3')!subcyc
      call history_add_field ('ADRAIN',    'Average rain effective Diameter',                         'lev','avg','m')!subcyc
      call history_add_field ('ADSNOW',    'Average snow effective Diameter',                         'lev','avg','m')!subcyc
      call history_add_field ('FREQR',     'Fractional occurrence of rain',                           'lev','avg','fraction')!subcyc
      call history_add_field ('FREQS',     'Fractional occurrence of snow',                           'lev','avg','fraction')!subcyc

      ! precipitation efficiency & other diagnostic fields
      call history_add_field('PE'    , 'Stratiform Precipitation Efficiency  (precip/cmeliq)', horiz_only, 'avg', '1') !subcyc
      call history_add_field('APRL'  , &
            'Average Stratiform Precip Rate over efficiency calculation',                      horiz_only, 'avg', 'm s-1') !subcyc

      call history_add_field('PEFRAC', 'Fraction of timesteps precip efficiency reported',     horiz_only, 'avg', '1') !subcyc

      call history_add_field('VPRCO' , 'Vertical average of autoconversion rate',             horiz_only,'avg','kg kg-1 s-1')!subcyc
      call history_add_field('VPRAO' , 'Vertical average of accretion rate',                  horiz_only,'avg','kg kg-1 s-1')!subcyc
      call history_add_field('RACAU', 'Accretion/autoconversion ratio from vertical average', horiz_only,'avg','kg kg-1 s-1')!subcyc

      call history_add_field('UMR', 'Mass-weighted rain  fallspeed',                          'trop_cld_lev','avg','m s-1') !subcyc
      call history_add_field('UMS', 'Mass-weighted snow fallspeed',                           'trop_cld_lev','avg','m s-1') !subcyc

      call history_add_field('UMG',         'Mass-weighted graupel/hail  fallspeed',       'trop_cld_lev','avg','m s-1') !subcyc
      call history_add_field('FREQG',       'Fractional occurrence of Graupel',            'lev',         'avg','fraction') !subcyc
      call history_add_field('LS_REFFGRAU', 'ls stratiform graupel/hail effective radius', 'lev',         'avg','um') !subcyc
      call history_add_field('AQGRAU',      'Average graupel/hail mixing ratio',           'lev',         'avg','kg kg-1') !subcyc
      call history_add_field('ANGRAU',      'Average graupel/hail number conc',            'lev',         'avg','m-3') !subcyc


      ! qc limiter (only output in versions 1.5 and later)
      call history_add_field('QCRAT',       'Qc Limiter: Fraction of qc tendency applied', 'lev', 'avg', 'fraction') !subcyc

   end subroutine micro_pumas_diagnostics_init

   !> \section arg_table_micro_pumas_pre_diagnostics_run  Argument Table
   !! \htmlinclude micro_pumas_pre_diagnostics_run.html
   subroutine micro_pumas_pre_diagnostics_run(ncol, pver, top_lev, gravit, ast, cldliq, cldice, pdel, &
                                              micro_mg_adjust_cpt, temp, zm, tropLev_cpp, tropT_cpp,  &
                                              tropZ_cpp, errmsg, errflg)

      use cam_history,       only: history_out_field
      use micro_pumas_utils, only: mincld

      integer,            intent(in)  :: ncol
      integer,            intent(in)  :: pver
      integer,            intent(in)  :: top_lev      ! top model level for cloud physics (index)
      real(kind_phys),    intent(in)  :: gravit       ! gravitational acceleration (m s-2)
      real(kind_phys),    intent(in)  :: ast(:,:)     ! stratiform cloud area fraction (fraction)
      real(kind_phys),    intent(in)  :: cldliq(:,:)  ! cloud liquid water mixing ratio before microphysics (kg kg-1)
      real(kind_phys),    intent(in)  :: cldice(:,:)  ! cloud ice mixing ratio before microphysics (kg kg-1)
      real(kind_phys),    intent(in)  :: pdel(:,:)    ! air pressure thickness (Pa)
      logical,            intent(in)  :: micro_mg_adjust_cpt ! flag for the cold point tropopause temperature adjustment
      real(kind_phys),    intent(in)  :: temp(:,:)    ! air temperature (K)
      real(kind_phys),    intent(in)  :: zm(:,:)      ! geopotential height wrt surface (m)
      integer,            intent(in)  :: tropLev_cpp(:) ! cold point tropopause vertical layer index (index)
      real(kind_phys),    intent(in)  :: tropT_cpp(:) ! cold point tropopause air temperature (K)
      real(kind_phys),    intent(in)  :: tropZ_cpp(:) ! cold point tropopause geopotential height wrt surface (m)

      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables
      real(kind_phys) :: iclwpi(ncol) ! vertically-integrated in-cloud liquid water path before microphysics (kg m-2)
      real(kind_phys) :: iciwpi(ncol) ! vertically-integrated in-cloud ice water path before microphysics (kg m-2)
      real(kind_phys) :: cp_dt(ncol)  ! cold point temperature adjustment (K)
      real(kind_phys) :: cp_dz(ncol)  ! distance of the cold point from the coldest model level (m)
      real(kind_phys) :: tadj(ncol,pver) ! air temperature after the cold point adjustment (K)
      integer         :: i, k

      errmsg = ''
      errflg = 0

      iclwpi(:) = 0._kind_phys
      iciwpi(:) = 0._kind_phys

      do k = top_lev, pver
         iclwpi(:) = iclwpi(:) + &
              min(cldliq(:ncol,k) / max(real(mincld, kind_phys), ast(:ncol,k)), 0.005_kind_phys) &
              * pdel(:ncol,k) / gravit
         iciwpi(:) = iciwpi(:) + &
              min(cldice(:ncol,k) / max(real(mincld, kind_phys), ast(:ncol,k)), 0.005_kind_phys) &
              * pdel(:ncol,k) / gravit
      end do

      ! Because of the of limited vertical resolution, there can be a signifcant
      ! warm bias at the cold point tropopause, which can create a wet bias in the
      ! stratosphere. For the microphysics only, update the cold point temperature, with
      ! an estimate of the coldest point between the model layers.
      ! (The cold point tropopause is found by tropopause_find, the adjusted temperature is
      !  passed to PUMAS by micro_pumas_ccpp_dimensions_pre, and it is included in the
      !  post-microphysics air temperature in micro_pumas_diagnostics_run.)
      if (micro_mg_adjust_cpt) then
         cp_dt(:ncol)         = 0._kind_phys
         cp_dz(:ncol)         = 0._kind_phys
         tadj(:ncol,:)        = temp(:ncol,:)

         do i = 1, ncol

            ! Update statistics and output results.
            if (tropLev_cpp(i) /= NOTFOUND) then
               cp_dt(i) = tropT_cpp(i) - tadj(i,tropLev_cpp(i))
               cp_dz(i) = tropZ_cpp(i) - zm(i,tropLev_cpp(i))

               ! NOTE: This change in temperature is just for the microphysics
               ! and should not be added to any tendencies or used to update
               ! any states
               tadj(i,tropLev_cpp(i)) = tadj(i,tropLev_cpp(i)) + cp_dt(i)
            end if
         end do

         ! Output all of the statistics related to the cold point
         ! tropopause adjustment. Th cold point information itself is
         ! output in tropopause.F90.
         call history_out_field('TROPF_TADJ', tadj)
         call history_out_field('TROPF_CDT',  cp_dt)
         call history_out_field('TROPF_CDZ',  cp_dz)
      end if

      call history_out_field('MPICLWPI', iclwpi)
      call history_out_field('MPICIWPI', iciwpi)

   end subroutine micro_pumas_pre_diagnostics_run

   !> \section arg_table_micro_pumas_diagnostics_run  Argument Table
   !! \htmlinclude micro_pumas_diagnostics_run.html
   subroutine micro_pumas_diagnostics_run(ncol, pver, top_lev, col_start, col_end, dtime,  &
                   gravit, rair, tmelt, pi, rhoh2o, cpairv, micro_mg_warm_rain,             &
                   pmid, pdel, temp, q, cldliq, cldice, numliq, numice, qr, nr, qs, ns, qg, &
                   ng, ast, cld, concld, cmeliq, prec_dp, prec_sh, naai, naai_hom,          &
                   calday, lat, lon, eccen, obliqr, lambm0, mvelpp,                         &
                   micro_mg_adjust_cpt, tropLev_cpp, tropT_cpp,                             &
                   proc_rates, tlat, qvlat, qcten, qiten, ncten, niten,                     &
                   effi, rel, des, sadice, sadsnow, nevapr, prain, prect,                   &
                   lflx, iflx, gflx, rflx, sflx, qrout, qsout, qgout, nrout, nsout,         &
                   refl, arefl, areflz, frefl, csrfl, acsrfl, fcsrfl, refl10cm, reflz10cm,  &
                   rercld, ncai, ncal, qrout2, qsout2, nrout2, nsout2, qgout2, ngout2,      &
                   freqg, freqs, freqr, nfice, qcrat, errmsg, errflg)

      use cam_history,       only: history_out_field
      use micro_pumas_diags, only: proc_rates_type
      use micro_pumas_utils, only: size_dist_param_liq, mg_liq_props, avg_diameter
      use micro_pumas_utils, only: rhow, rhosn, rhog, qsmall, mincld
      use shr_orb_mod,       only: shr_orb_decl, shr_orb_cosz
      use wv_saturation,     only: qsat

      !------------------------------------------------
      !   Input / output parameters
      !------------------------------------------------
      integer,          intent(in) :: ncol
      integer,          intent(in) :: pver
      integer,          intent(in) :: top_lev            ! top model level for cloud physics (index)
      integer,          intent(in) :: col_start          ! starting column index of this chunk
      integer,          intent(in) :: col_end            ! ending column index of this chunk
      real(kind_phys),  intent(in) :: dtime              ! physics timestep (s)
      real(kind_phys),  intent(in) :: gravit             ! gravitational acceleration (m s-2)
      real(kind_phys),  intent(in) :: rair               ! gas constant of dry air (J kg-1 K-1)
      real(kind_phys),  intent(in) :: tmelt              ! freezing point of water (K)
      real(kind_phys),  intent(in) :: pi                 ! pi
      real(kind_phys),  intent(in) :: rhoh2o             ! density of liquid water (kg m-3)
      real(kind_phys),  intent(in) :: cpairv(:,:)        ! specific heat of dry air (J kg-1 K-1)
      character(len=*), intent(in) :: micro_mg_warm_rain ! warm rain method (kk2000, sb2001, tau, emulated)

      ! State variables (constituents have been updated by microphysics, air temperature has not)
      real(kind_phys),  intent(in) :: pmid(:,:)          ! air pressure (Pa)
      real(kind_phys),  intent(in) :: pdel(:,:)          ! air pressure thickness (Pa)
      real(kind_phys),  intent(in) :: temp(:,:)          ! air temperature before microphysics heating (K)
      real(kind_phys),  intent(in) :: q(:,:)             ! water vapor mixing ratio (kg kg-1)
      real(kind_phys),  intent(in) :: cldliq(:,:)        ! cloud liquid water mixing ratio (kg kg-1)
      real(kind_phys),  intent(in) :: cldice(:,:)        ! cloud ice mixing ratio (kg kg-1)
      real(kind_phys),  intent(in) :: numliq(:,:)        ! cloud liquid number concentration (kg-1)
      real(kind_phys),  intent(in) :: numice(:,:)        ! cloud ice number concentration (kg-1)
      real(kind_phys),  intent(in) :: qr(:,:)            ! rain mixing ratio (kg kg-1)
      real(kind_phys),  intent(in) :: nr(:,:)            ! rain number concentration (kg-1)
      real(kind_phys),  intent(in) :: qs(:,:)            ! snow mixing ratio (kg kg-1)
      real(kind_phys),  intent(in) :: ns(:,:)            ! snow number concentration (kg-1)
      real(kind_phys),  intent(in) :: qg(:,:)            ! graupel mixing ratio (kg kg-1)
      real(kind_phys),  intent(in) :: ng(:,:)            ! graupel number concentration (kg-1)

      ! Inputs from other schemes
      real(kind_phys),  intent(in) :: ast(:,:)           ! stratiform cloud area fraction (fraction)
      real(kind_phys),  intent(in) :: cld(:,:)           ! cloud area fraction (fraction)
      real(kind_phys),  intent(in) :: concld(:,:)        ! convective cloud area fraction (fraction)
      real(kind_phys),  intent(in) :: cmeliq(:,:)        ! rate of cond-evap of cloud liquid from macrophysics (kg kg-1 s-1)
      real(kind_phys),  intent(in) :: prec_dp(:)         ! surface precip rate from deep convection (m s-1)
      real(kind_phys),  intent(in) :: prec_sh(:)         ! surface precip rate from shallow convection (m s-1)
      real(kind_phys),  intent(in) :: naai(:,:)          ! number of activated ice nuclei (kg-1 s-1)
      real(kind_phys),  intent(in) :: naai_hom(:,:)      ! number of activated ice nuclei from homogeneous freezing (kg-1 s-1)

      ! Inputs for the solar zenith angle
      real(kind_phys),  intent(in) :: calday             ! fractional calendar day (days)
      real(kind_phys),  intent(in) :: lat(:)             ! latitude (rad)
      real(kind_phys),  intent(in) :: lon(:)             ! longitude (rad)
      real(kind_phys),  intent(in) :: eccen              ! Earth's orbital eccentricity factor (1)
      real(kind_phys),  intent(in) :: obliqr             ! Earth's obliquity (rad)
      real(kind_phys),  intent(in) :: lambm0             ! mean longitude of perihelion at the vernal equinox (rad)
      real(kind_phys),  intent(in) :: mvelpp             ! moving vernal equinox longitude of perihelion plus pi (rad)
      logical,          intent(in) :: micro_mg_adjust_cpt ! flag for the cold point tropopause temperature adjustment
      integer,          intent(in) :: tropLev_cpp(:)     ! cold point tropopause vertical layer index (index)
      real(kind_phys),  intent(in) :: tropT_cpp(:)       ! cold point tropopause air temperature (K)

      ! PUMAS outputs
      type(proc_rates_type), intent(in) :: proc_rates    ! PUMAS process rates (on PUMAS vertical levels)
      real(kind_phys),  intent(in) :: tlat(:,:)          ! latent heating rate (J kg-1 s-1)
      ! The microphysical tendencies are the PUMAS outputs on the PUMAS vertical levels: the
      ! constituent tendencies on the host model levels have already been applied to the state and
      ! reset to zero (by apply_constituent_tendencies) by the time this scheme runs
      real(pumas_r8),   intent(in) :: qvlat(:,:)         ! microphysical tendency qv (kg kg-1 s-1)
      real(pumas_r8),   intent(in) :: qcten(:,:)         ! microphysical tendency qc (kg kg-1 s-1)
      real(pumas_r8),   intent(in) :: qiten(:,:)         ! microphysical tendency qi (kg kg-1 s-1)
      real(pumas_r8),   intent(in) :: ncten(:,:)         ! microphysical tendency nc (kg-1 s-1)
      real(pumas_r8),   intent(in) :: niten(:,:)         ! microphysical tendency ni (kg-1 s-1)
      real(kind_phys),  intent(in) :: effi(:,:)          ! cloud ice effective radius, recomputed in pumas_post_main (um)
      real(kind_phys),  intent(in) :: rel(:,:)           ! cloud liquid effective radius, recomputed in pumas_post_main (um)
      real(kind_phys),  intent(in) :: des(:,:)           ! snow effective diameter, recomputed in pumas_post_main (um)
      real(kind_phys),  intent(in) :: sadice(:,:)        ! cloud ice surface area density (cm2 cm-3)
      real(kind_phys),  intent(in) :: sadsnow(:,:)       ! cloud snow surface area density (cm2 cm-3)
      real(kind_phys),  intent(in) :: nevapr(:,:)        ! evaporation rate of rain + snow (kg kg-1 s-1)
      real(kind_phys),  intent(in) :: prain(:,:)         ! production of rain + snow (kg kg-1 s-1)
      real(kind_phys),  intent(in) :: prect(:)           ! surface precip rate (m s-1)
      real(kind_phys),  intent(in) :: lflx(:,:)          ! grid-box average liquid condensate flux (kg m-2 s-1)
      real(kind_phys),  intent(in) :: iflx(:,:)          ! grid-box average ice condensate flux (kg m-2 s-1)
      real(kind_phys),  intent(in) :: gflx(:,:)          ! grid-box average graupel/hail flux (kg m-2 s-1)
      real(kind_phys),  intent(in) :: rflx(:,:)          ! grid-box average rain flux (kg m-2 s-1)
      real(kind_phys),  intent(in) :: sflx(:,:)          ! grid-box average snow flux (kg m-2 s-1)
      real(kind_phys),  intent(in) :: qrout(:,:)         ! grid-box average rain mixing ratio (kg kg-1)
      real(kind_phys),  intent(in) :: qsout(:,:)         ! grid-box average snow mixing ratio (kg kg-1)
      real(kind_phys),  intent(in) :: qgout(:,:)         ! grid-box average graupel/hail mixing ratio (kg kg-1)
      real(kind_phys),  intent(in) :: nrout(:,:)         ! rain number concentration (m-3)
      real(kind_phys),  intent(in) :: nsout(:,:)         ! snow number concentration (m-3)
      real(kind_phys),  intent(in) :: refl(:,:)          ! analytic radar reflectivity (94GHZ, cloud radar) (dBZ)
      real(kind_phys),  intent(in) :: arefl(:,:)         ! average reflectivity will zero points inside valid range (dBZ)
      real(kind_phys),  intent(in) :: areflz(:,:)        ! average reflectivity in z. (mm6 m-3)
      real(kind_phys),  intent(in) :: frefl(:,:)         ! fractional occurrence of radar reflectivity (fraction)
      real(kind_phys),  intent(in) :: csrfl(:,:)         ! cloudsat reflectivity (dBZ)
      real(kind_phys),  intent(in) :: acsrfl(:,:)        ! cloudsat average (dBZ)
      real(kind_phys),  intent(in) :: fcsrfl(:,:)        ! cloudsat fractional occurrence of radar reflectivity (fraction)
      real(kind_phys),  intent(in) :: refl10cm(:,:)      ! 10cm (rain) analytic radar reflectivity (dBZ)
      real(kind_phys),  intent(in) :: reflz10cm(:,:)     ! 10cm (rain) analytic radar reflectivity (mm6 m-3)
      real(kind_phys),  intent(in) :: rercld(:,:)        ! effective radius calculation for rain + cloud (m)
      real(kind_phys),  intent(in) :: ncai(:,:)          ! number conc of ice nuclei available (m-3)
      real(kind_phys),  intent(in) :: ncal(:,:)          ! number conc of CCN (m-3)
      real(kind_phys),  intent(in) :: qrout2(:,:)        ! average rain mixing ratio (kg kg-1)
      real(kind_phys),  intent(in) :: qsout2(:,:)        ! average snow mixing ratio (kg kg-1)
      real(kind_phys),  intent(in) :: nrout2(:,:)        ! average rain number concentration (m-3)
      real(kind_phys),  intent(in) :: nsout2(:,:)        ! average snow number concentration (m-3)
      real(kind_phys),  intent(in) :: qgout2(:,:)        ! average graupel/hail mixing ratio (kg kg-1)
      real(kind_phys),  intent(in) :: ngout2(:,:)        ! average graupel/hail number concentration (m-3)
      real(kind_phys),  intent(in) :: freqg(:,:)         ! fractional occurrence of graupel (fraction)
      real(kind_phys),  intent(in) :: freqs(:,:)         ! fractional occurrence of snow (fraction)
      real(kind_phys),  intent(in) :: freqr(:,:)         ! fractional occurrence of rain (fraction)
      real(kind_phys),  intent(in) :: nfice(:,:)         ! fraction of frozen water to total condensed water (fraction)
      real(kind_phys),  intent(in) :: qcrat(:,:)         ! limiter for qc process rates (1=no limit --> 0. no qc) (fraction)

      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      !------------------------------------------------
      !   Local variables
      !------------------------------------------------
      ! Work array for outputting fields which are on the PUMAS vertical levels
      real(kind_phys) :: work(ncol,pver)

      ! Process rate combinations (on PUMAS vertical levels)
      real(pumas_r8)  :: ftem(ncol,pver-top_lev+1)

      real(kind_phys) :: t_post(ncol,pver)            ! air temperature after microphysics heating (K)
      real(kind_phys) :: cp_dt                        ! cold point temperature adjustment (K)
      real(kind_phys) :: cp_rh(ncol,pver)             ! relative humidity after cold point adjustment (percent)
      real(kind_phys) :: es                           ! saturation vapor pressure (Pa)
      real(kind_phys) :: qsat_out                     ! saturation specific humidity (kg kg-1)
      real(kind_phys) :: mnuccdohet(ncol,pver)      ! heterogeneous nucleation from vapor (kg kg-1 s-1)
      real(kind_phys) :: mgflxprc(ncol,pver+1)        ! stratiform rain+snow flux (kg m-2 s-1)
      real(kind_phys) :: mgflxsnw(ncol,pver+1)        ! stratiform snow flux (kg m-2 s-1)
      real(kind_phys) :: cvreffliq(ncol,pver)         ! convective cloud liquid effective radius (um)
      real(kind_phys) :: cvreffice(ncol,pver)         ! convective cloud ice effective radius (um)
      real(kind_phys) :: refl10cm_out(ncol,pver)      ! 10cm radar reflectivity with fill above top_lev (dBZ)

      real(kind_phys) :: icimrst(ncol,pver)           ! in-stratus ice mixing ratio (kg kg-1)
      real(kind_phys) :: icwmrst(ncol,pver)           ! in-stratus water mixing ratio (kg kg-1)
      real(kind_phys) :: icimrst_out(ncol,pver)       ! in-stratus ice mixing ratio for output (kg kg-1)
      real(kind_phys) :: icwmrst_out(ncol,pver)       ! in-stratus water mixing ratio for output (kg kg-1)
      real(kind_phys) :: icinc(ncol,pver)             ! in-cloud ice number concentration (m-3)
      real(kind_phys) :: icwnc(ncol,pver)             ! in-cloud water number concentration (m-3)

      ! Size distribution parameters (PUMAS kind, as these are passed to the PUMAS routines)
      real(pumas_r8)  :: rho(ncol,pver)               ! air density (kg m-3)
      real(pumas_r8)  :: icwmrst_r8(ncol,pver)        ! in-stratus water mixing ratio (kg kg-1)
      real(pumas_r8)  :: ncic(ncol,pver)              ! in-cloud droplet number concentration (kg-1)
      real(pumas_r8)  :: mu(ncol,pver)                ! droplet size distribution shape parameter (1)
      real(pumas_r8)  :: lambdac(ncol,pver)           ! droplet size distribution slope parameter (m-1)

      real(kind_phys) :: rel_fn(ncol,pver)            ! liquid effective radius at fixed number (um)
      real(kind_phys) :: drout2(ncol,pver)            ! mean rain particle diameter (m)
      real(kind_phys) :: dsout2(ncol,pver)            ! mean snow particle diameter (m)
      real(kind_phys) :: dgout2(ncol,pver)            ! mean graupel particle diameter (m)
      real(kind_phys) :: reff_rain(ncol,pver)         ! rain effective radius (um)
      real(kind_phys) :: reff_snow(ncol,pver)         ! snow effective radius (um)
      real(kind_phys) :: reff_grau(ncol,pver)         ! graupel effective radius (um)

      ! Vertical average of accretion, autoconversion and ratio
      real(kind_phys) :: vprao(ncol)
      real(kind_phys) :: vprco(ncol)
      real(kind_phys) :: racau(ncol)
      integer         :: cnt(ncol)

      real(kind_phys) :: zen_angle(ncol)              ! daytime solar zenith angles (radians)
      real(kind_phys) :: sza(ncol)                    ! solar zenith angle (degrees)
      real(kind_phys) :: rad2deg                      ! radians to degrees conversion factor
      real(kind_phys) :: delta                        ! solar declination angle (radians)
      real(kind_phys) :: eccf                         ! Earth orbit eccentricity factor

      ! History output fields
      real(kind_phys) :: cdnumc(ncol)
      real(kind_phys) :: efcout(ncol,pver)
      real(kind_phys) :: efiout(ncol,pver)
      real(kind_phys) :: ncout(ncol,pver)
      real(kind_phys) :: niout(ncol,pver)
      real(kind_phys) :: freql(ncol,pver)
      real(kind_phys) :: freqi(ncol,pver)
      real(kind_phys) :: freqm(ncol,pver)
      real(kind_phys) :: freqsl(ncol,pver)
      real(kind_phys) :: freqslm(ncol,pver)
      real(kind_phys) :: ctrel(ncol)
      real(kind_phys) :: ctrei(ncol)
      real(kind_phys) :: ctnl(ncol)
      real(kind_phys) :: ctni(ncol)
      real(kind_phys) :: fcti(ncol)
      real(kind_phys) :: fctl(ncol)
      real(kind_phys) :: fctm(ncol)
      real(kind_phys) :: fctsl(ncol)
      real(kind_phys) :: fctslm(ncol)

      ! Macrophysics dependent diagnostics (CME, CLDFSNOW, CLDFGRAU, precipitation efficiency and rainbow)
      real(kind_phys), parameter :: minlwp = 0.01_kind_phys ! minimum lwp threshold (kg/m3)

      real(kind_phys) :: qme(ncol,pver)               ! net microphysics condensation rate (kg kg-1 s-1)
      real(kind_phys) :: iclwpst(ncol,pver)           ! in-stratus liquid water path (kg m-2)
      real(kind_phys) :: cldfsnow(ncol,pver)          ! cloud fraction adjusted for snow (1)
      real(kind_phys) :: cldfgrau(ncol,pver)          ! cloud fraction adjusted for graupel (1)

      ! Precipitation efficiency
      real(kind_phys) :: tgliqwp(ncol)                ! column liquid water path (kg m-2)
      real(kind_phys) :: tgcmeliq(ncol)               ! column condensation rate (m s-1)
      real(kind_phys) :: pe(ncol)                     ! precipitation efficiency (1)
      real(kind_phys) :: tpr(ncol)                    ! average precip rate over efficiency calculation (m s-1)
      real(kind_phys) :: pefrac(ncol)                 ! fraction of timesteps precip efficiency reported (1)

      ! Rainbow diagnostic
      real(kind_phys) :: rbfrac(ncol)                 ! fraction of sky covered by a potential rainbow
      real(kind_phys) :: rbfreq(ncol)                 ! potential rainbow frequency
      real(kind_phys) :: precc(ncol)                  ! convective precip rate (m s-1)
      real(kind_phys) :: cldtot, cldmx, convmx, frlow, rmax
      logical         :: rval
      integer         :: top_idx

      integer :: i, k, icol

      errmsg = ''
      errflg = 0

      if (col_end - col_start + 1 /= ncol) then
         errflg = 1
         errmsg = 'micro_pumas_diagnostics_run: column indices are inconsistent with the number of columns'
         return
      end if

      ! Merge Precipitation rates (multi-process)
      precc(:) = prec_dp(:ncol) + prec_sh(:ncol)

      ! The air temperature after the microphysics heating has been applied
      ! (this was state_loc%t in CAM)
      !
      ! FOR_REVIEWER_FEEDBACK: CAM used both temperatures: the temperature BEFORE the microphysics
      !    heating (state%t) for the air density in the size distribution calculations, and the
      !    temperature AFTER the heating (state_loc%t) for the in-cloud number concentrations and the
      !    supercooled liquid (tmelt) tests.  This scheme gets the pre-heating temperature, so it
      !    must be placed after micro_pumas_ccpp_optics_limiter and BEFORE apply_heating_rate in the
      !    SDF (where pumas_post_main also requires the pre-heating temperature), and the post-heating
      !    temperature is recomputed here the same way apply_heating_rate computes it.  If this scheme
      !    is placed after apply_heating_rate, the results will be wrong without any error.  Options:
      !     1. Keep this approach (documented in the SDF files)
      !     2. Place this scheme after apply_heating_rate and back out the pre-heating temperature
      !        instead (equally placement dependent)
      !     3. Have an earlier scheme save the pre-heating air temperature as its own variable, so
      !        that this scheme can be placed anywhere after apply_heating_rate
      !    Which should be used?
      !    With the cold point tropopause adjustment (micro_mg_adjust_cpt), CAM's state_loc%t also
      !    included the adjustment of the cold point tropopause level (see
      !    micro_pumas_pre_diagnostics_run), so it is added here before the heating.
      t_post(:,:) = temp(:ncol,:)
      if (micro_mg_adjust_cpt) then
         do i = 1, ncol
            if (tropLev_cpp(i) /= NOTFOUND) then
               cp_dt = tropT_cpp(i) - temp(i,tropLev_cpp(i))
               t_post(i,tropLev_cpp(i)) = t_post(i,tropLev_cpp(i)) + cp_dt
            end if
         end do
      end if
      t_post(:,:) = t_post(:,:) + tlat(:ncol,:) * dtime / cpairv(:ncol,:)

      !-----------------------------------------------------------------------
      ! Get the solar zenith angle in degrees (also used for the rainbow diagnostics below)
      ! (this is CAM's call zenith(calday, rlats, rlons, zen_angle, ncol), with the zenith
      !  routine from CAM's orbit.F90 written out)
      !
      ! FOR_REVIEWER_FEEDBACK: The CAM-SIMA host already calculates this angle every timestep
      !    (orbital_data_advance in src/physics/utils/orbital_data.F90 sets the host variable
      !    solar_zenith_angle (rad) to acos(shr_orb_cosz(calday, lat, lon, solar_declination)) with
      !    the same calendar day and orbital parameters).  We suggest replacing the calculation below,
      !    and the calday, lat, lon, eccen, obliqr, lambm0 and mvelpp arguments, with that host
      !    variable (standard name solar_zenith_angle), i.e. sza(:) = solar_zenith_angle(:ncol) * rad2deg.
      !    This removes the duplicate orbital calculation and the shr_orb_mod dependency.  (The
      !    cosine_of_solar_zenith_angle_for_radiation host variable should NOT be used: it can be
      !    averaged over the radiation timestep, and it is not computed in snapshot runs.)
      !    Should the host variable be used?
      !-----------------------------------------------------------------------
      rad2deg = 180._kind_phys/pi

      zen_angle(:) = 0.0_kind_phys
      call shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, delta, eccf)
      do i = 1, ncol
         zen_angle(i) = shr_orb_cosz(calday, lat(i), lon(i), delta)
      end do
      where (zen_angle(:) <= 1.0_kind_phys .and. zen_angle(:) >= -1.0_kind_phys)
         zen_angle(:) = acos(zen_angle(:))
      elsewhere
         zen_angle(:) = 0.0_kind_phys
      end where

      sza(:) = zen_angle(:) * rad2deg
      call history_out_field('rbSZA', sza)

      ! Zero out diagnostic rainbow arrays
      rbfreq(:) = 0._kind_phys
      rbfrac(:) = 0._kind_phys

      !-----------------------------------------------------------------------
      ! Heterogeneous nucleation from vapor
      !-----------------------------------------------------------------------
      mnuccdohet(:,:) = 0._kind_phys
      do k = top_lev, pver
         do i = 1, ncol
            if (naai(i,k) > 0._kind_phys) then
               mnuccdohet(i,k) = real(proc_rates%mnuccdtot(i,k-top_lev+1), kind_phys) - &
                    (naai_hom(i,k)/naai(i,k)) * real(proc_rates%mnuccdtot(i,k-top_lev+1), kind_phys)
            end if
         end do
      end do

      !-----------------------------------------------------------------------
      ! Stratiform precipitation fluxes (MG2 and MG3 add the condensate and graupel fluxes)
      !-----------------------------------------------------------------------
      mgflxprc(:,:) = 0._kind_phys
      mgflxsnw(:,:) = 0._kind_phys
      mgflxprc(:,top_lev:) = rflx(:ncol,top_lev:) + sflx(:ncol,top_lev:) + iflx(:ncol,top_lev:) + &
                             lflx(:ncol,top_lev:) + gflx(:ncol,top_lev:)
      mgflxsnw(:,top_lev:) = sflx(:ncol,top_lev:) + iflx(:ncol,top_lev:) + gflx(:ncol,top_lev:)

      ! Effective radius of convective liquid and ice
      ! (hard-coded as average of hard-coded values used for deep/shallow convective detrainment)
      cvreffliq(:,:) = 0._kind_phys
      cvreffice(:,:) = 0._kind_phys
      cvreffliq(:,top_lev:) = 9.0_kind_phys
      cvreffice(:,top_lev:) = 37.0_kind_phys

      ! Net microphysics condensation rate
      qme(:,:) = 0._kind_phys
      qme(:,top_lev:) = cmeliq(:ncol,top_lev:) + real(proc_rates%cmeitot(:ncol,:), kind_phys)

      ! ------------------------------------------------------------ !
      ! Compute in cloud ice and liquid mixing ratios                !
      ! PUMAS uses the stratiform cloud fraction for both liquid and !
      ! ice cloud (liqcldf = icecldf = ast)                          !
      ! ------------------------------------------------------------ !
      icimrst(:,:)  = 0._kind_phys
      icwmrst(:,:)  = 0._kind_phys
      icinc(:,:)    = 0._kind_phys
      icwnc(:,:)    = 0._kind_phys
      iclwpst(:,:)  = 0._kind_phys
      cldfsnow(:,:) = 0._kind_phys
      cldfgrau(:,:) = 0._kind_phys

      do k = top_lev, pver
         do i = 1, ncol
            ! Limits for in-cloud mixing ratios consistent with MG microphysics
            ! in-cloud mixing ratio maximum limit of 0.005 kg/kg
            icimrst(i,k) = min(cldice(i,k) / max(real(mincld, kind_phys), ast(i,k)), 0.005_kind_phys)
            icwmrst(i,k) = min(cldliq(i,k) / max(real(mincld, kind_phys), ast(i,k)), 0.005_kind_phys)
            icinc(i,k)   = numice(i,k) / max(real(mincld, kind_phys), ast(i,k)) * &
                 pmid(i,k) / (rair*t_post(i,k))
            icwnc(i,k)   = numliq(i,k) / max(real(mincld, kind_phys), ast(i,k)) * &
                 pmid(i,k) / (rair*t_post(i,k))
            ! In-stratus liquid water path
            ! Note: uses stratiform cloud fraction!
            iclwpst(i,k) = min(cldliq(i,k) / max(real(mincld, kind_phys), ast(i,k)), 0.005_kind_phys) * &
                 pdel(i,k) / gravit

            ! ------------------------------ !
            ! Adjust cloud fraction for snow !
            ! ------------------------------ !
            cldfsnow(i,k) = cld(i,k)
            ! If cloud and only ice ( no convective cloud or ice ), then set to 0.
            if ((cldfsnow(i,k) > 1.e-4_kind_phys) .and. &
                (concld(i,k)   < 1.e-4_kind_phys) .and. &
                (cldliq(i,k)   < 1.e-10_kind_phys)) then
               cldfsnow(i,k) = 0._kind_phys
            end if
            ! If no cloud and snow, then set to 0.25
            if ((cldfsnow(i,k) <= 1.e-4_kind_phys) .and. (qsout(i,k) > 1.e-6_kind_phys)) then
               cldfsnow(i,k) = 0.25_kind_phys
            end if

            ! --------------------------------- !
            ! Adjust cloud fraction for graupel !
            ! --------------------------------- !
            cldfgrau(i,k) = cld(i,k)
            ! If cloud and only ice ( no convective cloud or ice ), then set to 0.
            if ((cldfgrau(i,k) > 1.e-4_kind_phys) .and. &
                (concld(i,k)   < 1.e-4_kind_phys) .and. &
                (cldliq(i,k)   < 1.e-10_kind_phys)) then
               cldfgrau(i,k) = 0._kind_phys
            end if
            ! If no cloud and graupel, then set to 0.25
            if ((cldfgrau(i,k) <= 1.e-4_kind_phys) .and. (qgout(i,k) > 1.e-9_kind_phys)) then
               cldfgrau(i,k) = 0.25_kind_phys
            end if
         end do
      end do

      ! ------------------------------------- !
      ! Size distribution calculation         !
      ! ------------------------------------- !

      ! Air density uses the temperature BEFORE the microphysics heating is applied
      ! (CAM used "state instead of state_loc to preserve answers for MG1")
      rho(:,:) = 0._pumas_r8
      rho(:,top_lev:) = real(pmid(:ncol,top_lev:), pumas_r8) / &
                        (real(rair, pumas_r8) * real(temp(:ncol,top_lev:), pumas_r8))
      icwmrst_r8(:,:) = real(icwmrst(:,:), pumas_r8)

      ! Effective radius for cloud liquid, fixed number.
      mu(:,:)      = 0._pumas_r8
      lambdac(:,:) = 0._pumas_r8
      rel_fn(:,:)  = 10._kind_phys
      ncic(:,:)    = 1.e8_pumas_r8

      do k = top_lev, pver
         call size_dist_param_liq(mg_liq_props, icwmrst_r8(:,k), ncic(:,k), rho(:,k), &
                                  mu(:,k), lambdac(:,k), ncol)
      end do

      where (icwmrst_r8(:,top_lev:) > qsmall)
         rel_fn(:,top_lev:) = real((mu(:,top_lev:) + 3._pumas_r8) / lambdac(:,top_lev:) / &
                                   2._pumas_r8 * 1.e6_pumas_r8, kind_phys)
      end where

      ! (The cloud liquid effective radius, rel, is recomputed in pumas_post_main.)

      ! Rain/snow/graupel effective diameter from the updated precipitation
      drout2(:,:)    = 0._kind_phys
      reff_rain(:,:) = 0._kind_phys
      dsout2(:,:)    = 0._kind_phys
      reff_snow(:,:) = 0._kind_phys
      reff_grau(:,:) = 0._kind_phys
      dgout2(:,:)    = 0._kind_phys

      where (qr(:ncol,top_lev:) >= 1.e-7_kind_phys)
         drout2(:,top_lev:) = real(avg_diameter( &
              real(qr(:ncol,top_lev:), pumas_r8), &
              real(nr(:ncol,top_lev:), pumas_r8) * rho(:,top_lev:), &
              rho(:,top_lev:), rhow), kind_phys)
         ! FOR_REVIEWER_FEEDBACK: 1.5 is the correct value (PUMAS uses 1.5 in micro_pumas_v1.F90, and
         !    CAM uses 1.5 for graupel below).  CAM's micro_pumas_cam.F90 uses shapeparam = 1.e5_r8
         !    here, which makes CAM's LS_REFFRAIN (and the mgreffrain physics buffer field used by COSP)
         !    about 66,667 times larger.  This appears to be a bug in CAM, so LS_REFFRAIN will not
         !    match CAM.  Is 1.5 correct?
         reff_rain(:,top_lev:) = drout2(:,top_lev:) * 1.5_kind_phys * 1.e6_kind_phys
      end where

      where (qs(:ncol,top_lev:) >= 1.e-7_kind_phys)
         dsout2(:,top_lev:) = real(avg_diameter( &
              real(qs(:ncol,top_lev:), pumas_r8), &
              real(ns(:ncol,top_lev:), pumas_r8) * rho(:,top_lev:), &
              rho(:,top_lev:), rhosn), kind_phys)
         ! FOR_REVIEWER_FEEDBACK: 1.5 is the correct value; CAM uses shapeparam = 1.e5_r8 here (see
         !    the comment for reff_rain above), so LS_REFFSNOW will not match CAM.
         reff_snow(:,top_lev:) = dsout2(:,top_lev:) * 1.5_kind_phys * 1.e6_kind_phys
      end where

      where (qg(:ncol,top_lev:) >= 1.e-7_kind_phys)
         dgout2(:,top_lev:) = real(avg_diameter( &
              real(qg(:ncol,top_lev:), pumas_r8), &
              real(ng(:ncol,top_lev:), pumas_r8) * rho(:,top_lev:), &
              rho(:,top_lev:), rhog), kind_phys)
         reff_grau(:,top_lev:) = dgout2(:,top_lev:) * 1.5_kind_phys * 1.e6_kind_phys
      end where

      ! ------------------------------------- !
      ! Precipitation efficiency Calculation  !
      ! ------------------------------------- !

      ! Compute liquid water paths, and column condensation
      tgliqwp(:)  = 0._kind_phys
      tgcmeliq(:) = 0._kind_phys
      do k = top_lev, pver
         do i = 1, ncol
            tgliqwp(i) = tgliqwp(i) + iclwpst(i,k)*cld(i,k)

            if (cmeliq(i,k) > 1.e-12_kind_phys) then
               !convert cmeliq to right units:  kgh2o/kgair/s  *  kgair/m2  / kgh2o/m3  = m/s
               tgcmeliq(i) = tgcmeliq(i) + cmeliq(i,k) * (pdel(i,k) / gravit) / rhoh2o
            end if
         end do
      end do

      ! zero out precip efficiency and total averaged precip
      pe(:)     = 0._kind_phys
      tpr(:)    = 0._kind_phys
      pefrac(:) = 0._kind_phys

      ! accumulate precip and condensation
      ! (the stratiform precip is the PUMAS surface precip as there is no sedimentation precip)
      ! FOR_REVIEWER_FEEDBACK: When restarts are implemented, the accumulators will be passed in as
      !    arguments on the horizontal_loop_extent dimension, so icol will be replaced by i (see the
      !    comment at the acprecl, acgcme and acnum declarations at the top of this module).
      do i = 1, ncol
         icol = col_start + i - 1

         acgcme(icol)  = acgcme(icol) + tgcmeliq(i)
         acprecl(icol) = acprecl(icol) + prect(i)
         acnum(icol)   = acnum(icol) + 1

         ! if LWP is zero, then 'end of cloud': calculate precip efficiency
         if (tgliqwp(i) < minlwp) then
            if (acprecl(icol) > 5.e-8_kind_phys) then
               tpr(i) = max(acprecl(icol)/acnum(icol), 1.e-15_kind_phys)
               if (acgcme(icol) > 1.e-10_kind_phys) then
                  pe(i) = min(max(acprecl(icol)/acgcme(icol), 1.e-15_kind_phys), 1.e5_kind_phys)
                  pefrac(i) = 1._kind_phys
               end if
            end if

            ! reset counters
            acprecl(icol) = 0._kind_phys
            acgcme(icol)  = 0._kind_phys
            acnum(icol)   = 0
         end if

         ! if never find any rain....(after 10^3 timesteps...)
         if (acnum(icol) > 1000) then
            acnum(icol)   = 0
            acprecl(icol) = 0._kind_phys
            acgcme(icol)  = 0._kind_phys
         end if
      end do

      !-----------------------------------------------------------------------
      ! vertical average of non-zero accretion, autoconversion and ratio.
      !-----------------------------------------------------------------------
      vprao(:) = 0._kind_phys
      cnt(:) = 0
      do k = 1, pver-top_lev+1
         vprao(:) = vprao(:) + real(proc_rates%pratot(:ncol,k), kind_phys)
         where (proc_rates%pratot(:ncol,k) /= 0._pumas_r8) cnt(:) = cnt(:) + 1
      end do

      where (cnt(:) > 0) vprao(:) = vprao(:)/cnt(:)

      vprco(:) = 0._kind_phys
      cnt(:) = 0
      do k = 1, pver-top_lev+1
         vprco(:) = vprco(:) + real(proc_rates%prctot(:ncol,k), kind_phys)
         where (proc_rates%prctot(:ncol,k) /= 0._pumas_r8) cnt(:) = cnt(:) + 1
      end do

      where (cnt(:) > 0)
         vprco(:) = vprco(:)/cnt(:)
         racau(:) = vprao(:)/vprco(:)
      elsewhere
         racau(:) = 0._kind_phys
      end where

      racau(:) = min(racau(:), 1.e10_kind_phys)

      !-----------------------------------------------------------------------
      ! Diagnostic Rainbow Calculation. Seriously.
      !-----------------------------------------------------------------------

      ! (the solar zenith angle, sza, was calculated above)

      do i = 1, ncol

         top_idx = pver
         convmx  = 0._kind_phys
         frlow   = 0._kind_phys
         cldmx   = 0._kind_phys
         cldtot  = maxval(ast(i,top_lev:))

         ! Find levels in surface layer
         do k = top_lev, pver
            if (pmid(i,k) > rb_pmin) then
               top_idx = min(k,top_idx)
            end if
         end do

         ! For all fractional precip calculated below, use maximum in surface layer.
         ! For convective precip, base on convective cloud area
         convmx = maxval(concld(i,top_idx:))
         ! For stratiform precip, base on precip fraction
         cldmx = maxval(freqr(i,top_idx:))
         ! Combine and use maximum of strat or conv fraction
         frlow = max(cldmx,convmx)

         ! max precip
         rmax = maxval(qrout(i,top_idx:))

         ! Stratiform precip mixing ratio OR some convective precip
         ! (rval = true  if any sig precip)
         rval = ((precc(i) > rb_rcmin) .or. (rmax > rb_rmin))

         ! Now can find conditions for a rainbow:
         ! Maximum cloud cover (CLDTOT) < 0.5
         ! 48 < SZA < 90
         ! freqr (below rb_pmin) > 0.25
         ! Some rain (liquid > 1.e-6 kg/kg, convective precip > 1.e-7 m/s
         if ((cldtot < 0.5_kind_phys) .and. (sza(i) > 48._kind_phys) .and. (sza(i) < 90._kind_phys) .and. rval) then

            ! Rainbow 'probability' (area) derived from solid angle theory
            ! as the fraction of the hemisphere for a spherical cap with angle phi=sza-48.
            ! This is only valid between 48 < sza < 90 (controlled for above).
            rbfrac(i) = max(0._kind_phys, (1._kind_phys-cos((sza(i)-48._kind_phys)*(pi/180._kind_phys)))/2._kind_phys) &
                        * frlow
            rbfreq(i) = 1.0_kind_phys
         end if

      end do

      call history_out_field('RBFRAC', rbfrac)
      call history_out_field('RBFREQ', rbfreq)

      ! --------------------- !
      ! History Output Fields !
      ! --------------------- !

      ! Column droplet concentration
      cdnumc(:) = sum(numliq(:ncol,top_lev:pver) * pdel(:ncol,top_lev:pver)/gravit, dim=2)

      ! Averaging for new output fields
      efcout(:,:)      = 0._kind_phys
      efiout(:,:)      = 0._kind_phys
      ncout(:,:)       = 0._kind_phys
      niout(:,:)       = 0._kind_phys
      freql(:,:)       = 0._kind_phys
      freqi(:,:)       = 0._kind_phys
      icwmrst_out(:,:) = 0._kind_phys
      icimrst_out(:,:) = 0._kind_phys
      freqm(:,:)       = 0._kind_phys
      freqsl(:,:)      = 0._kind_phys
      freqslm(:,:)     = 0._kind_phys

      do k = top_lev, pver
         do i = 1, ncol
            if (ast(i,k) > 0.01_kind_phys .and. icwmrst(i,k) > 5.e-5_kind_phys) then
               efcout(i,k) = rel(i,k) * ast(i,k)
               ncout(i,k)  = icwnc(i,k) * ast(i,k)
               freql(i,k)  = ast(i,k)
               icwmrst_out(i,k) = icwmrst(i,k)
            end if
            if (ast(i,k) > 0.01_kind_phys .and. icimrst(i,k) > 1.e-6_kind_phys) then
               efiout(i,k) = effi(i,k) * ast(i,k)
               niout(i,k)  = icinc(i,k) * ast(i,k)
               freqi(i,k)  = ast(i,k)
               icimrst_out(i,k) = icimrst(i,k)
            end if

            ! Supercooled liquid
            if (freql(i,k) > 0.01_kind_phys .and. freqi(i,k) > 0.01_kind_phys) then
               freqm(i,k) = ast(i,k)
            end if
            if (freql(i,k) > 0.01_kind_phys .and. freqi(i,k) < 0.01_kind_phys .and. t_post(i,k) < tmelt) then
               freqsl(i,k) = ast(i,k)
            end if
            if (freql(i,k) > 0.01_kind_phys .and. freqi(i,k) > 0.01_kind_phys .and. t_post(i,k) < tmelt) then
               freqslm(i,k) = ast(i,k)
            end if
         end do
      end do

      ! Cloud top effective radius and number.
      fcti(:)   = 0._kind_phys
      fctl(:)   = 0._kind_phys
      ctrel(:)  = 0._kind_phys
      ctrei(:)  = 0._kind_phys
      ctnl(:)   = 0._kind_phys
      ctni(:)   = 0._kind_phys
      fctm(:)   = 0._kind_phys
      fctsl(:)  = 0._kind_phys
      fctslm(:) = 0._kind_phys

      do i = 1, ncol
         do k = top_lev, pver
            if (ast(i,k) > 0.01_kind_phys .and. icwmrst(i,k) > 1.e-7_kind_phys) then
               ctrel(i) = rel(i,k) * ast(i,k)
               ctnl(i)  = icwnc(i,k) * ast(i,k)
               fctl(i)  = ast(i,k)

               ! Cloud Top Mixed phase, supercooled liquid only and supercooled liquid mixed
               if (freqi(i,k) > 0.01_kind_phys) then
                  fctm(i) = ast(i,k)
               end if
               if (freqi(i,k) < 0.01_kind_phys .and. t_post(i,k) < tmelt) then
                  fctsl(i) = ast(i,k)
               end if
               if (freqi(i,k) > 0.01_kind_phys .and. t_post(i,k) < tmelt) then
                  fctslm(i) = ast(i,k)
               end if

               exit
            end if

            if (ast(i,k) > 0.01_kind_phys .and. icimrst(i,k) > 1.e-7_kind_phys) then
               ctrei(i) = effi(i,k) * ast(i,k)
               ctni(i)  = icinc(i,k) * ast(i,k)
               fcti(i)  = ast(i,k)
               exit
            end if
         end do
      end do

      ! --------------------------------------------- !
      ! General outfield calls for microphysics       !
      ! --------------------------------------------- !

      ! Output a handle of variables which are calculated on the fly
      ftem(:,:) = proc_rates%qcrestot(:ncol,:)
      call out_pumas_lev_field('MPDW2V', ftem, top_lev, work)

      ftem(:,:) = proc_rates%melttot(:ncol,:) - proc_rates%mnuccctot(:ncol,:) - proc_rates%mnuccttot(:ncol,:) &
                - proc_rates%bergtot(:ncol,:) - proc_rates%homotot(:ncol,:) - proc_rates%msacwitot(:ncol,:)
      call out_pumas_lev_field('MPDW2I', ftem, top_lev, work)

      ! The graupel terms are zero if graupel is not turned on
      ftem(:,:) = -proc_rates%pratot(:ncol,:) - proc_rates%prctot(:ncol,:) - proc_rates%psacwstot(:ncol,:) &
                - proc_rates%bergstot(:ncol,:) - proc_rates%psacwgtot(:ncol,:) - proc_rates%pgsacwtot(:ncol,:)
      call out_pumas_lev_field('MPDW2P', ftem, top_lev, work)

      ftem(:,:) = proc_rates%cmeitot(:ncol,:) + proc_rates%qirestot(:ncol,:)
      call out_pumas_lev_field('MPDI2V', ftem, top_lev, work)

      ftem(:,:) = -proc_rates%melttot(:ncol,:) + proc_rates%mnuccctot(:ncol,:) + proc_rates%mnuccttot(:ncol,:) &
                + proc_rates%bergtot(:ncol,:) + proc_rates%homotot(:ncol,:) + proc_rates%msacwitot(:ncol,:)  &
                - proc_rates%qmultgtot(:ncol,:)
      call out_pumas_lev_field('MPDI2W', ftem, top_lev, work)

      ftem(:,:) = -proc_rates%prcitot(:ncol,:) - proc_rates%praitot(:ncol,:)
      call out_pumas_lev_field('MPDI2P', ftem, top_lev, work)

      ! TAU / emulated stochastic collection (machine learning) warm rain diagnostics
      ! (the binned fields need a 3-D version of out_pumas_trop_cld_field)
      if (trim(micro_mg_warm_rain) == 'tau' .or. trim(micro_mg_warm_rain) == 'emulated') then
         call out_pumas_trop_cld_field('scale_qc',        proc_rates%scale_qc,        ncol)
         call out_pumas_trop_cld_field('scale_nc',        proc_rates%scale_nc,        ncol)
         call out_pumas_trop_cld_field('scale_qr',        proc_rates%scale_qr,        ncol)
         call out_pumas_trop_cld_field('scale_nr',        proc_rates%scale_nr,        ncol)
!        call out_pumas_trop_cld_field('amk_c',           proc_rates%amk_c,           ncol) !UNCOMMENT WHEN 3D IMPLEMENTED
!        call out_pumas_trop_cld_field('ank_c',           proc_rates%ank_c,           ncol) !UNCOMMENT WHEN 3D IMPLEMENTED
!        call out_pumas_trop_cld_field('amk_r',           proc_rates%amk_r,           ncol) !UNCOMMENT WHEN 3D IMPLEMENTED
!        call out_pumas_trop_cld_field('ank_r',           proc_rates%ank_r,           ncol) !UNCOMMENT WHEN 3D IMPLEMENTED
!        call out_pumas_trop_cld_field('amk',             proc_rates%amk,             ncol) !UNCOMMENT WHEN 3D IMPLEMENTED
!        call out_pumas_trop_cld_field('ank',             proc_rates%ank,             ncol) !UNCOMMENT WHEN 3D IMPLEMENTED
!        call out_pumas_trop_cld_field('amk_out',         proc_rates%amk_out,         ncol) !UNCOMMENT WHEN 3D IMPLEMENTED
!        call out_pumas_trop_cld_field('ank_out',         proc_rates%ank_out,         ncol) !UNCOMMENT WHEN 3D IMPLEMENTED
         call out_pumas_trop_cld_field('QC_TAU_out',      proc_rates%qc_out_TAU,      ncol)
         call out_pumas_trop_cld_field('NC_TAU_out',      proc_rates%nc_out_TAU,      ncol)
         call out_pumas_trop_cld_field('QR_TAU_out',      proc_rates%qr_out_TAU,      ncol)
         call out_pumas_trop_cld_field('NR_TAU_out',      proc_rates%nr_out_TAU,      ncol)
         call out_pumas_trop_cld_field('qctend_TAU',      proc_rates%qctend_TAU,      ncol)
         call out_pumas_trop_cld_field('nctend_TAU',      proc_rates%nctend_TAU,      ncol)
         call out_pumas_trop_cld_field('qrtend_TAU',      proc_rates%qrtend_TAU,      ncol)
         call out_pumas_trop_cld_field('nrtend_TAU',      proc_rates%nrtend_TAU,      ncol)
         call out_pumas_trop_cld_field('gmnnn_lmnnn_TAU', proc_rates%gmnnn_lmnnn_TAU, ncol)
         call out_pumas_trop_cld_field('ML_fixer',        proc_rates%ML_fixer,        ncol)
         call out_pumas_trop_cld_field('qc_fixer',        proc_rates%qc_fixer,        ncol)
         call out_pumas_trop_cld_field('nc_fixer',        proc_rates%nc_fixer,        ncol)
         call out_pumas_trop_cld_field('qr_fixer',        proc_rates%qr_fixer,        ncol)
         call out_pumas_trop_cld_field('nr_fixer',        proc_rates%nr_fixer,        ncol)
         call out_pumas_trop_cld_field('QC_TAU_in',       proc_rates%qc_in_TAU,       ncol)
         call out_pumas_trop_cld_field('NC_TAU_in',       proc_rates%nc_in_TAU,       ncol)
         call out_pumas_trop_cld_field('QR_TAU_in',       proc_rates%qr_in_TAU,       ncol)
         call out_pumas_trop_cld_field('NR_TAU_in',       proc_rates%nr_in_TAU,       ncol)
      end if
      if (trim(micro_mg_warm_rain) == 'sb2001') then
         call out_pumas_trop_cld_field('qctend_SB2001', proc_rates%qctend_SB2001, ncol)
         call out_pumas_trop_cld_field('nctend_SB2001', proc_rates%nctend_SB2001, ncol)
         call out_pumas_trop_cld_field('qrtend_SB2001', proc_rates%qrtend_SB2001, ncol)
         call out_pumas_trop_cld_field('nrtend_SB2001', proc_rates%nrtend_SB2001, ncol)
      end if
      if (trim(micro_mg_warm_rain) == 'kk2000') then
         call out_pumas_trop_cld_field('qctend_KK2000', proc_rates%qctend_KK2000, ncol)
         call out_pumas_trop_cld_field('nctend_KK2000', proc_rates%nctend_KK2000, ncol)
         call out_pumas_trop_cld_field('qrtend_KK2000', proc_rates%qrtend_KK2000, ncol)
         call out_pumas_trop_cld_field('nrtend_KK2000', proc_rates%nrtend_KK2000, ncol)
      end if

      call out_pumas_trop_cld_field('LAMC', proc_rates%lamc_out, ncol)
      call out_pumas_trop_cld_field('LAMR', proc_rates%lamr_out, ncol)
      call out_pumas_trop_cld_field('PGAM', proc_rates%pgam_out, ncol)
      call out_pumas_trop_cld_field('N0R',  proc_rates%n0r_out,  ncol)

      call history_out_field('REFL',      refl)
      call history_out_field('AREFL',     arefl)
      call history_out_field('AREFLZ',    areflz)
      call history_out_field('FREFL',     frefl)
      call history_out_field('CSRFL',     csrfl)
      call history_out_field('ACSRFL',    acsrfl)
      call history_out_field('FCSRFL',    fcsrfl)
      ! CAM set the 10 cm reflectivity to -9999 above the top level for cloud physics
      refl10cm_out(:,:) = refl10cm(:ncol,:)
      refl10cm_out(:,:top_lev-1) = -9999._kind_phys
      call history_out_field('REFL10CM',  refl10cm_out)
      call history_out_field('REFLZ10CM', reflz10cm)
      call history_out_field('RERCLD',    rercld)
      call history_out_field('NCAL',      ncal)
      call history_out_field('NCAI',      ncai)
      call history_out_field('AQRAIN',    qrout2)
      call history_out_field('AQSNOW',    qsout2)
      call history_out_field('ANRAIN',    nrout2)
      call history_out_field('ANSNOW',    nsout2)
      call history_out_field('FREQR',     freqr)
      call history_out_field('FREQS',     freqs)
      call history_out_field('MPDT',      tlat)
      call out_pumas_lev_field('MPDQ',    qvlat, top_lev, work)
      call out_pumas_lev_field('MPDLIQ',  qcten, top_lev, work)
      call out_pumas_lev_field('MPDICE',  qiten, top_lev, work)
      call out_pumas_lev_field('MPDNLIQ', ncten, top_lev, work)
      call out_pumas_lev_field('MPDNICE', niten, top_lev, work)

      call out_pumas_trop_cld_field('EVAPSNOW', proc_rates%evapsnow,   ncol)
      call out_pumas_trop_cld_field('QCSEVAP',  proc_rates%qcsevap,    ncol)
      call out_pumas_trop_cld_field('QISEVAP',  proc_rates%qisevap,    ncol)
      call out_pumas_trop_cld_field('QVRES',    proc_rates%qvres,      ncol)
      call out_pumas_trop_cld_field('VTRMC',    proc_rates%vtrmc,      ncol)
      call out_pumas_trop_cld_field('VTRMI',    proc_rates%vtrmi,      ncol)
      call out_pumas_trop_cld_field('QCSEDTEN', proc_rates%qcsedten,   ncol)
      call out_pumas_trop_cld_field('QISEDTEN', proc_rates%qisedten,   ncol)
      call out_pumas_trop_cld_field('QRSEDTEN', proc_rates%qrsedten,   ncol)
      call out_pumas_trop_cld_field('QSSEDTEN', proc_rates%qssedten,   ncol)
      call out_pumas_trop_cld_field('MNUCCRIO', proc_rates%mnuccritot, ncol)
      call out_pumas_trop_cld_field('MNUDEPO',  proc_rates%mnudeptot,  ncol)
      call out_pumas_trop_cld_field('MELTSTOT', proc_rates%meltstot,   ncol)
      call out_pumas_trop_cld_field('MNUCCDO',  proc_rates%mnuccdtot,  ncol)
      call history_out_field('MNUCCDOhet', mnuccdohet)
      call out_pumas_trop_cld_field('MNUCCRO',  proc_rates%mnuccrtot,  ncol)
      call out_pumas_trop_cld_field('PRACSO',   proc_rates%pracstot,   ncol)
      call out_pumas_trop_cld_field('VAPDEPSO', proc_rates%vapdepstot, ncol)
      call out_pumas_trop_cld_field('MELTSDT',  proc_rates%meltsdttot, ncol)
      call out_pumas_trop_cld_field('FRZRDT',   proc_rates%frzrdttot,  ncol)
      call history_out_field('FICE',       nfice)
      call history_out_field('CLDFSNOW',   cldfsnow)

      call out_pumas_trop_cld_field('NNUCCCO',  proc_rates%nnuccctot,  ncol)
      call out_pumas_trop_cld_field('NNUCCTO',  proc_rates%nnuccttot,  ncol)
      call out_pumas_trop_cld_field('NNUCCDO',  proc_rates%nnuccdtot,  ncol)
      call out_pumas_trop_cld_field('NNUDEPO',  proc_rates%nnudeptot,  ncol)
      call out_pumas_trop_cld_field('NHOMO',    proc_rates%nhomotot,   ncol)
      call out_pumas_trop_cld_field('NNUCCRO',  proc_rates%nnuccrtot,  ncol)
      call out_pumas_trop_cld_field('NNUCCRIO', proc_rates%nnuccritot, ncol)
      call out_pumas_trop_cld_field('NSACWIO',  proc_rates%nsacwitot,  ncol)
      call out_pumas_trop_cld_field('NPRAO',    proc_rates%npratot,    ncol)
      call out_pumas_trop_cld_field('NPSACWSO', proc_rates%npsacwstot, ncol)
      call out_pumas_trop_cld_field('NPRAIO',   proc_rates%npraitot,   ncol)
      call out_pumas_trop_cld_field('NPRACSO',  proc_rates%npracstot,  ncol)
      call out_pumas_trop_cld_field('NPRCO',    proc_rates%nprctot,    ncol)
      call out_pumas_trop_cld_field('NPRCIO',   proc_rates%nprcitot,   ncol)
      call out_pumas_trop_cld_field('NCSEDTEN', proc_rates%ncsedten,   ncol)
      call out_pumas_trop_cld_field('NISEDTEN', proc_rates%nisedten,   ncol)
      call out_pumas_trop_cld_field('NRSEDTEN', proc_rates%nrsedten,   ncol)
      call out_pumas_trop_cld_field('NSSEDTEN', proc_rates%nssedten,   ncol)
      call out_pumas_trop_cld_field('NMELTO',   proc_rates%nmelttot,   ncol)
      call out_pumas_trop_cld_field('NMELTS',   proc_rates%nmeltstot,  ncol)

      call out_pumas_trop_cld_field('UMR',      proc_rates%umr,        ncol)
      call out_pumas_trop_cld_field('UMS',      proc_rates%ums,        ncol)

      call history_out_field('QCRAT',      qcrat)

      ! Graupel/hail fields
      call out_pumas_trop_cld_field('UMG',      proc_rates%umg,        ncol)
      call out_pumas_trop_cld_field('QGSEDTEN', proc_rates%qgsedten,   ncol)
      call history_out_field('FREQG',      freqg)
      call history_out_field('AQGRAU',     qgout2)
      call history_out_field('ANGRAU',     ngout2)
      call history_out_field('CLDFGRAU',   cldfgrau)
      call out_pumas_trop_cld_field('MELTGTOT', proc_rates%meltgtot,   ncol)
      call out_pumas_trop_cld_field('NMELTG',   proc_rates%nmeltgtot,  ncol)
      call out_pumas_trop_cld_field('NGSEDTEN', proc_rates%ngsedten,   ncol)

      ! Output fields which were on the grid in CAM
      call history_out_field('QRAIN',       qrout)
      call history_out_field('QSNOW',       qsout)
      call history_out_field('NRAIN',       nrout)
      call history_out_field('NSNOW',       nsout)
      call history_out_field('CV_REFFLIQ',  cvreffliq)
      call history_out_field('CV_REFFICE',  cvreffice)
      call history_out_field('LS_FLXPRC',   mgflxprc)
      call history_out_field('LS_FLXSNW',   mgflxsnw)
      call history_out_field('CME',         qme)
      call history_out_field('PRODPREC',    prain)
      call history_out_field('EVAPPREC',    nevapr)
      call out_pumas_lev_field('QCRESO',    proc_rates%qcrestot, top_lev, work)
      call history_out_field('LS_REFFRAIN', reff_rain)
      call history_out_field('LS_REFFSNOW', reff_snow)
      call history_out_field('DSNOW',       des)
      call history_out_field('ADRAIN',      drout2)
      call history_out_field('ADSNOW',      dsout2)
      call history_out_field('PE',          pe)
      call history_out_field('PEFRAC',      pefrac)
      call history_out_field('APRL',        tpr)
      call history_out_field('VPRAO',       vprao)
      call history_out_field('VPRCO',       vprco)
      call history_out_field('RACAU',       racau)
      call history_out_field('AREL',        efcout)
      call history_out_field('AREI',        efiout)
      call history_out_field('AWNC',        ncout)
      call history_out_field('AWNI',        niout)
      call history_out_field('FREQL',       freql)
      call history_out_field('FREQI',       freqi)
      call history_out_field('ACTREL',      ctrel)
      call history_out_field('ACTREI',      ctrei)
      call history_out_field('ACTNL',       ctnl)
      call history_out_field('ACTNI',       ctni)
      call history_out_field('FCTL',        fctl)
      call history_out_field('FCTI',        fcti)
      call history_out_field('ICINC',       icinc)
      call history_out_field('ICWNC',       icwnc)
      call history_out_field('EFFLIQ_IND',  rel_fn)
      call history_out_field('CDNUMC',      cdnumc)
      call history_out_field('REL',         rel)
      call history_out_field('REI',         effi)
      call history_out_field('MG_SADICE',   sadice)
      call history_out_field('MG_SADSNOW',  sadsnow)
      call history_out_field('ICIMRST',     icimrst_out)
      call history_out_field('ICWMRST',     icwmrst_out)
      call out_pumas_trop_cld_field('CMEIOUT',   proc_rates%cmeitot,   ncol)
      call out_pumas_lev_field('PRAO',      proc_rates%pratot,    top_lev, work)
      call out_pumas_lev_field('PRCO',      proc_rates%prctot,    top_lev, work)
      call out_pumas_lev_field('MNUCCCO',   proc_rates%mnuccctot, top_lev, work)
      call out_pumas_lev_field('MNUCCTO',   proc_rates%mnuccttot, top_lev, work)
      call out_pumas_lev_field('MSACWIO',   proc_rates%msacwitot, top_lev, work)
      call out_pumas_lev_field('PSACWSO',   proc_rates%psacwstot, top_lev, work)
      call out_pumas_lev_field('BERGSO',    proc_rates%bergstot,  top_lev, work)
      call out_pumas_lev_field('BERGO',     proc_rates%bergtot,   top_lev, work)
      call out_pumas_lev_field('MELTO',     proc_rates%melttot,   top_lev, work)
      call out_pumas_lev_field('HOMOO',     proc_rates%homotot,   top_lev, work)
      call out_pumas_lev_field('PRCIO',     proc_rates%prcitot,   top_lev, work)
      call out_pumas_lev_field('PRAIO',     proc_rates%praitot,   top_lev, work)
      call out_pumas_lev_field('QIRESO',    proc_rates%qirestot,  top_lev, work)
      call history_out_field('FREQM',       freqm)
      call history_out_field('FREQSL',      freqsl)
      call history_out_field('FREQSLM',     freqslm)
      call history_out_field('FCTM',        fctm)
      call history_out_field('FCTSL',       fctsl)
      call history_out_field('FCTSLM',      fctslm)

      call out_pumas_lev_field('PRACGO',    proc_rates%pracgtot,   top_lev, work)
      call out_pumas_lev_field('PSACRO',    proc_rates%psacrtot,   top_lev, work)
      call out_pumas_lev_field('PSACWGO',   proc_rates%psacwgtot,  top_lev, work)
      call out_pumas_lev_field('PGSACWO',   proc_rates%pgsacwtot,  top_lev, work)
      call out_pumas_lev_field('PGRACSO',   proc_rates%pgracstot,  top_lev, work)
      call out_pumas_lev_field('PRDGO',     proc_rates%prdgtot,    top_lev, work)
      call out_pumas_lev_field('QMULTGO',   proc_rates%qmultgtot,  top_lev, work)
      call out_pumas_lev_field('QMULTRGO',  proc_rates%qmultrgtot, top_lev, work)
      call history_out_field('LS_REFFGRAU', reff_grau)
      call out_pumas_lev_field('NPRACGO',   proc_rates%npracgtot,  top_lev, work)
      call out_pumas_lev_field('NSCNGO',    proc_rates%nscngtot,   top_lev, work)
      call out_pumas_lev_field('NGRACSO',   proc_rates%ngracstot,  top_lev, work)
      call out_pumas_lev_field('NMULTGO',   proc_rates%nmultgtot,  top_lev, work)
      call out_pumas_lev_field('NMULTRGO',  proc_rates%nmultrgtot, top_lev, work)
      call out_pumas_lev_field('NPSACWGO',  proc_rates%npsacwgtot, top_lev, work)

      if (micro_mg_adjust_cpt) then
         cp_rh(:ncol, :pver)  = 0._kind_phys

         do i = 1, ncol

            ! Calculate the RH including any T change that we make.
            do k = top_lev, pver
              call qsat(t_post(i,k), pmid(i,k), es, qsat_out)
              cp_rh(i,k) = q(i, k) / qsat_out * 100._kind_phys
            end do
         end do

         call history_out_field('TROPF_RHADJ', cp_rh)
      end if

   end subroutine micro_pumas_diagnostics_run

   !> \section arg_table_micro_pumas_diagnostics_final  Argument Table
   !! \htmlinclude micro_pumas_diagnostics_final.html
   subroutine micro_pumas_diagnostics_final(errmsg, errflg)

      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! FOR_REVIEWER_FEEDBACK: When restarts are implemented, this final phase can be removed
      !    (see the comment at the acprecl, acgcme and acnum declarations at the top of this module).
      if (allocated(acprecl)) deallocate(acprecl)
      if (allocated(acgcme))  deallocate(acgcme)
      if (allocated(acnum))   deallocate(acnum)

   end subroutine micro_pumas_diagnostics_final

   !=======================================================================

   ! Output a field which is on the PUMAS vertical levels (top_lev:pver) on the
   ! host model vertical levels ('lev'), with zeros above top_lev
   subroutine out_pumas_lev_field(fname, field, top_lev, work)

      use cam_history, only: history_out_field

      character(len=*), intent(in)  :: fname       ! history field name
      real(pumas_r8),   intent(in)  :: field(:,:)  ! field on the PUMAS vertical levels
      integer,          intent(in)  :: top_lev     ! top model level for cloud physics (index)
      real(kind_phys),  intent(out) :: work(:,:)   ! work array on the host model vertical levels

      work(:,:top_lev-1) = 0._kind_phys
      work(:,top_lev:)   = real(field(:size(work,1),:), kind_phys)

      call history_out_field(fname, work)

   end subroutine out_pumas_lev_field

   ! Output a field which is on the PUMAS vertical levels (top_lev:pver) on the
   ! 'trop_cld_lev' history coordinate
   subroutine out_pumas_trop_cld_field(fname, field, ncol)

      use cam_history, only: history_out_field

      character(len=*), intent(in)  :: fname       ! history field name
      real(pumas_r8),   intent(in)  :: field(:,:)  ! field on the PUMAS vertical levels
      integer,          intent(in)  :: ncol        ! number of columns

      call history_out_field(fname, real(field(:ncol,:), kind_phys))

   end subroutine out_pumas_trop_cld_field

end module micro_pumas_diagnostics
