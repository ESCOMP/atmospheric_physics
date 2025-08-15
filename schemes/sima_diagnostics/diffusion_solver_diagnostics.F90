! Diagnostic scheme for vertical diffusion solver
module diffusion_solver_diagnostics
  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  ! public CCPP-compliant subroutines
  public :: vertical_diffusion_tendencies_diagnostics_init
  public :: vertical_diffusion_tendencies_diagnostics_run

contains

!> \section arg_table_vertical_diffusion_tendencies_diagnostics_init Argument Table
!! \htmlinclude vertical_diffusion_tendencies_diagnostics_init.html
  subroutine vertical_diffusion_tendencies_diagnostics_init(const_props, errmsg, errflg)

      ! framework dependency for const_props
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

      ! dependency to get constituent index
      use ccpp_const_utils,          only: ccpp_const_get_idx

      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      type(ccpp_constituent_prop_ptr_t), intent(in) :: const_props(:)
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! History add field calls
      call history_add_field('QT',          'Total water mixing ratio',                           'lev',   'inst', 'kg/kg')
      call history_add_field('SL',          'Liquid water static energy',                         'lev',   'inst', 'J/kg')
      call history_add_field('SLV',         'Liquid water virtual static energy',                 'lev',   'inst', 'J/kg')
      call history_add_field('SLFLX',       'Liquid static energy flux',                          'ilev',  'inst', 'W/m2')
      call history_add_field('QTFLX',       'Total water flux',                                   'ilev',  'inst', 'kg/m2/s')
      call history_add_field('UFLX',        'Zonal momentum flux',                                'ilev',  'inst', 'N/m2')
      call history_add_field('VFLX',        'Meridional momentum flux',                           'ilev',  'inst', 'N/m2')

      ! Pre-PBL profiles
      call history_add_field('qt_pre_PBL',  'Total water mixing ratio before PBL',                'lev',   'inst', 'kg/kg')
      call history_add_field('sl_pre_PBL',  'Liquid water static energy before PBL',              'lev',   'inst', 'J/kg')
      call history_add_field('slv_pre_PBL', 'Liquid water virtual static energy before PBL',      'lev',   'inst', 'J/kg')
      call history_add_field('u_pre_PBL',   'Zonal wind before PBL',                              'lev',   'inst', 'm/s')
      call history_add_field('v_pre_PBL',   'Meridional wind before PBL',                         'lev',   'inst', 'm/s')
      call history_add_field('qv_pre_PBL',  'Water vapor mixing ratio before PBL',                'lev',   'inst', 'kg/kg')
      call history_add_field('ql_pre_PBL',  'Cloud liquid water mixing ratio before PBL',         'lev',   'inst', 'kg/kg')
      call history_add_field('qi_pre_PBL',  'Cloud ice water mixing ratio before PBL',            'lev',   'inst', 'kg/kg')
      call history_add_field('t_pre_PBL',   'Temperature before PBL',                             'lev',   'inst', 'K')
      call history_add_field('rh_pre_PBL',  'Relative humidity before PBL',                       'lev',   'inst', '%')

      ! Post-PBL profiles
      call history_add_field('qt_aft_PBL',  'Total water mixing ratio after PBL',                 'lev',   'inst', 'kg/kg')
      call history_add_field('sl_aft_PBL',  'Liquid water static energy after PBL',               'lev',   'inst', 'J/kg')
      call history_add_field('slv_aft_PBL', 'Liquid water virtual static energy after PBL',       'lev',   'inst', 'J/kg')
      call history_add_field('u_aft_PBL',   'Zonal wind after PBL',                               'lev',   'inst', 'm/s')
      call history_add_field('v_aft_PBL',   'Meridional wind after PBL',                          'lev',   'inst', 'm/s')
      call history_add_field('qv_aft_PBL',  'Water vapor mixing ratio after PBL',                 'lev',   'inst', 'kg/kg')
      call history_add_field('ql_aft_PBL',  'Cloud liquid water mixing ratio after PBL',          'lev',   'inst', 'kg/kg')
      call history_add_field('qi_aft_PBL',  'Cloud ice water mixing ratio after PBL',             'lev',   'inst', 'kg/kg')
      call history_add_field('t_aft_PBL',   'Temperature after PBL',                              'lev',   'inst', 'K')
      call history_add_field('rh_aft_PBL',  'Relative humidity after PBL',                        'lev',   'inst', '%')

      ! PBL fluxes
      call history_add_field('slflx_PBL',   'Liquid static energy flux by PBL',                   'ilev',  'inst', 'W/m2')
      call history_add_field('qtflx_PBL',   'Total water flux by PBL',                            'ilev',  'inst', 'kg/m2/s')
      call history_add_field('uflx_PBL',    'Zonal momentum flux by PBL',                         'ilev',  'inst', 'N/m2')
      call history_add_field('vflx_PBL',    'Meridional momentum flux by PBL',                    'ilev',  'inst', 'N/m2')

      ! Counter-gradient fluxes
      call history_add_field('slflx_cg_PBL','Counter-gradient liquid static energy flux by PBL',  'ilev',  'inst', 'W/m2')
      call history_add_field('qtflx_cg_PBL','Counter-gradient total water flux by PBL',           'ilev',  'inst', 'kg/m2/s')
      call history_add_field('uflx_cg_PBL', 'Counter-gradient zonal momentum flux by PBL',        'ilev',  'inst', 'N/m2')
      call history_add_field('vflx_cg_PBL', 'Counter-gradient meridional momentum flux by PBL',   'ilev',  'inst', 'N/m2')

      ! PBL tendencies
      call history_add_field('qtten_PBL',   'Total water mixing ratio tendency by PBL',           'lev',   'inst', 'kg/kg/s')
      call history_add_field('slten_PBL',   'Liquid static energy tendency by PBL',               'lev',   'inst', 'J/kg/s')
      call history_add_field('uten_PBL',    'Zonal wind tendency by PBL',                         'lev',   'inst', 'm/s2')
      call history_add_field('vten_PBL',    'Meridional wind tendency by PBL',                    'lev',   'inst', 'm/s2')
      call history_add_field('qvten_PBL',   'Water vapor tendency by PBL',                        'lev',   'inst', 'kg/kg/s')
      call history_add_field('qlten_PBL',   'Cloud liquid water tendency by PBL',                 'lev',   'inst', 'kg/kg/s')
      call history_add_field('qiten_PBL',   'Cloud ice water tendency by PBL',                    'lev',   'inst', 'kg/kg/s')
      call history_add_field('tten_PBL',    'Temperature tendency by PBL',                        'lev',   'inst', 'K/s')
      call history_add_field('rhten_PBL',   'Relative humidity tendency by PBL',                  'lev',   'inst', '%/s')

      ! Vertical diffusion tendencies
      call history_add_field('DTVKE',       'dT/dt vertical diffusion KE dissipation',            'lev',   'inst', 'K s-1')
      call history_add_field('DTV',         'Temperature vertical diffusion',                     'lev',   'inst', 'K s-1')
      call history_add_field('DUV',         'Zonal wind vertical diffusion',                      'lev',   'inst', 'm s-2')
      call history_add_field('DVV',         'Meridional wind vertical diffusion',                 'lev',   'inst', 'm s-2')

      ! These are constituent dependent and should maybe be initialized using the const props.
      ! call history_add_field('VD01',        'Vertical diffusion of water vapor',                  'lev',   'inst', 'kg/kg/s')
      ! call history_add_field('VDCLDLIQ',    'Vertical diffusion of cloud liquid water',           'lev',   'inst', 'kg/kg/s')
      ! call history_add_field('VDCLDICE',    'Vertical diffusion of cloud ice water',              'lev',   'inst', 'kg/kg/s')

  end subroutine vertical_diffusion_tendencies_diagnostics_init

!> \section arg_table_vertical_diffusion_tendencies_diagnostics_run Argument Table
!! \htmlinclude vertical_diffusion_tendencies_diagnostics_run.html
  subroutine vertical_diffusion_tendencies_diagnostics_run( &
             ncol, pver, pverp, ztodt, &
             const_props, &
             latvap, latice, zvir, cpair, gravit, rair, &
             pmid, pint, zi, zm, &
             kvh, kvm, cgh, cgs, &
             cam_in_shf, cam_in_cflx, &
             tautotx, tautoty, &
             t0, &
             q0, s0, u0, v0, &
             q1, s1, u1, v1, &
             dtk, &
             tend_s, tend_u, tend_v, tend_q, &
             errmsg, errflg)

      use cam_history, only: history_out_field
      use runtime_obj, only: wv_stdname

      ! framework dependency for const_props
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

      ! dependency to get constituent index
      use ccpp_const_utils,          only: ccpp_const_get_idx

      ! utility subroutine for calculating diagnostics based on outputs
      ! from the vertical diffusion solver.
      use vertical_diffusion_diagnostic_utils, only: vertical_diffusion_diagnostic_profiles

      ! Input arguments
      integer,          intent(in) :: ncol
      integer,          intent(in) :: pver
      integer,          intent(in) :: pverp
      real(kind_phys),  intent(in) :: ztodt
      type(ccpp_constituent_prop_ptr_t), intent(in) :: const_props(:)
      real(kind_phys),  intent(in) :: latvap                  ! Latent heat of vaporization [J kg-1]
      real(kind_phys),  intent(in) :: latice                  ! Latent heat of fusion [J kg-1]
      real(kind_phys),  intent(in) :: zvir                    ! rh2o/rair - 1 [1]
      real(kind_phys),  intent(in) :: cpair                   ! Specific heat of dry air at constant pressure [J kg-1 K-1]
      real(kind_phys),  intent(in) :: gravit                  ! Acceleration due to gravity [m s-2]
      real(kind_phys),  intent(in) :: rair                    ! Gas constant for dry air [J kg-1 K-1]
      real(kind_phys),  intent(in) :: pmid(:, :)              ! Midpoint pressure [Pa]
      real(kind_phys),  intent(in) :: pint(:, :)              ! Interface pressure [Pa]
      real(kind_phys),  intent(in) :: zi(:, :)                ! Geopotential height at interfaces [m]
      real(kind_phys),  intent(in) :: zm(:, :)                ! Geopotential height at midpoints [m]

      ! Diffusion coefficients and PBL diagnostics
      real(kind_phys),  intent(in) :: kvh(:, :)               ! Eddy diffusivity for heat at interfaces [m2 s-1]
      real(kind_phys),  intent(in) :: kvm(:, :)               ! Eddy diffusivity for momentum at interfaces [m2 s-1]
      real(kind_phys),  intent(in) :: cgh(:, :)               ! Counter-gradient term for heat [K m s-1]
      real(kind_phys),  intent(in) :: cgs(:, :)               ! Counter-gradient star [s m-2]

      ! Surface fluxes from coupler
      real(kind_phys),  intent(in) :: cam_in_shf(:)           ! Surface sensible heat flux [W m-2]
      real(kind_phys),  intent(in) :: cam_in_cflx(:,:)        ! Surface constituent fluxes [kg m-2 s-1]

      ! Surface stresses
      real(kind_phys),  intent(in) :: tautotx(:)              ! Total zonal surface stress [N m-2]
      real(kind_phys),  intent(in) :: tautoty(:)              ! Total meridional surface stress [N m-2]

      ! Initial state (before vertical diffusion)
      real(kind_phys),  intent(in) :: t0(:, :)                ! Temperature before diffusion [K]
      real(kind_phys),  intent(in) :: q0(:, :, :)             ! Constituent mixing ratios before diffusion [kg kg-1]
      real(kind_phys),  intent(in) :: s0(:, :)                ! Dry static energy before diffusion [J kg-1]
      real(kind_phys),  intent(in) :: u0(:, :)                ! Zonal wind before diffusion [m s-1]
      real(kind_phys),  intent(in) :: v0(:, :)                ! Meridional wind before diffusion [m s-1]

      ! "After state" (provisionally updated) quantities from diffusion solver
      real(kind_phys),  intent(in) :: q1(:, :, :)             ! Constituent mixing ratios after diffusion [J kg-1]
      real(kind_phys),  intent(in) :: s1(:, :)                ! Dry static energy after diffusion [J kg-1]
      real(kind_phys),  intent(in) :: u1(:, :)                ! Zonal wind after diffusion [m s-1]
      real(kind_phys),  intent(in) :: v1(:, :)                ! Meridional wind after diffusion [m s-1]

      ! Output from vertical diffusion kinetic energy (KE) dissipation
      ! (in vertical_diffusion_diffuse_horizontal_momentum_run)
      real(kind_phys),  intent(in)  :: dtk(:, :)              ! T tendency from KE dissipation [J kg-1]

      ! Tendencies from vertical diffusion
      real(kind_phys),  intent(in) :: tend_s(:, :)            ! Dry static energy tendency [J kg-1 s-1]
      real(kind_phys),  intent(in) :: tend_u(:, :)            ! Zonal wind tendency [m s-2]
      real(kind_phys),  intent(in) :: tend_v(:, :)            ! Meridional wind tendency [m s-2]
      real(kind_phys),  intent(in) :: tend_q(:, :, :)         ! Constituent tendency from vertical diffusion [kg kg-1 s-1]

      ! Error handling
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables for constituent indices
      integer :: const_wv_idx, const_cldliq_idx, const_cldice_idx

      ! Local variables for diagnostic profiles
      real(kind_phys) :: qt_pre_PBL(ncol, pver)        ! Total water mixing ratio before PBL [kg kg-1]
      real(kind_phys) :: sl_pre_PBL(ncol, pver)        ! Liquid water static energy before PBL [J kg-1]
      real(kind_phys) :: slv_pre_PBL(ncol, pver)       ! Virtual liquid water static energy before PBL [J kg-1]
      real(kind_phys) :: ftem_pre_PBL(ncol, pver)      ! Relative humidity before PBL [%]
      real(kind_phys) :: qt_aft_PBL(ncol, pver)        ! Total water mixing ratio after PBL [kg kg-1]
      real(kind_phys) :: sl_aft_PBL(ncol, pver)        ! Liquid water static energy after PBL [J kg-1]
      real(kind_phys) :: slv_aft_PBL(ncol, pver)       ! Virtual liquid water static energy after PBL [J kg-1]
      real(kind_phys) :: qv_aft_PBL(ncol, pver)        ! Water vapor mixing ratio after PBL [kg kg-1]
      real(kind_phys) :: ql_aft_PBL(ncol, pver)        ! Cloud liquid water mixing ratio after PBL [kg kg-1]
      real(kind_phys) :: qi_aft_PBL(ncol, pver)        ! Cloud ice water mixing ratio after PBL [kg kg-1]
      real(kind_phys) :: t_aft_PBL(ncol, pver)         ! Temperature after PBL [K]
      real(kind_phys) :: ftem_aft_PBL(ncol, pver)      ! Relative humidity after PBL [%]
      real(kind_phys) :: u_aft_PBL(ncol, pver)         ! Zonal wind after PBL [m s-1]
      real(kind_phys) :: v_aft_PBL(ncol, pver)         ! Meridional wind after PBL [m s-1]
      real(kind_phys) :: slflx(ncol, pverp)            ! Liquid static energy flux at interfaces [W m-2]
      real(kind_phys) :: qtflx(ncol, pverp)            ! Total water flux at interfaces [kg m-2 s-1]
      real(kind_phys) :: uflx(ncol, pverp)             ! Zonal momentum flux at interfaces [N m-2]
      real(kind_phys) :: vflx(ncol, pverp)             ! Meridional momentum flux at interfaces [N m-2]
      real(kind_phys) :: slflx_cg(ncol, pverp)         ! Counter-gradient liquid static energy flux [W m-2]
      real(kind_phys) :: qtflx_cg(ncol, pverp)         ! Counter-gradient total water flux [kg m-2 s-1]
      real(kind_phys) :: uflx_cg(ncol, pverp)          ! Counter-gradient zonal momentum flux [N m-2]
      real(kind_phys) :: vflx_cg(ncol, pverp)          ! Counter-gradient meridional momentum flux [N m-2]
      real(kind_phys) :: slten(ncol, pver)             ! Liquid water static energy tendency [J kg-1 s-1]
      real(kind_phys) :: qtten(ncol, pver)             ! Total water mixing ratio tendency [kg kg-1 s-1]
      real(kind_phys) :: tten(ncol, pver)              ! Temperature tendency [K s-1]
      real(kind_phys) :: rhten(ncol, pver)             ! Relative humidity tendency [% s-1]
      real(kind_phys) :: dtvke(ncol, pver)             ! Normalized KE dissipation heating for history [K s-1]
      real(kind_phys) :: dtv(ncol, pver)               ! Normalized temperature tendency for history [K s-1]

      ! Subset constituent upward ccpp constituent fluxes from coupler
      real(kind_phys) :: cam_in_cflx_wv(ncol)
      real(kind_phys) :: cam_in_cflx_cldliq(ncol)
      real(kind_phys) :: cam_in_cflx_cldice(ncol)

      ! Subset q0, q1, tendencies for constituents from coupler
      real(kind_phys) :: q0_wv(ncol, pver)
      real(kind_phys) :: q0_cldliq(ncol, pver)
      real(kind_phys) :: q0_cldice(ncol, pver)
      real(kind_phys) :: q1_cldliq(ncol, pver)
      real(kind_phys) :: q1_cldice(ncol, pver)
      real(kind_phys) :: tend_q_wv(ncol, pver)
      real(kind_phys) :: tend_q_cldliq(ncol, pver)
      real(kind_phys) :: tend_q_cldice(ncol, pver)

      errmsg = ''
      errflg = 0

      ! Get constituent indices for wv, cldliq, cldice
      call ccpp_const_get_idx(const_props, &
           wv_stdname, &
           const_wv_idx, errmsg, errflg)
      if (errflg /= 0) return

      call ccpp_const_get_idx(const_props, &
           'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water', &
           const_cldliq_idx, errmsg, errflg)
      if (errflg /= 0) return

      call ccpp_const_get_idx(const_props, &
           'cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water', &
           const_cldice_idx, errmsg, errflg)
      if (errflg /= 0) return

      ! Extract constituent fluxes from coupler upward constituent fluxes
      if(const_wv_idx > 0) then
        cam_in_cflx_wv(:ncol)     = cam_in_cflx(:ncol, const_wv_idx)
        q0_wv(:ncol, :)           = q0(:ncol, :, const_wv_idx)
        tend_q_wv(:ncol, :)       = tend_q(:ncol, :, const_wv_idx)
      else
        cam_in_cflx_wv(:ncol)     = 0._kind_phys
        q0_wv(:ncol, :) = 0._kind_phys
        tend_q_wv(:ncol, :) = 0._kind_phys
      endif

      if(const_cldliq_idx > 0) then
        cam_in_cflx_cldliq(:ncol) = cam_in_cflx(:ncol, const_cldliq_idx)
        q0_cldliq(:ncol, :)       = q0(:ncol, :, const_cldliq_idx)
        q1_cldliq(:ncol, :)       = q1(:ncol, :, const_cldliq_idx)
        tend_q_cldliq(:ncol, :)   = tend_q(:ncol, :, const_cldliq_idx)
      else
        cam_in_cflx_cldliq(:ncol) = 0._kind_phys
        q0_cldliq(:ncol, :) = 0._kind_phys
        q1_cldliq(:ncol, :) = 0._kind_phys
        tend_q_cldliq(:ncol, :) = 0._kind_phys
      endif

      if(const_cldice_idx > 0) then
        cam_in_cflx_cldice(:ncol) = cam_in_cflx(:ncol, const_cldice_idx)
        q0_cldice(:ncol, :)       = q0(:ncol, :, const_cldice_idx)
        q1_cldice(:ncol, :)       = q1(:ncol, :, const_cldice_idx)
        tend_q_cldice(:ncol, :)   = tend_q(:ncol, :, const_cldice_idx)
      else
        cam_in_cflx_cldice(:ncol) = 0._kind_phys
        q0_cldice(:ncol, :) = 0._kind_phys
        q1_cldice(:ncol, :) = 0._kind_phys
        tend_q_cldice(:ncol, :) = 0._kind_phys
      endif

      ! Call the diagnostic profile calculation subroutine
      call vertical_diffusion_diagnostic_profiles( &
           ncol                 = ncol, &
           pver                 = pver, &
           pverp                = pverp, &
           ztodt                = ztodt, &
           latvap               = latvap, &
           latice               = latice, &
           zvir                 = zvir, &
           cpair                = cpair, &
           gravit               = gravit, &
           rair                 = rair, &
           pint                 = pint(:ncol, :pverp), &
           pmid                 = pmid(:ncol, :pver), &
           zi                   = zi(:ncol, :pverp), &
           zm                   = zm(:ncol, :pver), &
           kvh                  = kvh(:ncol, :pverp), &
           kvm                  = kvm(:ncol, :pverp), &
           cgh                  = cgh(:ncol, :pverp), &
           cgs                  = cgs(:ncol, :pverp), &
           cam_in_shf           = cam_in_shf(:ncol), &
           cam_in_cflx_wv       = cam_in_cflx_wv(:ncol), &
           cam_in_cflx_cldliq   = cam_in_cflx_cldliq(:ncol), &
           cam_in_cflx_cldice   = cam_in_cflx_cldice(:ncol), &
           tautotx              = tautotx(:ncol), & ! total surface stresses to PBL coefficient scheme, not from provisional
           tautoty              = tautoty(:ncol), & ! total surface stresses to PBL coefficient scheme, not from provisional
           t0                   = t0(:ncol, :pver), &
           q0_wv                = q0_wv(:ncol, :pver), &
           q0_cldliq            = q0_cldliq(:ncol, :pver), &
           q0_cldice            = q0_cldice(:ncol, :pver), &
           s0                   = s0(:ncol, :pver), &
           u0                   = u0(:ncol, :pver), &
           v0                   = v0(:ncol, :pver), &
           q1_cldliq            = q1_cldliq(:ncol, :pver), &
           q1_cldice            = q1_cldice(:ncol, :pver), &
           s1                   = s1(:ncol, :pver), &
           u1                   = u1(:ncol, :pver), &
           v1                   = v1(:ncol, :pver), &
           tend_s               = tend_s(:ncol, :pver), &
           tend_u               = tend_u(:ncol, :pver), &
           tend_v               = tend_v(:ncol, :pver), &
           tend_q_wv            = tend_q_wv(:ncol, :pver), &
           tend_q_cldliq        = tend_q_cldliq(:ncol, :pver), &
           tend_q_cldice        = tend_q_cldice(:ncol, :pver), &
           ! below output:
           qt_pre_PBL           = qt_pre_PBL(:ncol, :pver), &
           sl_pre_PBL           = sl_pre_PBL(:ncol, :pver), &
           slv_pre_PBL          = slv_pre_PBL(:ncol, :pver), &
           ftem_pre_PBL         = ftem_pre_PBL(:ncol, :pver), &
           qt_aft_PBL           = qt_aft_PBL(:ncol, :pver), &
           sl_aft_PBL           = sl_aft_PBL(:ncol, :pver), &
           slv_aft_PBL          = slv_aft_PBL(:ncol, :pver), &
           qv_aft_PBL           = qv_aft_PBL(:ncol, :pver), &
           ql_aft_PBL           = ql_aft_PBL(:ncol, :pver), &
           qi_aft_PBL           = qi_aft_PBL(:ncol, :pver), &
           t_aft_PBL            = t_aft_PBL(:ncol, :pver), &
           ftem_aft_PBL         = ftem_aft_PBL(:ncol, :pver), &
           u_aft_PBL            = u_aft_PBL(:ncol, :pver), &
           v_aft_PBL            = v_aft_PBL(:ncol, :pver), &
           slflx                = slflx(:ncol, :pverp), &
           qtflx                = qtflx(:ncol, :pverp), &
           uflx                 = uflx(:ncol, :pverp), &
           vflx                 = vflx(:ncol, :pverp), &
           slflx_cg             = slflx_cg(:ncol, :pverp), &
           qtflx_cg             = qtflx_cg(:ncol, :pverp), &
           uflx_cg              = uflx_cg(:ncol, :pverp), &
           vflx_cg              = vflx_cg(:ncol, :pverp), &
           slten                = slten(:ncol, :pver), &
           qtten                = qtten(:ncol, :pver), &
           tten                 = tten(:ncol, :pver), &
           rhten                = rhten(:ncol, :pver))

      ! Standard output variables
      call history_out_field('QT',          qt_aft_PBL)
      call history_out_field('SL',          sl_aft_PBL)
      call history_out_field('SLV',         slv_aft_PBL)
      call history_out_field('SLFLX',       slflx)
      call history_out_field('QTFLX',       qtflx)
      call history_out_field('UFLX',        uflx)
      call history_out_field('VFLX',        vflx)

      ! Post-PBL profiles
      call history_out_field('sl_aft_PBL',  sl_aft_PBL)
      call history_out_field('qt_aft_PBL',  qt_aft_PBL)
      call history_out_field('slv_aft_PBL', slv_aft_PBL)
      call history_out_field('u_aft_PBL',   u_aft_PBL)
      call history_out_field('v_aft_PBL',   v_aft_PBL)
      call history_out_field('qv_aft_PBL',  qv_aft_PBL)
      call history_out_field('ql_aft_PBL',  ql_aft_PBL)
      call history_out_field('qi_aft_PBL',  qi_aft_PBL)
      call history_out_field('t_aft_PBL',   t_aft_PBL)
      call history_out_field('rh_aft_PBL',  ftem_aft_PBL)

      ! PBL fluxes
      call history_out_field('slflx_PBL',   slflx)
      call history_out_field('qtflx_PBL',   qtflx)
      call history_out_field('uflx_PBL',    uflx)
      call history_out_field('vflx_PBL',    vflx)

      ! Counter-gradient fluxes
      call history_out_field('slflx_cg_PBL', slflx_cg)
      call history_out_field('qtflx_cg_PBL', qtflx_cg)
      call history_out_field('uflx_cg_PBL',  uflx_cg)
      call history_out_field('vflx_cg_PBL',  vflx_cg)

      ! PBL tendencies
      call history_out_field('slten_PBL',   slten)
      call history_out_field('qtten_PBL',   qtten)
      call history_out_field('uten_PBL',    tend_u)
      call history_out_field('vten_PBL',    tend_v)
      call history_out_field('qvten_PBL',   tend_q_wv)
      call history_out_field('qlten_PBL',   tend_q_cldliq)
      call history_out_field('qiten_PBL',   tend_q_cldice)
      call history_out_field('tten_PBL',    tten)
      call history_out_field('rhten_PBL',   rhten)

      ! Pre-PBL profiles
      call history_out_field('qt_pre_PBL',  qt_pre_PBL)
      call history_out_field('sl_pre_PBL',  sl_pre_PBL)
      call history_out_field('slv_pre_PBL', slv_pre_PBL)
      call history_out_field('u_pre_PBL',   u0)
      call history_out_field('v_pre_PBL',   v0)
      call history_out_field('qv_pre_PBL',  q0_wv)
      call history_out_field('ql_pre_PBL',  q0_cldliq)
      call history_out_field('qi_pre_PBL',  q0_cldice)
      call history_out_field('t_pre_PBL',   t0)
      call history_out_field('rh_pre_PBL',  ftem_pre_PBL)

      ! Vertical diffusion tendencies
      call history_out_field('DUV',         tend_u)
      call history_out_field('DVV',         tend_v)
      call history_out_field('VD01',        tend_q_wv)
      call history_out_field('VDCLDLIQ',    tend_q_cldliq)
      call history_out_field('VDCLDICE',    tend_q_cldice)

      ! Normalize kinetic energy dissipation heating for history output
      ! Convert from [J kg-1] to [K s-1]
      dtvke(:ncol, :pver) = dtk(:ncol, :pver) / cpair / ztodt

      ! Normalize dry static energy tendency for history output
      ! Convert from [J kg-1 s-1] to [K s-1]
      dtv(:ncol, :pver) = tend_s(:ncol, :pver) / cpair

      ! Output diagnostics
      call history_out_field('DTVKE', dtvke)
      call history_out_field('DTV',   dtv)

  end subroutine vertical_diffusion_tendencies_diagnostics_run


end module diffusion_solver_diagnostics
