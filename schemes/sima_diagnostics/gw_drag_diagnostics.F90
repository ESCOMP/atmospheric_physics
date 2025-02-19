()! Diagnostics for gravity wave parameterizations
module gw_drag_diagnostics
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: gw_drag_diagnostics_init

  ! The gravity wave diagnostics have several phases,
  ! which have to run before the tendency updaters in order to capture
  ! the appropriate tendencies.

  ! 1) after the gw_oro scheme runs
  public :: gw_drag_diagnostics_after_gw_oro_scheme_run
  ! 2) after the gw_front scheme runs
  public :: gw_drag_diagnostics_after_gw_front_scheme_run
  ! 3) after the gw_front_igw scheme runs
  public :: gw_drag_diagnostics_after_gw_front_igw_scheme_run
  ! 4) after deep convect gw scheme
  public :: gw_drag_diagnostics_after_gw_convect_dp_scheme_run
  ! 5) after shallow convect gw scheme
  public :: gw_drag_diagnostics_after_gw_convect_sh_scheme_run
  ! 6) after shallow convect gw scheme
  public :: gw_drag_diagnostics_after_gw_movmtn_pbl_scheme_run
  ! 7) common gw spectral init for all schemes
  public :: gw_spec_addfld
  ! 8) common gw spectral output for all schemes
  public :: gw_spec_outfld
  ! 9) after shallow convective quantities are merged with deep convective quantities
  !    this outputs the final convection diagnostics taking into account both deep and gravity wave
  public :: gw_drag_diagnostics_after_sum_all_gw_schemes_run

contains

!> \section arg_table_gw_drag_diagnostics_init  Argument Table
!! \htmlinclude gw_drag_diagnostics_init.html
  subroutine gw_drag_diagnostics_init( &
       errmsg, errflg)
    use cam_history,         only: history_add_field
    use cam_history_support, only: horiz_only

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    character(len=1) :: cn,nn
    integer          :: i, l, k

    errmsg = ''
    errflg = 0

    !=======================================================
    ! Gravity Wave diagnostics for shallow schemes (common)
    !=======================================================

    if (use_gw_oro .or. use_gw_rdg_beta .or. use_gw_rdg_gamma) then
       call history_add_field ('TAUAORO', 'Total stress from original OGW scheme', 'lev', 'inst','N m-2')
       call history_add_field ('TTGWORO', 'T tendency - orographic gravity wave drag', 'lev', 'avg','K s-1')
       call history_add_field ('TTGWSDFORO', 'T tendency - orographic gravity wave, diffusion.', 'lev', 'avg','K s-1')
       call history_add_field ('TTGWSKEORO', 'T tendency - orographic gravity wave, breaking KE.', 'lev', 'avg','K s-1')
       call history_add_field ('UTGWORO', 'U tendency - orographic gravity wave drag', 'lev', 'avg','m s-2')
       call history_add_field ('VTGWORO', 'V tendency - orographic gravity wave drag', 'lev', 'avg','m s-2')
       call history_add_field ('TAUGWX', 'Zonal gravity wave surface stress',     horiz_only,  'avg','N m-2')
       call history_add_field ('TAUGWY', 'Meridional gravity wave surface stress',     horiz_only,  'avg','N m-2')
    end if

    if (use_gw_rdg_beta) then
       call history_add_field ('WBR_HT1', 'Wave breaking height for DSW',     horiz_only,  'inst','m')
       call history_add_field ('TLB_HT1', 'Form drag layer height',     horiz_only,  'inst','m')
       call history_add_field ('BWV_HT1', 'Bottom of freely-propagating OGW regime',     horiz_only,  'inst','m')
       call history_add_field ('TAUDSW1', 'DSW enhanced drag',     horiz_only,  'inst','Nm-2')
       call history_add_field ('TAUORO1', 'lower BC on propagating wave stress',     horiz_only,  'inst','Nm-2')
       call history_add_field ('UBMSRC1', 'below-peak-level on-ridge wind',     horiz_only,  'inst','ms-1')
       call history_add_field ('USRC1', 'below-peak-level Zonal wind',     horiz_only,  'inst','ms-1')
       call history_add_field ('VSRC1', 'below-peak-level Meridional wind',     horiz_only,  'inst','ms-1')
       call history_add_field ('NSRC1', 'below-peak-level stratification',     horiz_only,  'inst','s-1')
       call history_add_field ('MXDIS1', 'Ridge/obstacle height',     horiz_only,  'inst','m')
       call history_add_field ('ANGLL1', 'orientation clockwise w/resp north-south',     horiz_only,  'inst','degrees')
       call history_add_field ('ANIXY1', 'Ridge quality',     horiz_only,  'inst','1')
       call history_add_field ('HWDTH1', 'Ridge width',     horiz_only,  'inst','km')
       call history_add_field ('CLNGT1', 'Ridge length',     horiz_only,  'inst','km')
       call history_add_field ('GBXAR1', 'grid box area',     horiz_only,  'inst','km+2')
       call history_add_field ('Fr1_DIAG', 'Critical Froude number for linear waves',     horiz_only,  'inst','1')
       call history_add_field ('Fr2_DIAG', 'Critical Froude number for blocked flow',     horiz_only,  'inst','1')
       call history_add_field ('Frx_DIAG', 'Obstacle Froude Number',     horiz_only,  'inst','1')
       call history_add_field('UEGW', 'Zonal wind profile-entry to GW ', 'lev', 'avg'  ,'s-1' )
       call history_add_field('VEGW', 'Merdional wind profile-entry to GW ', 'lev', 'avg'  ,'s-1' )
       !jt missing vector register?
       call history_add_field('TEGW', 'Temperature profile-entry to GW ', 'lev', 'avg'  ,'K' )
       call history_add_field('ZEGW', 'interface geopotential heights in GW code ', 'ilev', 'avg'  ,'m' )
       call history_add_field('ZMGW', 'midlayer geopotential heights in GW code ', 'lev', 'avg'  ,'m' )
       call history_add_field('NIEGW', 'interface BV freq in GW code ', 'ilev', 'inst'  ,'s-1' )
       call history_add_field('NMEGW', 'midlayer BV freq in GW code ', 'lev', 'inst'  ,'s-1' )
       call history_add_field('RHOIEGW', 'interface density in GW code ', 'ilev', 'inst'  ,'kg m-33' )
       call history_add_field('PINTEGW', 'interface air pressure in GW code ', 'ilev', 'inst'  ,'Pa' )
       call history_add_field('TAUM1_DIAG' , 'Ridge based momentum flux profile', 'ilev', 'inst'  ,'N m-2' )
       call history_add_field('TAU1RDGBETAM' , 'Ridge based momentum flux profile', 'ilev', 'inst'  ,'N m-2' )
       call history_add_field('UBM1BETA', 'On-ridge wind profile           ', 'lev', 'avg'  ,'m s-1' )
       call history_add_field('UBT1RDGBETA' , 'On-ridge wind tendency from ridge 1     ', 'lev', 'inst'  ,'m s-2' )
       call history_add_field('TAURESIDBETAM' , 'Ridge based momentum flux profile', 'ilev', 'inst'  ,'N m-2' )
       call history_add_field('UBMRESIDBETA', 'On-ridge wind profile           ', 'lev', 'inst'  ,'m s-1' )
       call history_add_field('UBIRESIDBETA', 'On-ridge wind profile (interface)          ', 'ilev', 'inst'  ,'m s-1' )
       call history_add_field('SRC_LEVEL_RESIDBETA', 'src level index for ridge residual         ',  horiz_only , 'inst'  ,'1' )
       call history_add_field('TAUORO_RESID', 'Surface momentum flux from ridge residual       ',  horiz_only , 'inst'  ,'N m-2' )
       call history_add_field('TAUDIAG_RESID' , 'Ridge based momentum flux profile', 'ilev', 'inst'  ,'N m-2' )

       do i = 1, 6
          write(cn, '(i1)') i
          call history_add_field('TAU'//cn//'RDGBETAY' , 'Ridge based momentum flux profile', 'ilev', 'inst', 'N m-2')
          call history_add_field('TAU'//cn//'RDGBETAX' , 'Ridge based momentum flux profile', 'ilev', 'inst', 'N m-2')
          !jt missing vector register?
          call history_add_field('UT'//cn//'RDGBETA', 'U wind tendency from ridge '//cn, 'lev',  'inst', 'm s-1')
          call history_add_field('VT'//cn//'RDGBETA', 'V wind tendency from ridge '//cn, 'lev',  'inst', 'm s-1')
          !jt missing vector register?
       end do

       call history_add_field('TAUARDGBETAY' , 'Ridge based momentum flux profile', 'ilev', 'inst'  ,'N m-2' )
       call history_add_field('TAUARDGBETAX' , 'Ridge based momentum flux profile', 'ilev', 'inst'  ,'N m-2' )
       !jt missing vector register?
       call history_add_field('TAURESIDBETAY' , 'Ridge based momentum flux profile', 'ilev', 'inst'  ,'N m-2' )
       call history_add_field('TAURESIDBETAX' , 'Ridge based momentum flux profile', 'ilev', 'inst'  ,'N m-2' )
       !jt missing vector register?
    end if
    if (use_gw_rdg_gamma) then
       call history_add_field ('TAU1RDGGAMMAM' , 'Ridge based momentum flux profile', 'ilev', 'inst'  ,'N m-2' )
       call history_add_field ('UBM1GAMMA', 'On-ridge wind profile           ', 'lev', 'avg'  ,'s-1' )
       call history_add_field ('UBT1RDGGAMMA' , 'On-ridge wind tendency from ridge 1     ', 'lev', 'inst'  ,'m s-1' )

       do i = 1, 6
          write(cn, '(i1)') i
          call history_add_field('TAU'//cn//'RDGGAMMAY', 'Ridge based momentum flux profile', 'ilev', 'inst', 'N m-2')
          call history_add_field('TAU'//cn//'RDGGAMMAX', 'Ridge based momentum flux profile', 'ilev', 'inst', 'N m-2')
          call history_add_field('UT'//cn//'RDGGAMMA' , 'U wind tendency from ridge '//cn, 'lev',  'inst', 'm s-1')
          call history_add_field('VT'//cn//'RDGGAMMA' , 'V wind tendency from ridge '//cn, 'lev',  'inst', 'm s-1')
          !jt  missing register vector?
       end do

       call history_add_field ('TAUARDGGAMMAY' , 'Ridge based momentum flux profile', 'ilev', 'inst'  ,'N m-2' )
       call history_add_field ('TAUARDGGAMMAX' , 'Ridge based momentum flux profile', 'ilev', 'inst'  ,'N m-2' )
       !jt  missing register vector?
       call history_add_field ('TAURDGGMX', 'Zonal gravity wave surface stress',     horiz_only,  'avg','N m-2')
       call history_add_field ('TAURDGGMY', 'Meridional gravity wave surface stress',     horiz_only,  'avg','N m-2')
       !jt  missing register vector?
       call history_add_field ('UTRDGGM' , 'U wind tendency from ridge 6     ', 'lev', 'inst'  ,'m s-1' )
       call history_add_field ('VTRDGGM' , 'V wind tendency from ridge 6     ', 'lev', 'inst'  ,'m s-1' )
       !jt  missing register vector?
    end if

    if (use_gw_front .or. use_gw_front_igw) then
       call history_add_field ('FRONTGF', 'Frontogenesis function at gws src level', 'lev', 'avg', 'K+2 M-2 S-1')
       call history_add_field ('FRONTGFA', 'Frontogenesis function at gws src level', 'lev', 'avg','K+2 M-2 S-1')
    end if

    if (use_gw_movmtn_pbl) then
       call history_add_field ('VORT4GW', 'Vorticity', 'lev', 'avg', 's-1')
    end if

    if (use_gw_movmtn_pbl) then
       call history_add_field ('GWUT_MOVMTN', 'Mov Mtn dragforce - ubm component', 'lev', 'inst','m s-2')
       call history_add_field ('UTGW_MOVMTN', 'Mov Mtn dragforce - u component', 'lev', 'inst','m s-2')
       call history_add_field ('VTGW_MOVMTN', 'Mov Mtn dragforce - v component', 'lev', 'inst','m s-2')
       call history_add_field('TAU_MOVMTN', 'Moving Mountain momentum flux profile', 'ilev', 'inst', 'N m-2')
       call history_add_field('U_MOVMTN_IN', 'Moving Mountain - midpoint zonal input wind', 'lev', 'inst', 'm s-1')
       call history_add_field('V_MOVMTN_IN', 'Moving Mountain - midpoint meridional input wind', 'lev', 'inst', 'm s-1')
       call history_add_field('UBI_MOVMTN', 'Moving Mountain - interface wind in direction of wave', 'ilev', 'inst', 'm s-1')
       call history_add_field('UBM_MOVMTN', 'Moving Mountain - midpoint wind in direction of wave', 'lev', 'inst', 'm s-1')
       call history_add_field ('HDEPTH_MOVMTN', 'Heating Depth',horiz_only,'inst','km')
       call history_add_field ('UCELL_MOVMTN', 'Gravity Wave Moving Mountain - Source-level X-wind',horiz_only,'inst','m s-1')
       call history_add_field ('VCELL_MOVMTN', 'Gravity Wave Moving Mountain - Source-level Y-wind',horiz_only,'inst','m s-1')
       call history_add_field ('CS_MOVMTN', 'Gravity Wave Moving Mountain - phase speed in direction of wave',horiz_only,'inst','m s-1')
       call history_add_field ('STEER_LEVEL_MOVMTN', 'Gravity Wave Moving Mountain - steering level for movmtn GW',horiz_only,'inst','1')
       call history_add_field ('SRC_LEVEL_MOVMTN', 'Gravity Wave Moving Mountain - launch level for movmtn GW',horiz_only,'inst','1')
       call history_add_field ('TND_LEVEL_MOVMTN', 'Gravity Wave Moving Mountain - tendency lowest level for movmtn GW',horiz_only,'inst','1')
       call history_add_field ('NETDT_MOVMTN', 'Gravity Wave Moving Mountain - Net heating rate', 'lev','inst','K s-1')
       call history_add_field ('TTEND_CLUBB', 'Gravity Wave Moving Mountain - CLUBB Net heating rate', 'lev','avg','K s-1')
       call history_add_field ('THLP2_CLUBB_GW', 'Gravity Wave Moving Mountain - THLP variance from CLUBB to GW', 'ilev','avg','K+2')
       call history_add_field ('WPTHLP_CLUBB_GW', 'Gravity Wave Moving Mountain - WPTHLP from CLUBB to GW', 'ilev','avg','Km s-2')
       call history_add_field ('UPWP_CLUBB_GW', 'Gravity Wave Moving Mountain - X-momflux from CLUBB to GW', 'ilev','avg','m+2 s-2')
       call history_add_field ('VPWP_CLUBB_GW', 'Gravity Wave Moving Mountain - Y-momflux from CLUBB to GW', 'ilev','avg','m+2 s-2')
       call history_add_field ('XPWP_SRC_MOVMTN', 'Gravity Wave Moving Mountain - flux source for moving mtn',horiz_only,'inst','m+2 s-2')
    end if
    if (use_gw_convect_dp) then
       call history_add_field ('NETDT', 'Net heating rate', 'lev', 'avg','K s-1')
       call history_add_field ('MAXQ0', 'Max column heating rate',horiz_only  ,  'avg','K day-1')
       call history_add_field ('HDEPTH', 'Heating Depth',horiz_only,    'avg','km')
    end if
    if (use_gw_convect_sh) then
       call history_add_field ('SNETDT', 'Net heating rate', 'lev', 'avg','K s-1')
       call history_add_field ('SMAXQ0', 'Max column heating rate',horiz_only  ,  'avg','K day-1')
       call history_add_field ('SHDEPTH', 'Heating Depth',horiz_only,    'avg','km')
    end if

    call history_add_field ('EKGW' , 'Effective Kzz due to diffusion by gravity waves', 'ilev', 'avg','m+2 s-1')
    call history_add_field ('UTGW_TOTAL', 'Total U tendency due to gravity wave drag', 'lev', 'avg','m s-2')
    call history_add_field ('VTGW_TOTAL', 'Total V tendency due to gravity wave drag', 'lev', 'avg','m s-2')
    !jt  missing register vector?
    call history_add_field ('TTGW', 'T tendency - gravity wave drag', 'lev', 'avg', 'K s-1')
    call history_add_field('QTGW', 'Q tendency - gravity wave drag', 'lev', 'avg','kg/kg/s')
    call history_add_field('CLDLIQTGW', 'CLDLIQ tendency - gravity wave drag', 'lev', 'avg','kg kg-1 s-1')
    call history_add_field('CLDICETGW', 'CLDICE tendency - gravity wave drag', 'lev', 'avg','kg kg-1 s-1')
  end subroutine gw_drag_diagnostics_init

!> \section arg_table_gw_spec_addflds  Argument Table
!! \htmlinclude gw_spec_addflds.html
  subroutine gw_spec_addflds(prefix, scheme, band, history_defaults, errmsg, errflg)
    ! Add all history fields for a gravity wave spectrum source.
    use cam_history,         only: history_add_field
    use cam_history_support, only: horiz_only

    ! One character prefix prepended to output fields.
    character(len=1), intent(in) :: prefix
    ! Gravity wave scheme name prepended to output field descriptions.
    character(len=*), intent(in) :: scheme
    ! Wave speeds.
    type(GWBand), intent(in) :: band
    ! Whether or not to call add_default for fields output by WACCM.
    logical, intent(in) :: history_defaults
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg
    !---------------------------Local storage-------------------------------

    integer :: l
    ! 7 chars is enough for "-100.00"
    character(len=7)  :: fnum
    ! 10 chars is enough for "BTAUXSn32"
    character(len=10) :: dumc1x, dumc1y
    ! Allow 80 chars for description
    character(len=80) dumc2

    !-----------------------------------------------------------------------

    ! Overall wind tendencies.

    call history_add_field (trim(prefix)//'UTGWSPEC', ' U tendency - gravity wave spectrum', 'lev', 'avg','m s-2')
    call history_add_field (trim(prefix)//'VTGWSPEC', ' V tendency - gravity wave spectrum', 'lev', 'avg','m s-2')
    !jt  missing register vector?
    call history_add_field (trim(prefix)//'TTGWSPEC', ' T tendency - gravity wave spectrum', 'lev', 'avg','K s-1')

    ! Wind tendencies broken across five spectral bins.
    call history_add_field (trim(prefix)//'UTEND1', ' U tendency   c < -40', 'lev', 'avg','m s-2')
    call history_add_field (trim(prefix)//'UTEND2', ' U tendency  -40 < c < -15', 'lev', 'avg','m s-2')
    call history_add_field (trim(prefix)//'UTEND3', ' U tendency  -15 < c <  15', 'lev', 'avg','m s-2')
    call history_add_field (trim(prefix)//'UTEND4', ' U tendency   15 < c <  40', 'lev', 'avg','m s-2')
    call history_add_field (trim(prefix)//'UTEND5', ' U tendency   40 < c ', 'lev', 'avg','m s-2')

    ! Reynold's stress toward each cardinal direction, and net zonal stress.
    call history_add_field (trim(prefix)//'TAUE' , ' Eastward Reynolds stress', 'lev', 'avg','Pa')
    call history_add_field (trim(prefix)//'TAUW' , ' Westward Reynolds stress', 'lev', 'avg','Pa')
    call history_add_field (trim(prefix)//'TAUNET' , ' E+W Reynolds stress', 'ilev', 'avg','Pa')
    call history_add_field (trim(prefix)//'TAUN' , ' Northward Reynolds stress', 'lev', 'avg','Pa')
    call history_add_field (trim(prefix)//'TAUS' , ' Southward Reynolds stress', 'lev', 'avg','Pa')

    ! Momentum flux in each direction.
    call history_add_field (trim(prefix)//'EMF', ' Eastward MF',    'lev', 'avg','Pa')
    call history_add_field (trim(prefix)//'WMF', ' Westward MF',    'lev', 'avg','Pa')
    call history_add_field (trim(prefix)//'NMF', ' Northward MF',    'lev', 'avg','Pa')
    call history_add_field (trim(prefix)//'SMF', ' Southward MF',    'lev', 'avg','Pa')

    ! Temperature tendency terms.
    call history_add_field (trim(prefix)//'TTGWSDF' , ' t tendency - diffusion term', 'lev', 'avg','K s-1')
    call history_add_field (trim(prefix)//'TTGWSKE' , ' t tendency - kinetic energy conversion term', 'lev', 'avg','K s-1')

    ! Gravity wave source spectra by wave number.
    do l=-band%ngwv,band%ngwv
       ! String containing reference speed.
       write (fnum,fmt='(f7.2)') band%cref(l)

       dumc1x = tau_fld_name(l, prefix, x_not_y=.true.)
       dumc1y = tau_fld_name(l, prefix, x_not_y=.false.)
       dumc2 = trim(scheme)//" tau at c= "//trim(fnum)//" m s-1"

       call history_add_field (trim(dumc1x), dumc2, 'lev', 'avg','Pa')
       call history_add_field (trim(dumc1y), dumc2, 'lev', 'avg','Pa')
    end do

  end subroutine gw_spec_addflds

  subroutine gw_drag_diagnostics_after_gw_movmtn_pbl_scheme_run( &
       ncol, pver, pcnst, &
       const_props, &
       cpair, &
       cmfdt, dq, cmfdqr, &
       qc_sh, icwmr, cmfsl, cmflq, cmfmc_sh, &
       errmsg, errflg)

    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    use cam_history, only: history_out_field

    ! Input arguments
    integer,                          intent(in)  :: ncol, pver, pcnst
    type(ccpp_constituent_prop_ptr_t),intent(in)  :: const_props(:)
    real(kind_phys),                  intent(in)  :: cpair
    real(kind_phys),                  intent(in)  :: cmfdt (:,:)
    real(kind_phys),                  intent(in)  :: dq    (:,:,:)
    real(kind_phys),                  intent(in)  :: cmfdqr(:,:)
    real(kind_phys),                  intent(in)  :: qc_sh (:,:)
    real(kind_phys),                  intent(in)  :: icwmr (:,:)
    real(kind_phys),                  intent(in)  :: cmfsl (:,:)
    real(kind_phys),                  intent(in)  :: cmflq (:,:)
    real(kind_phys),                  intent(in)  :: cmfmc_sh(:,:)

    ! Output arguments
    character(len=512),               intent(out) :: errmsg
    integer,                          intent(out) :: errflg

    ! Local variables
    real(kind_phys)                 :: freqsh(ncol)            ! Fractional occurrence of gravity wave [fraction]
    integer                         :: const_wv_idx, const_cldliq_idx, const_cldice_idx
    character(len=512)              :: const_standard_name
    integer                         :: i, m

    errmsg = ''
    errflg = 0

    ! Find constituent indices for water vapor, cloud liquid, cloud ice
    const_wv_idx     = -1
    const_cldliq_idx = -1
    const_cldice_idx = -1
    const_check_loop: do m = 1, pcnst
      call const_props(m)%standard_name(const_standard_name)
      if (const_standard_name == 'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water') then
        const_wv_idx = m
      endif

      if (const_standard_name == 'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water') then
        const_cldliq_idx = m
      endif

      if (const_standard_name == 'cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water') then
        const_cldice_idx = m
      endif
    enddo const_check_loop

    ! CMFDT - temperature tendency
    ! Calculate temperature tendency from specific heat flux [K s-1]
    call history_out_field('CMFDT', cmfdt(:,:)/cpair)
    if (use_gw_movmtn_pbl) then
       call history_out_field('U_MOVMTN_IN', u)
       call history_out_field('V_MOVMTN_IN', v)
       call history_out_field('SRC_LEVEL_MOVMTN', real(src_level,kind_phys))
       call history_out_field('TND_LEVEL_MOVMTN', real(tend_level,kind_phys))
       call history_out_field('UBI_MOVMTN', ubi)
       call history_out_field('UBM_MOVMTN', ubm)
       call history_out_field('TAU_MOVMTN', tau(:,0,:)
       call history_out_field('GWUT_MOVMTN', gwut(:,:,0)
       call history_out_field('VTGW_MOVMTN', vtgw)
       call history_out_field('UTGW_MOVMTN', utgw)
       call history_out_field('HDEPTH_MOVMTN', hdepth/1000._kind_phys)
       call history_out_field('NETDT_MOVMTN', ttend_dp)
       call history_out_field('TTEND_CLUBB', ttend_clubb)
       call history_out_field('THLP2_CLUBB_GW', thlp2_clubb_gw)
       call history_out_field('WPTHLP_CLUBB_GW', wpthlp_clubb_gw)
       call history_out_field('UPWP_CLUBB_GW', upwp_clubb_gw)
       call history_out_field('VPWP_CLUBB_GW', vpwp_clubb_gw)
       call history_out_field ('VORT4GW', vort4gw)
    end if
  end subroutine gw_drag_diagnostics_after_gw_movmtn_pbl_scheme_run
  subroutine gw_drag_diagnostics_after_gw_movmtn_src( &
       ncol, pver, pcnst, &
       const_props, &
       cpair, &
       cmfdt, dq, cmfdqr, &
       qc_sh, icwmr, cmfsl, cmflq, cmfmc_sh, &
       errmsg, errflg)

    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    use cam_history, only: history_out_field

    ! Input arguments
    integer,                          intent(in)  :: ncol, pver, pcnst
    type(ccpp_constituent_prop_ptr_t),intent(in)  :: const_props(:)
    real(kind_phys),                  intent(in)  :: cpair
    real(kind_phys),                  intent(in)  :: cmfdt (:,:)
    real(kind_phys),                  intent(in)  :: dq    (:,:,:)
    real(kind_phys),                  intent(in)  :: cmfdqr(:,:)
    real(kind_phys),                  intent(in)  :: qc_sh (:,:)
    real(kind_phys),                  intent(in)  :: icwmr (:,:)
    real(kind_phys),                  intent(in)  :: cmfsl (:,:)
    real(kind_phys),                  intent(in)  :: cmflq (:,:)
    real(kind_phys),                  intent(in)  :: cmfmc_sh(:,:)

    ! Output arguments
    character(len=512),               intent(out) :: errmsg
    integer,                          intent(out) :: errflg

    ! Local variables
    real(kind_phys)                 :: freqsh(ncol)            ! Fractional occurrence of gravity wave [fraction]
    integer                         :: const_wv_idx, const_cldliq_idx, const_cldice_idx
    character(len=512)              :: const_standard_name
    integer                         :: i, m

    errmsg = ''
    errflg = 0

    ! Find constituent indices for water vapor, cloud liquid, cloud ice
    const_wv_idx     = -1
    const_cldliq_idx = -1
    const_cldice_idx = -1
    const_check_loop: do m = 1, pcnst
      call const_props(m)%standard_name(const_standard_name)
      if (const_standard_name == 'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water') then
        const_wv_idx = m
      endif

      if (const_standard_name == 'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water') then
        const_cldliq_idx = m
      endif

      if (const_standard_name == 'cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water') then
        const_cldice_idx = m
      endif
    enddo const_check_loop

    if (use_gw_movmtn_pbl) then
       call history_out_field('UCELL_MOVMTN', usteer)
       call history_out_field('VCELL_MOVMTN', vsteer)
       call history_out_field('CS_MOVMTN', CS)
       call history_out_field('STEER_LEVEL_MOVMTN',steer_level)
       call history_out_field('XPWP_SRC_MOVMTN', xpwp_src )
    end if
  end subroutine gw_drag_diagnostics_after_gw_movmtn_src

  subroutine gw_drag_diagnostics_after_gw_convect_dp_scheme_run( &
       ncol, pver, pcnst, &
       const_props, &
       cpair, &
       cmfdt, dq, cmfdqr, &
       qc_sh, icwmr, cmfsl, cmflq, cmfmc_sh, &
       errmsg, errflg)

    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    use cam_history, only: history_out_field

    ! Input arguments
    integer,                          intent(in)  :: ncol, pver, pcnst
    type(ccpp_constituent_prop_ptr_t),intent(in)  :: const_props(:)
    real(kind_phys),                  intent(in)  :: cpair
    real(kind_phys),                  intent(in)  :: cmfdt (:,:)
    real(kind_phys),                  intent(in)  :: dq    (:,:,:)
    real(kind_phys),                  intent(in)  :: cmfdqr(:,:)
    real(kind_phys),                  intent(in)  :: qc_sh (:,:)
    real(kind_phys),                  intent(in)  :: icwmr (:,:)
    real(kind_phys),                  intent(in)  :: cmfsl (:,:)
    real(kind_phys),                  intent(in)  :: cmflq (:,:)
    real(kind_phys),                  intent(in)  :: cmfmc_sh(:,:)

    ! Output arguments
    character(len=512),               intent(out) :: errmsg
    integer,                          intent(out) :: errflg

    ! Local variables
    real(kind_phys)                 :: freqsh(ncol)            ! Fractional occurrence of gravity wave [fraction]
    integer                         :: const_wv_idx, const_cldliq_idx, const_cldice_idx
    character(len=512)              :: const_standard_name
    integer                         :: i, m

    errmsg = ''
    errflg = 0

    ! Find constituent indices for water vapor, cloud liquid, cloud ice
    const_wv_idx     = -1
    const_cldliq_idx = -1
    const_cldice_idx = -1
    const_check_loop: do m = 1, pcnst
      call const_props(m)%standard_name(const_standard_name)
      if (const_standard_name == 'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water') then
        const_wv_idx = m
      endif

      if (const_standard_name == 'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water') then
        const_cldliq_idx = m
      endif

      if (const_standard_name == 'cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water') then
        const_cldice_idx = m
      endif
    enddo const_check_loop

    if (use_gw_convect_dp) then
       ! Change ttgw to a temperature tendency before outputing it.
       ttgw = ttgw / cpair
       call gw_spec_outflds(beres_dp_pf, lchnk, ncol, band_mid, phase_speeds, u, v, &
            xv, yv, gwut, dttdf, dttke, tau(:,:,2:), utgw, vtgw, ttgw, &
            taucd)

       call history_out_field('NETDT', ttend_dp)
       call history_out_field('HDEPTH', hdepth/1000._kind_phys)
       call history_out_field('MAXQ0', maxq0)
    end if
  end subroutine gw_drag_diagnostics_after_gw_convect_dp_scheme_run

  subroutine gw_drag_diagnostics_after_gw_convect_sh_scheme_run( &
       ncol, pver, pcnst, &
       const_props, &
       cpair, &
       cmfdt, dq, cmfdqr, &
       qc_sh, icwmr, cmfsl, cmflq, cmfmc_sh, &
       errmsg, errflg)

    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    use cam_history, only: history_out_field

    ! Input arguments
    integer,                          intent(in)  :: ncol, pver, pcnst
    type(ccpp_constituent_prop_ptr_t),intent(in)  :: const_props(:)
    real(kind_phys),                  intent(in)  :: cpair
    real(kind_phys),                  intent(in)  :: cmfdt (:,:)
    real(kind_phys),                  intent(in)  :: dq    (:,:,:)
    real(kind_phys),                  intent(in)  :: cmfdqr(:,:)
    real(kind_phys),                  intent(in)  :: qc_sh (:,:)
    real(kind_phys),                  intent(in)  :: icwmr (:,:)
    real(kind_phys),                  intent(in)  :: cmfsl (:,:)
    real(kind_phys),                  intent(in)  :: cmflq (:,:)
    real(kind_phys),                  intent(in)  :: cmfmc_sh(:,:)

    ! Output arguments
    character(len=512),               intent(out) :: errmsg
    integer,                          intent(out) :: errflg

    ! Local variables
    real(kind_phys)                 :: freqsh(ncol)            ! Fractional occurrence of gravity wave [fraction]
    integer                         :: const_wv_idx, const_cldliq_idx, const_cldice_idx
    character(len=512)              :: const_standard_name
    integer                         :: i, m

    errmsg = ''
    errflg = 0

    ! Find constituent indices for water vapor, cloud liquid, cloud ice
    const_wv_idx     = -1
    const_cldliq_idx = -1
    const_cldice_idx = -1
    const_check_loop: do m = 1, pcnst
      call const_props(m)%standard_name(const_standard_name)
      if (const_standard_name == 'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water') then
        const_wv_idx = m
      endif

      if (const_standard_name == 'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water') then
        const_cldliq_idx = m
      endif

      if (const_standard_name == 'cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water') then
        const_cldice_idx = m
      endif
    enddo const_check_loop

    if (use_gw_convect_sh) then
       ! Change ttgw to a temperature tendency before outputing it.
       ttgw = ttgw / cpair
       call gw_spec_outflds(beres_sh_pf, lchnk, ncol, band_mid, phase_speeds, u, v, &
            xv, yv, gwut, dttdf, dttke, tau(:,:,2:), utgw, vtgw, ttgw, &
            taucd)
       call history_out_field ('SNETDT', ttend_sh)
       call history_out_field ('SHDEPTH', hdepth/1000._kind_phys)
       call history_out_field ('SMAXQ0', maxq0)
    end if
  end subroutine gw_drag_diagnostics_after_gw_convect_sh_scheme_run

  subroutine gw_drag_diagnostics_after_gw_front_scheme_run( &
       ncol, pver, pcnst, &
       const_props, &
       cpair, &
       cmfdt, dq, cmfdqr, &
       qc_sh, icwmr, cmfsl, cmflq, cmfmc_sh, &
       errmsg, errflg)

    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    use cam_history, only: history_out_field

    ! Input arguments
    integer,                          intent(in)  :: ncol, pver, pcnst
    type(ccpp_constituent_prop_ptr_t),intent(in)  :: const_props(:)
    real(kind_phys),                  intent(in)  :: cpair
    real(kind_phys),                  intent(in)  :: cmfdt (:,:)
    real(kind_phys),                  intent(in)  :: dq    (:,:,:)
    real(kind_phys),                  intent(in)  :: cmfdqr(:,:)
    real(kind_phys),                  intent(in)  :: qc_sh (:,:)
    real(kind_phys),                  intent(in)  :: icwmr (:,:)
    real(kind_phys),                  intent(in)  :: cmfsl (:,:)
    real(kind_phys),                  intent(in)  :: cmflq (:,:)
    real(kind_phys),                  intent(in)  :: cmfmc_sh(:,:)

    ! Output arguments
    character(len=512),               intent(out) :: errmsg
    integer,                          intent(out) :: errflg

    ! Local variables
    real(kind_phys)                 :: freqsh(ncol)            ! Fractional occurrence of gravity wave [fraction]
    integer                         :: const_wv_idx, const_cldliq_idx, const_cldice_idx
    character(len=512)              :: const_standard_name
    integer                         :: i, m

    errmsg = ''
    errflg = 0

    ! Find constituent indices for water vapor, cloud liquid, cloud ice
    const_wv_idx     = -1
    const_cldliq_idx = -1
    const_cldice_idx = -1
    const_check_loop: do m = 1, pcnst
      call const_props(m)%standard_name(const_standard_name)
      if (const_standard_name == 'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water') then
        const_wv_idx = m
      endif

      if (const_standard_name == 'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water') then
        const_cldliq_idx = m
      endif

      if (const_standard_name == 'cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water') then
        const_cldice_idx = m
      endif
    enddo const_check_loop

    if (use_gw_front .or. use_gw_front_igw) then
       call history_out_field ('FRONTGF', frontgf)
       call history_out_field ('FRONTGFA', frontga)
    end if
    if (use_gw_front) then
       ! Change ttgw to a temperature tendency before outputing it.
       ttgw = ttgw / cpair
       call gw_spec_outflds(cm_pf, lchnk, ncol, band_mid, phase_speeds, u, v, &
            xv, yv, gwut, dttdf, dttke, tau(:,:,2:), utgw, vtgw, ttgw, &
            taucd)
    end if

    if (use_gw_front_igw) then
       ! Change ttgw to a temperature tendency before outputing it.
       ttgw = ttgw / cpair
       call gw_spec_outflds(cm_igw_pf, lchnk, ncol, band_long, phase_speeds, u, v, &
            xv, yv, gwut, dttdf, dttke, tau(:,:,2:), utgw, vtgw, ttgw, &
            taucd)
    end if
  end subroutine gw_drag_diagnostics_after_gw_front_scheme_run

  subroutine gw_drag_diagnostics_after_gw_oro_scheme_run( &
       ncol, pver, pcnst, &
       const_props, &
       cpair, &
       cmfdt, dq, cmfdqr, &
       qc_sh, icwmr, cmfsl, cmflq, cmfmc_sh, &
       errmsg, errflg)

    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    use cam_history, only: history_out_field

    ! Input arguments
    integer,                          intent(in)  :: ncol, pver, pcnst
    type(ccpp_constituent_prop_ptr_t),intent(in)  :: const_props(:)
    real(kind_phys),                  intent(in)  :: cpair
    real(kind_phys),                  intent(in)  :: cmfdt (:,:)
    real(kind_phys),                  intent(in)  :: dq    (:,:,:)
    real(kind_phys),                  intent(in)  :: cmfdqr(:,:)
    real(kind_phys),                  intent(in)  :: qc_sh (:,:)
    real(kind_phys),                  intent(in)  :: icwmr (:,:)
    real(kind_phys),                  intent(in)  :: cmfsl (:,:)
    real(kind_phys),                  intent(in)  :: cmflq (:,:)
    real(kind_phys),                  intent(in)  :: cmfmc_sh(:,:)

    ! Output arguments
    character(len=512),               intent(out) :: errmsg
    integer,                          intent(out) :: errflg

    ! Local variables
    real(kind_phys)                 :: freqsh(ncol)            ! Fractional occurrence of gravity wave [fraction]
    integer                         :: const_wv_idx, const_cldliq_idx, const_cldice_idx
    character(len=512)              :: const_standard_name
    integer                         :: i, m

    errmsg = ''
    errflg = 0

    ! Find constituent indices for water vapor, cloud liquid, cloud ice
    const_wv_idx     = -1
    const_cldliq_idx = -1
    const_cldice_idx = -1
    const_check_loop: do m = 1, pcnst
      call const_props(m)%standard_name(const_standard_name)
      if (const_standard_name == 'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water') then
        const_wv_idx = m
      endif

      if (const_standard_name == 'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water') then
        const_cldliq_idx = m
      endif

      if (const_standard_name == 'cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water') then
        const_cldice_idx = m
      endif
    enddo const_check_loop

    if (use_gw_oro) then
       call history_out_field('TAUAORO', tau(:,0,:))
       call history_out_field('UTGWORO', utgw)
       call history_out_field('VTGWORO', vtgw)
       call history_out_field('TTGWORO', ttgw)
       call history_out_field('TTGWSDFORO', dttdf / cpair)
       call history_out_field('TTGWSKEORO', dttke / cpair)
       call history_out_field('TAUGWX', tau0x)
       call history_out_field('TAUGWY', tau0y)
    end if
  end subroutine gw_drag_diagnostics_after_gw_oro_scheme_run

!jt  subroutine gw_drag_diagnostics_after_gw_rdg_beta_scheme_run( &
  subroutine gw_drag_diagnostics_initial_state_output_run( &
       ncol, pver, pcnst, &
       const_props, &
       cpair, &
       cmfdt, dq, cmfdqr, &
       qc_sh, icwmr, cmfsl, cmflq, cmfmc_sh, &
       errmsg, errflg)

    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    use cam_history, only: history_out_field

    ! Input arguments
    integer,                          intent(in)  :: ncol, pver, pcnst
    type(ccpp_constituent_prop_ptr_t),intent(in)  :: const_props(:)
    real(kind_phys),                  intent(in)  :: cpair
    real(kind_phys),                  intent(in)  :: cmfdt (:,:)
    real(kind_phys),                  intent(in)  :: dq    (:,:,:)
    real(kind_phys),                  intent(in)  :: cmfdqr(:,:)
    real(kind_phys),                  intent(in)  :: qc_sh (:,:)
    real(kind_phys),                  intent(in)  :: icwmr (:,:)
    real(kind_phys),                  intent(in)  :: cmfsl (:,:)
    real(kind_phys),                  intent(in)  :: cmflq (:,:)
    real(kind_phys),                  intent(in)  :: cmfmc_sh(:,:)

    ! Output arguments
    character(len=512),               intent(out) :: errmsg
    integer,                          intent(out) :: errflg

    ! Local variables
    real(kind_phys)                 :: freqsh(ncol)            ! Fractional occurrence of gravity wave [fraction]
    integer                         :: const_wv_idx, const_cldliq_idx, const_cldice_idx
    character(len=512)              :: const_standard_name
    integer                         :: i, m

    errmsg = ''
    errflg = 0

    ! Find constituent indices for water vapor, cloud liquid, cloud ice
    const_wv_idx     = -1
    const_cldliq_idx = -1
    const_cldice_idx = -1
    const_check_loop: do m = 1, pcnst
       call const_props(m)%standard_name(const_standard_name)
       if (const_standard_name == 'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water') then
          const_wv_idx = m
       endif

       if (const_standard_name == 'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water') then
          const_cldliq_idx = m
       endif

       if (const_standard_name == 'cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water') then
          const_cldice_idx = m
       endif
    enddo const_check_loop

    if (use_gw_rdg_beta) then

       ! Save state at top of routine
       ! Useful for unit testing checks

       call history_out_field('UEGW', u )
       call history_out_field('VEGW', v )
       call history_out_field('TEGW', t )
       call history_out_field('ZEGW', zi )
       call history_out_field('ZMGW', zm )
    end if
  end subroutine gw_drag_diagnostics_initial_state_output_run
  subroutine gw_drag_diagnostics_after_sum_all_gw_schemes_run( &
       ncol, &
       cmfmc, cnt, cnb, p_cnt, p_cnb, &
       qc_total, qc_sh, &
       errmsg, errflg)

    use cam_history, only: history_out_field

    ! Input arguments
    integer,          intent(in)  :: ncol             ! Number of columns
    real(kind_phys),  intent(in)  :: cmfmc(:,:)       ! Total convective mass flux [kg m-2 s-1]
    integer,          intent(in)  :: cnt(:)           ! Cloud top level index [1]
    integer,          intent(in)  :: cnb(:)           ! Cloud base level index [1]
    real(kind_phys),  intent(in)  :: p_cnt(:)         ! Convective cloud top pressure [Pa]
    real(kind_phys),  intent(in)  :: p_cnb(:)         ! Convective cloud base pressure [Pa]
    real(kind_phys),  intent(in)  :: qc_total(:,:)    !  detrainment_of_cloud_liquid_water_wrt_moist_air_and_condensed_water_due_to_all_convection [kg kg-1 s-1]
    real(kind_phys),  intent(in)  :: qc_sh(:,:)       !  detrainment_of_cloud_liquid_water_wrt_moist_air_and_condensed_water_due_to_shallow_convection [kg kg-1 s-1]

    ! Output arguments
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    call history_out_field('EKGW', egwdffi_tot )
    call history_out_field('TTGW', ptend%s/cpairv(:,:))
    call history_out_field('UTGW_TOTAL', ptend%u)
    call history_out_field('VTGW_TOTAL', ptend%v)
    call history_out_field('QTGW', ptend%q(:,:,1))
    call history_out_field('CLDLIQTGW', ptend%q(:,:,ixcldliq))
    call history_out_field('CLDICETGW', ptend%q(:,:,ixcldice))
  end subroutine gw_drag_diagnostics_after_sum_all_gw_schemes_run
  subroutine gw_drag_diagnostics_after_gw_rdg_scheme_run( &
       ncol, pver, pcnst, type, &
       const_props, &
       cpair, &
       cmfdt, dq, cmfdqr, &
       qc_sh, icwmr, cmfsl, cmflq, cmfmc_sh, &
       errmsg, errflg)

    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    use cam_history, only: history_out_field

    ! Input arguments
    integer,                          intent(in)  :: ncol, pver, pcnst
    character(len=5),                 intent(in)  :: type             ! BETA or GAMMA
    type(ccpp_constituent_prop_ptr_t),intent(in)  :: const_props(:)
    real(kind_phys),                  intent(in)  :: cpair
    real(kind_phys),                  intent(in)  :: cmfdt (:,:)
    real(kind_phys),                  intent(in)  :: dq    (:,:,:)
    real(kind_phys),                  intent(in)  :: cmfdqr(:,:)
    real(kind_phys),                  intent(in)  :: qc_sh (:,:)
    real(kind_phys),                  intent(in)  :: icwmr (:,:)
    real(kind_phys),                  intent(in)  :: cmfsl (:,:)
    real(kind_phys),                  intent(in)  :: cmflq (:,:)
    real(kind_phys),                  intent(in)  :: cmfmc_sh(:,:)

    ! Output arguments
    character(len=512),               intent(out) :: errmsg
    integer,                          intent(out) :: errflg

    ! Local variables
    real(kind_phys)                 :: freqsh(ncol)            ! Fractional occurrence of gravity wave [fraction]
    integer                         :: const_wv_idx, const_cldliq_idx, const_cldice_idx
    character(len=512)              :: const_standard_name
    integer                         :: i, m

    errmsg = ''
    errflg = 0

    ! Find constituent indices for water vapor, cloud liquid, cloud ice
    const_wv_idx     = -1
    const_cldliq_idx = -1
    const_cldice_idx = -1
    const_check_loop: do m = 1, pcnst
       call history_con_siet_props(m)%standard_name(const_standard_name)
       if (const_standard_name == 'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water') then
          const_wv_idx = m
       endif

       if (const_standard_name == 'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water') then
          const_cldliq_idx = m
       endif

       if (const_standard_name == 'cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water') then
          const_cldice_idx = m
       endif
    enddo const_check_loop
    do nn = 1, n_rdg
       if (nn == 1) then
          call history_out_field('BWV_HT1', bwv)
          call history_out_field('TLB_HT1', tlb)
          call history_out_field('WBR_HT1', wbr)
          call history_out_field('TAUDSW1', taudsw)
          call history_out_field('TAUORO1', tauoro)
          call history_out_field('UBMSRC1', ubmsrc)
          call history_out_field('USRC1',   usrc)
          call history_out_field('VSRC1',   vsrc)
          call history_out_field('NSRC1'  , nsrc)
          ! Froude numbers
          call history_out_field('Fr1_DIAG' , Fr1)
          call history_out_field('Fr2_DIAG' , Fr2)
          call history_out_field('Frx_DIAG' , Frx)
          ! Ridge quantities - don't change.  Written for convenience
          call history_out_field('MXDIS1' , mxdis(:,nn))
          call history_out_field('ANGLL1' , angll(:,nn))
          call history_out_field('ANIXY1' , anixy(:,nn))
          call history_out_field('HWDTH1' , hwdth(:,nn))
          call history_out_field('CLNGT1' , clngt(:,nn))
          call history_out_field('GBXAR1' , gbxar )
          call history_out_field('TAUM1_DIAG' , tau_diag )
          call history_out_field('TAU1RDG'//trim(type)//'M', tau(:,0,:))
          call history_out_field('UBM1'//trim(type),         ubm)
          call history_out_field('UBT1RDG'//trim(type),      gwut)
       end if
       if (nn <= 6) then
          write(cn, '(i1)') nn
          call history_out_field('TAU'//cn//'RDG'//trim(type)//'X', taurx0)
          call history_out_field('TAU'//cn//'RDG'//trim(type)//'Y', taury0)
          call history_out_field('UT'//cn//'RDG'//trim(type),       utgw)
          call history_out_field('VT'//cn//'RDG'//trim(type),       vtgw)
       end if
       call history_out_field('TAUARDG'//trim(type)//'X', taurx)
       call history_out_field('TAUARDG'//trim(type)//'Y', taury)
    end do

    if (luse_gw_rdg_resid) then
       call history_out_field('TAUDIAG_RESID', tau_diag)
       call history_out_field('TAUORO_RESID', tauoro )
       call history_out_field('TAURESID'//trim(type)//'M', tau(:,0,:))
       call history_out_field('TAURESID'//trim(type)//'X', taurx)
       call history_out_field('TAURESID'//trim(type)//'Y', taury)

       call history_out_field('UBMRESID'//trim(type),      ubm)
       call history_out_field('UBIRESID'//trim(type),      ubi)
       call history_out_field('SRC_LEVEL_RESID'//trim(type),      real(src_level, kind_phys))
       ! end of residual variance calc
    end if
    if (trim(type) == 'BETA') then
       fname(1) = 'TAUGWX'
       fname(2) = 'TAUGWY'
       fname(3) = 'UTGWORO'
       fname(4) = 'VTGWORO'
    else if (trim(type) == 'GAMMA') then
       fname(1) = 'TAURDGGMX'
       fname(2) = 'TAURDGGMY'
       fname(3) = 'UTRDGGM'
       fname(4) = 'VTRDGGM'
    else
       call endrun('gw_rdg_calc: FATAL: type must be either BETA or GAMMA'&
            //' type= '//type)
    end if

    call history_out_field(fname(1), taurx(:,pver+1))
    call history_out_field(fname(2), taury(:,pver+1))
    call history_out_field(fname(3), utrdg)
    call history_out_field(fname(4), vtrdg)
    call history_out_field('TTGWORO', ttrdg / cpair)
  end subroutine gw_drag_diagnostics_after_gw_rdg_scheme_run
  subroutine gw_spec_outflds(prefix, lchnk, ncol, band, phase_speeds, u, v, xv, yv, &
       gwut, dttdf, dttke, tau, utgw, vtgw, ttgw, taucd)

    use gw_common, only: west, east, south, north

    ! One-character prefix prepended to output fields.
    character(len=1), intent(in) :: prefix
    ! Chunk and number of columns in the chunk.
    integer, intent(in) :: lchnk
    integer, intent(in) :: ncol
    ! Wave speeds.
    type(GWBand), intent(in) :: band
    ! Wave phase speeds for each column.
    real(kind_phys), intent(in) :: phase_speeds(ncol,-band%ngwv:band%ngwv)
    ! Winds at cell midpoints.
    real(kind_phys), intent(in) :: u(ncol,pver)
    real(kind_phys), intent(in) :: v(ncol,pver)
    ! Unit vector in the direction of wind at source level.
    real(kind_phys), intent(in) :: xv(ncol)
    real(kind_phys), intent(in) :: yv(ncol)
    ! Wind tendency for each wave.
    real(kind_phys), intent(in) :: gwut(ncol,pver,-band%ngwv:band%ngwv)
    ! Temperature tendencies from diffusion and kinetic energy.
    real(kind_phys) :: dttdf(ncol,pver)
    real(kind_phys) :: dttke(ncol,pver)
    ! Wave Reynolds stress.
    real(kind_phys), intent(in) :: tau(ncol,-band%ngwv:band%ngwv,pver)
    ! Zonal and meridional total wind tendencies.
    real(kind_phys), intent(in) :: utgw(ncol,pver)
    real(kind_phys), intent(in) :: vtgw(ncol,pver)
    ! Temperature tendencies.
    real(kind_phys), intent(in) :: ttgw(ncol,pver)
    ! Reynolds stress for waves propagating in each cardinal direction.
    real(kind_phys), intent(in) :: taucd(ncol,pver+1,4)

    ! Indices
    integer :: i, k, l
    integer :: ix(ncol, -band%ngwv:band%ngwv), iy(ncol, -band%ngwv:band%ngwv)
    integer :: iu(ncol), iv(ncol)

    ! Zonal wind tendency, broken up into five bins.
    real(kind_phys) :: utb(ncol, pver, 5)
    ! Definition of the bin boundaries.
    real(kind_phys), parameter :: bounds(4) = (/ -40._kind_phys, -15._kind_phys, &
         15._kind_phys, 40._kind_phys /)

    ! Momentum flux in the four cardinal directions.
    real(kind_phys) :: mf(ncol, pver, 4)

    ! Wave stress in zonal/meridional direction
    real(kind_phys) :: taux(ncol,-band%ngwv:band%ngwv,pver)
    real(kind_phys) :: tauy(ncol,-band%ngwv:band%ngwv,pver)

    ! Temporaries for output
    real(kind_phys) :: dummyx(ncol,pver)
    real(kind_phys) :: dummyy(ncol,pver)
    ! Variable names
    character(len=10) :: dumc1x, dumc1y


    ! Accumulate wind tendencies binned according to phase speed.

    utb = 0._kind_phys

    ! Find which output bin the phase speed corresponds to.
    ix = find_bin(phase_speeds)

    ! Put the wind tendency in that bin.
    do l = -band%ngwv, band%ngwv
       do k = 1, pver
          do i = 1, ncol
             utb(i,k,ix(i,l)) = utb(i,k,ix(i,l)) + gwut(i,k,l)
          end do
       end do
    end do

    ! Find just the zonal part.
    do l = 1, 5
       do k = 1, pver
          utb(:, k, l) = utb(:, k, l) * xv
       end do
    end do

    call history_out_field(trim(prefix)//'UTEND1', utb(:,:,1))
    call history_out_field(trim(prefix)//'UTEND2', utb(:,:,2))
    call history_out_field(trim(prefix)//'UTEND3', utb(:,:,3))
    call history_out_field(trim(prefix)//'UTEND4', utb(:,:,4))
    call history_out_field(trim(prefix)//'UTEND5', utb(:,:,5))

    call history_out_field(trim(prefix)//'TTGWSDF', dttdf / cpair)
    call history_out_field(trim(prefix)//'TTGWSKE', dttke / cpair)

    ! Output tau broken down into zonal and meridional components.

    taux = 0._kind_phys
    tauy = 0._kind_phys

    ! Project phase_speeds, and convert each component to a wavenumber index.
    ! These are mappings from the wavenumber index of tau to those of taux
    ! and tauy, respectively.
    do l=-band%ngwv,band%ngwv
       ix(:,l) = c_to_l(phase_speeds(:,l)*xv)
       iy(:,l) = c_to_l(phase_speeds(:,l)*yv)
    end do

    ! Find projection of tau.
    do k = 1, pver
       do l = -band%ngwv,band%ngwv
          do i = 1, ncol
             taux(i,ix(i,l),k) = taux(i,ix(i,l),k) &
                  + abs(tau(i,l,k)*xv(i))
             tauy(i,iy(i,l),k) = tauy(i,iy(i,l),k) &
                  + abs(tau(i,l,k)*yv(i))
          end do
       end do
    end do

    do l=-band%ngwv,band%ngwv

       dummyx = taux(:,l,:)
       dummyy = tauy(:,l,:)

       dumc1x = tau_fld_name(l, prefix, x_not_y=.true.)
       dumc1y = tau_fld_name(l, prefix, x_not_y=.false.)

       call history_out_field(dumc1x,dummyx)
       call history_out_field(dumc1y,dummyy)
    end do

    ! Output momentum flux in each cardinal direction.
    mf = 0._kind_phys

    do k = 1, pver

       ! Convert wind speed components to wavenumber indices.
       iu = c_to_l(u(:,k))
       iv = c_to_l(v(:,k))

       ! Sum tau components in each cardinal direction.
       ! Split west/east and north/south based on whether wave speed exceeds
       ! wind speed.
       do l = -band%ngwv, band%ngwv

          where (iu > l)
             mf(:,k,west) = mf(:,k,west) + taux(:,l,k)
          elsewhere
             mf(:,k,east) = mf(:,k,east) + taux(:,l,k)
          end where

          where (iv > l)
             mf(:,k,south) = mf(:,k,south) + tauy(:,l,k)
          elsewhere
             mf(:,k,north) = mf(:,k,north) + tauy(:,l,k)
          end where

       end do

    end do

    call history_out_field(trim(prefix)//'WMF',mf(:,:,west))
    call history_out_field(trim(prefix)//'EMF',mf(:,:,east))
    call history_out_field(trim(prefix)//'SMF',mf(:,:,south))
    call history_out_field(trim(prefix)//'NMF',mf(:,:,north))
    call history_out_field (trim(prefix)//'UTGWSPEC', utgw)
    call history_out_field (trim(prefix)//'VTGWSPEC', vtgw )
    call history_out_field (trim(prefix)//'TTGWSPEC', ttgw)
    call history_out_field (trim(prefix)//'TAUE', taucd(:,:,east))
    call history_out_field (trim(prefix)//'TAUW', taucd(:,:,west))
    call history_out_field (trim(prefix)//'TAUN', taucd(:,:,north))
    call history_out_field (trim(prefix)//'TAUS', taucd(:,:,south))
    call history_out_field (trim(prefix)//'TAUNET', taucd(:,:,east)+taucd(:,:,west))
    contains

      ! Given a value, finds which bin marked by "bounds" the value falls
      ! into.
      elemental function find_bin(val) result(idx)
        real(kind_phys), intent(in) :: val

        integer :: idx

        ! We just have to count how many bounds are exceeded.
        if (val >= 0._kind_phys) then
           idx = count(val > bounds) + 1
        else
           idx = count(val >= bounds) + 1
        end if

      end function find_bin

      ! Convert a speed to a wavenumber between -ngwv and ngwv.
      elemental function c_to_l(c) result(l)
        real(kind_phys), intent(in) :: c

        integer :: l

        l = min( max(int(c/band%dc),-band%ngwv), band%ngwv )

      end function c_to_l
      !==========================================================================

      ! Generates names for tau output across the wave spectrum (e.g.
      ! BTAUXSn01 or TAUYSp05).
      ! Probably this should use a wavenumber dimension on one field rather
      ! than creating a ton of numbered fields.
      character(len=9) pure function tau_fld_name(l, prefix, x_not_y)
        ! Wavenumber
        integer, intent(in) :: l
        ! Single-character prefix for output
        character(len=1), intent(in) :: prefix
        ! X or Y?
        logical, intent(in) :: x_not_y

        character(len=2) :: num_str

        tau_fld_name = trim(prefix)

        tau_fld_name = trim(tau_fld_name)//"TAU"

        if (x_not_y) then
           tau_fld_name = trim(tau_fld_name)//"XS"
        else
           tau_fld_name = trim(tau_fld_name)//"YS"
        end if

        if (l < 0) then
           tau_fld_name = trim(tau_fld_name)//"n"
        else
           tau_fld_name = trim(tau_fld_name)//"p"
        end if

        write(num_str,'(I2.2)') abs(l)

        tau_fld_name = trim(tau_fld_name)//num_str

      end function tau_fld_name
      !==========================================================================
    end subroutine gw_spec_outflds

    !> \section arg_table_gw_drag_diagnostics_after_shallow_scheme_run  Argument Table
    !! \htmlinclude gw_drag_diagnostics_after_shallow_scheme_run.html
    subroutine gw_drag_diagnostics_after_shallow_scheme_run( &
         ncol, pver, pcnst, &
         const_props, &
         cpair, &
         cmfdt, dq, cmfdqr, &
         qc_sh, icwmr, cmfsl, cmflq, cmfmc_sh, &
         errmsg, errflg)

      ! framework dependency for const_props
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

      use cam_history, only: history_out_field

      ! Input arguments
      integer,                          intent(in)  :: ncol, pver, pcnst
      type(ccpp_constituent_prop_ptr_t),intent(in)  :: const_props(:)
      real(kind_phys),                  intent(in)  :: cpair
      real(kind_phys),                  intent(in)  :: cmfdt (:,:)
      real(kind_phys),                  intent(in)  :: dq    (:,:,:)
      real(kind_phys),                  intent(in)  :: cmfdqr(:,:)
      real(kind_phys),                  intent(in)  :: qc_sh (:,:)
      real(kind_phys),                  intent(in)  :: icwmr (:,:)
      real(kind_phys),                  intent(in)  :: cmfsl (:,:)
      real(kind_phys),                  intent(in)  :: cmflq (:,:)
      real(kind_phys),                  intent(in)  :: cmfmc_sh(:,:)

      ! Output arguments
      character(len=512),               intent(out) :: errmsg
      integer,                          intent(out) :: errflg

      ! Local variables
      real(kind_phys)                 :: freqsh(ncol)            ! Fractional occurrence of gravity wave [fraction]
      integer                         :: const_wv_idx, const_cldliq_idx, const_cldice_idx
      character(len=512)              :: const_standard_name
      integer                         :: i, m

      errmsg = ''
      errflg = 0

      ! Find constituent indices for water vapor, cloud liquid, cloud ice
      const_wv_idx     = -1
      const_cldliq_idx = -1
      const_cldice_idx = -1
      const_check_loop: do m = 1, pcnst
         call const_props(m)%standard_name(const_standard_name)
         if (const_standard_name == 'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water') then
            const_wv_idx = m
         endif

         if (const_standard_name == 'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water') then
            const_cldliq_idx = m
         endif

         if (const_standard_name == 'cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water') then
            const_cldice_idx = m
         endif
      enddo const_check_loop

      ! CMFDT - temperature tendency
      ! Calculate temperature tendency from specific heat flux [K s-1]
      call history_out_field('CMFDT', cmfdt(:,:)/cpair)

      ! Constituent tendencies
      if (const_wv_idx > 0) then
         call history_out_field('CMFDQ',   dq(:,:,const_wv_idx))
      endif
      if (const_cldliq_idx > 0) then
         call history_out_field('CMFDLIQ', dq(:,:,const_cldliq_idx))
      endif
      if (const_cldice_idx > 0) then
         call history_out_field('CMFDICE', dq(:,:,const_cldice_idx))
      endif

      ! CMFDQR
      call history_out_field('CMFDQR', cmfdqr)

      ! Calculate fractional occurrence of gravity wave
      freqsh(:) = 0._kind_phys
      do i = 1, ncol
         if(maxval(cmfmc_sh(i,:)) > 0._kind_phys) then
            freqsh(i) = 1._kind_phys
         endif
      enddo
      call history_out_field('FREQSH', freqsh)

      call history_out_field('DQP',     qc_sh)
      call history_out_field('ICWMRSH', icwmr)
      call history_out_field('CMFSL',   cmfsl)
      call history_out_field('CMFLQ',   cmflq)
      call history_out_field('CMFMCSH', cmfmc_sh)

    end subroutine gw_drag_diagnostics_after_shallow_scheme_run

!> \section arg_table_gw_drag_diagnostics_after_sum_to_deep_run  Argument Table
!! \htmlinclude gw_drag_diagnostics_after_sum_to_deep_run.html
  subroutine gw_drag_diagnostics_after_sum_to_deep_run( &
    ncol, &
    cmfmc, cnt, cnb, p_cnt, p_cnb, &
    qc_total, qc_sh, &
    errmsg, errflg)

    use cam_history, only: history_out_field

    ! Input arguments
    integer,          intent(in)  :: ncol             ! Number of columns
    real(kind_phys),  intent(in)  :: cmfmc(:,:)       ! Total convective mass flux [kg m-2 s-1]
    integer,          intent(in)  :: cnt(:)           ! Cloud top level index [1]
    integer,          intent(in)  :: cnb(:)           ! Cloud base level index [1]
    real(kind_phys),  intent(in)  :: p_cnt(:)         ! Convective cloud top pressure [Pa]
    real(kind_phys),  intent(in)  :: p_cnb(:)         ! Convective cloud base pressure [Pa]
    real(kind_phys),  intent(in)  :: qc_total(:,:)    !  detrainment_of_cloud_liquid_water_wrt_moist_air_and_condensed_water_due_to_all_convection [kg kg-1 s-1]
    real(kind_phys),  intent(in)  :: qc_sh(:,:)       !  detrainment_of_cloud_liquid_water_wrt_moist_air_and_condensed_water_due_to_shallow_convection [kg kg-1 s-1]

    ! Output arguments
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    call history_out_field('CMFMC',   cmfmc)
    call history_out_field('CLDTOP',  real(cnt, kind_phys))
    call history_out_field('CLDBOT',  real(cnb, kind_phys))
    call history_out_field('PCLDTOP', p_cnt)
    call history_out_field('PCLDBOT', p_cnb)

    ! even though history notes this as ZMDLF
    ! (in rk_stratiform_tend)
    ! this appears to be dlf after shallow added
    call history_out_field('ZMDLF',   qc_total)
    call history_out_field('SHDLF',   qc_sh)

  end subroutine gw_drag_diagnostics_after_sum_to_deep_run

end module gw_drag_diagnostics
