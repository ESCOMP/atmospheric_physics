! Prognostic cloud water data and methods (cldwat)
! Original authors as marked in subroutines
! CCPP-ized: Haipeng Lin, January 2025
module prognostic_cloud_water
  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  ! public CCPP-compliant subroutines
  public :: prognostic_cloud_water_init
  public :: prognostic_cloud_water_run

  ! tuning parameters used by prognostic cloud water (RK stratiform)
  real(kind_phys), public :: icritc ! threshold for autoconversion of cold ice [kg kg-1]
  ! REMOVECAM: icritc does not have to be public after CAM is retired (namelist option; CCPP framework can provide it)
  real(kind_phys) :: icritw         ! threshold for autoconversion of warm ice [kg kg-1]
  real(kind_phys) :: conke          ! tunable constant for evaporation of precipitation [kg-0.5 m s-0.5]
  real(kind_phys) :: r3lcrit        ! critical radius where liq conversion begins [m]

  logical         :: do_psrhmin     ! set rh in stratosphere poleward of 50 degrees [flag]
  real(kind_phys) :: psrhmin        ! rh set in stratosphere poleward of 50 degrees [%]

  ! module private variables
  real(kind_phys) :: rhonot         ! air density at surface [g cm-3]
  real(kind_phys) :: rhos           ! assumed snow density [g cm-3]
  real(kind_phys) :: rhow           ! water density [g cm-3]
  real(kind_phys) :: rhoi           ! ice density [g cm-3]
  real(kind_phys) :: esi            ! Collection efficiency for ice by snow [1]
  real(kind_phys) :: esw            ! Collection efficiency for water by snow [1]
  real(kind_phys) :: t0             ! Approx. freezing temperature [K]
  real(kind_phys) :: cldmin         ! Assumed minimum cloud amount [1]
  real(kind_phys) :: small          ! Small number compared to unity [1]
  real(kind_phys) :: c              ! Constant for graupel-like snow [cm^(1-d)/s]

  real(kind_phys) :: d              ! Constant for graupel-like snow [1]
  real(kind_phys) :: thrpd          ! 3+d
  real(kind_phys) :: gam3pd         ! gamma(3+d)
  real(kind_phys) :: gam4pd         ! gamma(4+d)

  real(kind_phys) :: nos            ! Particles snow concentration [cm^-4]
  real(kind_phys) :: prhonos        ! pi * rhos * nos

  ! cloud microphysics constants
  real(kind_phys) :: mcon01
  real(kind_phys) :: mcon02
  real(kind_phys) :: mcon03
  real(kind_phys) :: mcon04
  real(kind_phys) :: mcon05
  real(kind_phys) :: mcon06
  real(kind_phys) :: mcon07
  real(kind_phys) :: mcon08

  ! Parameters for findmcnew
  real(kind_phys) :: capnw          ! Warm continental cloud particle number concentration [cm-3]
  real(kind_phys) :: capnc          ! Cold/oceanic cloud particle number concentration [cm-3]
  real(kind_phys) :: capnsi         ! Sea ice cloud particle number concentration [cm-3]
  real(kind_phys) :: kconst         ! Terminal velocity constant (Stokes regime) [1]
  real(kind_phys) :: effc           ! Autoconv collection efficiency [1]
  real(kind_phys) :: alpha          ! Ratio of 3rd moment radius to 2nd [1]
  real(kind_phys) :: capc           ! Autoconversion constant [1]
  real(kind_phys) :: convfw         ! Fall velocity calculation constant [1]
  real(kind_phys) :: cracw          ! Rain accreting water constant [1]
  real(kind_phys) :: critpr         ! Critical precip rate for collection efficiency [mm day-1]
  real(kind_phys) :: ciautb         ! Ice autoconversion coefficient [s-1]

contains

  ! Initialize prognostic cloud water module and calculate microphysical constants
!> \section arg_table_prognostic_cloud_water_init Argument Table
!! \htmlinclude arg_table_prognostic_cloud_water_init.html
  subroutine prognostic_cloud_water_init( &
    amIRoot, iulog, &
    tmelt, rhodair, pi, &
    icritc_in, icritw_in, &
    conke_in, r3lcrit_in, &
    do_psrhmin_in, psrhmin_in, &
    errmsg, errflg)

    ! Input arguments
    logical,            intent(in)    :: amIRoot       ! are we on the MPI root task?
    integer,            intent(in)    :: iulog         ! log output unit

    real(kind_phys),    intent(in)    :: tmelt         ! freezing_point_of_water [K]
    real(kind_phys),    intent(in)    :: rhodair       ! density_of_dry_air_at_stp [kg m-3]
    real(kind_phys),    intent(in)    :: pi

    real(kind_phys),    intent(in)    :: icritc_in
    real(kind_phys),    intent(in)    :: icritw_in
    real(kind_phys),    intent(in)    :: conke_in
    real(kind_phys),    intent(in)    :: r3lcrit_in
    logical,            intent(in)    :: do_psrhmin_in
    real(kind_phys),    intent(in)    :: psrhmin_in

    ! Output arguments
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    errmsg = ''
    errflg = 0

    ! First populate tuning parameters in-module
    icritc = icritc_in
    icritw = icritw_in
    conke  = conke_in
    r3lcrit = r3lcrit_in
    do_psrhmin = do_psrhmin_in
    psrhmin = psrhmin_in

    if(amIRoot) then
      write(iulog,*) 'tuning parameters prognostic_cloud_water_init: icritw ',icritw,' icritc ',icritc,' conke ',conke,' r3lcrit ',r3lcrit
      write(iulog,*) 'prognostic_cloud_water_init: do_psrhmin = ', do_psrhmin
    endif

    !--------------------------------------------------
    ! Initialize constants used for prognostic condensate (inimc)
    ! Original author: P. Rasch, April 1997
    !--------------------------------------------------

    rhonot = rhodair/1000.0_kind_phys     ! convert from kg m-3 to g cm-3

    ! assumed densities of snow, water, ice [g cm-3]
    rhos   = 0.1_kind_phys
    rhow   = 1._kind_phys
    rhoi   = 1._kind_phys ! unused.

    esi    = 1._kind_phys
    esw    = 0.1_kind_phys
    t0     = tmelt
    cldmin = 0.02_kind_phys
    small  = 1.e-22_kind_phys
    c      = 152.93_kind_phys
    d      = 0.25_kind_phys               ! constant for graupel like snow
                                          ! values other than 0.25 are not supported.
    if(d == 0.25_kind_phys) then
      gam3pd = 2.549256966718531_kind_phys
      gam4pd = 8.285085141835282_kind_phys

      ! on cray machines this can be calculated using gamma function:
      ! call gamma(3._kind_phys+d, signgam, gam3pd)
      ! gam3pd = sign(exp(gam3pd),signgam)
      ! call gamma(4._kind_phys+d, signgam, gam4pd)
      ! gam4pd = sign(exp(gam4pd),signgam)
      ! write(iulog,*) ' d, gamma(3+d), gamma(4+d) =', gam3pd, gam4pd
    else
      errflg = 1
      errmsg = 'prognostic_cloud_water_init: cannot use d /= 0.25'
    endif

    nos    = 3.e-2_kind_phys

    ! note: pi was originally calculated here as 4._kind_phys * atan(1.0_kind_phys)
    ! changed in ccpp-ization to use physconst value.
    prhonos = pi*rhos*nos
    thrpd  = 3._kind_phys + d

    mcon01 = pi*nos*c*gam3pd/4._kind_phys
    mcon02 = 1._kind_phys/(c*gam4pd*sqrt(rhonot)/(6*prhonos**(d/4._kind_phys)))
    mcon03 = -(0.5_kind_phys+d/4._kind_phys)
    mcon04 = 4._kind_phys/(4._kind_phys+d)
    mcon05 = (3+d)/(4+d)
    mcon06 = (3+d)/4._kind_phys
    mcon07 = mcon01*sqrt(rhonot)*mcon02**mcon05/prhonos**mcon06
    mcon08 = -0.5_kind_phys/(4._kind_phys+d)

    ! Initialize parameters used by findmcnew
    ! Cloud particle densities [cm-3] for ...
    capnw  = 400._kind_phys  ! ...warm continental
    capnc  = 150._kind_phys  ! ...cold and oceanic
    capnsi = 75._kind_phys   ! ...sea ice

    kconst = 1.18e6_kind_phys

    ! Autoconv collection efficiencies.
    ! Tripoli and Cotton (default) = 0.55
    ! Boucher 96 = 1.
    ! Baker 93 = 0.55*0.05
    ! Turn off = 0
    effc   = 0.55_kind_phys

    alpha  = 1.1_kind_phys ** 4
    ! constant for autoconversion
    capc = pi**(-.333_kind_phys)*kconst*effc *(0.75_kind_phys)**(1.333_kind_phys)*alpha

    ! critical precip rate at which we assume the collector drops can change the
    ! drop size enough to enhance the auto-conversion process (mm/day)
    critpr = 0.5_kind_phys
    convfw = 1.94_kind_phys*2.13_kind_phys*sqrt(rhow*1000._kind_phys*9.81_kind_phys*2.7e-4_kind_phys)

    ! liquid microphysics - rain accreting water constant
    ! Tripoli and Cotton = default
    ! Beheng = 6.
    cracw = .884_kind_phys*sqrt(9.81_kind_phys/(rhow*1000._kind_phys*2.7e-4_kind_phys)) ! tripoli and cotton

    ! ice microphysics autoconversion coefficient
    ciautb = 5.e-4_kind_phys

    if(amIRoot) then
      write(iulog,*) 'tuning parameters prognostic_cloud_water_init: capnw ',capnw,' capnc ',capnc,' capnsi ',capnsi,' kconst ',kconst
      write(iulog,*) 'tuning parameters prognostic_cloud_water_init: effc ',effc,' alpha ',alpha,' capc ',capc
      write(iulog,*) 'tuning parameters prognostic_cloud_water_init: critpr ',critpr,' convfw ',convfw,' cracw ',cracw,' ciautb ',ciautb
    endif

  end subroutine prognostic_cloud_water_init

  ! Calculate prognostic condensate.
  ! Cloud water parameterization, returns tendencies to water vapor, temperature,
  ! and cloud water variables.
  !
  ! Rasch, P. J., and J. E. Kristjánsson, 1998: A Comparison of the CCM3 Model Climate Using Diagnosed and Predicted Condensate Parameterizations. J. Climate, 11, 1587–1614
  ! https://doi.org/10.1175/1520-0442(1998)011<1587:ACOTCM>2.0.CO;2
  !
  ! with modifications to improve the method of determining condensation/evaporation
  ! Zhang, M., W. Lin, C. Bretherton, J. Hack, and P. J. Rasch, 2003, A modified formulation of fractional stratiform condensation rate in the NCAR Community Atmospheric Model (CAM2), J. Geophys. Res., 108(D1), 4035
  ! https://doi.org/10.1029/2002JD002523
  !
  ! Original authors: M. Zhang, W. Lin, P. Rasch and J.E. Kristjansson
  !                   B. A. Boville (latent heat of fusion)
!> \section arg_table_prognostic_cloud_water_run Argument Table
!! \htmlinclude arg_table_prognostic_cloud_water_run.html
  subroutine prognostic_cloud_water_run( &
    ncol, pver, top_lev, deltat, &
    iulog, &
    pi, gravit, rh2o, epsilo, latvap, latice, cpair, &
    dlat, &
    pmid, pdel, &
    zi, &
    troplev, &
    ttend, tn, &
    qtend, qn, &
    ltend, &
    cldice, cldliq, &
    omega, &
    cldn, &
    fice, fsnow, &
    rhdfda, rhu00, &
    landm, seaicef, snowh, &
    qme, &
    prodprec, prodsnow, &
    evapprec, evapsnow, &
    evapheat, prfzheat, meltheat, &
    precip, snowab, &
    ice2pr, liq2pr, liq2snow, &
    lsflxprc, &
    lsflxsnw, &
    pracwo, psacwo, psacio, &
    fwaut, fsaut, fracw, fsacw, fsaci, &
    errmsg, errflg)
    ! Dependencies to-be-ccppized.
    use wv_saturation,  only: qsat, estblf, svp_to_qsat, findsp_vc

    ! Input arguments
    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    integer,            intent(in)    :: top_lev        ! vertical_layer_index_of_troposphere_cloud_physics_top [index]
    real(kind_phys),    intent(in)    :: deltat         ! timestep [s]
    integer,            intent(in)    :: iulog          ! log output unit [1]
    real(kind_phys),    intent(in)    :: pi             ! pi_constant [1]
    real(kind_phys),    intent(in)    :: gravit         ! gravitational acceleration [m s-2]
    real(kind_phys),    intent(in)    :: rh2o           ! gas_constant_of_water_vapor [J K-1 kg-1]
    real(kind_phys),    intent(in)    :: epsilo         ! ratio_of_water_vapor_to_dry_air_molecular_weights [1]
    real(kind_phys),    intent(in)    :: latvap         ! latent_heat_of_vaporization_of_water_at_0c [J kg-1]
    real(kind_phys),    intent(in)    :: latice         ! latent_heat_of_fusion_of_water_at_0c [J kg-1]
    real(kind_phys),    intent(in)    :: cpair          ! specific_heat_of_dry_air_at_constant_pressure [J K-1 kg-1]
    real(kind_phys),    intent(in)    :: dlat(:)        ! latitude_degrees_north [degrees]
    real(kind_phys),    intent(in)    :: pmid(:,:)      ! air_pressure [Pa]
    real(kind_phys),    intent(in)    :: pdel(:,:)      ! air_pressure_thickness [Pa]
    real(kind_phys),    intent(in)    :: zi(:,:)        ! geopotential_height_wrt_surface_at_interface [m]
    integer,            intent(in)    :: tropLev(:)     ! tropopause_vertical_layer_index [index]

    real(kind_phys),    intent(in)    :: ttend(:,:)     ! Temperature tendency [K s-1] -- from non-micro/macrophysics
    real(kind_phys),    intent(in)    :: tn(:,:)        ! New Temperature [K]
    real(kind_phys),    intent(in)    :: qtend(:,:)     ! Water vapor tendency [kg kg-1 s-1] -- from non-micro/macrophysics
    real(kind_phys),    intent(in)    :: qn(:,:)        ! New Water vapor mixing ratio [kg kg-1]
    real(kind_phys),    intent(in)    :: ltend(:,:)     ! Cloud liquid water tendency [kg kg-1 s-1] -- from non-micro/macrophysics
    real(kind_phys),    intent(in)    :: cldice(:,:)    ! adv: cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
    real(kind_phys),    intent(in)    :: cldliq(:,:)    ! adv: cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]

    real(kind_phys),    intent(in)    :: omega(:,:)     ! lagrangian_tendency_of_air_pressure [Pa s-1]
    real(kind_phys),    intent(in)    :: cldn(:,:)      ! New Cloud fraction [fraction]
    real(kind_phys),    intent(in)    :: fice(:,:)      ! Ice fraction within cwat [fraction]
    real(kind_phys),    intent(in)    :: fsnow(:,:)     ! Fraction of rain that freezes to snow [fraction]

    real(kind_phys),    intent(in)    :: rhdfda(:,:)    ! dG(a)/da, rh=G(a), when rh>u00 [??]
    real(kind_phys),    intent(in)    :: rhu00(:,:)     ! Rhlim for cloud [percent]
    real(kind_phys),    intent(in)    :: landm(:)       ! smoothed_land_area_fraction [fraction]
    real(kind_phys),    intent(in)    :: seaicef(:)     ! Sea ice fraction [1 (fraction)]
    real(kind_phys),    intent(in)    :: snowh(:)       ! Snow depth over land, water equivalent [m]

    ! Output arguments
    real(kind_phys),    intent(out)   :: qme(:,:)       ! Rate of condensation-evaporation of condensate (net_condensation_rate_due_to_microphysics) [kg kg-1 s-1]
    real(kind_phys),    intent(out)   :: prodprec(:,:)  ! Conversion rate of condensate to precip (precipitation_production_due_to_microphysics) [kg kg-1 s-1]
    real(kind_phys),    intent(out)   :: prodsnow(:,:)  ! Snow production rate (ignored in RK?) [kg kg-1 s-1]
    real(kind_phys),    intent(out)   :: evapprec(:,:)  ! Falling precipitation evaporation rate (precipitation_evaporation_due_to_microphysics) [kg kg-1 s-1] -- & combined to apply q(wv) tendency
    real(kind_phys),    intent(out)   :: evapsnow(:,:)  ! Falling snow evaporation rate [kg kg-1 s-1]
    real(kind_phys),    intent(out)   :: evapheat(:,:)  ! heating rate due to evaporation of precipitation [J kg-1 s-1]
    real(kind_phys),    intent(out)   :: prfzheat(:,:)  ! heating rate due to freezing of precipitation [J kg-1 s-1]
    real(kind_phys),    intent(out)   :: meltheat(:,:)  ! heating rate due to snow melt [J kg-1 s-1]
    ! note -- these are precip units for atmosphere-surface exchange.
    real(kind_phys),    intent(out)   :: precip(:)      ! Precipitation rate (lwe_stratiform_precipitation_rate_at_surface) [m s-1]
    real(kind_phys),    intent(out)   :: snowab(:)      ! Snow rate (lwe_snow_precipitation_rate_at_surface_due_to_microphysics) [m s-1]
    ! / note
    real(kind_phys),    intent(out)   :: ice2pr(:,:)    ! Conversion rate of ice to precip [kg kg-1 s-1] -- apply q(cldice) tendency
    real(kind_phys),    intent(out)   :: liq2pr(:,:)    ! Conversion rate of liquid to precip [kg kg-1 s-1] -- apply q(cldliq) tendency
    real(kind_phys),    intent(out)   :: liq2snow(:,:)  ! Conversion rate of liquid to snow [kg kg-1 s-1] -- discard??
    real(kind_phys),    intent(out)   :: lsflxprc(:,:)  ! Grid-box mean, flux large scale cloud rain+snow at interfaces (stratiform_rain_and_snow_flux_at_interface) [kg m-2 s-1]
    real(kind_phys),    intent(out)   :: lsflxsnw(:,:)  ! Grid-box mean, flux large scale cloud snow at interfaces (stratiform_snow_flux_at_interface) [kg m-2 s-1]
    real(kind_phys),    intent(out)   :: pracwo(:,:)    ! accretion of cloud water by rain [s-1]
    real(kind_phys),    intent(out)   :: psacwo(:,:)    ! accretion of cloud water by snow [s-1]
    real(kind_phys),    intent(out)   :: psacio(:,:)    ! accretion of cloud ice by snow [s-1]
    real(kind_phys),    intent(out)   :: fwaut(ncol,pver)                 ! Relative importance of liquid autoconversion [fraction]
    real(kind_phys),    intent(out)   :: fsaut(ncol,pver)                 ! Relative importance of ice autoconversion [fraction]
    real(kind_phys),    intent(out)   :: fracw(ncol,pver)                 ! Relative importance of liquid collection by rain [fraction]
    real(kind_phys),    intent(out)   :: fsacw(ncol,pver)                 ! Relative importance of liquid collection by snow [fraction]
    real(kind_phys),    intent(out)   :: fsaci(ncol,pver)                 ! Relative importance of ice collection by snow [fraction]

    ! Output arguments
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    ! Local variables
    integer :: i, k, l                                  ! Iteration index [1]
    integer :: iter                                     ! # of iterations for precipitation calculation [1]
    logical :: error_found                              ! Flag for error detection [flag]

    ! Total cloud water (from cldice+cldliq)
    real(kind_phys) :: cwat(ncol,pver)                  ! Cloud water mixing ratio [kg kg-1]

    ! Precipitation and conversion rates
    real(kind_phys) :: nice2pr                          ! Rate of conversion from ice to snow [kg kg-1 s-1]
    real(kind_phys) :: nliq2pr                          ! Rate of conversion from liquid to precipitation [kg kg-1 s-1]
    real(kind_phys) :: nliq2snow                        ! Rate of conversion from liquid to snow [kg kg-1 s-1]
    real(kind_phys) :: precab(ncol)                     ! Rate of precipitation entering layer [kg m-2 s-1]
    real(kind_phys) :: prprov(ncol)                     ! Provisional precipitation at bottom of layer [kg m-2 s-1]

    ! Cloud properties
    real(kind_phys) :: cldm(ncol)                       ! Mean cloud fraction over timestep [1]
    real(kind_phys) :: cldmax(ncol)                     ! Maximum cloud fraction above current level [1]
    real(kind_phys) :: icwc(ncol)                       ! In-cloud water content [kg kg-1]
    real(kind_phys) :: cwm(ncol)                        ! Cloud water mixing ratio at midpoint of timestep [kg kg-1]
    real(kind_phys) :: cwn(ncol)                        ! Cloud water mixing ratio at end of timestep [kg kg-1]
    real(kind_phys) :: coef(ncol)                       ! Conversion timescale for condensate to rain [s-1]

    ! Thermodynamic variables
    real(kind_phys) :: t(ncol,pver)                     ! Temperature before timestep [K]
    real(kind_phys) :: tsp(ncol,pver)                   ! Saturation point temperature [K]
    real(kind_phys) :: es(ncol)                         ! Saturation vapor pressure [Pa]
    real(kind_phys) :: q(ncol,pver)                     ! Water vapor mixing ratio [kg kg-1]
    real(kind_phys) :: qs(ncol)                         ! Saturation specific humidity [kg kg-1]
    real(kind_phys) :: qsp(ncol,pver)                   ! Saturation point mixing ratio [kg kg-1]
    real(kind_phys) :: relhum(ncol)                     ! Relative humidity [1]
    real(kind_phys) :: relhum1(ncol)                    ! Updated relative humidity [1]
    real(kind_phys) :: rhu_adj(ncol,pver)               ! Adjusted relative humidity limit for strat dehydration [1]

    ! Cloud-Evaporation scheme variables
    real(kind_phys) :: calpha(ncol)                     ! Alpha term in C-E formulation [1]
    real(kind_phys) :: cbeta(ncol)                      ! Beta term in C-E formulation [1]
    real(kind_phys) :: cbetah(ncol)                     ! Beta-hat at saturation portion [1]
    real(kind_phys) :: cgamma(ncol)                     ! Gamma term in C-E formulation [1]
    real(kind_phys) :: cgamah(ncol)                     ! Gamma-hat at saturation portion [1]
    real(kind_phys) :: rcgama(ncol)                     ! Ratio of gamma to gamma-hat [1]
    real(kind_phys) :: csigma(ncol)                     ! Sigma term in C-E formulation [1]
    real(kind_phys) :: qmeres(ncol)                     ! Residual condensation after C-E and evapprec [kg kg-1 s-1]
    real(kind_phys) :: qmec1(ncol)                      ! Cloud condensation coefficient 1 C-E formulation [1]
    real(kind_phys) :: qmec2(ncol)                      ! Cloud condensation coefficient 2 C-E formulation [1]
    real(kind_phys) :: qmec3(ncol)                      ! Cloud condensation coefficient 3 C-E formulation [1]
    real(kind_phys) :: qmec4(ncol)                      ! Cloud condensation coefficient 4 C-E formulation [1]

    ! Diagnostic arrays for cloud water budget
    ! Hardcoded 2 here = number of iterations (iter below)
    real(kind_phys) :: rcwn(ncol,2,pver)               ! Cloud water ratio evolution [kg kg-1]
    real(kind_phys) :: rliq(ncol,2,pver)               ! Liquid water ratio evolution [kg kg-1]
    real(kind_phys) :: rice(ncol,2,pver)               ! Ice ratio evolution [kg kg-1]

    ! Constants and parameters
    real(kind_phys), parameter :: omsm = 0.99999_kind_phys  ! unity minus small number for rounding
    real(kind_phys) :: mincld                           ! Minimum cloud fraction [1]
    real(kind_phys) :: cpohl                            ! Ratio of specific heat to latent heat [K-1]
    real(kind_phys) :: hlocp                            ! Ratio of latent heat to specific heat [K]
    real(kind_phys) :: dto2                             ! Half timestep [s]

    ! Work variables
    real(kind_phys) :: dqsdt                            ! Change in saturation specific humidity with temperature [kg kg-1 K-1]
    real(kind_phys) :: gamma(ncol)                     ! Temperature derivative of saturation specific humidity [kg kg-1 K-1]
    real(kind_phys) :: qtl(ncol)                       ! Saturation tendency [kg kg-1 s-1]
    real(kind_phys) :: qtmp                             ! Temporary mixing ratio [kg kg-1]
    real(kind_phys) :: ttmp                             ! Temporary temperature [K]
    real(kind_phys) :: qsn                              ! Updated saturation specific humidity [kg kg-1]
    real(kind_phys) :: esn                              ! Updated saturation vapor pressure [Pa]
    real(kind_phys) :: prtmp                            ! Temporary precipitation rate [kg m-2 s-1]
    real(kind_phys) :: ctmp                             ! Temporary condensation rate [kg kg-1 s-1]
    real(kind_phys) :: cdt                              ! Timestep factor [1]
    real(kind_phys) :: wtthick                          ! Layer thickness weight [1]
    real(kind_phys) :: pol                              ! Production over loss ratio [1]

    errmsg = ''
    errflg = 0
    error_found = .false.

    cpohl  = cpair/latvap
    hlocp  = latvap/cpair
    dto2   = 0.5_kind_phys * deltat

#ifdef PERGRO
    mincld = 1.e-4_kind_phys
    iter = 1           ! number of times to iterate the precipitation calculation
#else
    mincld = 1.e-4_kind_phys
    iter = 2
#endif

    ! initialize total cloud water [kg kg-1]
    cwat(:ncol,:pver) = cldice(:ncol,:pver) + cldliq(:ncol,:pver)

    ! initialize single level and multi level fields
    do i = 1, ncol
      precip(i) = 0.0_kind_phys
      precab(i) = 0.0_kind_phys
      snowab(i) = 0.0_kind_phys
      cldmax(i) = 0.0_kind_phys
    end do

    do k = 1, pver
      do i = 1, ncol
        q(i,k) = qn(i,k)
        t(i,k) = tn(i,k)
        ! q(i,k)=qn(i,k)-qtend(i,k)*deltat
        ! t(i,k)=tn(i,k)-ttend(i,k)*deltat
      end do
    end do

    qme     (:ncol,:) = 0._kind_phys
    evapprec(:ncol,:) = 0._kind_phys
    prodprec(:ncol,:) = 0._kind_phys
    evapsnow(:ncol,:) = 0._kind_phys
    prodsnow(:ncol,:) = 0._kind_phys
    evapheat(:ncol,:) = 0._kind_phys
    meltheat(:ncol,:) = 0._kind_phys
    prfzheat(:ncol,:) = 0._kind_phys
    ice2pr(:ncol,:)   = 0._kind_phys
    liq2pr(:ncol,:)   = 0._kind_phys
    liq2snow(:ncol,:) = 0._kind_phys
    fwaut(:ncol,:)    = 0._kind_phys
    fsaut(:ncol,:)    = 0._kind_phys
    fracw(:ncol,:)    = 0._kind_phys
    fsacw(:ncol,:)    = 0._kind_phys
    fsaci(:ncol,:)    = 0._kind_phys
    lsflxprc(:ncol,:) = 0._kind_phys
    lsflxsnw(:ncol,:) = 0._kind_phys
    pracwo(:ncol,:)   = 0._kind_phys
    psacwo(:ncol,:)   = 0._kind_phys
    psacio(:ncol,:)   = 0._kind_phys

    ! reset diagnostic arrays
    rcwn(:,:,:) = 0._kind_phys
    rliq(:,:,:) = 0._kind_phys
    rice(:,:,:) = 0._kind_phys

    ! find the wet bulb temp and saturation value
    ! for the provisional t and q without condensation
    mp_level_loop: do k = top_lev, pver
      ! "True" means that ice will be taken into account.
      call findsp_vc(qn(:ncol,k), tn(:ncol,k), pmid(:ncol,k), .true., &
           tsp(:ncol,k), qsp(:ncol,k))

      call qsat(t(1:ncol,k), pmid(1:ncol,k), es(1:ncol), qs(1:ncol), ncol, gam=gamma(1:ncol))

      do i = 1, ncol
        relhum(i) = q(i,k)/qs(i)
        cldm(i) = max(cldn(i,k),mincld)

        ! the max cloud fraction above this level
        cldmax(i) = max(cldmax(i), cldm(i))

        ! define the coefficients for C - E calculation
        calpha(i) = 1.0_kind_phys/qs(i)
        cbeta (i) = q(i,k)/qs(i)**2*gamma(i)*cpohl
        cbetah(i) = 1.0_kind_phys/qs(i)*gamma(i)*cpohl
        cgamma(i) = calpha(i)+latvap*cbeta(i)/cpair
        cgamah(i) = calpha(i)+latvap*cbetah(i)/cpair
        rcgama(i) = cgamma(i)/cgamah(i)

        if(cldm(i) > mincld) then
          icwc(i) = max(0._kind_phys,cwat(i,k)/cldm(i))
        else
          icwc(i) = 0.0_kind_phys
        endif
        ! PJR the above logic give zero icwc with nonzero cwat, dont like it!
        ! PJR generates problems with csigma
        ! PJR set the icwc to a very small number, so we can start from zero cloud cover and make some clouds
        ! icwc(i) = max(1.e-8_kind_phys,cwat(i,k)/cldm(i))

        ! initial guess of evaporation, will be updated within iteration
        evapprec(i,k) = conke*(1._kind_phys - cldm(i))*sqrt(precab(i)) &
                        *(1._kind_phys - min(relhum(i),1._kind_phys))

        ! zero qmeres before iteration for each level
        qmeres(i) = 0.0_kind_phys
      end do

      do i = 1, ncol
        ! calculate the cooling due to a phase change of the rainwater from above
        if (t(i,k) >= t0) then
          meltheat(i,k) = -latice * snowab(i) * gravit/pdel(i,k)
          snowab(i) = 0._kind_phys
        else
          meltheat(i,k) = 0._kind_phys
        endif
      end do

      ! calculate qme and formation of precip.
      !
      ! The cloud microphysics is highly nonlinear and coupled with qme
      ! Both rain processes and qme are calculated iteratively.
      qme_iter_loop: do l = 1, iter
        qme_update_loop: do i = 1, ncol
          ! calculation of qme has 4 scenarios
          call relhum_min_adj(ncol, pver, tropLev, dlat, rhu00, rhu_adj)

          if(relhum(i) > rhu_adj(i,k)) then
            ! 1. whole grid saturation
            if(relhum(i) >= 0.999_kind_phys .or. cldm(i) >= 0.999_kind_phys) then
              qme(i,k) = (calpha(i)*qtend(i,k)-cbetah(i)*ttend(i,k))/cgamah(i)
            ! 2. fractional saturation
            else
              if (rhdfda(i,k) .eq. 0._kind_phys .and. icwc(i) .eq. 0._kind_phys) then
                write(iulog,*) 'prognostic_cloud_water: empty rh cloud @ ', i, k
                write(iulog,*) 'relhum, iter ', relhum(i), l, rhu_adj(i,k), cldm(i), cldn(i,k)

                errflg = 1
                errmsg = 'prognostic_cloud_water: empty RH cloud'
                return
              endif

              csigma(i) = 1.0_kind_phys/(rhdfda(i,k)+cgamma(i)*icwc(i))
              qmec1(i) = (1.0_kind_phys-cldm(i))*csigma(i)*rhdfda(i,k)
              qmec2(i) = cldm(i)*calpha(i)/cgamah(i)+(1.0_kind_phys-rcgama(i)*cldm(i))*   &
                         csigma(i)*calpha(i)*icwc(i)
              qmec3(i) = cldm(i)*cbetah(i)/cgamah(i) +  &
                       (cbeta(i)-rcgama(i)*cldm(i)*cbetah(i))*csigma(i)*icwc(i)
              qmec4(i) = csigma(i)*cgamma(i)*icwc(i)

              ! Q = C-E = -C1*Al + C2*Aq - C3*At + C4*Er
              qme(i,k) = -qmec1(i)*ltend(i,k) + qmec2(i)*qtend(i,k)  &
                         -qmec3(i)*ttend(i,k) + qmec4(i)*evapprec(i,k)
            endif
          ! 3. when rh < rhu00, evaporate existing cloud water
          else if(cwat(i,k) > 0.0_kind_phys) then
            ! liquid water should be evaporated but not to exceed
            ! saturation point. if qn > qsp, not to evaporate cwat
            qme(i,k) = -min(max(0._kind_phys,qsp(i,k)-qn(i,k)),cwat(i,k))/deltat
          ! 4. no condensation nor evaporation
          else
            qme(i,k) = 0.0_kind_phys
          endif
        end do qme_update_loop

        ! Because of the finite time step, place a bound here not to exceed wet bulb point
        ! and not to evaporate more than available water
        qme_bound_loop: do i = 1, ncol
          qtmp = qn(i,k) - qme(i,k)*deltat

          ! possibilities to have qtmp > qsp
          !
          !   1. if qn > qs(tn), it condenses;
          !      if after applying qme,  qtmp > qsp,  more condensation is applied.
          !
          !   2. if qn < qs, evaporation should not exceed qsp,
          if(qtmp > qsp(i,k)) then
            qme(i,k) = qme(i,k) + (qtmp-qsp(i,k))/deltat
          endif

          ! if net evaporation, it should not exceed available cwat
          if(qme(i,k) < -cwat(i,k)/deltat) then
             qme(i,k) = -cwat(i,k)/deltat
          endif

          ! addition of residual condensation from previous step of iteration
          qme(i,k) = qme(i,k) + qmeres(i)

          ! limit qme for roundoff errors (multiply by slightly less than unity)
          qme(i,k) = qme(i,k) * omsm
        end do qme_bound_loop

        do i = 1, ncol
          ! as a safe limit, condensation should not reduce grid mean rh below rhu00
          if(qme(i,k) > 0.0_kind_phys .and. relhum(i) > rhu_adj(i,k)) then
            qme(i,k) = min(qme(i,k), (qn(i,k)-qs(i)*rhu_adj(i,k))/deltat)
          endif

          ! initial guess for cwm (mean cloud water over time step) if 1st iteration
          if(l < 2) then
            cwm(i) = max(cwat(i,k)+qme(i,k)*dto2, 0._kind_phys)
          endif
        enddo

        ! provisional precipitation falling through model layer
        prprov_update_loop: do i = 1, ncol
           ! prprov(i) =  precab(i) + prodprec(i,k)*pdel(i,k)/gravit
           ! rain produced in this layer not too effective in collection process
           wtthick = max(0._kind_phys,min(0.5_kind_phys,((zi(i,k)-zi(i,k+1))/1000._kind_phys)**2))
           prprov(i) = precab(i) + wtthick*prodprec(i,k)*pdel(i,k)/gravit
        end do prprov_update_loop

        ! calculate conversion of condensate to precipitation by cloud microphysics
        call findmcnew( &
          ncol    = ncol, &
          pver    = pver, &
          pi      = pi,   &
          k       = k,    &
          precab  = prprov(:ncol), &
          snowab  = snowab(:ncol), &
          t       = t(:ncol,:), &
          p       = pmid(:ncol,:), &
          cwm     = cwm(:ncol), &
          cldm    = cldm(:ncol), &
          cldmax  = cldmax(:ncol), &
          fice    = fice(:ncol,k),  &
          landm   = landm(:ncol),   &
          seaicef = seaicef(:ncol), &
          snowh   = snowh(:ncol),   &
          ! below output
          coef    = coef(:ncol), &
          fwaut   = fwaut(:ncol,k), &
          fsaut   = fsaut(:ncol,k), &
          fracw   = fracw(:ncol,k), &
          fsacw   = fsacw(:ncol,k), &
          fsaci   = fsaci(:ncol,k), &
          pracwo  = pracwo(:ncol,k),&
          psacwo  = psacwo(:ncol,k),&
          psacio  = psacio(:ncol,k))

        ! calculate the precip rate
        error_found = .false.
        precip_update_loop: do i = 1, ncol
          if (cldm(i) > 0) then
            ! first predict the cloud water
            cdt = coef(i)*deltat
            if(cdt > 0.01_kind_phys) then
              pol = qme(i,k)/coef(i) ! production over loss
              cwn(i) = max(0._kind_phys,(cwat(i,k)-pol)*exp(-cdt)+ pol)
            else
              cwn(i) = max(0._kind_phys,(cwat(i,k) + qme(i,k)*deltat)/(1+cdt))
            endif

            ! now back out the tendency of net rain production
            prodprec(i,k) = max(0._kind_phys,qme(i,k)-(cwn(i)-cwat(i,k))/deltat)
          else
            prodprec(i,k) = 0.0_kind_phys
            cwn(i) = 0._kind_phys
          endif

          ! provisional calculation of conversion terms
          ice2pr(i,k) = prodprec(i,k)*(fsaut(i,k)+fsaci(i,k))
          liq2pr(i,k) = prodprec(i,k)*(fwaut(i,k)+fsacw(i,k)+fracw(i,k))

          ! revision suggested by Jim McCaa
          ! it controls the amount of snow hitting the sfc
          ! by forcing a lot of conversion of cloud liquid to snow phase
          ! it might be better done later by an explicit representation of
          ! rain accreting ice (and freezing), or by an explicit freezing of raindrops
          !
          ! old:
          ! liq2snow(i,k) = prodprec(i,k)*fsacw(i,k)
          liq2snow(i,k) = max(prodprec(i,k)*fsacw(i,k), fsnow(i,k)*liq2pr(i,k))

          ! bounds
          nice2pr = min(ice2pr(i,k),(cwat(i,k)+qme(i,k)*deltat)*fice(i,k)/deltat)
          nliq2pr = min(liq2pr(i,k),(cwat(i,k)+qme(i,k)*deltat)*(1._kind_phys-fice(i,k))/deltat)

          if (liq2pr(i,k) .ne. 0._kind_phys) then
            nliq2snow = liq2snow(i,k)*nliq2pr/liq2pr(i,k)   ! correction
          else
            nliq2snow = liq2snow(i,k)
          endif

          ! avoid roundoff problems generating negatives
          nliq2snow = nliq2snow*omsm
          nliq2pr   = nliq2pr*omsm
          nice2pr   = nice2pr*omsm

          ! final estimates of conversion to precip and snow
          prodprec(i,k) = (nliq2pr + nice2pr)
          prodsnow(i,k) = (nice2pr + nliq2snow)

          ! compute diagnostics and sanity checks
          rcwn(i,l,k) =  cwat(i,k) + (qme(i,k) - prodprec(i,k))*deltat
          rliq(i,l,k) = (cwat(i,k) + qme(i,k)*deltat)*(1._kind_phys-fice(i,k)) - nliq2pr*deltat
          rice(i,l,k) = (cwat(i,k) + qme(i,k)*deltat) * fice(i,k) - nice2pr*deltat

          ! Sanity checks
          if(abs(rcwn(i,l,k)) < 1.e-300_kind_phys) rcwn(i,l,k) = 0._kind_phys
          if(abs(rliq(i,l,k)) < 1.e-300_kind_phys) rliq(i,l,k) = 0._kind_phys
          if(abs(rice(i,l,k)) < 1.e-300_kind_phys) rice(i,l,k) = 0._kind_phys
          if(rcwn(i,l,k) < 0._kind_phys) error_found = .true.
          if(rliq(i,l,k) < 0._kind_phys) error_found = .true.
          if(rice(i,l,k) < 0._kind_phys) error_found = .true.
          if(error_found) then
            if(rcwn(i,l,k) < 0._kind_phys) then
              write(iulog,*) 'prognostic_cloud_water: prob with neg rcwn1 ', rcwn(i,l,k), cwn(i)
              write(iulog,*) ' cwat, qme*deltat, prodprec*deltat ', &
                 cwat(i,k), qme(i,k)*deltat,               &
                 prodprec(i,k)*deltat,                     &
                 (qme(i,k)-prodprec(i,k))*deltat

              errflg = 1
              errmsg = 'prognostic_cloud_water: negative rcwn1'
              return
            endif

            if (rliq(i,l,k) < 0._kind_phys) then
              write(iulog,*) 'prognostic_cloud_water: prob with neg rliq1 ', rliq(i,l,k)
              errflg = 1
              errmsg = 'prognostic_cloud_water: negative rliq1'
              return
            endif

            if (rice(i,l,k) < 0._kind_phys) then
              write(iulog,*) 'prognostic_cloud_water: prob with neg rice ', rice(i,l,k)
              errflg = 1
              errmsg = 'prognostic_cloud_water: negative rice'
              return
            endif
          endif

          ! Final version of condensate to precip terms
          liq2pr(i,k) = nliq2pr
          liq2snow(i,k) = nliq2snow
          ice2pr(i,k) = nice2pr
          cwn(i) = rcwn(i,l,k)

          ! update any remaining provisional values
          cwm(i) = (cwn(i) + cwat(i,k))*0.5_kind_phys

          ! update in cloud water
          if(cldm(i) > mincld) then
            icwc(i) = cwm(i)/cldm(i)
          else
            icwc(i) = 0.0_kind_phys
          endif
          ! PJR the above logic give zero icwc with nonzero cwat, dont like it!
          ! PJR generates problems with csigma
          ! PJR set the icwc to a very small number, so we can start from zero cloud cover and make some clouds
          ! icwc(i) = max(1.e-8_kind_phys,cwm(i)/cldm(i))
        end do precip_update_loop

        ! calculate provisional value of cloud water for
        ! evaporation of precipitate (evapprec) calculation
        qtl_update_loop: do i = 1, ncol
          qtmp = qn(i,k) - qme(i,k)*deltat
          ttmp = tn(i,k) + deltat/cpair * ( meltheat(i,k) + (latvap + latice*fice(i,k)) * qme(i,k) )
          esn = estblf(ttmp)
          qsn = svp_to_qsat(esn, pmid(i,k))
          qtl(i) = max((qsn - qtmp)/deltat,0._kind_phys)
          relhum1(i) = qtmp/qsn
        end do qtl_update_loop

        evap_update_loop: do i = 1, ncol
#ifdef PERGRO
          evapprec(i,k) = conke*(1._kind_phys - max(cldm(i),mincld))* &
                          sqrt(precab(i))*(1._kind_phys - min(relhum1(i),1._kind_phys))
#else
          evapprec(i,k) = conke*(1._kind_phys - cldm(i))*sqrt(precab(i)) &
                          *(1._kind_phys - min(relhum1(i),1._kind_phys))
#endif

          ! limit the evaporation to the amount which is entering the box
          ! or saturates the box
          prtmp = precab(i)*gravit/pdel(i,k)
          evapprec(i,k) = min(evapprec(i,k), prtmp, qtl(i))*omsm
#ifdef PERGRO
          ! zeroing needed for pert growth
          evapprec(i,k) = 0._kind_phys
#endif

          ! Partition evaporation of precipitate between rain and snow using
          ! the fraction of snow falling into the box. Determine the heating
          ! due to evaporation. Note that evaporation is positive (loss of precip,
          ! gain of vapor) and that heating is negative.
          if (evapprec(i,k) > 0._kind_phys) then
            evapsnow(i,k) = evapprec(i,k) * snowab(i) / precab(i)
            evapheat(i,k) = -latvap * evapprec(i,k) - latice * evapsnow(i,k)
          else
            evapsnow(i,k) = 0._kind_phys
            evapheat(i,k) = 0._kind_phys
          end if

          ! Account for the latent heat of fusion for liquid drops collected by falling snow
          prfzheat(i,k) = latice * liq2snow(i,k)
        end do evap_update_loop

        ! now remove the residual of any over-saturation. Normally,
        ! the oversaturated water vapor should have been removed by
        ! qme formulation plus constraints by wet bulb tsp/qsp
        ! as computed above. However, because of non-linearity,
        ! addition of (qme-evapprec) to update t and q may still cause
        ! a very small amount of over saturation. It is called a
        ! residual of over-saturation because theoretically, qme
        ! should have taken care of all of large scale condensation.
        qmeres_update_loop: do i = 1, ncol
          qtmp = qn(i,k)-(qme(i,k)-evapprec(i,k))*deltat
          ttmp = tn(i,k) + deltat/cpair * ( meltheat(i,k) + evapheat(i,k) + prfzheat(i,k)      &
                 + (latvap + latice*fice(i,k)) * qme(i,k) )

          call qsat(ttmp, pmid(i,k), esn, qsn, dqsdt=dqsdt)

          if(qtmp > qsn) then
            ! now extra condensation to bring air to just saturation
            ctmp = (qtmp-qsn)/(1._kind_phys+hlocp*dqsdt)/deltat
            qme(i,k) = qme(i,k)+ctmp
            ! save residual on qmeres to addtion to qme on entering next iteration
            ! qme exit here contain the residual but overrided if back to iteration
            qmeres(i) = ctmp
          else
            qmeres(i) = 0.0_kind_phys
          endif
       end do qmeres_update_loop
      end do qme_iter_loop  ! loop over l (iteration)

      ! precipitation
      precip_snow_loop: do i = 1, ncol
        precip(i) = precip(i) + pdel(i,k)/gravit * (prodprec(i,k) - evapprec(i,k))
        precab(i) = precab(i) + pdel(i,k)/gravit * (prodprec(i,k) - evapprec(i,k))
        if(precab(i) < 0._kind_phys) then
          precab(i) = 0._kind_phys
        endif

        ! snowab(i) = snowab(i) + pdel(i,k)/gravit * (prodprec(i,k)*fice(i,k) - evapsnow(i,k))
        snowab(i) = snowab(i) + pdel(i,k)/gravit * (prodsnow(i,k) - evapsnow(i,k))

        ! If temperature above freezing, all precip is rain flux. if temperature below freezing, all precip is snow flux.
        lsflxprc(i,k+1) = precab(i)   !! making this consistent with other precip fluxes.  prc = rain + snow
        ! lsflxprc(i,k+1) = precab(i) - snowab(i)
        lsflxsnw(i,k+1) = snowab(i)
      end do precip_snow_loop
    end do mp_level_loop    ! loop over k (level)

    ! Convert precip (lwe_stratiform_precipitation_rate_at_surface)
    ! and snowab (lwe_snow_precipitation_rate_at_surface_due_to_microphysics)
    ! from kg m-2 s-1 to m s-1 (precipitation units) for surface exchange.
    !
    ! If this conversion is removed in the future, the metadata needs to
    ! be updated.
    precip(:ncol) = precip(:ncol)/1000._kind_phys
    snowab(:ncol) = snowab(:ncol)/1000._kind_phys
  end subroutine prognostic_cloud_water_run

  ! Calculate the conversion of condensate to precipitate
  ! Rasch, P. J., and J. E. Kristjánsson, 1998: A Comparison of the CCM3 Model Climate Using Diagnosed and Predicted Condensate Parameterizations. J. Climate, 11, 1587–1614
  ! https://doi.org/10.1175/1520-0442(1998)011<1587:ACOTCM>2.0.CO;2
  !
  ! Original author: P. Rasch
  subroutine findmcnew( &
    ncol, pver, &
    pi, &
    k, &
    precab, snowab, t, p, cwm, cldm, cldmax, fice, &
    landm, seaicef, snowh, &
    coef, fwaut, fsaut, &
    fracw, fsacw, fsaci, &
    pracwo, psacwo, psacio)

    ! input arguments
    integer,         intent(in) :: ncol             ! number of atmospheric columns
    integer,         intent(in) :: pver             ! number of levels
    real(kind_phys), intent(in) :: pi               ! pi_constant [1]
    integer,         intent(in) :: k                ! level index

    real(kind_phys), intent(in) :: precab(:)        ! rate of precipitation from above [kg m-2 s-1]
    real(kind_phys), intent(in) :: t(:,:)           ! temperature [K]
    real(kind_phys), intent(in) :: p(:,:)           ! air_pressure [Pa]
    real(kind_phys), intent(in) :: cwm(:)           ! condensate mixing ratio [kg kg-1]
    real(kind_phys), intent(in) :: cldm(:)          ! cloud fraction
    real(kind_phys), intent(in) :: cldmax(:)        ! max cloud fraction above this level
    real(kind_phys), intent(in) :: fice(:)          ! fraction of cwat that is ice
    real(kind_phys), intent(in) :: landm(:)         ! Land fraction ramped over water
    real(kind_phys), intent(in) :: seaicef(:)       ! sea ice fraction
    real(kind_phys), intent(in) :: snowab(:)        ! rate of snow from above [kg m-2 s-1]
    real(kind_phys), intent(in) :: snowh(:)         ! Snow depth over land, water equivalent [m]

    ! output arguments
    real(kind_phys), intent(out) :: coef(:)         ! conversion rate [s-1]
    real(kind_phys), intent(out) :: fwaut(:)        ! relative importance of liquid autoconversion (a diagnostic)
    real(kind_phys), intent(out) :: fsaut(:)        ! relative importance of ice autoconversion    (a diagnostic)
    real(kind_phys), intent(out) :: fracw(:)        ! relative importance of rain accreting liquid (a diagnostic)
    real(kind_phys), intent(out) :: fsacw(:)        ! relative importance of snow accreting liquid (a diagnostic)
    real(kind_phys), intent(out) :: fsaci(:)        ! relative importance of snow accreting ice    (a diagnostic)
    real(kind_phys), intent(out) :: pracwo(:)       ! accretion of cloud water by rain [s-1]
    real(kind_phys), intent(out) :: psacwo(:)       ! accretion of cloud water by snow [s-1]
    real(kind_phys), intent(out) :: psacio(:)       ! accretion of cloud ice by snow [s-1]

    ! local variables
    integer :: i, ii                                ! Loop index [index]
    integer :: ncols                                ! Number of active columns for microphysics (different from ncol!!) [count]
    integer :: ind(ncol)                            ! Active column indices [index]
    real(kind_phys) :: capn                         ! Local cloud particle number concentration [cm-3]
    real(kind_phys) :: cldloc(ncol)                ! Non-zero cloud fraction [1]
    real(kind_phys) :: cldpr(ncol)                 ! Cloud fraction for precipitation [1]
    real(kind_phys) :: totmr(ncol)                 ! In-cloud total water mixing ratio [kg kg-1]
    real(kind_phys) :: icemr(ncol)                 ! In-cloud ice mixing ratio [kg kg-1]
    real(kind_phys) :: liqmr(ncol)                 ! In-cloud liquid water mixing ratio [kg kg-1]
    real(kind_phys) :: rainmr(ncol)                ! In-cloud rain mixing ratio [kg kg-1]
    real(kind_phys) :: ciaut                        ! Ice autoconversion coefficient [s-1]
    real(kind_phys) :: pracw                        ! Rate of rain collecting cloud water [s-1]
    real(kind_phys) :: psaci                        ! Rate of snow collecting cloud ice [s-1]
    real(kind_phys) :: psacw                        ! Rate of snow collecting cloud water [s-1]
    real(kind_phys) :: psaut                        ! Rate of ice autoconversion [s-1]
    real(kind_phys) :: pwaut                        ! Rate of liquid water autoconversion [s-1]
    real(kind_phys) :: ptot                         ! Total conversion rate [s-1]
    real(kind_phys) :: prlloc(ncol)                ! Local rain flux [mm day-1]
    real(kind_phys) :: prscgs(ncol)                ! Local snow amount [g cm-2]
    real(kind_phys) :: snowfr                       ! Snow fraction of precipitation [1]
    real(kind_phys) :: vfallw                       ! Fall speed of liquid precipitation [m s-1]
    real(kind_phys) :: rho(ncol)                   ! Air density [kg m-3]
    real(kind_phys) :: rhocgs                       ! Air density in CGS units [g cm-3]
    real(kind_phys) :: r3l                          ! Cloud droplet volume radius [m]
    real(kind_phys) :: icrit                        ! Ice autoconversion threshold [kg kg-1]
    real(kind_phys) :: wt                           ! Ice fraction weight [1]
    real(kind_phys) :: con1                         ! Work constant for radius calculation [m]
    real(kind_phys) :: con2                         ! Work constant for density ratios [1]
    real(kind_phys) :: csacx                        ! Constant used for snow accreting liquid or ice [??]
    real(kind_phys) :: rat1                         ! Density ratio work variable [1]
    real(kind_phys) :: rat2                         ! Mixing ratio ratio work variable [1]

    ! find all the points where we need to do the microphysics
    ! and set the output variables to zero
    ncols = 0
    do i = 1, ncol
      coef(i) = 0._kind_phys
      fwaut(i) = 0._kind_phys
      fsaut(i) = 0._kind_phys
      fracw(i) = 0._kind_phys
      fsacw(i) = 0._kind_phys
      fsaci(i) = 0._kind_phys
      liqmr(i) = 0._kind_phys
      rainmr(i) = 0._kind_phys
      if (cwm(i) > 1.e-20_kind_phys) then
        ncols = ncols + 1
        ind(ncols) = i
      endif
    end do

    do ii = 1, ncols
      ! map from active microphysics columns to physics columns
      i = ind(ii)

      ! the local cloudiness at this level
      cldloc(i) = max(cldmin,cldm(i))

      ! a weighted mean between max cloudiness above, and this layer
      cldpr(i) = max(cldmin,(cldmax(i)+cldm(i))*0.5_kind_phys)

      ! decompose the suspended condensate into
      ! an incloud liquid and ice phase component
      totmr(i) = cwm(i)/cldloc(i)
      icemr(i) = totmr(i)*fice(i)
      liqmr(i) = totmr(i)*(1._kind_phys-fice(i))

      ! density
      rho(i) = p(i,k)/(287._kind_phys*t(i,k))
      rhocgs = rho(i)*1.e-3_kind_phys     ! density in cgs units

      ! decompose the precipitate into a liquid and ice phase
      if (t(i,k) > t0) then
         vfallw = convfw/sqrt(rho(i))
         rainmr(i) = precab(i)/(rho(i)*vfallw*cldpr(i))
         snowfr = 0
      else
         snowfr = 1
         rainmr(i) = 0._kind_phys
      endif

      ! local snow amount in cgs units
      prscgs(i) = precab(i)/cldpr(i)*0.1_kind_phys*snowfr

      ! local rain amount in mm/day
      prlloc(i) = precab(i)*86400._kind_phys/cldpr(i)
    end do

    con1 = 1._kind_phys/(1.333_kind_phys*pi)**0.333_kind_phys * 0.01_kind_phys ! meters

    ! calculate the conversion terms
    do ii = 1, ncols
      i = ind(ii)
      rhocgs = rho(i)*1.e-3_kind_phys     ! density in cgs units

      ! some temperature dependent constants
      wt = fice(i)
      icrit = icritc*wt + icritw*(1-wt)

      ! jrm Reworked droplet number concentration algorithm
      ! Start with pressure-dependent value appropriate for continental air
      ! Note: reltab has a temperature dependence here
      capn = capnw + (capnc-capnw) * min(1._kind_phys,max(0._kind_phys,1.0_kind_phys-(p(i,k)-0.8_kind_phys*p(i,pver))/(0.2_kind_phys*p(i,pver))))
      ! Modify for snow depth over land
      capn = capn + (capnc-capn) * min(1.0_kind_phys,max(0.0_kind_phys,snowh(i)*10._kind_phys))
      ! Ramp between polluted value over land to clean value over ocean.
      capn = capn + (capnc-capn) * min(1.0_kind_phys,max(0.0_kind_phys,1.0_kind_phys-landm(i)))
      ! Ramp between the resultant value and a sea ice value in the presence of ice.
      capn = capn + (capnsi-capn) * min(1.0_kind_phys,max(0.0_kind_phys,seaicef(i)))
      ! end jrm

      ! useful terms in following calculations
      rat1 = rhocgs/rhow
      rat2 = liqmr(i)/capn
      con2 = (rat1*rat2)**0.333_kind_phys

      ! volume radius
      r3l = con1*con2

      ! critical threshold for autoconversion if modified for mixed phase
      ! clouds to mimic a bergeron findeisen process
      ! r3lc2 = r3lcrit*(1.-0.5*fice(i)*(1-fice(i)))
      !
      ! autoconversion of liquid
      !
      !        cwaut = 2.e-4
      !        cwaut = 1.e-3
      !        lcrit = 2.e-4
      !        lcrit = 5.e-4
      !        pwaut = max(0._kind_phys,liqmr(i)-lcrit)*cwaut
      !
      ! pwaut is following tripoli and cotton (and many others)
      ! we reduce the autoconversion below critpr, because these are regions where
      ! the drop size distribution is likely to imply much smaller collector drops than
      ! those relevant for a cloud distribution corresponding to the value of effc = 0.55
      ! suggested by cotton (see austin 1995 JAS, baker 1993)

      ! easy to follow form
      !        pwaut = capc*liqmr(i)**2*rhocgs/rhow
      !                *(liqmr(i)*rhocgs/(rhow*capn))**(.333)
      !                *heavy(r3l,r3lcrit)
      !                *max(0.10_kind_phys,min(1._kind_phys,prlloc(i)/critpr))
      ! somewhat faster form

      ! using modified Heaviside function:
      pwaut = capc*liqmr(i)**2*rat1*con2*heavymp(r3l,r3lcrit) * &
              max(0.10_kind_phys,min(1._kind_phys,prlloc(i)/critpr))

      ! not using modified Heaviside function:
      ! pwaut = capc*liqmr(i)**2*rat1*con2*heavym(r3l,r3lcrit)* &
      !         max(0.10_kind_phys,min(1._kind_phys,prlloc(i)/critpr))

      ! autoconversion of ice
      ciaut = ciautb

      ! autoconversion of ice condensate
#ifdef PERGRO
      psaut = heavyp(icemr(i),icrit)*icemr(i)*ciaut
#else
      psaut = max(0._kind_phys,icemr(i)-icrit)*ciaut
#endif

      ! collection of liquid by rain
      ! pracw = cracw*rho(i)*liqmr(i)*rainmr(i) !(beheng 1994)
      pracw = cracw*rho(i)*sqrt(rho(i))*liqmr(i)*rainmr(i) !(tripoli and cotton)
      pracwo(i) = pracw

      ! the following lines calculate the slope parameter and snow mixing ratio
      ! from the precip rate using the equations found in lin et al 83
      ! in the most natural form, but it is expensive, so after some tedious
      ! algebraic manipulation you can use the cheaper form found below
      !            vfalls = c*gam4pd/(6*lamdas**d)*sqrt(rhonot/rhocgs)
      !     $               *0.01   ! convert from cm/s to m/s
      !            snowmr(i) = snowfr*precab(i)/(rho(i)*vfalls*cldpr(i))
      !            snowmr(i) = ( prscgs(i)*mcon02 * (rhocgs**mcon03) )**mcon04
      !            lamdas = (prhonos/max(rhocgs*snowmr(i),small))**0.25
      !            csacw = mcon01*sqrt(rhonot/rhocgs)/(lamdas**thrpd)
      !
      ! coefficient for collection by snow independent of phase
      csacx = mcon07*rhocgs**mcon08*prscgs(i)**mcon05

      ! collection of liquid by snow (lin et al 1983)
      psacw = csacx*liqmr(i)*esw

#ifdef PERGRO
      ! this is necessary for pergro
      psacw = 0._kind_phys
#endif

      psacwo(i) = psacw

      ! collection of ice by snow (lin et al 1983)
      psaci = csacx*icemr(i)*esi
      psacio(i) = psaci

      ! total conversion of condensate to precipitate
      ptot = pwaut + psaut + pracw + psacw + psaci

      ! the reciprocal of cloud water amount (or zero if no cloud water)
      !  rcwm =  totmr(i)/(max(totmr(i),small)**2)

      ! turn the tendency back into a loss rate [s-1]
      if (totmr(i) > 0._kind_phys) then
         coef(i) = ptot/totmr(i)
      else
         coef(i) = 0._kind_phys
      endif

      if (ptot.gt.0._kind_phys) then
         fwaut(i) = pwaut/ptot
         fsaut(i) = psaut/ptot
         fracw(i) = pracw/ptot
         fsacw(i) = psacw/ptot
         fsaci(i) = psaci/ptot
      else
         fwaut(i) = 0._kind_phys
         fsaut(i) = 0._kind_phys
         fracw(i) = 0._kind_phys
         fsacw(i) = 0._kind_phys
         fsaci(i) = 0._kind_phys
      endif
    end do

  end subroutine findmcnew

  ! For special WACCM/CAM-chem cases with CAM4 physics:
  ! Sets rhu to a different value poleward of +/- 50 deg latitude and
  ! levels above the tropopause if polstrat_rhmin is specified
  subroutine relhum_min_adj(ncol, pver, tropLev, dlat, rhu, rhu_adj)
    integer,  intent(in)  :: ncol
    integer,  intent(in)  :: pver
    integer,  intent(in)  :: tropLev(:)
    real(kind_phys), intent(in)  :: dlat(:)        ! latitudes in degrees
    real(kind_phys), intent(in)  :: rhu(:,:)
    real(kind_phys), intent(out) :: rhu_adj(:,:)

    integer :: i, k

    rhu_adj(:,:) = rhu(:,:)
    if (.not. do_psrhmin) return

    do k = 1, pver
      do i = 1, ncol
        ! assumption that up is lower index
        if ((k .lt. tropLev(i)) .and. (abs(dlat(i)) .gt. 50._kind_phys)) then
          rhu_adj(i,k) = psrhmin
        endif
      enddo
    enddo
  end subroutine relhum_min_adj

  ! Heavyside step functions used in findmcnew
  ! Heavyside step function
  pure function heavy(a1, a2) result(h)
    ! h: 1 if a1 > a2, 0 otherwise
    real(kind_phys), intent(in) :: a1, a2
    real(kind_phys) :: h
    h = max(0.0_kind_phys, sign(1.0_kind_phys, a1-a2))
  end function heavy

  ! Modified Heavyside function with minimum value of 0.01
  pure function heavym(a1, a2) result(h)
    ! h: 1 if a1 > a2, 0.01 otherwise
    real(kind_phys), intent(in) :: a1, a2
    real(kind_phys) :: h
    h = max(0.01_kind_phys, sign(1.0_kind_phys, a1-a2))
  end function heavym

  ! Smoothed Heavyside function using ratio formulation
  pure function heavyp(a1, a2) result(h)
    ! h: Ratio of a1/(a1+a2+eps)
    real(kind_phys), intent(in) :: a1, a2
    real(kind_phys) :: h
    real(kind_phys), parameter :: eps = 1.0e-36_kind_phys
    h = a1/(a2 + a1 + eps)
  end function heavyp

  ! Modified smooth Heavyside with offset term
  pure function heavymp(a1, a2) result(h)
    ! h: Ratio (a1 + 0.01*a2)/(a1+a2+eps)
    real(kind_phys), intent(in) :: a1, a2
    real(kind_phys) :: h
    real(kind_phys), parameter :: eps = 1.0e-36_kind_phys
    h = (a1 + 0.01_kind_phys*a2)/(a2 + a1 + eps)
  end function heavymp

end module prognostic_cloud_water
