! Droplet activation by bulk aerosols.
module ndrop_bam
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: ndrop_bam_init
  public :: ndrop_bam_calc

  ! these are currently public parameters for use by CAM.
  public :: naer_all, aername, psat, ccn_name

  ! module parameters:

  ! # of supersaturations to calc ccn concentration
  integer, parameter :: psat = 6
  ! supersaturation (%) to determine ccn concentration
  real(kind_phys), parameter :: supersat(psat) = &
                                (/0.02_kind_phys, 0.05_kind_phys, 0.1_kind_phys, 0.2_kind_phys, 0.5_kind_phys, 1.0_kind_phys/)
  ! corresponding diagnostic names:
  character(len=4), parameter :: ccn_name(psat) = &
                                 (/'CCN1', 'CCN2', 'CCN3', 'CCN4', 'CCN5', 'CCN6'/)
  real(kind_phys)              :: super(psat)

  real(kind_phys), allocatable :: ccnfact(:, :)
  real(kind_phys), allocatable :: alogsig(:)       ! natl log of geometric standard dev of aerosol
  real(kind_phys), allocatable :: exp45logsig(:)
  real(kind_phys), allocatable :: argfactor(:)
  real(kind_phys), allocatable :: amcube(:)        ! cube of dry mode radius (m)
  real(kind_phys), allocatable :: smcrit(:)        ! critical supersaturation for activation
  real(kind_phys), allocatable :: lnsm(:)          ! ln(smcrit)
  real(kind_phys), allocatable :: amcubefactor(:)  ! factors for calculating mode radius
  real(kind_phys), allocatable :: smcritfactor(:)  ! factors for calculating critical supersaturation
  real(kind_phys), allocatable :: f1(:), f2(:)     ! abdul-razzak functions of width

  real(kind_phys) :: pi     ! pi
  real(kind_phys) :: aten
  real(kind_phys) :: third, sixth
  real(kind_phys) :: sq2, sqpi
  real(kind_phys) :: alogten, alog2, alog3, alogaten

  ! aerosol properties
  character(len=20), allocatable :: aername(:)
  real(kind_phys),   allocatable :: dryrad_aer(:)
  real(kind_phys),   allocatable :: density_aer(:)
  real(kind_phys),   allocatable :: hygro_aer(:)
  real(kind_phys),   allocatable :: dispersion_aer(:)
  real(kind_phys),   allocatable :: num_to_mass_aer(:)

  integer :: naer_all     ! number of aerosols affecting climate
  integer :: idxsul = -1  ! index in aerosol list for sulfate
  integer :: idxdst2 = -1 ! index in aerosol list for dust2
  integer :: idxdst3 = -1 ! index in aerosol list for dust3
  integer :: idxdst4 = -1 ! index in aerosol list for dust4

contains

  ! Initialize constants for droplet activation by bulk aerosols
  subroutine ndrop_bam_init(amIRoot, iulog, mwh2o, r_universal, tmelt, rhoh2o)

    use shr_spfn_mod, only: erf => shr_spfn_erf

    use aerosol_instances_mod, only: aerosol_instances_get_props, aerosol_instances_get_num_models
    use aerosol_properties_mod, only: aerosol_properties

    logical, intent(in) :: amIRoot
    integer, intent(in) :: iulog
    real(kind_phys), intent(in) :: mwh2o       ! molecular weight of water (kg/kmol)
    real(kind_phys), intent(in) :: r_universal  ! universal gas constant (J/K/kmol)
    real(kind_phys), intent(in) :: tmelt       ! freezing point of water (K)
    real(kind_phys), intent(in) :: rhoh2o      ! density of liquid water (kg/m3)

    integer  :: l, m, iaer, iaermod
    real(kind_phys) :: surften       ! surface tension of water w/respect to air (N/m)
    real(kind_phys) :: arg
    class(aerosol_properties), pointer :: aero_props_bam
    !-------------------------------------------------------------------------------

    ! Access the physical properties of the bulk aerosols that are affecting the climate
    ! by using the abstract aerosol properties interface.

    ! Find BAM properties object from factory
    aero_props_bam => null()
    do iaermod = 1, aerosol_instances_get_num_models()
      aero_props_bam => aerosol_instances_get_props(iaermod, 0)
      if (associated(aero_props_bam)) then
        if (aero_props_bam%model_is('BAM')) exit
      end if
      aero_props_bam => null()
    end do

    if (associated(aero_props_bam)) then
      naer_all = aero_props_bam%nbins()
    else
      naer_all = 0
    end if

    allocate ( &
      aername(naer_all), &
      dryrad_aer(naer_all), &
      density_aer(naer_all), &
      hygro_aer(naer_all), &
      dispersion_aer(naer_all), &
      num_to_mass_aer(naer_all))

    do iaer = 1, naer_all
      call aero_props_bam%get(iaer, 1, &
                              specname=aername(iaer), &
                              dryrad=dryrad_aer(iaer), &
                              density=density_aer(iaer), &
                              hygro=hygro_aer(iaer), &
                              num_to_mass_aer=num_to_mass_aer(iaer))
      dispersion_aer(iaer) = exp(aero_props_bam%alogsig(iaer))

      ! Look for sulfate and dust aerosols in this list (Bulk aerosol only)
      if (trim(aername(iaer)) == 'SULFATE') idxsul = iaer
      if (trim(aername(iaer)) == 'DUST2') idxdst2 = iaer
      if (trim(aername(iaer)) == 'DUST3') idxdst3 = iaer
      if (trim(aername(iaer)) == 'DUST4') idxdst4 = iaer

    end do

    if (amIRoot) then
      write (iulog, *) 'ndrop_bam_init: iaer, name, dryrad, density, hygro, dispersion, num_to_mass'
      do iaer = 1, naer_all
        write (iulog, *) iaer, aername(iaer), dryrad_aer(iaer), density_aer(iaer), hygro_aer(iaer), &
          dispersion_aer(iaer), num_to_mass_aer(iaer)
      end do
      if (idxsul < 1) then
        write (iulog, *) 'ndrop_bam_init: SULFATE aerosol properties NOT FOUND'
      else
        write (iulog, *) 'ndrop_bam_init: SULFATE aerosol properties FOUND at index ', idxsul
      end if
    end if

    ! set parameters for droplet activation,
    ! following abdul-razzak and ghan 2000, JGR
    third = 1._kind_phys/3._kind_phys
    sixth = 1._kind_phys/6._kind_phys
    sq2 = sqrt(2._kind_phys)
    pi = 4._kind_phys*atan(1.0_kind_phys)     ! note: this should probably use physconst pi but here for b4b
    sqpi = sqrt(pi)
    surften = 0.076_kind_phys
    aten = 2._kind_phys*mwh2o*surften/(r_universal*tmelt*rhoh2o)
    alogaten = log(aten)
    alog2 = log(2._kind_phys)
    alog3 = log(3._kind_phys)
    super(:) = 0.01_kind_phys*supersat(:)

    allocate ( &
      alogsig(naer_all), &
      exp45logsig(naer_all), &
      argfactor(naer_all), &
      f1(naer_all), &
      f2(naer_all), &
      amcubefactor(naer_all), &
      smcritfactor(naer_all), &
      amcube(naer_all), &
      smcrit(naer_all), &
      lnsm(naer_all), &
      ccnfact(psat, naer_all))

    do m = 1, naer_all

      ! Skip aerosols that don't have a dispersion defined.
      if (dispersion_aer(m) == 0._kind_phys) cycle

      alogsig(m) = log(dispersion_aer(m))
      exp45logsig(m) = exp(4.5_kind_phys*alogsig(m)*alogsig(m))
      argfactor(m) = 2._kind_phys/(3._kind_phys*sqrt(2._kind_phys)*alogsig(m))
      f1(m) = 0.5_kind_phys*exp(2.5_kind_phys*alogsig(m)*alogsig(m))
      f2(m) = 1._kind_phys + 0.25_kind_phys*alogsig(m)
      amcubefactor(m) = 3._kind_phys/(4._kind_phys*pi*exp45logsig(m)*density_aer(m))
      smcritfactor(m) = 2._kind_phys*aten*sqrt(aten/(27._kind_phys*max(1.e-10_kind_phys, hygro_aer(m))))
      amcube(m) = amcubefactor(m)/num_to_mass_aer(m)

      if (hygro_aer(m) .gt. 1.e-10_kind_phys) then
        smcrit(m) = smcritfactor(m)/sqrt(amcube(m))
      else
        smcrit(m) = 100._kind_phys
      end if
      lnsm(m) = log(smcrit(m))

      do l = 1, psat
        arg = argfactor(m)*log(smcrit(m)/super(l))
        if (arg < 2) then
          if (arg < -2) then
            ccnfact(l, m) = 1.e-6_kind_phys
          else
            ccnfact(l, m) = 1.e-6_kind_phys*0.5_kind_phys*erfc(arg)
          end if
        else
          ccnfact(l, m) = 0._kind_phys
        end if
      end do

    end do

  end subroutine ndrop_bam_init

  ! BAM droplet activation, contact freezing dust, and CCN diagnostics.
  ! Computes naer2/maerosol from aerosol state internally via get_bulk_num_and_mass.
  subroutine ndrop_bam_calc( &
    aero_state, aero_props, &
    ncol, pver, top_lev, &
    gravit, rair, tmelt, cpair, rh2o, rhoh2o, latvap, &
    rho, tair, wsub, qcld, qsmall_in, ast, numliq, deltatin, &
    npccn, nacon, ccn, naer2_diag, &
    errmsg, errflg)

    use aerosol_properties_mod, only: aerosol_properties
    use aerosol_state_mod, only: aerosol_state
    use bulk_aerosol_state_mod, only: bulk_aerosol_state

    class(aerosol_state), intent(in) :: aero_state
    class(aerosol_properties), intent(in) :: aero_props
    integer, intent(in)  :: ncol
    integer, intent(in)  :: pver
    integer, intent(in)  :: top_lev
    real(kind_phys), intent(in)  :: gravit          ! gravitational acceleration (m/s2)
    real(kind_phys), intent(in)  :: rair            ! dry air gas constant (J/K/kg)
    real(kind_phys), intent(in)  :: tmelt           ! freezing point of water (K)
    real(kind_phys), intent(in)  :: cpair           ! specific heat of dry air (J/K/kg)
    real(kind_phys), intent(in)  :: rh2o            ! water vapor gas constant (J/K/kg)
    real(kind_phys), intent(in)  :: rhoh2o          ! density of liquid water (kg/m3)
    real(kind_phys), intent(in)  :: latvap          ! latent heat of vaporization (J/kg)
    real(kind_phys), intent(in)  :: rho(:, :)          ! air density (kg/m3)
    real(kind_phys), intent(in)  :: tair(:, :)         ! temperature (K)
    real(kind_phys), intent(in)  :: wsub(:, :)         ! sub-grid vertical velocity (m/s)
    real(kind_phys), intent(in)  :: qcld(:, :)         ! cloud liquid mmr (kg/kg)
    real(kind_phys), intent(in)  :: qsmall_in          ! minimum cloud liquid threshold
    real(kind_phys), intent(in)  :: ast(:, :)          ! stratiform_cloud_area_fraction [fraction]
    real(kind_phys), intent(in)  :: numliq(:, :)       ! droplet number (#/kg)
    real(kind_phys), intent(in)  :: deltatin           ! timestep (s)
    real(kind_phys), intent(out) :: npccn(:, :)        ! droplet number tendency (#/kg/s)
    real(kind_phys), intent(out) :: nacon(:, :, :)      ! contact freezing dust (#/m3), (ncol,pver,4)
    real(kind_phys), intent(out) :: ccn(:, :, :)        ! CCN at 6 supersaturations (#/cm3), (ncol,pver,psat)
    real(kind_phys), intent(out) :: naer2_diag(:, :, :) ! aerosol number conc (#/m3), (ncol,pver,naer_all)

    character(len=512), intent(out) :: errmsg
    integer, intent(out) :: errflg

    ! minimum allowed cloud fraction
    real(kind_phys), parameter :: mincld = 0.0001_kind_phys

    ! local workspace
    real(kind_phys), allocatable :: naer2(:, :, :)      ! aerosol number concentration [m-3]
    real(kind_phys), allocatable :: maerosol(:, :, :)   ! aerosol mass conc [kg m-3]
    real(kind_phys) :: lcldm(ncol, pver)
    real(kind_phys) :: nact
    integer  :: i, k, m

    errmsg = ''
    errflg = 0

    allocate (naer2(ncol, pver, naer_all), stat=errflg, errmsg=errmsg)
    allocate (maerosol(ncol, pver, naer_all), stat=errflg, errmsg=errmsg)
    if (errflg /= 0) return

    ! Compute aerosol number and mass concentrations from aerosol state.
    ! b4b operation order: (mmr * rho) * ntm [* 2.0 for sulfate].
    select type (bam => aero_state)
    type is (bulk_aerosol_state)
      do m = 1, naer_all
        call bam%get_bulk_num_and_mass(m, ncol, rho, naer2(:, :, m), maerosol(:, :, m))
      end do
    class default
      errmsg = 'ndrop_bam_calc: aero_state must be bulk_aerosol_state'
      errflg = 1
      return
    end select

    do k = top_lev, pver
      do i = 1, ncol
        lcldm(i, k) = max(ast(i, k), mincld)
      end do
    end do

    ! Droplet activation
    npccn(:, :) = 0._kind_phys
    do k = top_lev, pver
      do i = 1, ncol
        if (naer_all > 0 .and. qcld(i, k) >= qsmall_in) then
          call activate(wsub(i, k), tair(i, k), rho(i, k), naer2(i, k, :), &
                        naer_all, naer_all, maerosol(i, k, :), &
                        gravit, rair, tmelt, cpair, rh2o, rhoh2o, latvap, &
                        nact, errmsg, errflg)
          if (errflg /= 0) return
        else
          nact = 0._kind_phys
        end if
        npccn(i, k) = (nact*lcldm(i, k) - numliq(i, k))/deltatin
      end do
    end do

    ! Contact freezing: dust number concentrations for bins 2-4
    nacon(:, :, :) = 0._kind_phys
    do k = top_lev, pver
      do i = 1, ncol
        if (tair(i, k) < 269.15_kind_phys) then
          if (idxdst2 > 0) nacon(i, k, 2) = naer2(i, k, idxdst2)
          if (idxdst3 > 0) nacon(i, k, 3) = naer2(i, k, idxdst3)
          if (idxdst4 > 0) nacon(i, k, 4) = naer2(i, k, idxdst4)
        end if
      end do
    end do

    ! CCN diagnostics
    call ccn_diag(ncol, pver, top_lev, maerosol, naer2, ccn)

    ! Copy naer2 for diagnostic output
    naer2_diag(:ncol, :, :) = naer2(:ncol, :, :)

    deallocate (naer2, maerosol)

  end subroutine ndrop_bam_calc

!===============================================================================

  subroutine activate( &
    wbar, tair, rhoair, na, pmode, &
    nmode, ma, &
    gravit, rair, tmelt, cpair, rh2o, rhoh2o, latvap, &
    nact, errmsg, errflg)

    use shr_spfn_mod, only: erf => shr_spfn_erf
    use wv_saturation, only: qsat

    ! calculates number fraction of aerosols activated as CCN
    ! assumes an internal mixture within each of up to pmode multiple aerosol modes
    ! a gaussian spectrum of updrafts can be treated.

    !      mks units

    !      Abdul-Razzak and Ghan, A parameterization of aerosol activation.
    !      2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.
    !      https://doi.org/10.1029/1999JD901161

    ! input
    integer, intent(in) :: pmode         ! dimension of modes
    integer, intent(in) :: nmode         ! number of aerosol modes
    real(kind_phys), intent(in) :: wbar          ! grid cell mean vertical velocity (m/s)
    real(kind_phys), intent(in) :: tair          ! air temperature (K)
    real(kind_phys), intent(in) :: rhoair        ! air density (kg/m3)
    real(kind_phys), intent(in) :: na(pmode)     ! aerosol number concentration (1/m3)
    real(kind_phys), intent(in) :: ma(pmode)     ! aerosol mass concentration (kg/m3)
    real(kind_phys), intent(in) :: gravit        ! gravitational acceleration (m/s2)
    real(kind_phys), intent(in) :: rair          ! dry air gas constant (J/K/kg)
    real(kind_phys), intent(in) :: tmelt         ! freezing point of water (K)
    real(kind_phys), intent(in) :: cpair         ! specific heat of dry air (J/K/kg)
    real(kind_phys), intent(in) :: rh2o          ! water vapor gas constant (J/K/kg)
    real(kind_phys), intent(in) :: rhoh2o        ! density of liquid water (kg/m3)
    real(kind_phys), intent(in) :: latvap        ! latent heat of vaporization (J/kg)

    ! output
    real(kind_phys), intent(out) :: nact         ! number fraction of aerosols activated
    character(len=512), intent(out) :: errmsg
    integer, intent(out) :: errflg

    ! local variables
    integer :: maxmodes

    real(kind_phys), allocatable :: volc(:) ! total aerosol volume  concentration (m3/m3)
    real(kind_phys), allocatable :: eta(:)
    real(kind_phys), allocatable :: smc(:)
    real(kind_phys), allocatable :: etafactor2(:)
    real(kind_phys), allocatable :: amcubeloc(:)
    real(kind_phys), allocatable :: lnsmloc(:)

    real(kind_phys), parameter   :: pref = 1013.25e2_kind_phys ! reference pressure [Pa]

    real(kind_phys) :: pres     ! pressure (Pa)
    real(kind_phys) :: diff0
    real(kind_phys) :: conduct0 ! thermal conductivity (Joule/m/sec/deg)
    real(kind_phys) :: qs       ! water vapor saturation mixing ratio
    real(kind_phys) :: dqsdt    ! change in qs with temperature
    real(kind_phys) :: gloc     ! thermodynamic function (m2/s)
    real(kind_phys) :: zeta
    real(kind_phys) :: lnsmax   ! ln(smax)
    real(kind_phys) :: alpha
    real(kind_phys) :: gammaloc
    real(kind_phys) :: beta
    real(kind_phys) :: sqrtg
    real(kind_phys) :: wnuc
    real(kind_phys) :: alw
    real(kind_phys) :: sqrtalw
    real(kind_phys) :: smax
    real(kind_phys) :: x
    real(kind_phys) :: etafactor1
    real(kind_phys) :: etafactor2max
    real(kind_phys) :: es
    integer  :: m

    errflg = 0
    errmsg = ''

    maxmodes = naer_all
    allocate ( &
      volc(maxmodes), &
      eta(maxmodes), &
      smc(maxmodes), &
      etafactor2(maxmodes), &
      amcubeloc(maxmodes), &
      lnsmloc(maxmodes), &
      stat=errflg, errmsg=errmsg)
    if (errflg /= 0) return

    if (maxmodes < pmode) then
      write (errmsg, *) 'ndrop_bam activate: maxmodes,pmode=', maxmodes, pmode
      errflg = 1
      return
    end if

    nact = 0._kind_phys

    if (nmode .eq. 1 .and. na(1) .lt. 1.e-20_kind_phys) return

    if (wbar .le. 0._kind_phys) return

    pres = rair*rhoair*tair
    diff0 = 0.211e-4_kind_phys*(pref/pres)*(tair/tmelt)**1.94_kind_phys
    conduct0 = (5.69_kind_phys + 0.017_kind_phys*(tair - tmelt))*4.186e2_kind_phys*1.e-5_kind_phys ! convert to J/m/s/deg
    call qsat(tair, pres, es, qs)
    dqsdt = latvap/(rh2o*tair*tair)*qs
    alpha = gravit*(latvap/(cpair*rh2o*tair*tair) - 1._kind_phys/(rair*tair))
    gammaloc = (1 + latvap/cpair*dqsdt)/(rhoair*qs)
    ! growth coefficent Abdul-Razzak & Ghan 1998 eqn 16
    ! should depend on mean radius of mode to account for gas kinetic effects
    gloc = 1._kind_phys/(rhoh2o/(diff0*rhoair*qs) &
                         + latvap*rhoh2o/(conduct0*tair)*(latvap/(rh2o*tair) - 1._kind_phys))
    sqrtg = sqrt(gloc)
    beta = 4._kind_phys*pi*rhoh2o*gloc*gammaloc
    etafactor2max = 1.e10_kind_phys/(alpha*wbar)**1.5_kind_phys ! this should make eta big if na is very small.

    do m = 1, nmode
      ! skip aerosols with no dispersion, since they aren't meant to be CCN
      if (dispersion_aer(m) == 0._kind_phys) then
        smc(m) = 100._kind_phys
        cycle
      end if
      ! internal mixture of aerosols
      volc(m) = ma(m)/(density_aer(m)) ! only if variable size dist
      if (volc(m) > 1.e-39_kind_phys .and. na(m) > 1.e-39_kind_phys) then
        etafactor2(m) = 1._kind_phys/(na(m)*beta*sqrtg)  !fixed or variable size dist
        ! number mode radius (m)
        amcubeloc(m) = (3._kind_phys*volc(m)/(4._kind_phys*pi*exp45logsig(m)*na(m)))  ! only if variable size dist
        smc(m) = smcrit(m) ! only for prescribed size dist

        if (hygro_aer(m) > 1.e-10_kind_phys) then   ! loop only if variable size dist
          smc(m) = 2._kind_phys*aten*sqrt(aten/(27._kind_phys*hygro_aer(m)*amcubeloc(m)))
        else
          smc(m) = 100._kind_phys
        end if
      else
        smc(m) = 1._kind_phys
        etafactor2(m) = etafactor2max ! this should make eta big if na is very small.
      end if
      lnsmloc(m) = log(smc(m)) ! only if variable size dist
    end do

    ! single  updraft
    wnuc = wbar
    alw = alpha*wnuc
    sqrtalw = sqrt(alw)
    zeta = 2._kind_phys*sqrtalw*aten/(3._kind_phys*sqrtg)
    etafactor1 = 2._kind_phys*alw*sqrtalw

    do m = 1, nmode
      ! skip aerosols with no dispersion, since they aren't meant to be CCN
      if (dispersion_aer(m) /= 0._kind_phys) eta(m) = etafactor1*etafactor2(m)
    end do

    call maxsat(zeta, eta, nmode, smc, smax)
    lnsmax = log(smax)

    nact = 0._kind_phys
    do m = 1, nmode
      ! skip aerosols with no dispersion, since they aren't meant to be CCN
      if (dispersion_aer(m) == 0._kind_phys) cycle
      x = 2*(lnsmloc(m) - lnsmax)/(3*sq2*alogsig(m))
      nact = nact + 0.5_kind_phys*(1._kind_phys - erf(x))*na(m)
    end do
    nact = nact/rhoair ! convert from #/m3 to #/kg

    deallocate ( &
      volc, &
      eta, &
      smc, &
      etafactor2, &
      amcubeloc, &
      lnsmloc)

  end subroutine activate

!===============================================================================

  subroutine ccn_diag(ncol, pver, top_lev, maerosol, naer2, ccn_out)
    use shr_spfn_mod, only: erfc => shr_spfn_erfc

    !-------------------------------------------------------------------------------
    !
    ! Compute diagnostic bulk aerosol ccn concentration
    !
    !-------------------------------------------------------------------------------

    ! Input arguments
    integer, intent(in)  :: ncol
    integer, intent(in)  :: pver
    integer, intent(in)  :: top_lev         ! top level of troposphere cloud physics
    real(kind_phys), intent(in)  :: naer2(:, :, :)    ! aerosol number concentration (1/m3)
    real(kind_phys), intent(in)  :: maerosol(:, :, :) ! aerosol mass conc (kg/m3)
    real(kind_phys), intent(out) :: ccn_out(:, :, :)  ! CCN at 6 supersaturations (#/cm3)

    ! Local variables
    integer :: i, k, l, m
    real(kind_phys) :: arg

    real(kind_phys) :: amcubesulfate(ncol)  ! cube of dry mode radius (m) of sulfate
    real(kind_phys) :: smcritsulfate(ncol)  ! critical supersatuation for activation of sulfate
    real(kind_phys) :: ccnfactsulfate
    !-------------------------------------------------------------------------------

    ccn_out(:ncol, :, :) = 0._kind_phys

    do k = top_lev, pver

      do m = 1, naer_all

        if (m == idxsul) then
          ! Lohmann treatment for sulfate has variable size distribution
          do i = 1, ncol
            if (naer2(i, k, m) > 0._kind_phys) then
              amcubesulfate(i) = amcubefactor(m)*maerosol(i, k, m)/(naer2(i, k, m))
              smcritsulfate(i) = smcritfactor(m)/sqrt(amcubesulfate(i))
            else
              smcritsulfate(i) = 1._kind_phys
            end if
          end do
        end if

        do l = 1, psat

          if (m == idxsul) then
            ! This code is modifying ccnfact for sulfate only.
            do i = 1, ncol
              arg = argfactor(m)*log(smcritsulfate(i)/super(l))
              if (arg < 2) then
                if (arg < -2) then
                  ccnfactsulfate = 1.0e-6_kind_phys
                else
                  ccnfactsulfate = 0.5e-6_kind_phys*erfc(arg)
                end if
              else
                ccnfactsulfate = 0.0_kind_phys
              end if
              ccn_out(i, k, l) = ccn_out(i, k, l) + naer2(i, k, m)*ccnfactsulfate
            end do
          else
            ! Non-sulfate species use ccnfact computed by the init routine
            ccn_out(:ncol, k, l) = ccn_out(:ncol, k, l) + naer2(:ncol, k, m)*ccnfact(l, m)
          end if

        end do   ! supersaturation
      end do      ! bulk aerosol
    end do         ! level

  end subroutine ccn_diag

!===============================================================================

  subroutine maxsat(zeta, eta, nmode, smc, smax)

    ! calculates maximum supersaturation for multiple
    ! competing aerosol modes.

    ! Abdul-Razzak and Ghan, A parameterization of aerosol activation.
    ! 2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.

    real(kind_phys), intent(in) :: zeta
    integer, intent(in) :: nmode ! number of modes
    real(kind_phys), intent(in) :: smc(:) ! critical supersaturation for number mode radius
    real(kind_phys), intent(in) :: eta(:)

    real(kind_phys), intent(out) :: smax ! maximum supersaturation

    integer :: m  ! mode index
    real(kind_phys) :: sum, g1, g2

    do m = 1, nmode
      if (zeta .gt. 1.e5_kind_phys*eta(m) .or. smc(m)*smc(m) .gt. 1.e5_kind_phys*eta(m)) then
        ! weak forcing. essentially none activated
        smax = 1.e-20_kind_phys
      else
        ! significant activation of this mode. calc activation all modes.
        go to 1
      end if
    end do

    return

1   continue

    sum = 0
    do m = 1, nmode
      if (eta(m) .gt. 1.e-20_kind_phys) then
        g1 = sqrt(zeta/eta(m))
        g1 = g1*g1*g1
        g2 = smc(m)/sqrt(eta(m) + 3*zeta)
        g2 = sqrt(g2)
        g2 = g2*g2*g2
        sum = sum + (f1(m)*g1 + f2(m)*g2)/(smc(m)*smc(m))
      else
        sum = 1.e20_kind_phys
      end if
    end do

    smax = 1._kind_phys/sqrt(sum)

  end subroutine maxsat

!===============================================================================

end module ndrop_bam
