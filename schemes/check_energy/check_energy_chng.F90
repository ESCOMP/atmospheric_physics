module check_energy_chng
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public  :: check_energy_chng_init
  public  :: check_energy_chng_timestep_init
  public  :: check_energy_chng_run

  ! Private module options.
  logical :: print_energy_errors = .false.    ! Turn on verbose output identifying columns that fail
                                              ! energy/water checks?

contains

!> \section arg_table_check_energy_chng_init Argument Table
!! \htmlinclude arg_table_check_energy_chng_init.html
  subroutine check_energy_chng_init(print_energy_errors_in)
    ! Input arguments
    logical,            intent(in)    :: print_energy_errors_in

    print_energy_errors = print_energy_errors_in
  end subroutine check_energy_chng_init

  ! Compute initial values of energy and water integrals,
  ! and zero out cumulative boundary tendencies.
!> \section arg_table_check_energy_chng_timestep_init Argument Table
!! \htmlinclude arg_table_check_energy_chng_timestep_init.html
  subroutine check_energy_chng_timestep_init( &
       ncol, pver, pcnst, &
       is_first_timestep, &
       q, pdel, &
       u, v, T, &
       pintdry, phis, zm, &
       cp_phys, &              ! cpairv generally, cpair fixed size for subcolumns code
       cp_or_cv_dycore, &
       te_ini_phys, te_ini_dyn, &
       tw_ini, &
       te_cur_phys, te_cur_dyn, &
       tw_cur, &
       tend_te_tnd, tend_tw_tnd, &
       temp_ini, z_ini, &
       count, &
       teout, &
       energy_formula_physics, energy_formula_dycore, &
       errmsg, errflg)

    ! This scheme is non-portable due to dependencies on cam_thermo
    ! for hydrostatic energy calculation (physics and dycore formulas)
    use cam_thermo,         only: get_hydrostatic_energy
    use cam_thermo_formula, only: ENERGY_FORMULA_DYCORE_SE, ENERGY_FORMULA_DYCORE_MPAS

    ! Input arguments
    integer,            intent(in)    :: ncol           ! number of atmospheric columns
    integer,            intent(in)    :: pver           ! number of vertical layers
    integer,            intent(in)    :: pcnst          ! number of ccpp constituents
    logical,            intent(in)    :: is_first_timestep ! is first step of initial run?
    real(kind_phys),    intent(in)    :: q(:,:,:)       ! constituent mass mixing ratios [kg kg-1]
    real(kind_phys),    intent(in)    :: pdel(:,:)      ! layer thickness [Pa]
    real(kind_phys),    intent(in)    :: u(:,:)         ! zonal wind [m s-1]
    real(kind_phys),    intent(in)    :: v(:,:)         ! meridional wind [m s-1]
    real(kind_phys),    intent(in)    :: T(:,:)         ! temperature [K]
    real(kind_phys),    intent(in)    :: pintdry(:,:)   ! interface pressure dry [Pa]
    real(kind_phys),    intent(in)    :: phis(:)        ! surface geopotential [m2 s-2]
    real(kind_phys),    intent(in)    :: zm(:,:)        ! geopotential height at layer midpoints [m]
    real(kind_phys),    intent(in)    :: cp_phys(:,:)   ! enthalpy (cpairv generally) [J kg-1 K-1]
    real(kind_phys),    intent(in)    :: cp_or_cv_dycore(:,:)  ! enthalpy or heat capacity, dycore dependent [J K-1 kg-1]
    integer,            intent(in)    :: energy_formula_physics! total energy formulation physics
    integer,            intent(in)    :: energy_formula_dycore ! total energy formulation dycore

    ! Output arguments
    real(kind_phys),    intent(out)   :: temp_ini(:,:)  ! initial temperature [K]
    real(kind_phys),    intent(out)   :: z_ini(:,:)     ! initial geopotential height [m]
    integer,            intent(out)   :: count          ! count of values with significant energy or water imbalances [1]
    real(kind_phys),    intent(out)   :: teout(:)       ! total energy for global fixer in next timestep [J m-2]
    real(kind_phys),    intent(out)   :: tend_te_tnd(:) ! total energy tendency [J m-2 s-1]
    real(kind_phys),    intent(out)   :: tend_tw_tnd(:) ! total water tendency [kg m-2 s-1]

    ! Input/Output arguments
    real(kind_phys),    intent(inout) :: te_ini_phys(:) ! physics formula: initial total energy [J m-2]
    real(kind_phys),    intent(inout) :: te_ini_dyn (:) ! dycore  formula: initial total energy [J m-2]
    real(kind_phys),    intent(inout) :: tw_ini     (:) ! initial total water [kg m-2]
    real(kind_phys),    intent(inout) :: te_cur_phys(:) ! physics formula: current total energy [J m-2]
    real(kind_phys),    intent(inout) :: te_cur_dyn (:) ! dycore  formula: current total energy [J m-2]
    real(kind_phys),    intent(inout) :: tw_cur     (:) ! current total water [kg m-2]

    ! Output arguments
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    !------------------------------------------------
    ! Physics total energy.
    !------------------------------------------------
    call get_hydrostatic_energy(                           &
        tracer             = q(1:ncol,1:pver,1:pcnst),     & ! moist mixing ratios
        moist_mixing_ratio = .true.,                       &
        pdel_in            = pdel       (1:ncol,1:pver),   &
        cp_or_cv           = cp_phys    (1:ncol,1:pver),   &
        U                  = u          (1:ncol,1:pver),   &
        V                  = v          (1:ncol,1:pver),   &
        T                  = T          (1:ncol,1:pver),   &
        vcoord             = energy_formula_physics,       & ! energy formula for physics
        ptop               = pintdry    (1:ncol,1),        &
        phis               = phis       (1:ncol),          &
        te                 = te_ini_phys(1:ncol),          & ! vertically integrated total energy
        H2O                = tw_ini     (1:ncol)           & ! v.i. total water
    )

    ! Save initial state temperature and geopotential height for use in run phase
    temp_ini(:ncol,:) = T (:ncol, :)
    z_ini   (:ncol,:) = zm(:ncol, :)

    !------------------------------------------------
    ! Dynamical core total energy.
    !------------------------------------------------
    if (energy_formula_dycore == ENERGY_FORMULA_DYCORE_SE) then
      ! SE dycore specific hydrostatic energy (enthalpy)
      call get_hydrostatic_energy(                               &
          tracer             = q(1:ncol,1:pver,1:pcnst),         & ! moist mixing ratios
          moist_mixing_ratio = .true.,                           &
          pdel_in            = pdel           (1:ncol,1:pver),   &
          cp_or_cv           = cp_or_cv_dycore(1:ncol,1:pver),   &
          U                  = u              (1:ncol,1:pver),   &
          V                  = v              (1:ncol,1:pver),   &
          T                  = T              (1:ncol,1:pver),   &
          vcoord             = energy_formula_dycore,            & ! energy formula for dycore
          ptop               = pintdry        (1:ncol,1),        &
          phis               = phis           (1:ncol),          &
          te                 = te_ini_dyn     (1:ncol)           & ! WRITE OPERATION - vertically integrated total energy
      )

    else if (energy_formula_dycore == ENERGY_FORMULA_DYCORE_MPAS) then
      ! MPAS dycore: compute cv if vertical coordinate is height: cv = cp - R (internal energy)
      call get_hydrostatic_energy(                               &
          tracer             = q(1:ncol,1:pver,1:pcnst),         & ! moist mixing ratios
          moist_mixing_ratio = .true.,                           &
          pdel_in            = pdel           (1:ncol,1:pver),   &
          cp_or_cv           = cp_or_cv_dycore(1:ncol,1:pver),   &
          U                  = u              (1:ncol,1:pver),   &
          V                  = v              (1:ncol,1:pver),   &
          T                  = T              (1:ncol,1:pver),   & ! enthalpy-scaled temperature for energy consistency
          vcoord             = energy_formula_dycore,            & ! energy formula for dycore
          ptop               = pintdry        (1:ncol,1),        &
          phis               = phis           (1:ncol),          &
          z_mid              = z_ini          (1:ncol,:),        & ! unique for MPAS
          te                 = te_ini_dyn     (1:ncol)           & ! WRITE OPERATION - vertically integrated total energy
      )
    else
      ! FV dycore: dycore energy is the same as physics
      te_ini_dyn(:ncol) = te_ini_phys(:ncol)
    endif

    ! Set current state to be the same as initial
    te_cur_phys(:ncol) = te_ini_phys(:ncol)
    te_cur_dyn (:ncol) = te_ini_dyn (:ncol)
    tw_cur     (:ncol) = tw_ini     (:ncol)

    ! Zero out current energy unbalances count
    count = 0

    ! Zero out cumulative boundary fluxes
    tend_te_tnd(:ncol) = 0._kind_phys
    tend_tw_tnd(:ncol) = 0._kind_phys

    ! If first timestep, initialize value of teout
    if(is_first_timestep) then
      teout(:ncol) = te_ini_dyn(:ncol)
    endif

  end subroutine check_energy_chng_timestep_init


  ! Check that energy and water change matches the boundary fluxes.
!> \section arg_table_check_energy_chng_run Argument Table
!! \htmlinclude arg_table_check_energy_chng_run.html
  subroutine check_energy_chng_run( &
       ncol, pver, pcnst, &
       iulog, &
       q, pdel, &
       u, v, T, &
       pintdry, phis, zm, &
       cp_phys, &              ! cpairv generally, cpair fixed size for subcolumns code
       cp_or_cv_dycore, &
       scaling_dycore, &       ! From check_energy_scaling to work around subcolumns code
       te_cur_phys, te_cur_dyn, &
       tw_cur, &
       tend_te_tnd, tend_tw_tnd, &
       temp_ini, z_ini, &
       count, ztodt, &
       latice, latvap, &
       energy_formula_physics, energy_formula_dycore, &
       name, flx_vap, flx_cnd, flx_ice, flx_sen, &
       errmsg, errflg)

    ! This scheme is non-portable due to dependencies on cam_thermo
    ! for hydrostatic energy calculation (physics and dycore formulas)
    use cam_thermo,         only: get_hydrostatic_energy

    ! Dependency for energy formula used by physics and dynamical cores
    use cam_thermo_formula, only: ENERGY_FORMULA_DYCORE_FV, ENERGY_FORMULA_DYCORE_SE, ENERGY_FORMULA_DYCORE_MPAS

    ! Input arguments
    integer,            intent(in)    :: ncol           ! number of atmospheric columns
    integer,            intent(in)    :: pver           ! number of vertical layers
    integer,            intent(in)    :: pcnst          ! number of ccpp constituents
    real(kind_phys),    intent(in)    :: q(:,:,:)       ! constituent mass mixing ratios [kg kg-1]
    real(kind_phys),    intent(in)    :: pdel(:,:)      ! layer thickness [Pa]
    real(kind_phys),    intent(in)    :: u(:,:)         ! zonal wind [m s-1]
    real(kind_phys),    intent(in)    :: v(:,:)         ! meridional wind [m s-1]
    real(kind_phys),    intent(in)    :: T(:,:)         ! temperature [K]
    real(kind_phys),    intent(in)    :: pintdry(:,:)   ! interface pressure dry [Pa]
    real(kind_phys),    intent(in)    :: phis(:)        ! surface geopotential [m2 s-2]
    real(kind_phys),    intent(in)    :: zm(:,:)        ! geopotential height at layer midpoints [m]
    real(kind_phys),    intent(in)    :: temp_ini(:,:)  ! initial temperature [K]
    real(kind_phys),    intent(in)    :: z_ini(:,:)     ! initial geopotential height [m]
    real(kind_phys),    intent(in)    :: cp_phys(:,:)   ! enthalpy (cpairv generally) [J kg-1 K-1]
    real(kind_phys),    intent(in)    :: cp_or_cv_dycore(:,:)  ! enthalpy or heat capacity, dycore dependent [J K-1 kg-1]
    real(kind_phys),    intent(in)    :: scaling_dycore(:,:)   ! scaling for conversion of temperature increment [1]
    real(kind_phys),    intent(in)    :: ztodt          ! physics timestep [s]
    real(kind_phys),    intent(in)    :: latice         ! constant, latent heat of fusion of water at 0 C [J kg-1]
    real(kind_phys),    intent(in)    :: latvap         ! constant, latent heat of vaporization of water at 0 C [J kg-1]
    integer,            intent(in)    :: energy_formula_physics! total energy formulation physics
    integer,            intent(in)    :: energy_formula_dycore ! total energy formulation dycore

    ! Input from CCPP-scheme being checked:
    ! parameterization name; surface fluxes of (1) vapor, (2) liquid+ice, (3) ice, (4) sensible heat
    ! to pass in the values to be checked, call check_energy_zero_input_fluxes to reset these values
    ! before a parameterization that is checked, then update these values as-needed
    ! (can be all zero; in fact, most parameterizations calling _chng call with zero arguments)
    !
    ! Original comment from BAB:
    ! Note that the precip and ice fluxes are in precip units (m/s).
    ! I would prefer to have kg/m2/s.
    ! I would also prefer liquid (not total) and ice fluxes
    character(len=*),   intent(in)    :: name           ! parameterization name for fluxes
    real(kind_phys),    intent(in)    :: flx_vap(:)     ! boundary flux of vapor [kg m-2 s-1]
    real(kind_phys),    intent(in)    :: flx_cnd(:)     ! boundary flux of liquid+ice (precip?) [m s-1]
    real(kind_phys),    intent(in)    :: flx_ice(:)     ! boundary flux of ice (snow?) [m s-1]
    real(kind_phys),    intent(in)    :: flx_sen(:)     ! boundary flux of sensible heat [W m-2]

    ! Input/Output arguments
    real(kind_phys),    intent(inout) :: te_cur_phys(:) ! physics formula: current total energy [J m-2]
    real(kind_phys),    intent(inout) :: te_cur_dyn (:) ! dycore  formula: current total energy [J m-2]
    real(kind_phys),    intent(inout) :: tw_cur     (:) ! current total water [kg m-2]
    integer,            intent(inout) :: count          ! count of columns with significant energy or water imbalances [1]
    real(kind_phys),    intent(inout) :: tend_te_tnd(:) ! total energy tendency [J m-2 s-1]
    real(kind_phys),    intent(inout) :: tend_tw_tnd(:) ! total water tendency [kg m-2 s-1]

    ! Output arguments
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    ! Local variables
    real(kind_phys) :: te_xpd(ncol)                     ! expected value (f0 + dt*boundary_flux)
    real(kind_phys) :: te_dif(ncol)                     ! energy of input state - original energy
    real(kind_phys) :: te_tnd(ncol)                     ! tendency from last process
    real(kind_phys) :: te_rer(ncol)                     ! relative error in energy column

    real(kind_phys) :: tw_xpd(ncol)                     ! expected value (w0 + dt*boundary_flux)
    real(kind_phys) :: tw_dif(ncol)                     ! tw_inp - original water
    real(kind_phys) :: tw_tnd(ncol)                     ! tendency from last process
    real(kind_phys) :: tw_rer(ncol)                     ! relative error in water column

    real(kind_phys) :: te(ncol)                         ! vertical integral of total energy
    real(kind_phys) :: tw(ncol)                         ! vertical integral of total water
    real(kind_phys) :: temp(ncol,pver)                  ! temperature

    real(kind_phys) :: se(ncol)                         ! enthalpy or internal energy (J/m2)
    real(kind_phys) :: po(ncol)                         ! surface potential or potential energy (J/m2)
    real(kind_phys) :: ke(ncol)                         ! kinetic energy    (J/m2)
    real(kind_phys) :: wv(ncol)                         ! column integrated vapor       (kg/m2)
    real(kind_phys) :: liq(ncol)                        ! column integrated liquid      (kg/m2)
    real(kind_phys) :: ice(ncol)                        ! column integrated ice         (kg/m2)

    integer :: i

    !------------------------------------------------
    ! Physics total energy.
    !------------------------------------------------
    call get_hydrostatic_energy(                       &
        tracer             = q(1:ncol,1:pver,1:pcnst), & ! moist mixing ratios
        moist_mixing_ratio = .true.,                   &
        pdel_in            = pdel   (1:ncol,1:pver),   &
        cp_or_cv           = cp_phys(1:ncol,1:pver),   &
        U                  = u      (1:ncol,1:pver),   &
        V                  = v      (1:ncol,1:pver),   &
        T                  = T      (1:ncol,1:pver),   &
        vcoord             = energy_formula_physics,   & ! energy formula for physics
        ptop               = pintdry(1:ncol,1),        &
        phis               = phis   (1:ncol),          &
        te                 = te     (1:ncol),          & ! vertically integrated total energy
        H2O                = tw     (1:ncol),          & ! v.i. total water
        se                 = se     (1:ncol),          & ! v.i. enthalpy
        po                 = po     (1:ncol),          & ! v.i. PHIS term
        ke                 = ke     (1:ncol),          & ! v.i. kinetic energy
        wv                 = wv     (1:ncol),          & ! v.i. water vapor
        liq                = liq    (1:ncol),          & ! v.i. liquid
        ice                = ice    (1:ncol)           & ! v.i. ice
    )

    ! compute expected values and tendencies
    do i = 1, ncol
       ! change in static energy and total water
       te_dif(i) = te(i) - te_cur_phys(i)
       tw_dif(i) = tw(i) - tw_cur     (i)

       ! expected tendencies from boundary fluxes for last process
       te_tnd(i) = flx_vap(i)*(latvap+latice) - (flx_cnd(i) - flx_ice(i))*1000._kind_phys*latice + flx_sen(i)
       tw_tnd(i) = flx_vap(i) - flx_cnd(i) *1000._kind_phys

       ! cummulative tendencies from boundary fluxes
       tend_te_tnd(i) = tend_te_tnd(i) + te_tnd(i)
       tend_tw_tnd(i) = tend_tw_tnd(i) + tw_tnd(i)

       ! expected new values from previous state plus boundary fluxes
       te_xpd(i) = te_cur_phys(i) + te_tnd(i)*ztodt
       tw_xpd(i) = tw_cur     (i) + tw_tnd(i)*ztodt

       ! relative error, expected value - input state / previous state
       te_rer(i) = (te_xpd(i) - te(i)) / te_cur_phys(i)
    end do

    ! relative error for total water (allow for dry atmosphere)
    tw_rer = 0._kind_phys
    where (tw_cur(:ncol) > 0._kind_phys)
       tw_rer(:ncol) = (tw_xpd(:ncol) - tw(:ncol)) / tw_cur(:ncol)
    end where

    ! error checking
    if (print_energy_errors) then
       if (any(abs(te_rer(1:ncol)) > 1.E-14_kind_phys .or. abs(tw_rer(1:ncol)) > 1.E-10_kind_phys)) then
          do i = 1, ncol
             ! the relative error threshold for the water budget has been reduced to 1.e-10
             ! to avoid messages generated by QNEG3 calls
             ! PJR- change to identify if error in energy or water
             if (abs(te_rer(i)) > 1.E-14_kind_phys ) then
                count = count + 1
                write(iulog,*) "significant energy conservation error after ", name,        &
                      " count", count, "col", i
                write(iulog,*) te(i),te_xpd(i),te_dif(i),tend_te_tnd(i)*ztodt,  &
                      te_tnd(i)*ztodt,te_rer(i)
             endif
             if ( abs(tw_rer(i)) > 1.E-10_kind_phys) then
                count = count + 1
                write(iulog,*) "significant water conservation error after ", name,        &
                      " count", count, "col", i
                write(iulog,*) tw(i),tw_xpd(i),tw_dif(i),tend_tw_tnd(i)*ztodt,  &
                      tw_tnd(i)*ztodt,tw_rer(i)
             end if
          end do
       end if
    end if

    ! WRITE OPERATION - copy new value to state, including total water.
    ! the total water operations are consistent regardless of energy formula,
    ! so it only has to be written once.
    do i = 1, ncol
      te_cur_phys(i) = te(i)
      tw_cur(i)      = tw(i)
    end do

    !------------------------------------------------
    ! Dynamical core total energy.
    !------------------------------------------------
    if (energy_formula_dycore == ENERGY_FORMULA_DYCORE_SE) then
      ! SE dycore specific hydrostatic energy

      ! enthalpy scaling for energy consistency
      temp(1:ncol,:)   = temp_ini(1:ncol,:)+scaling_dycore(1:ncol,:)*(T(1:ncol,:)-temp_ini(1:ncol,:))

      call get_hydrostatic_energy(                               &
          tracer             = q(1:ncol,1:pver,1:pcnst),         & ! moist mixing ratios
          moist_mixing_ratio = .true.,                           &
          pdel_in            = pdel           (1:ncol,1:pver),   &
          cp_or_cv           = cp_or_cv_dycore(1:ncol,1:pver),   &
          U                  = u              (1:ncol,1:pver),   &
          V                  = v              (1:ncol,1:pver),   &
          T                  = temp           (1:ncol,1:pver),   & ! enthalpy-scaled temperature for energy consistency
          vcoord             = energy_formula_dycore,            & ! energy formula for dycore
          ptop               = pintdry        (1:ncol,1),        &
          phis               = phis           (1:ncol),          &
          te                 = te_cur_dyn     (1:ncol)           & ! WRITE OPERATION - vertically integrated total energy
      )

    else if (energy_formula_dycore == ENERGY_FORMULA_DYCORE_MPAS) then
      ! MPAS dycore: compute cv if vertical coordinate is height: cv = cp - R

      ! REMOVECAM: note this scaling is different with subcols off/on which is why it was put into separate scheme (hplin, 9/5/24)
      temp(1:ncol,:) = temp_ini(1:ncol,:)+scaling_dycore(1:ncol,:)*(T(1:ncol,:)-temp_ini(1:ncol,:))

      call get_hydrostatic_energy(                               &
          tracer             = q(1:ncol,1:pver,1:pcnst),         & ! moist mixing ratios
          moist_mixing_ratio = .true.,                           &
          pdel_in            = pdel           (1:ncol,1:pver),   &
          cp_or_cv           = cp_or_cv_dycore(1:ncol,1:pver),   &
          U                  = u              (1:ncol,1:pver),   &
          V                  = v              (1:ncol,1:pver),   &
          T                  = temp           (1:ncol,1:pver),   & ! enthalpy-scaled temperature for energy consistency
          vcoord             = energy_formula_dycore,            & ! energy formula for dycore
          ptop               = pintdry        (1:ncol,1),        &
          phis               = phis           (1:ncol),          &
          z_mid              = z_ini          (1:ncol,:),        & ! unique for MPAS
          te                 = te_cur_dyn     (1:ncol)           & ! WRITE OPERATION - vertically integrated total energy
      )

    else
      ! FV dycore
      te_cur_dyn(1:ncol) = te(1:ncol)
    end if
  end subroutine check_energy_chng_run

end module check_energy_chng
