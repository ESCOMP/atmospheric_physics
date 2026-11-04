module state_converters

  implicit none
  private

  ! Convert temperature to potential temperature and back
  public :: temp_to_potential_temp_run
  public :: potential_temp_to_temp_run

  ! Convert temperature to virtual temperature and back
  public :: temp_to_virtual_temp_run
  public :: virtual_temp_to_temp_run

  ! Calculate dry air density by equation of state/ideal gas law
  public :: calc_dry_air_ideal_gas_density_run

  ! Calculate air density by equation of state/ideal gas law
  public :: calc_air_ideal_gas_density_run

  ! Calculate atmosphere layer thickness
  public :: calc_atmosphere_layer_thickness_run

  ! Calculate exner function
  public :: calc_exner_run

  ! Convert between wet and dry mass mixing ratios
  public :: wet_to_dry_water_vapor_run
  public :: wet_to_dry_cloud_liquid_water_run
  public :: wet_to_dry_cloud_ice_run
  public :: wet_to_dry_rain_run
  public :: wet_to_dry_snow_run
  public :: wet_to_dry_graupel_run
  public :: dry_to_wet_water_vapor_run
  public :: dry_to_wet_cloud_liquid_water_run
  public :: dry_to_wet_cloud_ice_run
  public :: dry_to_wet_rain_run
  public :: dry_to_wet_snow_run
  public :: dry_to_wet_graupel_run

  ! Convert between wet and dry mass number concentrations
  public :: wet_to_dry_cloud_liquid_water_number_concentration_run
  public :: wet_to_dry_cloud_ice_number_concentration_run
  public :: wet_to_dry_rain_number_concentration_run
  public :: wet_to_dry_snow_number_concentration_run
  public :: wet_to_dry_graupel_number_concentration_run
  public :: dry_to_wet_cloud_liquid_water_number_concentration_run
  public :: dry_to_wet_cloud_ice_number_concentration_run
  public :: dry_to_wet_rain_number_concentration_run
  public :: dry_to_wet_snow_number_concentration_run
  public :: dry_to_wet_graupel_number_concentration_run

contains

  !> \section arg_table_temp_to_potential_temp_run Argument Table
  !! \htmlinclude temp_to_potential_temp_run.html
  pure subroutine temp_to_potential_temp_run(ncol, nz, temp, exner, theta, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol        ! Number of columns
    integer,          intent(in)  :: nz          ! Number of vertical layers
    real(kind_phys),  intent(in)  :: temp(:, :)  ! Temperature (K)
    real(kind_phys),  intent(in)  :: exner(:, :) ! Exner function (1)
    real(kind_phys),  intent(out) :: theta(:, :) ! Potential temperature (K)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: k

    do k = 1, nz
      theta(:ncol, k) = temp(:ncol, k) / exner(:ncol, k)
    end do

    errmsg = ''
    errflg = 0
  end subroutine temp_to_potential_temp_run

  !> \section arg_table_potential_temp_to_temp_run Argument Table
  !! \htmlinclude potential_temp_to_temp_run.html
  pure subroutine potential_temp_to_temp_run(ncol, nz, theta, exner, temp, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol        ! Number of columns
    integer,          intent(in)  :: nz          ! Number of vertical layers
    real(kind_phys),  intent(in)  :: theta(:, :) ! Potential temperature (K)
    real(kind_phys),  intent(in)  :: exner(:, :) ! Exner function (1)
    real(kind_phys),  intent(out) :: temp(:, :)  ! Temperature (K)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: k

    do k = 1, nz
      temp(:ncol, k) = theta(:ncol, k) * exner(:ncol, k)
    end do

    errmsg = ''
    errflg = 0
  end subroutine potential_temp_to_temp_run

  !> \section arg_table_temp_to_virtual_temp_run Argument Table
  !! \htmlinclude temp_to_virtual_temp_run.html
  pure subroutine temp_to_virtual_temp_run(ncol, nz, temp, zvirv, qv, virtual_temp, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol               ! Number of columns
    integer,          intent(in)  :: nz                 ! Number of vertical layers
    real(kind_phys),  intent(in)  :: temp(:, :)         ! Temperature (K)
    real(kind_phys),  intent(in)  :: zvirv(:, :)        ! Ratio of water vapor gas constant to composition-dependent
                                                        ! dry air gas constant minus one (1)
    real(kind_phys),  intent(in)  :: qv(:, :)           ! Water vapor mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys),  intent(out) :: virtual_temp(:, :) ! Virtual temperature (K)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: k

    do k = 1, nz
      virtual_temp(:ncol, k) = temp(:ncol, k) * (1.0_kind_phys + zvirv(:ncol, k) * qv(:ncol, k))
    end do

    errmsg = ''
    errflg = 0
  end subroutine temp_to_virtual_temp_run

  !> \section arg_table_virtual_temp_to_temp_run Argument Table
  !! \htmlinclude virtual_temp_to_temp_run.html
  pure subroutine virtual_temp_to_temp_run(ncol, nz, virtual_temp, zvirv, qv, temp, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol               ! Number of columns
    integer,          intent(in)  :: nz                 ! Number of vertical layers
    real(kind_phys),  intent(in)  :: virtual_temp(:, :) ! Virtual temperature (K)
    real(kind_phys),  intent(in)  :: zvirv(:, :)        ! Ratio of water vapor gas constant to composition-dependent
                                                        ! dry air gas constant minus one (1)
    real(kind_phys),  intent(in)  :: qv(:, :)           ! Water vapor mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys),  intent(out) :: temp(:, :)         ! Temperature (K)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: k

    do k = 1, nz
      temp(:ncol, k) = virtual_temp(:ncol, k) / (1.0_kind_phys + zvirv(:ncol, k) * qv(:ncol, k))
    end do

    errmsg = ''
    errflg = 0
  end subroutine virtual_temp_to_temp_run

  !> \section arg_table_calc_dry_air_ideal_gas_density_run Argument Table
  !! \htmlinclude calc_dry_air_ideal_gas_density_run.html
  pure subroutine calc_dry_air_ideal_gas_density_run(ncol, nz, rairv, pmiddry, temp, rho, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol          ! Number of columns
    integer,          intent(in)  :: nz            ! Number of vertical layers
    real(kind_phys),  intent(in)  :: rairv(:, :)   ! Composition-dependent gas constant of dry air (J kg-1 K-1)
    real(kind_phys),  intent(in)  :: pmiddry(:, :) ! Air pressure of dry air (Pa)
    real(kind_phys),  intent(in)  :: temp(:, :)    ! Air temperature (K)
    real(kind_phys),  intent(out) :: rho(:, :)     ! Dry air density (kg m-3)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: k

    do k = 1, nz
      rho(:ncol, k) = pmiddry(:ncol, k) / (rairv(:ncol, k) * temp(:ncol, k))
    end do

    errmsg = ''
    errflg = 0
  end subroutine calc_dry_air_ideal_gas_density_run

  !> \section arg_table_calc_air_ideal_gas_density_run Argument Table
  !! \htmlinclude calc_air_ideal_gas_density_run.html
  pure subroutine calc_air_ideal_gas_density_run(ncol, nz, rairv, pmid, virtual_temp, rho, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol               ! Number of columns
    integer,          intent(in)  :: nz                 ! Number of vertical layers
    real(kind_phys),  intent(in)  :: rairv(:, :)        ! Composition-dependent gas constant of dry air (J kg-1 K-1)
    real(kind_phys),  intent(in)  :: pmid(:, :)         ! Air pressure (Pa)
    real(kind_phys),  intent(in)  :: virtual_temp(:, :) ! Virtual temperature (K)
    real(kind_phys),  intent(out) :: rho(:, :)          ! Air density (kg m-3)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: k

    do k = 1, nz
      rho(:ncol, k) = pmid(:ncol, k) / (rairv(:ncol, k) * virtual_temp(:ncol, k))
    end do

    errmsg = ''
    errflg = 0
  end subroutine calc_air_ideal_gas_density_run

  !> \section arg_table_calc_atmosphere_layer_thickness_run Argument Table
  !! \htmlinclude calc_atmosphere_layer_thickness_run.html
  pure subroutine calc_atmosphere_layer_thickness_run(ncol, zi, dz, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    real(kind_phys),  intent(in)  :: zi(:, :)
    real(kind_phys),  intent(out) :: dz(:, :)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: i

    ! In CAM-SIMA, the first vertical index is at top of atmosphere.
    ! The last one is at bottom of atmosphere. The resulting `dz` is positive.
    do i = 1, ncol
      dz(i, :) = zi(i, 1:size(zi, 2) - 1) - zi(i, 2:size(zi, 2))
    end do

    errmsg = ''
    errflg = 0
  end subroutine calc_atmosphere_layer_thickness_run

  !> \section arg_table_calc_exner_run Argument Table
  !! \htmlinclude calc_exner_run.html
  pure subroutine calc_exner_run(ncol, nz, cpairv, rairv, pref, pmid, exner, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol         ! Number of columns
    integer,          intent(in)  :: nz           ! Number of vertical layers
    real(kind_phys),  intent(in)  :: cpairv(:, :) ! Composition-dependent specific heat of dry air at
                                                  ! constant pressure (J kg-1 K-1)
    real(kind_phys),  intent(in)  :: rairv(:, :)  ! Composition-dependent gas constant of dry air (J kg-1 K-1)
    real(kind_phys),  intent(in)  :: pref         ! Reference pressure (Pa)
    real(kind_phys),  intent(in)  :: pmid(:, :)   ! Mid-point air pressure (Pa)
    real(kind_phys),  intent(out) :: exner(:, :)  ! Exner function
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: k

    do k = 1, nz
      exner(:ncol, k) = (pmid(:ncol, k) / pref) ** (rairv(:ncol, k) / cpairv(:ncol, k))
    end do

    errmsg = ''
    errflg = 0
  end subroutine calc_exner_run

  elemental subroutine generic_wet_to_dry_mass_mixing_ratio_run(pdel, pdeldry, mmr, mmrdry)
    use ccpp_kinds, only: kind_phys

    real(kind_phys), intent(in)  :: pdel    ! Air pressure thickness (Pa)
    real(kind_phys), intent(in)  :: pdeldry ! Air pressure thickness of dry air (Pa)
    real(kind_phys), intent(in)  :: mmr     ! Constituent mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys), intent(out) :: mmrdry  ! Constituent mixing ratio wrt dry air (kg kg-1)

    mmrdry = mmr * (pdel / pdeldry)
  end subroutine generic_wet_to_dry_mass_mixing_ratio_run

  elemental subroutine generic_dry_to_wet_mass_mixing_ratio_run(pdel, pdeldry, mmrdry, mmr)
    use ccpp_kinds, only: kind_phys

    real(kind_phys), intent(in)  :: pdel    ! Air pressure thickness (Pa)
    real(kind_phys), intent(in)  :: pdeldry ! Air pressure thickness of dry air (Pa)
    real(kind_phys), intent(in)  :: mmrdry  ! Constituent mixing ratio wrt dry air (kg kg-1)
    real(kind_phys), intent(out) :: mmr     ! Constituent mixing ratio wrt moist air and condensed water (kg kg-1)

    mmr = mmrdry * (pdeldry / pdel)
  end subroutine generic_dry_to_wet_mass_mixing_ratio_run

  elemental subroutine generic_wet_to_dry_mass_number_concentration_run(pdel, pdeldry, num, numdry)
    use ccpp_kinds, only: kind_phys

    real(kind_phys), intent(in)  :: pdel    ! Air pressure thickness (Pa)
    real(kind_phys), intent(in)  :: pdeldry ! Air pressure thickness of dry air (Pa)
    real(kind_phys), intent(in)  :: num     ! Constituent mass number concentration in moist air and condensed water (kg-1)
    real(kind_phys), intent(out) :: numdry  ! Constituent mass number concentration in dry air (kg-1)

    numdry = num * (pdel / pdeldry)
  end subroutine generic_wet_to_dry_mass_number_concentration_run

  elemental subroutine generic_dry_to_wet_mass_number_concentration_run(pdel, pdeldry, numdry, num)
    use ccpp_kinds, only: kind_phys

    real(kind_phys), intent(in)  :: pdel    ! Air pressure thickness (Pa)
    real(kind_phys), intent(in)  :: pdeldry ! Air pressure thickness of dry air (Pa)
    real(kind_phys), intent(in)  :: numdry  ! Constituent mass number concentration in dry air (kg-1)
    real(kind_phys), intent(out) :: num     ! Constituent mass number concentration in moist air and condensed water (kg-1)

    num = numdry * (pdeldry / pdel)
  end subroutine generic_dry_to_wet_mass_number_concentration_run

  !> \section arg_table_wet_to_dry_water_vapor_run Argument Table
  !! \htmlinclude wet_to_dry_water_vapor_run.html
  pure subroutine wet_to_dry_water_vapor_run(ncol, nz, pdel, pdeldry, qv, qv_dry, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: nz
    real(kind_phys),  intent(in)  :: pdel(:, :)    ! Air pressure thickness (Pa)
    real(kind_phys),  intent(in)  :: pdeldry(:, :) ! Air pressure thickness of dry air (Pa)
    real(kind_phys),  intent(in)  :: qv(:, :)      ! Water vapor mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys),  intent(out) :: qv_dry(:, :)  ! Water vapor mixing ratio wrt dry air (kg kg-1)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call generic_wet_to_dry_mass_mixing_ratio_run( &
      pdel(:ncol, :nz), pdeldry(:ncol, :nz), qv(:ncol, :nz), qv_dry(:ncol, :nz))

    errmsg = ''
    errflg = 0
  end subroutine wet_to_dry_water_vapor_run

  !> \section arg_table_wet_to_dry_cloud_liquid_water_run Argument Table
  !! \htmlinclude wet_to_dry_cloud_liquid_water_run.html
  pure subroutine wet_to_dry_cloud_liquid_water_run(ncol, nz, pdel, pdeldry, qc, qc_dry, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: nz
    real(kind_phys),  intent(in)  :: pdel(:, :)    ! Air pressure thickness (Pa)
    real(kind_phys),  intent(in)  :: pdeldry(:, :) ! Air pressure thickness of dry air (Pa)
    real(kind_phys),  intent(in)  :: qc(:, :)      ! Cloud liquid water mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys),  intent(out) :: qc_dry(:, :)  ! Cloud liquid water mixing ratio wrt dry air (kg kg-1)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call generic_wet_to_dry_mass_mixing_ratio_run( &
      pdel(:ncol, :nz), pdeldry(:ncol, :nz), qc(:ncol, :nz), qc_dry(:ncol, :nz))

    errmsg = ''
    errflg = 0
  end subroutine wet_to_dry_cloud_liquid_water_run

  !> \section arg_table_wet_to_dry_cloud_ice_run Argument Table
  !! \htmlinclude wet_to_dry_cloud_ice_run.html
  pure subroutine wet_to_dry_cloud_ice_run(ncol, nz, pdel, pdeldry, qi, qi_dry, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: nz
    real(kind_phys),  intent(in)  :: pdel(:, :)    ! Air pressure thickness (Pa)
    real(kind_phys),  intent(in)  :: pdeldry(:, :) ! Air pressure thickness of dry air (Pa)
    real(kind_phys),  intent(in)  :: qi(:, :)      ! Cloud ice mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys),  intent(out) :: qi_dry(:, :)  ! Cloud ice mixing ratio wrt dry air (kg kg-1)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call generic_wet_to_dry_mass_mixing_ratio_run( &
      pdel(:ncol, :nz), pdeldry(:ncol, :nz), qi(:ncol, :nz), qi_dry(:ncol, :nz))

    errmsg = ''
    errflg = 0
  end subroutine wet_to_dry_cloud_ice_run

  !> \section arg_table_wet_to_dry_rain_run Argument Table
  !! \htmlinclude wet_to_dry_rain_run.html
  pure subroutine wet_to_dry_rain_run(ncol, nz, pdel, pdeldry, qr, qr_dry, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: nz
    real(kind_phys),  intent(in)  :: pdel(:, :)    ! Air pressure thickness (Pa)
    real(kind_phys),  intent(in)  :: pdeldry(:, :) ! Air pressure thickness of dry air (Pa)
    real(kind_phys),  intent(in)  :: qr(:, :)      ! Rain mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys),  intent(out) :: qr_dry(:, :)  ! Rain mixing ratio wrt dry air (kg kg-1)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call generic_wet_to_dry_mass_mixing_ratio_run( &
      pdel(:ncol, :nz), pdeldry(:ncol, :nz), qr(:ncol, :nz), qr_dry(:ncol, :nz))

    errmsg = ''
    errflg = 0
  end subroutine wet_to_dry_rain_run

  !> \section arg_table_wet_to_dry_snow_run Argument Table
  !! \htmlinclude wet_to_dry_snow_run.html
  pure subroutine wet_to_dry_snow_run(ncol, nz, pdel, pdeldry, qs, qs_dry, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: nz
    real(kind_phys),  intent(in)  :: pdel(:, :)    ! Air pressure thickness (Pa)
    real(kind_phys),  intent(in)  :: pdeldry(:, :) ! Air pressure thickness of dry air (Pa)
    real(kind_phys),  intent(in)  :: qs(:, :)      ! Snow mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys),  intent(out) :: qs_dry(:, :)  ! Snow mixing ratio wrt dry air (kg kg-1)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call generic_wet_to_dry_mass_mixing_ratio_run( &
      pdel(:ncol, :nz), pdeldry(:ncol, :nz), qs(:ncol, :nz), qs_dry(:ncol, :nz))

    errmsg = ''
    errflg = 0
  end subroutine wet_to_dry_snow_run

  !> \section arg_table_wet_to_dry_graupel_run Argument Table
  !! \htmlinclude wet_to_dry_graupel_run.html
  pure subroutine wet_to_dry_graupel_run(ncol, nz, pdel, pdeldry, qg, qg_dry, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: nz
    real(kind_phys),  intent(in)  :: pdel(:, :)    ! Air pressure thickness (Pa)
    real(kind_phys),  intent(in)  :: pdeldry(:, :) ! Air pressure thickness of dry air (Pa)
    real(kind_phys),  intent(in)  :: qg(:, :)      ! Graupel mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys),  intent(out) :: qg_dry(:, :)  ! Graupel mixing ratio wrt dry air (kg kg-1)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call generic_wet_to_dry_mass_mixing_ratio_run( &
      pdel(:ncol, :nz), pdeldry(:ncol, :nz), qg(:ncol, :nz), qg_dry(:ncol, :nz))

    errmsg = ''
    errflg = 0
  end subroutine wet_to_dry_graupel_run

  !> \section arg_table_dry_to_wet_water_vapor_run Argument Table
  !! \htmlinclude dry_to_wet_water_vapor_run.html
  pure subroutine dry_to_wet_water_vapor_run(ncol, nz, pdel, pdeldry, qv_dry, qv, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: nz
    real(kind_phys),  intent(in)  :: pdel(:, :)    ! Air pressure thickness (Pa)
    real(kind_phys),  intent(in)  :: pdeldry(:, :) ! Air pressure thickness of dry air (Pa)
    real(kind_phys),  intent(in)  :: qv_dry(:, :)  ! Water vapor mixing ratio wrt dry air (kg kg-1)
    real(kind_phys),  intent(out) :: qv(:, :)      ! Water vapor mixing ratio wrt moist air and condensed water (kg kg-1)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call generic_dry_to_wet_mass_mixing_ratio_run( &
      pdel(:ncol, :nz), pdeldry(:ncol, :nz), qv_dry(:ncol, :nz), qv(:ncol, :nz))

    errmsg = ''
    errflg = 0
  end subroutine dry_to_wet_water_vapor_run

  !> \section arg_table_dry_to_wet_cloud_liquid_water_run Argument Table
  !! \htmlinclude dry_to_wet_cloud_liquid_water_run.html
  pure subroutine dry_to_wet_cloud_liquid_water_run(ncol, nz, pdel, pdeldry, qc_dry, qc, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: nz
    real(kind_phys),  intent(in)  :: pdel(:, :)    ! Air pressure thickness (Pa)
    real(kind_phys),  intent(in)  :: pdeldry(:, :) ! Air pressure thickness of dry air (Pa)
    real(kind_phys),  intent(in)  :: qc_dry(:, :)  ! Cloud liquid water mixing ratio wrt dry air (kg kg-1)
    real(kind_phys),  intent(out) :: qc(:, :)      ! Cloud liquid water mixing ratio wrt moist air and condensed water (kg kg-1)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call generic_dry_to_wet_mass_mixing_ratio_run( &
      pdel(:ncol, :nz), pdeldry(:ncol, :nz), qc_dry(:ncol, :nz), qc(:ncol, :nz))

    errmsg = ''
    errflg = 0
  end subroutine dry_to_wet_cloud_liquid_water_run

  !> \section arg_table_dry_to_wet_cloud_ice_run Argument Table
  !! \htmlinclude dry_to_wet_cloud_ice_run.html
  pure subroutine dry_to_wet_cloud_ice_run(ncol, nz, pdel, pdeldry, qi_dry, qi, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: nz
    real(kind_phys),  intent(in)  :: pdel(:, :)    ! Air pressure thickness (Pa)
    real(kind_phys),  intent(in)  :: pdeldry(:, :) ! Air pressure thickness of dry air (Pa)
    real(kind_phys),  intent(in)  :: qi_dry(:, :)  ! Cloud ice mixing ratio wrt dry air (kg kg-1)
    real(kind_phys),  intent(out) :: qi(:, :)      ! Cloud ice mixing ratio wrt moist air and condensed water (kg kg-1)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call generic_dry_to_wet_mass_mixing_ratio_run( &
      pdel(:ncol, :nz), pdeldry(:ncol, :nz), qi_dry(:ncol, :nz), qi(:ncol, :nz))

    errmsg = ''
    errflg = 0
  end subroutine dry_to_wet_cloud_ice_run

  !> \section arg_table_dry_to_wet_rain_run Argument Table
  !! \htmlinclude dry_to_wet_rain_run.html
  pure subroutine dry_to_wet_rain_run(ncol, nz, pdel, pdeldry, qr_dry, qr, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: nz
    real(kind_phys),  intent(in)  :: pdel(:, :)    ! Air pressure thickness (Pa)
    real(kind_phys),  intent(in)  :: pdeldry(:, :) ! Air pressure thickness of dry air (Pa)
    real(kind_phys),  intent(in)  :: qr_dry(:, :)  ! Rain mixing ratio wrt dry air (kg kg-1)
    real(kind_phys),  intent(out) :: qr(:, :)      ! Rain mixing ratio wrt moist air and condensed water (kg kg-1)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call generic_dry_to_wet_mass_mixing_ratio_run( &
      pdel(:ncol, :nz), pdeldry(:ncol, :nz), qr_dry(:ncol, :nz), qr(:ncol, :nz))

    errmsg = ''
    errflg = 0
  end subroutine dry_to_wet_rain_run

  !> \section arg_table_dry_to_wet_snow_run Argument Table
  !! \htmlinclude dry_to_wet_snow_run.html
  pure subroutine dry_to_wet_snow_run(ncol, nz, pdel, pdeldry, qs_dry, qs, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: nz
    real(kind_phys),  intent(in)  :: pdel(:, :)    ! Air pressure thickness (Pa)
    real(kind_phys),  intent(in)  :: pdeldry(:, :) ! Air pressure thickness of dry air (Pa)
    real(kind_phys),  intent(in)  :: qs_dry(:, :)  ! Snow mixing ratio wrt dry air (kg kg-1)
    real(kind_phys),  intent(out) :: qs(:, :)      ! Snow mixing ratio wrt moist air and condensed water (kg kg-1)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call generic_dry_to_wet_mass_mixing_ratio_run( &
      pdel(:ncol, :nz), pdeldry(:ncol, :nz), qs_dry(:ncol, :nz), qs(:ncol, :nz))

    errmsg = ''
    errflg = 0
  end subroutine dry_to_wet_snow_run

  !> \section arg_table_dry_to_wet_graupel_run Argument Table
  !! \htmlinclude dry_to_wet_graupel_run.html
  pure subroutine dry_to_wet_graupel_run(ncol, nz, pdel, pdeldry, qg_dry, qg, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: nz
    real(kind_phys),  intent(in)  :: pdel(:, :)    ! Air pressure thickness (Pa)
    real(kind_phys),  intent(in)  :: pdeldry(:, :) ! Air pressure thickness of dry air (Pa)
    real(kind_phys),  intent(in)  :: qg_dry(:, :)  ! Graupel mixing ratio wrt dry air (kg kg-1)
    real(kind_phys),  intent(out) :: qg(:, :)      ! Graupel mixing ratio wrt moist air and condensed water (kg kg-1)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call generic_dry_to_wet_mass_mixing_ratio_run( &
      pdel(:ncol, :nz), pdeldry(:ncol, :nz), qg_dry(:ncol, :nz), qg(:ncol, :nz))

    errmsg = ''
    errflg = 0
  end subroutine dry_to_wet_graupel_run

  !> \section arg_table_wet_to_dry_cloud_liquid_water_number_concentration_run Argument Table
  !! \htmlinclude wet_to_dry_cloud_liquid_water_number_concentration_run.html
  pure subroutine wet_to_dry_cloud_liquid_water_number_concentration_run(ncol, nz, pdel, pdeldry, nc, nc_dry, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: nz
    real(kind_phys),  intent(in)  :: pdel(:, :)    ! Air pressure thickness (Pa)
    real(kind_phys),  intent(in)  :: pdeldry(:, :) ! Air pressure thickness of dry air (Pa)
    real(kind_phys),  intent(in)  :: nc(:, :)      ! Cloud liquid water mass number concentration
                                                   ! in moist air and condensed water (kg-1)
    real(kind_phys),  intent(out) :: nc_dry(:, :)  ! Cloud liquid water mass number concentration
                                                   ! in dry air (kg-1)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call generic_wet_to_dry_mass_number_concentration_run( &
      pdel(:ncol, :nz), pdeldry(:ncol, :nz), nc(:ncol, :nz), nc_dry(:ncol, :nz))

    errmsg = ''
    errflg = 0
  end subroutine wet_to_dry_cloud_liquid_water_number_concentration_run

  !> \section arg_table_wet_to_dry_cloud_ice_number_concentration_run Argument Table
  !! \htmlinclude wet_to_dry_cloud_ice_number_concentration_run.html
  pure subroutine wet_to_dry_cloud_ice_number_concentration_run(ncol, nz, pdel, pdeldry, ni, ni_dry, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: nz
    real(kind_phys),  intent(in)  :: pdel(:, :)    ! Air pressure thickness (Pa)
    real(kind_phys),  intent(in)  :: pdeldry(:, :) ! Air pressure thickness of dry air (Pa)
    real(kind_phys),  intent(in)  :: ni(:, :)      ! Cloud ice mass number concentration
                                                   ! in moist air and condensed water (kg-1)
    real(kind_phys),  intent(out) :: ni_dry(:, :)  ! Cloud ice mass number concentration
                                                   ! in dry air (kg-1)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call generic_wet_to_dry_mass_number_concentration_run( &
      pdel(:ncol, :nz), pdeldry(:ncol, :nz), ni(:ncol, :nz), ni_dry(:ncol, :nz))

    errmsg = ''
    errflg = 0
  end subroutine wet_to_dry_cloud_ice_number_concentration_run

  !> \section arg_table_wet_to_dry_rain_number_concentration_run Argument Table
  !! \htmlinclude wet_to_dry_rain_number_concentration_run.html
  pure subroutine wet_to_dry_rain_number_concentration_run(ncol, nz, pdel, pdeldry, nr, nr_dry, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: nz
    real(kind_phys),  intent(in)  :: pdel(:, :)    ! Air pressure thickness (Pa)
    real(kind_phys),  intent(in)  :: pdeldry(:, :) ! Air pressure thickness of dry air (Pa)
    real(kind_phys),  intent(in)  :: nr(:, :)      ! Rain mass number concentration
                                                   ! in moist air and condensed water (kg-1)
    real(kind_phys),  intent(out) :: nr_dry(:, :)  ! Rain mass number concentration
                                                   ! in dry air (kg-1)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call generic_wet_to_dry_mass_number_concentration_run( &
      pdel(:ncol, :nz), pdeldry(:ncol, :nz), nr(:ncol, :nz), nr_dry(:ncol, :nz))

    errmsg = ''
    errflg = 0
  end subroutine wet_to_dry_rain_number_concentration_run

  !> \section arg_table_wet_to_dry_snow_number_concentration_run Argument Table
  !! \htmlinclude wet_to_dry_snow_number_concentration_run.html
  pure subroutine wet_to_dry_snow_number_concentration_run(ncol, nz, pdel, pdeldry, ns, ns_dry, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: nz
    real(kind_phys),  intent(in)  :: pdel(:, :)    ! Air pressure thickness (Pa)
    real(kind_phys),  intent(in)  :: pdeldry(:, :) ! Air pressure thickness of dry air (Pa)
    real(kind_phys),  intent(in)  :: ns(:, :)      ! Snow mass number concentration
                                                   ! in moist air and condensed water (kg-1)
    real(kind_phys),  intent(out) :: ns_dry(:, :)  ! Snow mass number concentration
                                                   ! in dry air (kg-1)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call generic_wet_to_dry_mass_number_concentration_run( &
      pdel(:ncol, :nz), pdeldry(:ncol, :nz), ns(:ncol, :nz), ns_dry(:ncol, :nz))

    errmsg = ''
    errflg = 0
  end subroutine wet_to_dry_snow_number_concentration_run

  !> \section arg_table_wet_to_dry_graupel_number_concentration_run Argument Table
  !! \htmlinclude wet_to_dry_graupel_number_concentration_run.html
  pure subroutine wet_to_dry_graupel_number_concentration_run(ncol, nz, pdel, pdeldry, ng, ng_dry, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: nz
    real(kind_phys),  intent(in)  :: pdel(:, :)    ! Air pressure thickness (Pa)
    real(kind_phys),  intent(in)  :: pdeldry(:, :) ! Air pressure thickness of dry air (Pa)
    real(kind_phys),  intent(in)  :: ng(:, :)      ! Graupel mass number concentration
                                                   ! in moist air and condensed water (kg-1)
    real(kind_phys),  intent(out) :: ng_dry(:, :)  ! Graupel mass number concentration
                                                   ! in dry air (kg-1)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call generic_wet_to_dry_mass_number_concentration_run( &
      pdel(:ncol, :nz), pdeldry(:ncol, :nz), ng(:ncol, :nz), ng_dry(:ncol, :nz))

    errmsg = ''
    errflg = 0
  end subroutine wet_to_dry_graupel_number_concentration_run

  !> \section arg_table_dry_to_wet_cloud_liquid_water_number_concentration_run Argument Table
  !! \htmlinclude dry_to_wet_cloud_liquid_water_number_concentration_run.html
  pure subroutine dry_to_wet_cloud_liquid_water_number_concentration_run(ncol, nz, pdel, pdeldry, nc_dry, nc, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: nz
    real(kind_phys),  intent(in)  :: pdel(:, :)    ! Air pressure thickness (Pa)
    real(kind_phys),  intent(in)  :: pdeldry(:, :) ! Air pressure thickness of dry air (Pa)
    real(kind_phys),  intent(in)  :: nc_dry(:, :)  ! Cloud liquid water mass number concentration
                                                   ! in dry air (kg-1)
    real(kind_phys),  intent(out) :: nc(:, :)      ! Cloud liquid water mass number concentration
                                                   ! in moist air and condensed water (kg-1)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call generic_dry_to_wet_mass_number_concentration_run( &
      pdel(:ncol, :nz), pdeldry(:ncol, :nz), nc_dry(:ncol, :nz), nc(:ncol, :nz))

    errmsg = ''
    errflg = 0
  end subroutine dry_to_wet_cloud_liquid_water_number_concentration_run

  !> \section arg_table_dry_to_wet_cloud_ice_number_concentration_run Argument Table
  !! \htmlinclude dry_to_wet_cloud_ice_number_concentration_run.html
  pure subroutine dry_to_wet_cloud_ice_number_concentration_run(ncol, nz, pdel, pdeldry, ni_dry, ni, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: nz
    real(kind_phys),  intent(in)  :: pdel(:, :)    ! Air pressure thickness (Pa)
    real(kind_phys),  intent(in)  :: pdeldry(:, :) ! Air pressure thickness of dry air (Pa)
    real(kind_phys),  intent(in)  :: ni_dry(:, :)  ! Cloud ice mass number concentration
                                                   ! in dry air (kg-1)
    real(kind_phys),  intent(out) :: ni(:, :)      ! Cloud ice mass number concentration
                                                   ! in moist air and condensed water (kg-1)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call generic_dry_to_wet_mass_number_concentration_run( &
      pdel(:ncol, :nz), pdeldry(:ncol, :nz), ni_dry(:ncol, :nz), ni(:ncol, :nz))

    errmsg = ''
    errflg = 0
  end subroutine dry_to_wet_cloud_ice_number_concentration_run

  !> \section arg_table_dry_to_wet_rain_number_concentration_run Argument Table
  !! \htmlinclude dry_to_wet_rain_number_concentration_run.html
  pure subroutine dry_to_wet_rain_number_concentration_run(ncol, nz, pdel, pdeldry, nr_dry, nr, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: nz
    real(kind_phys),  intent(in)  :: pdel(:, :)    ! Air pressure thickness (Pa)
    real(kind_phys),  intent(in)  :: pdeldry(:, :) ! Air pressure thickness of dry air (Pa)
    real(kind_phys),  intent(in)  :: nr_dry(:, :)  ! Rain mass number concentration
                                                   ! in dry air (kg-1)
    real(kind_phys),  intent(out) :: nr(:, :)      ! Rain mass number concentration
                                                   ! in moist air and condensed water (kg-1)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call generic_dry_to_wet_mass_number_concentration_run( &
      pdel(:ncol, :nz), pdeldry(:ncol, :nz), nr_dry(:ncol, :nz), nr(:ncol, :nz))

    errmsg = ''
    errflg = 0
  end subroutine dry_to_wet_rain_number_concentration_run

  !> \section arg_table_dry_to_wet_snow_number_concentration_run Argument Table
  !! \htmlinclude dry_to_wet_snow_number_concentration_run.html
  pure subroutine dry_to_wet_snow_number_concentration_run(ncol, nz, pdel, pdeldry, ns_dry, ns, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: nz
    real(kind_phys),  intent(in)  :: pdel(:, :)    ! Air pressure thickness (Pa)
    real(kind_phys),  intent(in)  :: pdeldry(:, :) ! Air pressure thickness of dry air (Pa)
    real(kind_phys),  intent(in)  :: ns_dry(:, :)  ! Snow mass number concentration
                                                   ! in dry air (kg-1)
    real(kind_phys),  intent(out) :: ns(:, :)      ! Snow mass number concentration
                                                   ! in moist air and condensed water (kg-1)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call generic_dry_to_wet_mass_number_concentration_run( &
      pdel(:ncol, :nz), pdeldry(:ncol, :nz), ns_dry(:ncol, :nz), ns(:ncol, :nz))

    errmsg = ''
    errflg = 0
  end subroutine dry_to_wet_snow_number_concentration_run

  !> \section arg_table_dry_to_wet_graupel_number_concentration_run Argument Table
  !! \htmlinclude dry_to_wet_graupel_number_concentration_run.html
  pure subroutine dry_to_wet_graupel_number_concentration_run(ncol, nz, pdel, pdeldry, ng_dry, ng, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: nz
    real(kind_phys),  intent(in)  :: pdel(:, :)    ! Air pressure thickness (Pa)
    real(kind_phys),  intent(in)  :: pdeldry(:, :) ! Air pressure thickness of dry air (Pa)
    real(kind_phys),  intent(in)  :: ng_dry(:, :)  ! Graupel mass number concentration
                                                   ! in dry air (kg-1)
    real(kind_phys),  intent(out) :: ng(:, :)      ! Graupel mass number concentration
                                                   ! in moist air and condensed water (kg-1)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call generic_dry_to_wet_mass_number_concentration_run( &
      pdel(:ncol, :nz), pdeldry(:ncol, :nz), ng_dry(:ncol, :nz), ng(:ncol, :nz))

    errmsg = ''
    errflg = 0
  end subroutine dry_to_wet_graupel_number_concentration_run
end module state_converters
