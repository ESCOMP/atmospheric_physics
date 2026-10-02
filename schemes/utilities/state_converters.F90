module state_converters

  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

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

  ! Calculate exner
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

CONTAINS

!> \section arg_table_temp_to_potential_temp_run  Argument Table
!! \htmlinclude temp_to_potential_temp_run.html
  subroutine temp_to_potential_temp_run(ncol, nz, temp, exner, theta, errmsg, errflg)
    ! Dummy arguments
    integer,          intent(in)  :: ncol              ! Number of columns
    integer,          intent(in)  :: nz                ! Number of vertical levels
    real(kind_phys),         intent(in)  :: temp(:,:)  ! temperature (K)
    real(kind_phys),         intent(in)  :: exner(:,:) ! exner function
    real(kind_phys),         intent(out) :: theta(:,:) ! potential temperature (K)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg
    ! Local variable
    integer                       :: col

    do col = 1, nz
      theta(:ncol, col) = temp(:ncol, col) / exner(:ncol, col)
    end do
    errflg = 0
    errmsg = ''
  end subroutine temp_to_potential_temp_run

!> \section arg_table_potential_temp_to_temp_run  Argument Table
!! \htmlinclude potential_temp_to_temp_run.html
  subroutine potential_temp_to_temp_run(ncol, nz, theta, exner, temp, errmsg, errflg)
    ! Dummy arguments
    integer,          intent(in)  :: ncol               ! Number of columns
    integer,          intent(in)  :: nz                 ! Number of vertical levels
    real(kind_phys),         intent(in)  :: theta(:,:)  ! potential temperature (K)
    real(kind_phys),         intent(in)  :: exner(:,:)  ! exner function
    real(kind_phys),         intent(inout) :: temp(:,:) ! temperature (K)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg
    ! Local variable
    integer                       :: col

    do col = 1, nz
      temp(:ncol, col) = theta(:ncol, col) * exner(:ncol, col)
    end do
    errflg = 0
    errmsg = ''
  end subroutine potential_temp_to_temp_run

  !> \section arg_table_temp_to_virtual_temp_run Argument Table
  !! \htmlinclude temp_to_virtual_temp_run.html
  pure subroutine temp_to_virtual_temp_run(temp, zvirv, qv, virtual_temp, errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    real(kind_phys), intent(in) :: temp(:, :)          ! temperature (K)
    real(kind_phys), intent(in) :: zvirv(:, :)         ! ratio of water vapor gas constant to composition-dependent
                                                       ! dry air gas constant minus one (1)
    real(kind_phys), intent(in) :: qv(:, :)            ! water vapor mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys), intent(out) :: virtual_temp(:, :) ! virtual temperature (K)
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    virtual_temp(:, :) = temp(:, :) * (1.0_kind_phys + zvirv(:, :) * qv(:, :))

    errmsg = ''
    errflg = 0
  end subroutine temp_to_virtual_temp_run

  !> \section arg_table_virtual_temp_to_temp_run Argument Table
  !! \htmlinclude virtual_temp_to_temp_run.html
  pure subroutine virtual_temp_to_temp_run(virtual_temp, zvirv, qv, temp, errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    real(kind_phys), intent(in) :: virtual_temp(:, :) ! virtual temperature (K)
    real(kind_phys), intent(in) :: zvirv(:, :)        ! ratio of water vapor gas constant to composition-dependent
                                                      ! dry air gas constant minus one (1)
    real(kind_phys), intent(in) :: qv(:, :)           ! water vapor mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys), intent(out) :: temp(:, :)        ! temperature (K)
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    temp(:, :) = virtual_temp(:, :) / (1.0_kind_phys + zvirv(:, :) * qv(:, :))

    errmsg = ''
    errflg = 0
  end subroutine virtual_temp_to_temp_run

!> \section arg_table_calc_dry_air_ideal_gas_density_run  Argument Table
!! \htmlinclude calc_dry_air_ideal_gas_density_run.html
  subroutine calc_dry_air_ideal_gas_density_run(ncol, nz, rair, pmiddry, temp, rho, errmsg, errflg)
    integer,          intent(in)    :: ncol         ! Number of columns
    integer,          intent(in)    :: nz           ! Number of vertical levels
    real(kind_phys),  intent(in)    :: rair(:,:)    ! Gas constant of dry air (J kg-1 K-1)
    real(kind_phys),  intent(in)    :: pmiddry(:,:) ! Air pressure of dry air (Pa)
    real(kind_phys),  intent(in)    :: temp(:,:)    ! Air temperature (K)
    real(kind_phys),  intent(out)   :: rho(:,:)     ! Dry air density (kg m-3)
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    integer :: k

    do k = 1, nz
      rho(:ncol,k) = pmiddry(:ncol,k)/(rair(:ncol,k)*temp(:ncol,k))
    end do

    errmsg = ''
    errflg = 0

  end subroutine calc_dry_air_ideal_gas_density_run

  !> \section arg_table_calc_air_ideal_gas_density_run Argument Table
  !! \htmlinclude calc_air_ideal_gas_density_run.html
  pure subroutine calc_air_ideal_gas_density_run(pmid, rairv, virtual_temp, rho, errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    real(kind_phys), intent(in) :: pmid(:, :)         ! air pressure (Pa)
    real(kind_phys), intent(in) :: rairv(:, :)        ! composition-dependent gas constant of dry air (J kg-1 K-1)
    real(kind_phys), intent(in) :: virtual_temp(:, :) ! virtual temperature (K)
    real(kind_phys), intent(out) :: rho(:, :)         ! air density (kg m-3)
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    rho(:, :) = pmid(:, :) / (rairv(:, :) * virtual_temp(:, :))

    errmsg = ''
    errflg = 0
  end subroutine calc_air_ideal_gas_density_run

  !> \section arg_table_calc_atmosphere_layer_thickness_run Argument Table
  !! \htmlinclude calc_atmosphere_layer_thickness_run.html
  pure subroutine calc_atmosphere_layer_thickness_run( &
      ncol, &
      zisfc, &
      dz, &
      errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer, intent(in) :: ncol
    real(kind_phys), intent(in) :: zisfc(:, :)
    real(kind_phys), intent(out) :: dz(:, :)
    character(*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    integer :: i

    ! In CAM-SIMA, the first vertical index is at top of atmosphere.
    ! The last one is at bottom of atmosphere. The resulting `dz` is positive.
    do i = 1, ncol
        dz(i, :) = zisfc(i, 1:size(zisfc, 2) - 1) - zisfc(i, 2:size(zisfc, 2))
    end do

    errmsg = ''
    errflg = 0
  end subroutine calc_atmosphere_layer_thickness_run

!> \section arg_table_calc_exner_run  Argument Table
!! \htmlinclude calc_exner_run.html
  subroutine calc_exner_run(ncol, nz, cpair, rair, ref_pres, pmid, exner,     &
       errmsg, errflg)

    integer,          intent(in)  :: ncol       ! Number of columns
    integer,          intent(in)  :: nz         ! Number of vertical levels
    real(kind_phys),  intent(in)  :: rair(:,:)  ! Gas constant for dry air (J kg-1 K-1)
    real(kind_phys),  intent(in)  :: cpair(:,:) ! Heat capacity at constant pressure (J kg-1 K-1)
    real(kind_phys),  intent(in)  :: ref_pres   ! Reference pressure (Pa)
    real(kind_phys),  intent(in)  :: pmid(:,:)  ! Mid-point air pressure (Pa)
    real(kind_phys),  intent(out) :: exner(:,:) ! Exner function
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: i

    do i=1,nz
      exner(:ncol,i) = (pmid(:ncol,i)/ref_pres)**(rair(:ncol,i)/cpair(:ncol,i))
    end do

    errflg = 0
    errmsg = ''

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
end module state_converters
