module state_converters

  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  ! Convert temperature to potential temperature and back
  public :: temp_to_potential_temp_run
  public :: potential_temp_to_temp_run

  ! Calculate density from equation of state/ideal gas law
  public :: calc_dry_air_ideal_gas_density_run

  ! Calculate exner
  public :: calc_exner_run

  ! Convert between wet and dry
  public :: wet_to_dry_water_vapor_run
  public :: wet_to_dry_cloud_liquid_water_run
  public :: wet_to_dry_cloud_ice_run
  public :: wet_to_dry_rain_run
  public :: dry_to_wet_water_vapor_run
  public :: dry_to_wet_cloud_liquid_water_run
  public :: dry_to_wet_cloud_ice_run
  public :: dry_to_wet_rain_run

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

!> \section arg_table_calc_dry_air_ideal_gas_density_run  Argument Table
!! \htmlinclude calc_dry_air_ideal_gas_density_run.html
  subroutine calc_dry_air_ideal_gas_density_run(ncol, nz, rair, pmiddry, temp, rho, errmsg, errflg)
    integer,          intent(in)    :: ncol         ! Number of columns
    integer,          intent(in)    :: nz           ! Number of vertical levels
    real(kind_phys),  intent(in)    :: rair(:,:)   ! gas constant for dry air (J kg-1)
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

!> \section arg_table_wet_to_dry_water_vapor_run  Argument Table
!! \htmlinclude wet_to_dry_water_vapor_run.html
  subroutine wet_to_dry_water_vapor_run(ncol, nz, pdel, pdeldry, qv, qv_dry,  &
       errmsg, errflg)

     integer,          intent(in)  :: ncol
     integer,          intent(in)  :: nz
     real(kind_phys),  intent(in)  :: pdel(:,:)    ! pressure thickness of layer (Pa)
     real(kind_phys),  intent(in)  :: pdeldry(:,:) ! dry air pressure thickness of layer (Pa)
     real(kind_phys),  intent(in)  :: qv(:,:)      ! water vapor mixing ratio wrt moist air + condensates (kg/kg)
     real(kind_phys),  intent(out) :: qv_dry(:,:)  ! water vapor mixing ratio wrt dry air (kg/kg)
     character(len=*), intent(out) :: errmsg
     integer,          intent(out) :: errflg

     integer :: k

     errflg = 0
     errmsg = ''

     do k = 1, nz
       qv_dry(:ncol,k) = qv(:ncol,k) * (pdel(:ncol,k) / pdeldry(:ncol,k))
     end do

  end subroutine wet_to_dry_water_vapor_run

!> \section arg_table_wet_to_dry_cloud_liquid_water_run  Argument Table
!! \htmlinclude wet_to_dry_cloud_liquid_water_run.html
  subroutine wet_to_dry_cloud_liquid_water_run(ncol, nz, pdel, pdeldry,       &
       qc, qc_dry, errmsg, errflg)

     integer,          intent(in)  :: ncol
     integer,          intent(in)  :: nz
     real(kind_phys),  intent(in)  :: pdel(:,:)    ! pressure thickness of layer (Pa)
     real(kind_phys),  intent(in)  :: pdeldry(:,:) ! dry air pressure thickness of layer (Pa)
     real(kind_phys),  intent(in)  :: qc(:,:)      ! cloud liquid mixing ratio wrt moist air (kg/kg)
     real(kind_phys),  intent(out) :: qc_dry(:,:)  ! cloud liquid mixing ratio wrt dry air (kg/kg)
     character(len=*), intent(out) :: errmsg
     integer,          intent(out) :: errflg

     integer :: k

     errflg = 0
     errmsg = ''

     do k = 1, nz
       qc_dry(:ncol,k) = qc(:ncol,k) * (pdel(:ncol,k) / pdeldry(:ncol,k))
     end do

  end subroutine wet_to_dry_cloud_liquid_water_run

!> \section arg_table_wet_to_dry_cloud_ice_run Argument Table
!! \htmlinclude wet_to_dry_cloud_ice_run.html
  subroutine wet_to_dry_cloud_ice_run(ncol, nz, pdel, pdeldry, &
       qi, qi_dry, errmsg, errflg)

     integer,          intent(in)  :: ncol
     integer,          intent(in)  :: nz
     real(kind_phys),  intent(in)  :: pdel(:,:)    ! pressure thickness of layer (Pa)
     real(kind_phys),  intent(in)  :: pdeldry(:,:) ! dry air pressure thickness of layer (Pa)
     real(kind_phys),  intent(in)  :: qi(:,:)      ! cloud ice mixing ratio wrt moist air (kg/kg)
     real(kind_phys),  intent(out) :: qi_dry(:,:)  ! cloud ice mixing ratio wrt dry air (kg/kg)
     character(len=*), intent(out) :: errmsg
     integer,          intent(out) :: errflg

     integer :: k

     errflg = 0
     errmsg = ''

     do k = 1, nz
       qi_dry(:ncol,k) = qi(:ncol,k) * (pdel(:ncol,k) / pdeldry(:ncol,k))
     end do

  end subroutine wet_to_dry_cloud_ice_run

!> \section arg_table_wet_to_dry_rain_run  Argument Table
!! \htmlinclude wet_to_dry_rain_run.html
  subroutine wet_to_dry_rain_run(ncol, nz, pdel, pdeldry, qr, qr_dry,         &
       errmsg, errflg)

     integer,          intent(in)  :: ncol
     integer,          intent(in)  :: nz
     real(kind_phys),  intent(in)  :: pdel(:,:)     ! pressure thickness of layer (Pa)
     real(kind_phys),  intent(in)  :: pdeldry(:,:)  ! dry air pressure thickness of layer (Pa)
     real(kind_phys),  intent(in)  :: qr(:,:)       ! rain mixing ratio wrt moist air (kg/kg)
     real(kind_phys),  intent(out) :: qr_dry(:,:)   ! rain mixing ratio wrt dry air (kg/kg)
     character(len=*), intent(out) :: errmsg
     integer,          intent(out) :: errflg

     integer :: k

     errflg = 0
     errmsg = ''

     do k = 1, nz
       qr_dry(:ncol,k) = qr(:ncol,k) * (pdel(:ncol,k) / pdeldry(:ncol,k))
     end do

  end subroutine wet_to_dry_rain_run

!> \section arg_table_dry_to_wet_water_vapor_run  Argument Table
!! \htmlinclude dry_to_wet_water_vapor_run.html
  subroutine dry_to_wet_water_vapor_run(ncol, nz, pdel, pdeldry, qv_dry, qv,  &
       errmsg, errflg)

     integer,          intent(in)  :: ncol
     integer,          intent(in)  :: nz
     real(kind_phys),  intent(in)  :: pdel(:,:)     ! pressure thickness of layer (Pa)
     real(kind_phys),  intent(in)  :: pdeldry(:,:)  ! dry air pressure thickness of layer (Pa)
     real(kind_phys),  intent(in)  :: qv_dry(:,:)   ! water vapor mixing ratio wrt dry air (kg/kg)
     real(kind_phys),  intent(out) :: qv(:,:)       ! water vapor mixing ratio wrt moist air + condensates (kg/kg)
     character(len=*), intent(out) :: errmsg
     integer,          intent(out) :: errflg

     integer :: k

     errflg = 0
     errmsg = ''

     do k = 1, nz
       qv(:ncol,k) = qv_dry(:ncol,k) * (pdeldry(:ncol,k) / pdel(:ncol,k))
     end do

  end subroutine dry_to_wet_water_vapor_run

!> \section arg_table_dry_to_wet_cloud_liquid_water_run  Argument Table
!! \htmlinclude dry_to_wet_cloud_liquid_water_run.html
  subroutine dry_to_wet_cloud_liquid_water_run(ncol, nz, pdel, pdeldry,       &
       qc_dry, qc, errmsg, errflg)

     integer,          intent(in)  :: ncol
     integer,          intent(in)  :: nz
     real(kind_phys),  intent(in)  :: pdel(:,:)    ! pressure thickness of layer (Pa)
     real(kind_phys),  intent(in)  :: pdeldry(:,:) ! dry air pressure thickness of layer (Pa)
     real(kind_phys),  intent(in)  :: qc_dry(:,:)  ! cloud liquid mixing ratio wrt dry air (kg/kg)
     real(kind_phys),  intent(out) :: qc(:,:)      ! cloud liquid mixing ratio wrt moist air (kg/kg)
     character(len=*), intent(out) :: errmsg
     integer,          intent(out) :: errflg

     integer  :: k

     errflg = 0
     errmsg = ''

     do k = 1, nz
       qc(:ncol,k) = qc_dry(:ncol,k) * (pdeldry(:ncol,k) / pdel(:ncol,k))
     end do

  end subroutine dry_to_wet_cloud_liquid_water_run

!> \section arg_table_dry_to_wet_cloud_ice_run Argument Table
!! \htmlinclude dry_to_wet_cloud_ice_run.html
  subroutine dry_to_wet_cloud_ice_run(ncol, nz, pdel, pdeldry, &
       qi_dry, qi, errmsg, errflg)

     integer,          intent(in)  :: ncol
     integer,          intent(in)  :: nz
     real(kind_phys),  intent(in)  :: pdel(:,:)    ! pressure thickness of layer (Pa)
     real(kind_phys),  intent(in)  :: pdeldry(:,:) ! dry air pressure thickness of layer (Pa)
     real(kind_phys),  intent(in)  :: qi_dry(:,:)  ! cloud ice mixing ratio wrt dry air (kg/kg)
     real(kind_phys),  intent(out) :: qi(:,:)      ! cloud ice mixing ratio wrt moist air (kg/kg)
     character(len=*), intent(out) :: errmsg
     integer,          intent(out) :: errflg

     integer :: k

     errflg = 0
     errmsg = ''

     do k = 1, nz
       qi(:ncol,k) = qi_dry(:ncol,k) * (pdeldry(:ncol,k) / pdel(:ncol,k))
     end do

  end subroutine dry_to_wet_cloud_ice_run

!> \section arg_table_dry_to_wet_rain_run  Argument Table
!! \htmlinclude dry_to_wet_rain_run.html
  subroutine dry_to_wet_rain_run(ncol, nz, pdel, pdeldry, qr_dry, qr,         &
       errmsg, errflg)

     integer,          intent(in)  :: ncol
     integer,          intent(in)  :: nz
     real(kind_phys),  intent(in)  :: pdel(:,:)    ! pressure thickness of layer (Pa)
     real(kind_phys),  intent(in)  :: pdeldry(:,:) ! dry air pressure thickness of layer (Pa)
     real(kind_phys),  intent(in)  :: qr_dry(:,:)  ! rain mixing ratio wrt dry air (kg/kg)
     real(kind_phys),  intent(out) :: qr(:,:)      ! rain mixing ratio wrt moist air (kg/kg)
     character(len=*), intent(out) :: errmsg
     integer,          intent(out) :: errflg

     integer  :: k

     errflg = 0
     errmsg = ''

     do k = 1, nz
       qr(:ncol,k) = qr_dry(:ncol,k) * (pdeldry(:ncol,k) / pdel(:ncol,k))
     end do

  end subroutine dry_to_wet_rain_run

end module state_converters
