module state_converters

  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  ! Convert temperature to potential temperature and back
  public :: temp_to_potential_temp_run
  public :: potential_temp_to_temp_run

  ! Convert dry pressure to dry air density
  public :: pres_to_density_dry_init
  public :: pres_to_density_dry_run

  ! Calculate exner
  public :: calc_exner_init
  public :: calc_exner_run

  ! Convert between wet and dry
  public :: wet_to_dry_water_vapor_run
  public :: wet_to_dry_cloud_liquid_water_run
  public :: wet_to_dry_rain_run
  public :: dry_to_wet_water_vapor_run
  public :: dry_to_wet_cloud_liquid_water_run
  public :: dry_to_wet_rain_run

  ! Private module data (constants set at initialization)
  real(kind_phys), parameter :: unset = 98989.8e99_kind_phys
  real(kind_phys) :: rd = unset    ! gas constant for dry air, J/(kgK)
  real(kind_phys) :: cp = unset    ! heat capacity at constant pressure, J/(kgK)

  ! Private interfaces
  private :: safe_set ! Set constants checking for consistency

CONTAINS

  subroutine safe_set(var, set_val, var_name, errmsg, errflg)
    ! Dummy arguments
    real(kind_phys),  intent(inout):: var     ! variable to set
    real(kind_phys),  intent(in)  :: set_val ! value to set
    character(len=*), intent(in)  :: var_name
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    if (var == unset) then
      ! var has not been set, just set it
      var = set_val
      errflg = 0
      errmsg = ''
    else if (var /= set_val) then
      errflg = 1
      errmsg = 'attempt to set '//trim(var_name)//' to inconsistent value'
    else
      ! var is already set to correct value, no error
      errflg = 0
      errmsg = ''
    end if

  end subroutine safe_set

!> \section arg_table_temp_to_potential_temp_run  Argument Table
!! \htmlinclude temp_to_potential_temp_run.html
  subroutine temp_to_potential_temp_run(ncol, nz, temp, exner, theta, errmsg, errflg)
    ! Dummy arguments
    integer,          intent(in)  :: ncol       ! Number of columns
    integer,          intent(in)  :: nz         ! Number of vertical levels
    real(kind_phys),         intent(in)  :: temp(:,:)  ! temperature (K)
    real(kind_phys),         intent(in)  :: exner(:,:) ! inverse exner function
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
    integer,          intent(in)  :: ncol       ! Number of columns
    integer,          intent(in)  :: nz         ! Number of vertical levels
    real(kind_phys),         intent(in)  :: theta(:,:) ! potential temperature (K)
    real(kind_phys),         intent(in)  :: exner(:,:) ! inverse exner function
    real(kind_phys),         intent(inout) :: temp(:,:)  ! temperature (K)
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

!> \section arg_table_pres_to_density_dry_init  Argument Table
!! \htmlinclude res_to_density_dry_init.html
  subroutine pres_to_density_dry_init(cpair, rair, errmsg, errflg)
    real(kind_phys),  intent(in)  :: rair  ! gas constant for dry air
    real(kind_phys),  intent(in)  :: cpair ! heat capacity at constant pressure
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call safe_set(cp, cpair, 'cpair', errmsg, errflg)
    if (errflg /= 0) then
      errmsg = 'pres_to_density_dry_init: '//trim(errmsg)
    else
      call safe_set(rd, rair, 'rair', errmsg, errflg)
      if (errflg /= 0) then
        errmsg = 'pres_to_density_dry_init: '//trim(errmsg)
      end if
    end if

  end subroutine pres_to_density_dry_init

!> \section arg_table_pres_to_density_dry_run  Argument Table
!! \htmlinclude pres_to_density_dry_run.html
  subroutine pres_to_density_dry_run(ncol, nz, pmiddry, temp, rho, errmsg, errflg)
    integer,          intent(in)    :: ncol      ! Number of columns
    integer,          intent(in)    :: nz        ! Number of vertical levels
    real(kind_phys),  intent(in)    :: pmiddry(:,:)
    real(kind_phys),  intent(in)    :: temp(:,:)
    real(kind_phys),         intent(out)   :: rho(:,:)  ! Dry air density (kg/m^3)
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    integer :: k

    do k = 1, nz
      rho(:ncol,k) = pmiddry(:ncol,k)/(rd*temp(:ncol,k))
    end do

    errmsg = ''
    errflg = 0

  end subroutine pres_to_density_dry_run

!> \section arg_table_calc_exner_init  Argument Table
!! \htmlinclude calc_exner_init.html
  subroutine calc_exner_init(errmsg, errflg)

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errflg = 0
    errmsg = ''

  end subroutine calc_exner_init

!> \section arg_table_calc_exner_run  Argument Table
!! \htmlinclude calc_exner_run.html
  subroutine calc_exner_run(ncol, nz, cpair, rair, pmid, exner, errmsg, errflg)

    integer,          intent(in)  :: ncol       ! Number of columns
    integer,          intent(in)  :: nz         ! Number of vertical levels
    real(kind_phys),  intent(in)  :: rair  ! gas constant for dry air
    real(kind_phys),  intent(in)  :: cpair ! heat capacity at constant pressure
    real(kind_phys),  intent(in)  :: pmid(:,:)
    real(kind_phys),  intent(out) :: exner(:,:)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: i

    do i=1,nz
      exner(:ncol,i) = (pmid(:ncol,i)/1.e5_kind_phys)**(rair/cpair)
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
     real(kind_phys),  intent(in)  :: pdel(:,:)
     real(kind_phys),  intent(in)  :: pdeldry(:,:)
     real(kind_phys),  intent(in)  :: qv(:,:)
     real(kind_phys),  intent(out) :: qv_dry(:,:)
     character(len=*), intent(out) :: errmsg
     integer,          intent(out) :: errflg

     integer         :: k

     errflg = 0
     errmsg = ''

     do k = 1, nz
       qv_dry(:ncol,k) = qv(:ncol,k) * pdel(:ncol,k) / pdeldry(:ncol,k)
     end do

  end subroutine wet_to_dry_water_vapor_run

!> \section arg_table_wet_to_dry_cloud_liquid_water_run  Argument Table
!! \htmlinclude wet_to_dry_cloud_liquid_water_run.html
  subroutine wet_to_dry_cloud_liquid_water_run(ncol, nz, pdel, pdeldry,       &
       qc, qc_dry, errmsg, errflg)

     integer,          intent(in)  :: ncol
     integer,          intent(in)  :: nz
     real(kind_phys),  intent(in)  :: pdel(:,:)
     real(kind_phys),  intent(in)  :: pdeldry(:,:)
     real(kind_phys),  intent(in)  :: qc(:,:)
     real(kind_phys),  intent(out) :: qc_dry(:,:)
     character(len=*), intent(out) :: errmsg
     integer,          intent(out) :: errflg

     integer         :: k

     errflg = 0
     errmsg = ''

     do k = 1, nz
       qc_dry(:ncol,k) = qc(:ncol,k) * pdel(:ncol,k) / pdeldry(:ncol,k)
     end do

  end subroutine wet_to_dry_cloud_liquid_water_run

!> \section arg_table_wet_to_dry_rain_run  Argument Table
!! \htmlinclude wet_to_dry_rain_run.html
  subroutine wet_to_dry_rain_run(ncol, nz, pdel, pdeldry, qr, qr_dry,         &
       errmsg, errflg)

     integer,          intent(in)  :: ncol
     integer,          intent(in)  :: nz
     real(kind_phys),  intent(in)  :: pdel(:,:)
     real(kind_phys),  intent(in)  :: pdeldry(:,:)
     real(kind_phys),  intent(in)  :: qr(:,:)
     real(kind_phys),  intent(out) :: qr_dry(:,:)
     character(len=*), intent(out) :: errmsg
     integer,          intent(out) :: errflg

     integer         :: k

     errflg = 0
     errmsg = ''

     do k = 1, nz
       qr_dry(:ncol,k) = qr(:ncol,k) * pdel(:ncol,k) / pdeldry(:ncol,k)
     end do

  end subroutine wet_to_dry_rain_run

!> \section arg_table_dry_to_wet_water_vapor_run  Argument Table
!! \htmlinclude dry_to_wet_water_vapor_run.html
  subroutine dry_to_wet_water_vapor_run(ncol, nz, pdel, pdeldry, qv_dry, qv,  &
       errmsg, errflg)

     integer,          intent(in)  :: ncol
     integer,          intent(in)  :: nz
     real(kind_phys),  intent(in)  :: pdel(:,:)
     real(kind_phys),  intent(in)  :: pdeldry(:,:)
     real(kind_phys),  intent(in)  :: qv_dry(:,:)
     real(kind_phys),  intent(out) :: qv(:,:)
     character(len=*), intent(out) :: errmsg
     integer,          intent(out) :: errflg

     integer  :: k

     errflg = 0
     errmsg = ''

     do k = 1, nz
       qv(:ncol,k) = qv_dry(:ncol,k) * pdeldry(:ncol,k) / pdel(:ncol,k)
     end do

  end subroutine dry_to_wet_water_vapor_run

!> \section arg_table_dry_to_wet_cloud_liquid_water_run  Argument Table
!! \htmlinclude dry_to_wet_cloud_liquid_water_run.html
  subroutine dry_to_wet_cloud_liquid_water_run(ncol, nz, pdel, pdeldry,       &
       qc_dry, qc, errmsg, errflg)

     integer,          intent(in)  :: ncol
     integer,          intent(in)  :: nz
     real(kind_phys),  intent(in)  :: pdel(:,:)
     real(kind_phys),  intent(in)  :: pdeldry(:,:)
     real(kind_phys),  intent(in)  :: qc_dry(:,:)
     real(kind_phys),  intent(out) :: qc(:,:)
     character(len=*), intent(out) :: errmsg
     integer,          intent(out) :: errflg

     integer  :: k

     errflg = 0
     errmsg = ''

     do k = 1, nz
       qc(:ncol,k) = qc_dry(:ncol,k) * pdeldry(:ncol,k) / pdel(:ncol,k)
     end do

  end subroutine dry_to_wet_cloud_liquid_water_run

!> \section arg_table_dry_to_wet_rain_run  Argument Table
!! \htmlinclude dry_to_wet_rain_run.html
  subroutine dry_to_wet_rain_run(ncol, nz, pdel, pdeldry, qr_dry, qr,         &
       errmsg, errflg)

     integer,          intent(in)  :: ncol
     integer,          intent(in)  :: nz
     real(kind_phys),  intent(in)  :: pdel(:,:)
     real(kind_phys),  intent(in)  :: pdeldry(:,:)
     real(kind_phys),  intent(in)  :: qr_dry(:,:)
     real(kind_phys),  intent(out) :: qr(:,:)
     character(len=*), intent(out) :: errmsg
     integer,          intent(out) :: errflg

     integer  :: k

     errflg = 0
     errmsg = ''

     do k = 1, nz
       qr(:ncol,k) = qr_dry(:ncol,k) * pdeldry(:ncol,k) / pdel(:ncol,k)
     end do

  end subroutine dry_to_wet_rain_run

end module state_converters
