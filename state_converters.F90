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
  public :: wet_to_dry_run
  public :: dry_to_wet_run
 
  ! Private module data (constants set at initialization)
  real(kind_phys), parameter :: unset = 98989.8e99_kind_phys
  real(kind_phys) :: rd = unset    ! gas constant for dry air, J/(kgK)
  real(kind_phys) :: cp = unset    ! heat capacity at constant pressure, J/(kgK)
  real(kind_phys) :: lv = unset    ! latent heat of vaporization, J/kg
  real(kind_phys) :: psl = unset   ! reference pressure at sea level, mb
  real(kind_phys) :: rhoqr = unset ! density of liquid water, kg/m^3

  ! Private interfaces
  private :: safe_set ! Set constants checking for consistency

CONTAINS

  subroutine safe_set(var, set_val, var_name, errmsg, errflg)
    ! Dummy arguments
    real(kind_phys),  intent(out)  :: var     ! variable to set
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

    integer :: k, rk

    do k = 1, nz
      rk = nz - k + 1
      rho(:ncol,rk) = pmiddry(:ncol,k)/(rd*temp(:ncol,k))
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

!> \section arg_table_wet_to_dry_run  Argument Table
!! \htmlinclude wet_to_dry_run.html
  subroutine wet_to_dry_run(ncol, nz, pdel, pdeldry, qv, qc, qr, errmsg, errflg)

     integer, intent(in)     :: ncol
     integer, intent(in)     :: nz
     real(kind_phys), intent(in)    :: pdel(:,:)
     real(kind_phys), intent(in)    :: pdeldry(:,:)
     real(kind_phys), intent(inout) :: qv(:,:)
     real(kind_phys), intent(inout),optional :: qc(:,:)
     real(kind_phys), intent(inout),optional :: qr(:,:)
     character(len=*), intent(out) :: errmsg
     integer,          intent(out) :: errflg

     integer  :: k
     real(kind_phys) :: w_to_d(ncol)

     errflg = 0
     errmsg = ''

     do k=1,nz
       w_to_d(:ncol) = pdel(:ncol,k)/pdeldry(:ncol,k)
       qv(:ncol,k) = qv(:ncol,k)*w_to_d(:ncol)
       if (present(qc)) qc(:ncol,k) = qc(:ncol,k)*w_to_d(:ncol)
       if (present(qr)) qr(:ncol,k) = qr(:ncol,k)*w_to_d(:ncol)
     end do


  end subroutine wet_to_dry_run

!> \section arg_table_dry_to_wet_run  Argument Table
!! \htmlinclude dry_to_wet_run.html
  subroutine dry_to_wet_run(ncol, nz, pdel, pdeldry, qv, qc, qr, errmsg, errflg)

     integer, intent(in)     :: ncol
     integer, intent(in)     :: nz
     real(kind_phys), intent(in)    :: pdel(:,:)
     real(kind_phys), intent(in)    :: pdeldry(:,:)
     real(kind_phys), intent(inout) :: qv(:,:)
     real(kind_phys), intent(inout),optional :: qc(:,:)
     real(kind_phys), intent(inout),optional :: qr(:,:)
     character(len=*), intent(out) :: errmsg
     integer,          intent(out) :: errflg

     integer  :: k
     real(kind_phys) :: d_to_w(ncol)

     errflg = 0
     errmsg = ''

     do k=1,nz
       d_to_w(:ncol) = pdeldry(:ncol,k)/pdel(:ncol,k)
       qv(:ncol,k) = qv(:ncol,k)*d_to_w(:ncol)
       if (present(qc)) qc(:ncol,k) = qc(:ncol,k)*d_to_w(:ncol)
       if (present(qr)) qr(:ncol,k) = qr(:ncol,k)*d_to_w(:ncol)
     end do


  end subroutine dry_to_wet_run


end module state_converters
