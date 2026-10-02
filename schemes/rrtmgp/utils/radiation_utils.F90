module radiation_utils
  use ccpp_kinds,       only: kind_phys

  implicit none
  private

  public :: radiation_utils_init
  public :: get_sw_spectral_boundaries_ccpp
  public :: get_lw_spectral_boundaries_ccpp
  public :: get_mu_lambda_weights_ccpp
  public :: get_molar_mass_ratio
  public :: parse_rad_climate_entry

  real(kind_phys), allocatable :: wavenumber_low_shortwave(:)
  real(kind_phys), allocatable :: wavenumber_high_shortwave(:)
  real(kind_phys), allocatable :: wavenumber_low_longwave(:)
  real(kind_phys), allocatable :: wavenumber_high_longwave(:)
  integer :: nswbands
  integer :: nlwbands
  logical :: wavenumber_boundaries_set = .false.

contains

  subroutine radiation_utils_init(nswbands_in, nlwbands_in, low_shortwave, high_shortwave, &
                  low_longwave, high_longwave, errmsg, errflg)
    integer,          intent(in) :: nswbands_in         ! Number of shortwave bands
    integer,          intent(in) :: nlwbands_in         ! Number of longwave bands
    real(kind_phys),  intent(in) :: low_shortwave(:)    ! Low range values for shortwave bands  (cm-1)
    real(kind_phys),  intent(in) :: high_shortwave(:)   ! High range values for shortwave bands (cm-1)
    real(kind_phys),  intent(in) :: low_longwave(:)     ! Low range values for longwave bands   (cm-1)
    real(kind_phys),  intent(in) :: high_longwave(:)    ! High range values for longwave bands  (cm-1)
    integer,         intent(out) :: errflg
    character(len=*),intent(out) :: errmsg
    ! Local variables
    character(len=256) :: alloc_errmsg

    errflg = 0
    errmsg = ''
    nswbands = nswbands_in
    nlwbands = nlwbands_in
    allocate(wavenumber_low_shortwave(nswbands), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg,'(a,a)') 'radiation_utils_init: failed to allocate wavenumber_low_shortwave, message: ', &
          alloc_errmsg
    end if
    allocate(wavenumber_high_shortwave(nswbands), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg,'(a,a)') 'radiation_utils_init: failed to allocate wavenumber_high_shortwave, message: ', &
          alloc_errmsg
    end if
    allocate(wavenumber_low_longwave(nlwbands), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg,'(a,a)') 'radiation_utils_init: failed to allocate wavenumber_low_longwave, message: ', &
          alloc_errmsg
    end if
    allocate(wavenumber_high_longwave(nlwbands), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg,'(a,a)') 'radiation_utils_init: failed to allocate wavenumber_high_longwave, message: ', &
          alloc_errmsg
    end if

    wavenumber_low_shortwave = low_shortwave
    wavenumber_high_shortwave = high_shortwave
    wavenumber_low_longwave = low_longwave
    wavenumber_high_longwave = high_longwave

    wavenumber_boundaries_set = .true.

  end subroutine radiation_utils_init

!=========================================================================================

 subroutine get_sw_spectral_boundaries_ccpp(low_boundaries, high_boundaries, units, errmsg, errflg)

   ! provide spectral boundaries of each shortwave band in the units requested

   character(len=*),                intent(in) :: units               ! requested units
   real(kind_phys),  dimension(:), intent(out) :: low_boundaries      ! low range bounds for shortwave bands in requested units
   real(kind_phys),  dimension(:), intent(out) :: high_boundaries     ! high range bounds for shortwave bands in requested units
   character(len=*),               intent(out) :: errmsg
   integer,                        intent(out) :: errflg

   character(len=*), parameter :: sub = 'get_sw_spectral_boundaries_ccpp'
   !----------------------------------------------------------------------------

   ! Set error variables
   errmsg = ''
   errflg = 0

   if (.not. wavenumber_boundaries_set) then
      write(errmsg,'(a,a)') sub, ': ERROR, wavenumber boundaries not set.'
   end if

   select case (units)
   case ('inv_cm','cm^-1','cm-1')
      low_boundaries = wavenumber_low_shortwave
      high_boundaries = wavenumber_high_shortwave
   case('m','meter','meters')
      low_boundaries = 1.e-2_kind_phys/wavenumber_high_shortwave
      high_boundaries = 1.e-2_kind_phys/wavenumber_low_shortwave
   case('nm','nanometer','nanometers')
      low_boundaries = 1.e7_kind_phys/wavenumber_high_shortwave
      high_boundaries = 1.e7_kind_phys/wavenumber_low_shortwave
   case('um','micrometer','micrometers','micron','microns')
      low_boundaries = 1.e4_kind_phys/wavenumber_high_shortwave
      high_boundaries = 1.e4_kind_phys/wavenumber_low_shortwave
   case('cm','centimeter','centimeters')
      low_boundaries  = 1._kind_phys/wavenumber_high_shortwave
      high_boundaries = 1._kind_phys/wavenumber_low_shortwave
   case default
      write(errmsg, '(a,a,a)') sub, ': ERROR, requested spectral units not recognized: ', units
      errflg = 1
   end select

 end subroutine get_sw_spectral_boundaries_ccpp

!=========================================================================================

subroutine get_lw_spectral_boundaries_ccpp(low_boundaries, high_boundaries, units, errmsg, errflg)

   ! provide spectral boundaries of each longwave band in the units requested

   character(len=*),  intent(in) :: units                       ! requested units
   real(kind_phys),  intent(out) :: low_boundaries(nlwbands)    ! low range bounds for longwave bands in requested units
   real(kind_phys),  intent(out) :: high_boundaries(nlwbands)   ! high range bounds for longwave bands in requested units
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   character(len=*), parameter :: sub = 'get_lw_spectral_boundaries_ccpp'
   !----------------------------------------------------------------------------

   ! Set error variables
   errmsg = ''
   errflg = 0

   if (.not. wavenumber_boundaries_set) then
      write(errmsg,'(a,a)') sub, ': ERROR, wavenumber boundaries not set.'
   end if

   select case (units)
   case ('inv_cm','cm^-1','cm-1')
      low_boundaries  = wavenumber_low_longwave
      high_boundaries = wavenumber_high_longwave
   case('m','meter','meters')
      low_boundaries  = 1.e-2_kind_phys/wavenumber_high_longwave
      high_boundaries = 1.e-2_kind_phys/wavenumber_low_longwave
   case('nm','nanometer','nanometers')
      low_boundaries  = 1.e7_kind_phys/wavenumber_high_longwave
      high_boundaries = 1.e7_kind_phys/wavenumber_low_longwave
   case('um','micrometer','micrometers','micron','microns')
      low_boundaries  = 1.e4_kind_phys/wavenumber_high_longwave
      high_boundaries = 1.e4_kind_phys/wavenumber_low_longwave
   case('cm','centimeter','centimeters')
      low_boundaries  = 1._kind_phys/wavenumber_high_longwave
      high_boundaries = 1._kind_phys/wavenumber_low_longwave
   case default
      write(errmsg, '(a,a,a)') sub, ': ERROR, requested spectral units not recognized: ', units
      errflg = 1
   end select

end subroutine get_lw_spectral_boundaries_ccpp

!=========================================================================================

subroutine get_mu_lambda_weights_ccpp(nmu, nlambda, g_mu, g_lambda, lamc, pgam, &
                mu_wgts, lambda_wgts, errmsg, errflg)
  use interpolate_data, only: interp_type, lininterp_init, lininterp, &
                              extrap_method_bndry
  ! Get mu and lambda interpolation weights
  integer,            intent(in) :: nmu            ! number of mu values
  integer,            intent(in) :: nlambda        ! number of lambda values
  real(kind_phys),    intent(in) :: g_mu(:)        ! mu values
  real(kind_phys),    intent(in) :: g_lambda(:,:)  ! lambda table
  real(kind_phys),    intent(in) :: lamc           ! prognosed value of lambda for cloud
  real(kind_phys),    intent(in) :: pgam           ! prognosed value of mu for cloud
  ! Output interpolation weights. Caller is responsible for freeing these.
  type(interp_type), intent(out) :: mu_wgts        ! mu interpolation weights
  type(interp_type), intent(out) :: lambda_wgts    ! lambda interpolation weights
  character(len=*),  intent(out) :: errmsg
  integer,           intent(out) :: errflg

  integer :: ilambda
  real(kind_phys) :: g_lambda_interp(nlambda)

  ! Set error variables
  errmsg = ''
  errflg = 0

  ! Make interpolation weights for mu.
  ! (Put pgam in a temporary array for this purpose.)
  call lininterp_init(g_mu, nmu, [pgam], 1, extrap_method_bndry, mu_wgts)

  ! Use mu weights to interpolate to a row in the lambda table.
  do ilambda = 1, nlambda
     call lininterp(g_lambda(:,ilambda), nmu, &
          g_lambda_interp(ilambda:ilambda), 1, mu_wgts)
  end do

  ! Make interpolation weights for lambda.
  call lininterp_init(g_lambda_interp, nlambda, [lamc], 1, &
       extrap_method_bndry, lambda_wgts)

end subroutine get_mu_lambda_weights_ccpp

!=========================================================================================

subroutine get_molar_mass_ratio(gas_name, massratio, errmsg, errflg)
  use ccpp_kinds,              only: kind_phys

  ! return the molar mass ratio of dry air to gas based on gas_name

  character(len=*), intent(in)  :: gas_name
  real(kind_phys),  intent(out) :: massratio
  character(len=*), intent(out) :: errmsg
  integer,          intent(out) :: errflg

  ! local variables
  real(kind_phys), parameter :: amdw = 1.607793_kind_phys    ! Molecular weight of dry air / water vapor
  real(kind_phys), parameter :: amdc = 0.658114_kind_phys    ! Molecular weight of dry air / carbon dioxide
  real(kind_phys), parameter :: amdo = 0.603428_kind_phys    ! Molecular weight of dry air / ozone
  real(kind_phys), parameter :: amdm = 1.805423_kind_phys    ! Molecular weight of dry air / methane
  real(kind_phys), parameter :: amdn = 0.658090_kind_phys    ! Molecular weight of dry air / nitrous oxide
  real(kind_phys), parameter :: amdo2 = 0.905140_kind_phys   ! Molecular weight of dry air / oxygen
  real(kind_phys), parameter :: amdc1 = 0.210852_kind_phys   ! Molecular weight of dry air / CFC11
  real(kind_phys), parameter :: amdc2 = 0.239546_kind_phys   ! Molecular weight of dry air / CFC12

  character(len=*), parameter :: sub='get_molar_mass_ratio'
  !----------------------------------------------------------------------------
  ! Set error variables
  errmsg = ''
  errflg = 0

  select case (trim(gas_name)) 
     case ('H2O') 
        massratio = amdw
     case ('CO2')
        massratio = amdc
     case ('O3')
        massratio = amdo
     case ('CH4')
        massratio = amdm
     case ('N2O')
        massratio = amdn
     case ('O2')
        massratio = amdo2
     case ('CFC11')
        massratio = amdc1
     case ('CFC12')
        massratio = amdc2
     case default
        write(errmsg, '(a,a,a)') sub, ': Invalid gas: ', trim(gas_name)
        errflg = 1
  end select

end subroutine get_molar_mass_ratio

subroutine parse_rad_climate_entry(rad_climate_entry, source, identifier, gas_name, errmsg, errflg)

  ! Parse one entry of the rad_climate namelist variable, of the form
  ! "flag:identifier:gas_name", where flag is the gas source ('A' advected,
  ! 'N' non-advected, 'Z' zero), identifier is the constituent identifier
  ! (matched against the constituents' diagnostic names), and gas_name is
  ! the radiation gas ("standard") name.

  character(len=*), intent(in)  :: rad_climate_entry
  character(len=*), intent(out) :: source
  character(len=*), intent(out) :: identifier
  character(len=*), intent(out) :: gas_name
  character(len=*), intent(out) :: errmsg
  integer,          intent(out) :: errflg

  ! local variables
  character(len=len(rad_climate_entry)) :: tmpstr
  integer :: strlen, ipos, idx

  errmsg = ''
  errflg = 0

  ! There are no fields in the input strings in which a blank character is allowed.
  ! To simplify the parsing go through the input strings and remove blanks.
  tmpstr = adjustl(rad_climate_entry)
  do
     strlen = len_trim(tmpstr)
     ipos = index(tmpstr, ' ')
     if (ipos == 0 .or. ipos > strlen) exit
     tmpstr = tmpstr(:ipos-1) // tmpstr(ipos+1:strlen)
  end do

  ! Locate the ':' separating source from the constituent identifier.
  idx = index(tmpstr, ':')
  if (idx == 0) then
     errmsg = 'rad_climate namelist variable error: all entries must be of the format "flag:identifier:gas_name". Failed to parse "'//trim(tmpstr)//'"'
     errflg = 1
     return
  end if
  source = tmpstr(:idx-1)
  tmpstr = tmpstr(idx+1:)

  ! Locate the ':' separating the constituent identifier from the rad gas ("standard") name.
  idx = scan(tmpstr, ':')
  if (idx == 0) then
     errmsg = 'rad_climate namelist variable error: all entries must be of the format "flag:identifier:gas_name". Failed to parse "'//trim(tmpstr)//'"'
     errflg = 1
     return
  end if

  identifier = tmpstr(:idx-1)
  gas_name = tmpstr(idx+1:)

end subroutine parse_rad_climate_entry

end module radiation_utils
