! Reshape prescribed (non-advected) radiation GHGs in the stratosphere with a vertical profile.
! The mixing ratio of GHGs is uniform below a climatological tropopause.
! Above it they decay with gas- and latitude-dependent exponents.
! CO2 and O2 are vertically uniform.
!
! This scheme should run after GHGs are prescribed (but never when all GHGs are prognostic)
!
! This is the equivalent of CAM ghg_data.
module prescribe_radiative_gas_vertical_profiles
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: prescribe_radiative_gas_vertical_profiles_init
  public :: prescribe_radiative_gas_vertical_profiles_timestep_init

  ! Gases whose stratospheric profile is prescribed:
  integer, parameter :: num_gases = 4
  character(len=5), parameter :: gas_names(num_gases) = &
       [character(len=5) :: 'CH4', 'N2O', 'CFC11', 'CFC12']

  ! Stratospheric decay exponent coefficients per gas:
  real(kind_phys), parameter :: exponent_coefs(4, num_gases) = reshape( &
       [0.2353_kind_phys, 0.0_kind_phys,     0.2353_kind_phys, 0.0225489_kind_phys, & ! CH4
        0.3478_kind_phys, 0.00116_kind_phys, 0.4000_kind_phys, 0.013333_kind_phys,  & ! N2O
        0.7273_kind_phys, 0.00606_kind_phys, 1.00_kind_phys,   0.013333_kind_phys,  & ! CFC11
        0.4000_kind_phys, 0.00222_kind_phys, 0.50_kind_phys,   0.024444_kind_phys], & ! CFC12
       [4, num_gases])

  ! Constituent indices resolved at init
  ! Unresolved or advected gases (i.e., owned by prognostic chemistry) are not reshaped.
  integer :: gas_const_idx(num_gases)
  logical :: gas_handled(num_gases)

contains

!> \section arg_table_prescribe_radiative_gas_vertical_profiles_init Argument Table
!! \htmlinclude arg_table_prescribe_radiative_gas_vertical_profiles_init.html
  subroutine prescribe_radiative_gas_vertical_profiles_init(const_props, errmsg, errflg)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t, int_unassigned
    use ccpp_scheme_utils,         only: ccpp_constituent_index

    type(ccpp_constituent_prop_ptr_t), intent(in) :: const_props(:)

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: m
    logical :: is_advected
    character(len=*), parameter :: subname = 'prescribe_radiative_gas_vertical_profiles_init'

    errmsg = ''
    errflg = 0

    gas_handled(:) = .false.
    gas_const_idx(:) = int_unassigned

    do m = 1, num_gases
      call ccpp_constituent_index(trim(gas_names(m)), gas_const_idx(m), errflg, errmsg)
      if (errflg /= 0) then
        return
      end if
      if (gas_const_idx(m) == int_unassigned) then
        cycle
      end if

      ! An advected gas belongs to prognostic chemistry; its vertical structure
      ! must not be overwritten by this prescription.
      call const_props(gas_const_idx(m))%is_advected(is_advected, errflg, errmsg)
      if (errflg /= 0) then
        return
      end if
      gas_handled(m) = .not. is_advected
    end do

    if (.not. any(gas_handled)) then
      errmsg = subname // ': none of CH4, N2O, CFC11, CFC12 is available as a' // &
               ' non-advected constituent. Provide them (e.g. through' // &
               ' prescribe_radiative_gas_concentrations) or remove this scheme from the suite'
      errflg = 1
      return
    end if

  end subroutine prescribe_radiative_gas_vertical_profiles_init

!> \section arg_table_prescribe_radiative_gas_vertical_profiles_timestep_init Argument Table
!! \htmlinclude arg_table_prescribe_radiative_gas_vertical_profiles_timestep_init.html
  subroutine prescribe_radiative_gas_vertical_profiles_timestep_init(&
             ncol, pver, lat, pi, &
             pmid, constituents, errmsg, errflg)

    use atmos_phys_rad_utils, only: stratospheric_trace_gas_profile

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    real(kind_phys),  intent(in)    :: lat(:)              ! latitude [rad]
    real(kind_phys),  intent(in)    :: pi                  ! pi constant [1]
    real(kind_phys),  intent(in)    :: pmid(:,:)           ! air pressure at layer midpoints [Pa]
    real(kind_phys),  intent(inout) :: constituents(:,:,:) ! constituent array [kg kg-1]

    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    integer         :: m, i
    real(kind_phys) :: degrees_per_radian ! 180/pi [deg rad-1]
    real(kind_phys) :: trop_mmr           ! tropospheric mixing ratio of the gas [kg kg-1]

    errmsg = ''
    errflg = 0

    degrees_per_radian = 180._kind_phys / pi

    do m = 1, num_gases
      if (.not. gas_handled(m)) then
        cycle
      end if
      do i = 1, ncol
        ! The bottom layer is below the tropopause everywhere.
        ! It holds the tropospheric value set by the uniform fills and is preserved by the profile rewrite:
        trop_mmr = constituents(i, pver, gas_const_idx(m))
        call stratospheric_trace_gas_profile(trop_mmr, lat(i), degrees_per_radian, &
             exponent_coefs(:, m), pmid(i, :), constituents(i, :, gas_const_idx(m)))
      end do
    end do

  end subroutine prescribe_radiative_gas_vertical_profiles_timestep_init

end module prescribe_radiative_gas_vertical_profiles
