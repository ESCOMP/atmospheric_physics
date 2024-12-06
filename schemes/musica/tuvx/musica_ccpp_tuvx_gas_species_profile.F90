module musica_ccpp_tuvx_gas_species_profile
  use ccpp_kinds, only: kind_phys
  use musica_tuvx_profile, only: profile_t

  implicit none

  private
  public :: create_gas_species, configure_gas_species, create_gas_species_profile_group, &
            set_gas_species_values, deallocate_gas_species, deallocate_gas_species_profile_group

  !> Conversion factor from km to cm
  real(kind_phys), parameter, public :: km_to_cm = 1.0e5
  !> Conversion factor from m3 to cm3
  real(kind_phys), parameter, public :: m_3_to_cm_3 = 1.0e6
  !> Molar mass value of dry air is obtained from 'CAM-SIMA/src/utils/std_atm_profile.F90'
  ! TODO(jiwon)
  real(kind_phys), parameter, public :: MOLAR_MASS_DRY_AIR = 0.0289644_kind_phys ! kg mol-1

  !> Definition of the gas species type
  type, public :: gas_species_t
    character(len=50) :: label
    character(len=20) :: unit
    real(kind_phys)   :: molar_mass              ! kg mol-1
    real(kind_phys)   :: scale_height            ! km
    integer           :: index_constituent_props ! index of the constituents
  end type gas_species_t

  !> Defines a type to contaion a pointer to 'profile_t'
  type, public :: profile_group_t
    type(profile_t), pointer :: profile
  end type profile_group_t

contains
  !> Creates 'gas_species_t' object
  function create_gas_species(label, unit, molar_mass, scale_height, &
      index_constituent_props) result( gas_species )

    character(len=*), intent(in) :: label
    character(len=*), intent(in) :: unit
    real(kind_phys),  intent(in) :: molar_mass ! kg mol-1
    real(kind_phys),  intent(in) :: scale_height ! km
    integer,          intent(in) :: index_constituent_props
    type(gas_species_t)          :: gas_species

    gas_species%label = label
    gas_species%unit = unit
    gas_species%molar_mass = molar_mass
    gas_species%scale_height = scale_height
    gas_species%index_constituent_props = index_constituent_props

  end function create_gas_species

  !> Configures gas species.
  ! Note: The user can customize this function to configure gas species in the system
  !       This function allocates memory for the gas species that is needed in tuvx_init and tuvx_run. 
  !       The user is responsible for deallocating the memory.
  subroutine configure_gas_species(constituent_props, gas_species_group, errmsg, errcode)
    use ccpp_const_utils,          only: ccpp_const_get_idx
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    type(ccpp_constituent_prop_ptr_t), intent(in)  :: constituent_props(:)
    type(gas_species_t), allocatable,  intent(out) :: gas_species_group(:)
    character(len=512),                intent(out) :: errmsg
    integer,                           intent(out) :: errcode

    ! local variable
    integer                    :: num_gas_species = 3
    real(kind_phys), parameter :: SCALE_HEIGHT_DRY_AIR = 8.01_kind_phys    ! km
    real(kind_phys), parameter :: SCALE_HEIGHT_O2 = 7.0_kind_phys          ! km
    real(kind_phys), parameter :: SCALE_HEIGHT_O3 = 7.0_kind_phys          ! km
    real(kind_phys)            :: molar_mass_O2                            ! kg mol-1
    real(kind_phys)            :: molar_mass_O3                            ! kg mol-1
    integer,         parameter :: INDEX_NOT_KNOWN = -9999
    integer                    :: index_air, index_O2, index_O3

    allocate(gas_species_group(num_gas_species)) ! need to deallocate 
  
    ! TODO(jiwon) - I commented out this block of code that searches for the molar mass
    ! of air and instead hard-coded it using the value imported from CAM-SIMA.
    ! Should we register air in the register phase like cloud liquid water content,
    ! or should we include it in the configuration file, or is it good as it is?
  
    ! Air
    ! Find the index of the species and create a 'gas_species_t' for it
    ! call ccpp_const_get_idx(constituent_props, "air", index_air, errmsg, errcode)
    ! write(*,*) "air index: ", index_air
    ! if (errcode /= 0) return

    ! call constituent_props(index_air)%molar_mass(molar_mass_air, errcode, errmsg)
    ! if (errcode /= 0) return

    index_air = INDEX_NOT_KNOWN
    gas_species_group(1) = create_gas_species("air", "molecule cm-3", &
                        MOLAR_MASS_DRY_AIR, SCALE_HEIGHT_DRY_AIR, index_air)
    ! O2
    call ccpp_const_get_idx(constituent_props, "O2", index_O2, errmsg, errcode)
    if (errcode /= 0) return

    call constituent_props(index_O2)%molar_mass(molar_mass_O2, errcode, errmsg)
    if (errcode /= 0) return

    gas_species_group(2) = create_gas_species("O2", "molecule cm-3", &
                            molar_mass_O2, SCALE_HEIGHT_O2, index_O2)

    ! O3
    call ccpp_const_get_idx(constituent_props, "O3", index_O3, errmsg, errcode)
    if (errcode /= 0) return

    call constituent_props(index_O3)%molar_mass(molar_mass_O3, errcode, errmsg)
    if (errcode /= 0) return

    gas_species_group(3) = create_gas_species("O3", "molecule cm-3", &
                            molar_mass_O3, SCALE_HEIGHT_O3, index_O3)

  end subroutine configure_gas_species

  !> Creates an array of profile group that contains a 'profile_t' for each gas species
  subroutine create_gas_species_profile_group(height_grid, gas_species_group, &
      profile_group, errmsg, errcode)
    use musica_ccpp_util,    only: has_error_occurred
    use musica_tuvx_grid,    only: grid_t
    use musica_tuvx_profile, only: profile_t
    use musica_util,         only: error_t

    type(grid_t),          intent(in)    :: height_grid
    type(gas_species_t),   intent(in)    :: gas_species_group(:)
    type(profile_group_t), intent(inout) :: profile_group(:)
    character(len=512),    intent(out)   :: errmsg
    integer,               intent(out)   :: errcode

    ! local variables
    type(error_t) :: error
    integer       :: i_species

    do i_species = 1, size(gas_species_group)
      profile_group(i_species)%profile => profile_t( gas_species_group(i_species)%label, &
                            gas_species_group(i_species)%unit, height_grid, error )
      if ( has_error_occurred( error, errmsg, errcode ) ) return
    end do

  end subroutine create_gas_species_profile_group

  !> Sets interface, density, and above-column density values for gas species
  subroutine set_gas_species_values(profile, gas_species, constituent, &
      height_deltas, errmsg, errcode)
    use musica_ccpp_util,    only: has_error_occurred
    use musica_tuvx_profile, only: profile_t
    use musica_util,         only: error_t

    type(profile_t),     intent(inout) :: profile
    type(gas_species_t), intent(in)    :: gas_species
    real(kind_phys),     intent(in)    :: constituent(:)  ! mol m-3
    real(kind_phys),     intent(in)    :: height_deltas(:) ! km, change in height in each vertical layer
    character(len=*),    intent(out)   :: errmsg
    integer,             intent(out)   :: errcode

    ! local variables
    type(error_t)   :: error
    integer         :: num_constituents
    real(kind_phys) :: constituent_mol_per_cm_3(size(constituent)) ! mol cm-3
    real(kind_phys) :: interfaces(size(constituent) + 2)
    real(kind_phys) :: densities(size(constituent) + 1)

    num_constituents = size(constituent)
    constituent_mol_per_cm_3(:) = constituent(:) / m_3_to_cm_3


    interfaces(1) = constituent_mol_per_cm_3(num_constituents)
    interfaces(2:num_constituents+1) = constituent_mol_per_cm_3(num_constituents:1:-1)
    interfaces(num_constituents+2) = constituent_mol_per_cm_3(1)

    densities(:) = height_deltas(:) * km_to_cm             &
                  * sqrt(interfaces(1:num_constituents+1)) &
                  * sqrt(interfaces(2:num_constituents+2))

    call profile%set_edge_values( interfaces, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    call profile%set_layer_densities( densities, error)
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    call profile%calculate_exo_layer_density( gas_species%scale_height, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end subroutine set_gas_species_values

  !> Deallocates gas species object
  subroutine deallocate_gas_species(gas_species_group)
    type(gas_species_t), allocatable, intent(inout) :: gas_species_group(:)

    if (allocated( gas_species_group )) deallocate( gas_species_group )

  end subroutine

  !> Deallocates profiles of each gas species
  subroutine deallocate_gas_species_profile_group(profile_group)
    type(profile_group_t), allocatable, intent(inout) :: profile_group(:)

    ! local variables
    integer :: i_species

    if (allocated( profile_group )) then
      do i_species = 1, size(profile_group)
        if (associated( profile_group(i_species)%profile) ) then
          deallocate( profile_group(i_species)%profile)
        end if
      end do

      deallocate( profile_group )
    end if

  end subroutine deallocate_gas_species_profile_group

end module musica_ccpp_tuvx_gas_species_profile