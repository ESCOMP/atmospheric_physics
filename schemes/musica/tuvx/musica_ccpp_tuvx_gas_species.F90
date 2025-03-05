module musica_ccpp_tuvx_gas_species
  use ccpp_kinds, only: kind_phys

  implicit none

  private
  public :: create_dry_air_profile, create_O2_profile, create_O3_profile, &
            set_gas_species_values

  !> Conversion factor from km to cm
  real(kind_phys), parameter, public :: km_to_cm = 1.0e5_kind_phys
  !> Conversion factor from m3 to cm3
  real(kind_phys), parameter, public :: m_3_to_cm_3 = 1.0e6_kind_phys

contains

  !> Creates TUV-x dry air profile
  function create_dry_air_profile(height_grid, errmsg, errcode) &
      result(profile)
    use musica_ccpp_util,              only: has_error_occurred
    use musica_tuvx_grid,              only: grid_t
    use musica_tuvx_profile,           only: profile_t
    use musica_util,                   only: error_t
    use musica_ccpp_species,           only: tuvx_species_set
    use musica_ccpp_tuvx_load_species, only: index_dry_air

    type(grid_t),     intent(in)  :: height_grid
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errcode
    type(profile_t),  pointer     :: profile

    ! local variables
    type(error_t) :: error

    profile => profile_t( tuvx_species_set(index_dry_air)%name, &
                          tuvx_species_set(index_dry_air)%unit, &
                          height_grid, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end function create_dry_air_profile

  !> Creates TUV-x O2 profile
  function create_O2_profile(height_grid, errmsg, errcode) &
      result(profile)
    use musica_ccpp_util,              only: has_error_occurred
    use musica_tuvx_grid,              only: grid_t
    use musica_tuvx_profile,           only: profile_t
    use musica_util,                   only: error_t
    use musica_ccpp_species,           only: tuvx_species_set
    use musica_ccpp_tuvx_load_species, only: index_O2

    type(grid_t),     intent(in)  :: height_grid
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errcode
    type(profile_t),  pointer     :: profile

    ! local variables
    type(error_t) :: error

    profile => profile_t( tuvx_species_set(index_O2)%name, &
                          tuvx_species_set(index_O2)%unit, &
                          height_grid, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end function create_O2_profile

  !> Creates TUV-x O3 profile
  function create_O3_profile(height_grid, errmsg, errcode) &
      result(profile)
    use musica_ccpp_util,              only: has_error_occurred
    use musica_tuvx_grid,              only: grid_t
    use musica_tuvx_profile,           only: profile_t
    use musica_util,                   only: error_t
    use musica_ccpp_species,           only: tuvx_species_set
    use musica_ccpp_tuvx_load_species, only: index_O3

    type(grid_t),     intent(in)  :: height_grid
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errcode
    type(profile_t),  pointer     :: profile

    ! local variables
    type(error_t) :: error

    profile => profile_t( tuvx_species_set(index_O3)%name, &
                          tuvx_species_set(index_O3)%unit, &
                          height_grid, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end function create_O3_profile

  !> Sets the species constituents in the vertical layer
  subroutine set_gas_species_values(profile, dry_air_density, constituents, &
                                height_deltas, index_species, errmsg, errcode)
    use musica_ccpp_util,              only: has_error_occurred
    use musica_ccpp_species,           only: tuvx_species_set
    use musica_ccpp_tuvx_load_species, only: O3_LABEL
    use musica_tuvx_profile,           only: profile_t
    use musica_util,                   only: error_t

    type(profile_t),     intent(inout) :: profile
    real(kind_phys),     intent(in)    :: dry_air_density(:) 
    real(kind_phys),     intent(in)    :: constituents(:)  ! kg kg-1
    real(kind_phys),     intent(in)    :: height_deltas(:) ! km, change in height in each vertical layer
    integer,             intent(in)    :: index_species
    character(len=*),    intent(out)   :: errmsg
    integer,             intent(out)   :: errcode

    ! local variables
    type(error_t)   :: error
    integer         :: num_vertical_levels
    real(kind_phys) :: constituent_mol_per_cm_3(size(constituents)) ! mol cm-3
    real(kind_phys) :: interfaces(size(constituents) + 2)
    real(kind_phys) :: densities(size(constituents) + 1)
    real(kind_phys) :: molar_mass

    molar_mass = tuvx_species_set(index_species)%molar_mass
    constituent_mol_per_cm_3(:) = constituents(:) * dry_air_density(:) / molar_mass / m_3_to_cm_3

    num_vertical_levels = size(constituents)
    interfaces(1) = constituent_mol_per_cm_3(num_vertical_levels)
    interfaces(2:num_vertical_levels+1) = constituent_mol_per_cm_3(num_vertical_levels:1:-1)
    interfaces(num_vertical_levels+2) = constituent_mol_per_cm_3(1)

    if (tuvx_species_set(index_species)%name == O3_LABEL) then
      densities(:) = height_deltas(:) * km_to_cm * 0.5_kind_phys &
          * ( interfaces(1:num_vertical_levels+1) + interfaces(2:num_vertical_levels+2) )
    else
      densities(:) = height_deltas(:) * km_to_cm &
          * sqrt(interfaces(1:num_vertical_levels+1)) * sqrt(interfaces(2:num_vertical_levels+2))
    end if

    call profile%set_edge_values( interfaces, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    call profile%set_layer_densities( densities, error)
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    call profile%calculate_exo_layer_density( &
          tuvx_species_set(index_species)%scale_height, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end subroutine set_gas_species_values

end module musica_ccpp_tuvx_gas_species