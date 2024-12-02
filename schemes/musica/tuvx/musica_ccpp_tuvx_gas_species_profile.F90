module musica_ccpp_tuvx_gas_species_profile
  use ccpp_kinds, only: kind_phys

  implicit none

  private
  public :: create_gas_species, configure_gas_species, create_gas_species_profile_group

  !> Conversion factor from km to cm
  real(kind_phys), parameter, public :: km_to_cm = 1.0e5
  !> Conversion factor from m3 to cm3
  real(kind_phys), parameter, public :: m_3_to_cm_3 = 1.0e6
  !> Default value of number of wavelength grid bins
  integer, parameter :: DEFAULT_NUM_GRID_SECTIONS = 0
  !> Number of wavelength grid bins
  integer, protected :: num_midpoint_heights_ = DEFAULT_NUM_GRID_SECTIONS
  integer, protected :: num_interfaces_ = DEFAULT_NUM_GRID_SECTIONS

  type, public :: gas_species_t
    character(len=50) :: name
    character(len=20) :: unit
    real(kind_phys)   :: molar_mass
    real(kind_phys)   :: scale_height
    integer           :: index_constituent_props
  end type gas_species_t

contains

  function create_gas_species(name, unit, molar_mass, scale_height, &
      index_constituent_props) result( gas_species )
    use ccpp_kinds, only: kind_phys

    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: unit
    real(kind_phys),  intent(in) :: molar_mass
    real(kind_phys),  intent(in) :: scale_height
    integer,          intent(in) :: index_constituent_props
    type(gas_species_t)          :: gas_species

    gas_species%name = name
    gas_species%unit = unit
    gas_species%molar_mass = molar_mass
    gas_species%scale_height = scale_height
    gas_species%index_constituent_props = index_constituent_props

  end function create_gas_species

  !> TODO(jiwon) - write function description, error message
  subroutine configure_gas_species(constituent_props, gas_species_group, errmsg, errcode)
    use ccpp_kinds,                only: kind_phys
    use ccpp_const_utils,          only: ccpp_const_get_idx
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    type(ccpp_constituent_prop_ptr_t), intent(in)  :: constituent_props(:)
    type(gas_species_t), allocatable,  intent(out) :: gas_species_group(:)
    character(len=512),                intent(out) :: errmsg
    integer,                           intent(out) :: errcode

    ! local variable
    integer                    :: number_of_gas_species = 3
    ! The value of molar mass of dry air is from 'CAM-SIMA/src/utils/std_atm_profile.F90'
    real(kind_phys), parameter :: SCALE_HEIGHT_DRY_AIR = 8.01_kind_phys    ! km
    real(kind_phys), parameter :: SCALE_HEIGHT_O2 = 7.0_kind_phys          ! km
    real(kind_phys), parameter :: SCALE_HEIGHT_O3 = 7.0_kind_phys          ! km
    real(kind_phys), parameter :: MOLAR_MASS_DRY_AIR = 0.0289644_kind_phys ! kg mol-1
    real(kind_phys)            :: molar_mass_O2                            ! kg mol-1
    real(kind_phys)            :: molar_mass_O3                            ! kg mol-1

    integer :: index_air, index_O2, index_O3
    allocate(gas_species_group(number_of_gas_species)) ! TODO(jiwon) - when deallocate this?
  
    ! Air
    ! Find the index of the species and create a 'gas_species_t' for it
    ! call ccpp_const_get_idx(constituent_props, "air", index_air, errmsg, errcode)
    ! write(*,*) "air index: ", index_air
    ! if (errcode /= 0) return

    ! TODO(jiwon) air may not be registered
    ! call constituent_props(index_air)%molar_mass(molar_mass_air, errcode, errmsg)
    ! if (errcode /= 0) return

    gas_species_group(1) = create_gas_species("air", "molecule cm-3", &
                        MOLAR_MASS_DRY_AIR, SCALE_HEIGHT_DRY_AIR, index_air)
    ! O2
    call ccpp_const_get_idx(constituent_props, "O2", index_O2, errmsg, errcode)
    write(*,*) "air index: ", index_O2
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

  !> Creates
  subroutine create_gas_species_profile_group(height_grid, gas_species_group, &
      profile_group, errmsg, errcode)
    use musica_ccpp_util,    only: has_error_occurred
    use musica_tuvx_grid,    only: grid_t
    use musica_tuvx_profile, only: profile_t
    use musica_util,         only: error_t

    type(grid_t),              intent(in)  :: height_grid
    type(gas_species_t),       intent(in)  :: gas_species_group(:)
    type(profile_t), pointer,  intent(out) :: profile_group(:)
    character(len=512),        intent(out) :: errmsg
    integer,                   intent(out) :: errcode

    ! local variables
    type(error_t) :: error
    integer       :: number_of_gas_species
    integer       :: i_species

    number_of_gas_species = size(gas_species_group)
    allocate( profile_group( number_of_gas_species ) )

    do i_species = 1, number_of_gas_species
      allocate(profile_group(i_species))
      write(*,*) "gas_species_group(i_species)%name: ", gas_species_group(i_species)%name
      profile_group(i_species) = profile_t( gas_species_group(i_species)%name, &
                            gas_species_group(i_species)%unit, height_grid, error )
      if ( has_error_occurred( error, errmsg, errcode ) ) return
      ! TODO(jiwon) deallocate the previously allocated one?
    end do

  end subroutine create_gas_species_profile_group

  ! function get_height_delta(height_grid, num_interfaces, errmsg, errcode) result( height_delta )
  !   use musica_ccpp_util,    only: has_error_occurred
  !   use musica_tuvx_grid,    only: grid_t
  !   use musica_util,         only: error_t

  !   type(grid_t),     intent(in)  :: height_grid
  !   integer,          intent(in)  :: num_interfaces
  !   character(len=*), intent(out) :: errmsg
  !   integer,          intent(out) :: errcode

  !   ! local variables
  !   type(error_t)   :: error
  !   real(kind_phys) :: interfaces(num_interfaces)
  !   real(kind_phys) :: height_delta(num_interfaces-1)
  !   integer         :: i_elem

  !   call height_grid%get_edges(interfaces, error)
  !   if ( has_error_occurred( error, errmsg, errcode ) ) return

  !   do i_elem = 1, size(height_delta)
  !     height_delta[i] = interfaces[i+1] - interfaces[i]
  !   end do

  ! end function get_height_delta

  ! !> Creates a TUV-x surface albedo profile from the host-model wavelength grid
  ! function create_air_profile( height_grid, errmsg, errcode ) & 
  !     result( profile )
  !   use musica_ccpp_util,    only: has_error_occurred
  !   use musica_tuvx_grid,    only: grid_t
  !   use musica_tuvx_profile, only: profile_t
  !   use musica_util,         only: error_t

  !   type(grid_t),     intent(in)  :: height_grid
  !   character(len=*), intent(out) :: errmsg
  !   integer,          intent(out) :: errcode
  !   type(profile_t),  pointer     :: profile

  !   ! local variables
  !   type(error_t) :: error

  !   profile => profile_t( air_label, air_unit, height_grid, error )
  !   if ( has_error_occurred( error, errmsg, errcode ) ) return

  ! end function create_air_profile

  ! !> Creates a TUV-x surface albedo profile from the host-model wavelength grid
  ! function create_O2_profile( height_grid, errmsg, errcode ) & 
  !     result( profile )
  !   use musica_ccpp_util,    only: has_error_occurred
  !   use musica_tuvx_grid,    only: grid_t
  !   use musica_tuvx_profile, only: profile_t
  !   use musica_util,         only: error_t

  !   type(grid_t),     intent(in)  :: height_grid
  !   character(len=*), intent(out) :: errmsg
  !   integer,          intent(out) :: errcode
  !   type(profile_t),  pointer     :: profile

  !   ! local variables
  !   type(error_t) :: error

  !   profile => profile_t( O2_label, O2_unit, height_grid, error )
  !   if ( has_error_occurred( error, errmsg, errcode ) ) return

  ! end function create_O2_profile

  ! !> Creates a TUV-x surface albedo profile from the host-model wavelength grid
  ! function create_O3_profile( height_grid, errmsg, errcode ) & 
  !     result( profile )
  !   use musica_ccpp_util,    only: has_error_occurred
  !   use musica_tuvx_grid,    only: grid_t
  !   use musica_tuvx_profile, only: profile_t
  !   use musica_util,         only: error_t

  !   type(grid_t),     intent(in)  :: height_grid
  !   character(len=*), intent(out) :: errmsg
  !   integer,          intent(out) :: errcode
  !   type(profile_t),  pointer     :: profile

  !   ! local variables
  !   type(error_t) :: error

  !   profile => profile_t( O3_label, O3_unit, height_grid, error )
  !   if ( has_error_occurred( error, errmsg, errcode ) ) return

  ! end function create_O3_profile

  ! subroutine create_gas_species_profiles( air_profile, O2_profile, O3_profile, &
  !                                         height_grid, errmsg, errcode )
  !   use musica_ccpp_util,    only: has_error_occurred
  !   use musica_tuvx_grid,    only: grid_t
  !   use musica_tuvx_profile, only: profile_t
  !   use musica_util,         only: error_t

  !   type(profile_t),  intent(inout) :: air_profile
  !   type(profile_t),  intent(inout) :: O2_profile
  !   type(profile_t),  intent(inout) :: O3_profile
  !   type(grid_t),     intent(in)  :: height_grid
  !   character(len=*), intent(out) :: errmsg
  !   integer,          intent(out) :: errcode

  !   ! local variables
  !   type(error_t) :: error

  !   air_profile = create_air_profile( height_grid, errmsg, errcode )
  !   if (errcode /= 0) return

  !   O2_profile = create_O2_profile( height_grid, errmsg, errcode )
  !   if (errcode /= 0) return

  !   O3_profile = create_O3_profile( height_grid, errmsg, errcode )
  !   if (errcode /= 0) return

  !   num_midpoint_heights_ = height_grid%number_of_sections( error )
  !   if ( has_error_occurred( error, errmsg, errcode ) ) return

  !   num_interfaces_ = num_midpoint_heights_ + 1

  !   ! Set height delta values
  !   height_delta_ = get_height_delta( height_grid, num_interfaces_, errmsg, errcode )
  !   if ( has_error_occurred( error, errmsg, errcode ) ) return

  ! end subroutine create_gas_species_profiles

  ! !> Sets TUV-x surface albedo values
  ! subroutine set_dry_air_density_values( profile, dry_air_density, errmsg, errcode )
  !   use ccpp_kinds,          only: kind_phys
  !   use musica_ccpp_util,    only: has_error_occurred
  !   use musica_tuvx_profile, only: profile_t
  !   use musica_util,         only: error_t

  !   type(profile_t),  intent(inout) :: profile
  !   real(kind_phys),  intent(in)    :: dry_air_density(:) ! TODO(jiwon) : molecule cm-3 for CAM, kg m-3 as input
  !   integer,          intent(in)    :: index_air
  !   character(len=*), intent(out)   :: errmsg
  !   integer,          intent(out)   :: errcode

  !   ! local variables
  !   type(error_t)   :: error
  !   real(kind_phys) :: scale_height = 8.01_kind_phys ! km
  !   integer         :: num_midpoints = size(dry_air_density)
  !   real(kind_phys) :: interface(num_midpoints + 2)
  !   real(kind_phys) :: densities(num_midpoints + 1)

  !   interfaces(1) = dry_air_density(num_midpoints)
  !   interfaces(2:num_midpoints+1) = dry_air_density(num_midpoints:1:-1)
  !   interfaces(num_midpoints+2) = dry_air_density(1) ! use upper mid-point value for top edge
    
  !   densities(1:num_midpoints+1) = &
  !       height_delta_(1:num_midpoints+1) * km_to_cm & ! TODO(jiwon) : molecule cm-3 for CAM, kg m-3 as input
  !       * sqrt(interfaces(1:num_midpoints+1)) &
  !       * sqrt(interfaces(2:num_midpoints+2))

  !   call profile%set_edge_values( interfaces, error )
  !   if ( has_error_occurred( error, errmsg, errcode ) ) return

  !   call profile%set_layer_densities( densities, error)
  !   if ( has_error_occurred( error, errmsg, errcode ) ) return

  !   call profile%calculate_exo_layer_density( scale_height, error )
  !   if ( has_error_occurred( error, errmsg, errcode ) ) return

  ! end subroutine set_air_density_values

  ! !> Sets TUV-x surface albedo values
  ! subroutine set_O2_values(profile, O2_constituent, height_delta, errmsg, errcode)
  !   use ccpp_kinds,          only: kind_phys
  !   use musica_ccpp_util,    only: has_error_occurred
  !   use musica_tuvx_profile, only: profile_t
  !   use musica_util,         only: error_t

  !   type(profile_t),  intent(inout) :: profile
  !   real(kind_phys),  intent(in)    :: O2_constituent(:) ! mol m-3
  !   real(kind_phys),  intent(in)    :: height_delta(:)   ! km
  !   character(len=*), intent(out)   :: errmsg
  !   integer,          intent(out)   :: errcode

  !   ! local variables
  !   type(error_t)              :: error
  !   real(kind_phys), parameter :: SCALE_HEIGHT = 7.0_kind_phys                  ! km
  !   real(kind_phys)            :: O2_constituent_mol_cm-3(size(O2_constituent)) ! mol cm-3
  !   integer                    :: num_midpoints = size(O2_constituent)
  !   real(kind_phys)            :: interface(size(O2_constituent) + 2)
  !   real(kind_phys)            :: densities(size(O2_constituent) + 1)

  !   O2_constituent_mol_cm-3(:) = O2_constituent(:) * m_3_to_cm_3

  !   ! TODO(jiwon) if ( is_parameterized ) then

  !   interfaces(1) = O2_constituent_mol_cm-3(num_midpoints)
  !   interfaces(2:num_midpoints+1) = O2_constituent_mol_cm-3(num_midpoints:1:-1)
  !   interfaces(num_midpoints+2) = O2_constituent_mol_cm-3(1)

  !   densities(1:num_midpoints+1) = &
  !       height_delta(1:num_midpoints+1) * km_to_cm &
  !       * sqrt(interfaces(1:num_midpoints+1)) &
  !       * sqrt(interfaces(2:num_midpoints+2))

  !   call profile%set_edge_values( interfaces, error )
  !   if ( has_error_occurred( error, errmsg, errcode ) ) return

  !   call profile%set_layer_densities( densities, error)
  !   if ( has_error_occurred( error, errmsg, errcode ) ) return

  !   call profile%calculate_exo_layer_density( SCALE_HEIGHT, error )
  !   if ( has_error_occurred( error, errmsg, errcode ) ) return

  ! end subroutine set_O2_values

  ! !> Sets TUV-x surface albedo values
  ! subroutine set_O3_values(profile, O3_constituent, height_delta, above_column_density, &
  !                          errmsg, errcode)
  !   use ccpp_kinds,          only: kind_phys
  !   use musica_ccpp_util,    only: has_error_occurred
  !   use musica_tuvx_profile, only: profile_t
  !   use musica_util,         only: error_t

  !   type(profile_t),  intent(inout) :: profile
  !   real(kind_phys),  intent(in)    :: O2_constituent(:) ! mol m-3
  !   real(kind_phys),  intent(in)    :: height_delta(:)   ! km
  !   real(kind_phys),  intent(in)    :: above_column_density(:,:)        ! molecule cm-2
  !   character(len=*), intent(out)   :: errmsg
  !   integer,          intent(out)   :: errcode

  !   ! local variables
  !   type(error_t)              :: error
  !   real(kind_phys), parameter :: SCALE_HEIGHT = 7.0_kind_phys                  ! km
  !   real(kind_phys)            :: O2_constituent_mol_cm-3(size(O3_constituent)) ! mol cm-3
  !   integer                    :: num_midpoints = size(O3_constituent)
  !   real(kind_phys)            :: interface(size(O3_constituent) + 2)
  !   real(kind_phys)            :: densities(size(O3_constituent) + 1)

  !   O3_constituent_mol_cm-3(:) = O3_constituent(:) * m_3_to_cm_3  !TODO(devided maybe)

  !   interfaces(1) = O3_constituent_mol_cm-3(num_midpoints)
  !   interfaces(2:num_midpoints+1) = O3_constituent_mol_cm-3(num_midpoints:1:-1)
  !   interfaces(num_midpoints+2) = O3_constituent_mol_cm-3(1)

  !   call profile%set_edge_values( interfaces, error )
  !   if ( has_error_occurred( error, errmsg, errcode ) ) return

  !   densities(1:num_midpoints+1) = &
  !       height_delta(1:num_midpoints+1) * km_to_cm &
  !       * sqrt(interfaces(1:num_midpoints+1)) &
  !       * sqrt(interfaces(2:num_midpoints+2))

  !   if ( num_absorbing_species >= 1 ) then ! TODO - nabscol = 2 in chem_mods
  !     densities(1) = 0.5_kind_phys * above_column_density(num_midpoints, 1)
  !     densities(2:num_midpoints) = 0.5_kind_phys * ( above_column_density(num_midpoints-1:1:-1,1) &
  !         + above_column_density(num_midpoints:2:-1,1) )
  !     densities(num_midpoints+1) = above_column_density(0,1) &
  !         + 0.5_kind_phys * above_column_density(1,1)

  !     call profile%set_layer_densities( densities, error)
  !     if ( has_error_occurred( error, errmsg, errcode ) ) return

  !     call profile%set_exo_layer_density( above_column_density(0, 1), error) ! TODO: should be incorrect? exo_density = above_column_density(0,1) )
  !     if ( has_error_occurred( error, errmsg, errcode ) ) return 
  !   else
  !     densities(1:num_midpoints+1) = height_delta(1:num_midpoints+1) * km_to_cm &
  !       * ( interfaces(1:num_midpoints+1) + interfaces(2:num_midpoints+2) ) * 0.5_kind_phys

  !     call profile%set_layer_densities( densities, error)
  !     if ( has_error_occurred( error, errmsg, errcode ) ) return

  !     call profile%calculate_exo_layer_density( SCALE_HEIGHT, error )
  !     if ( has_error_occurred( error, errmsg, errcode ) ) return
  !   end if

  ! end subroutine set_O3_values

  ! subroutine set_gas_species_values( air_profile, O2_profile, O3_profile        &
  !     fixed_species_density, species_volume_mixing_ratio, above_column_density, &
  !     index_air, index_O2, index_O3, errmsg, errcode)
  !   use ccpp_kinds,          only: kind_phys
  !   use musica_tuvx_profile, only: profile_t
  !   use musica_util,         only: error_t

  !   type(profile_t),  intent(inout) :: air_profile
  !   type(profile_t),  intent(inout) :: O2_profile
  !   type(profile_t),  intent(inout) :: O3_profile
  !   real(kind_phys),  intent(in)    :: fixed_species_density(:,:)       ! molecule cm-3
  !   real(kind_phys),  intent(in)    :: species_volume_mixing_ratio(:,:) ! mol mol-1
  !   real(kind_phys),  intent(in)    :: above_column_density(:,:)        ! molecule cm-2
  !   integer,          intent(in)    :: index_air
  !   integer,          intent(in)    :: index_O2
  !   integer,          intent(in)    :: index_O3
  !   character(len=*), intent(out)   :: errmsg
  !   integer,          intent(out)   :: errcode

  !   ! local variables
  !   real(kind_phys) :: air_density = fixed_species_density(:, index_air)


  !   ! do the conversion to mol / cm3
  !   call set_air_density_values( air_profile, &
  !     fixed_species_density(:, index_air), errmsg, errcode )
  !   if (errcode /= 0) return

  !   call set_O2_values( O2_profile, fixed_species_density(:, index_air), &
  !                       fixed_species_density(:, index_O2),              &
  !                       species_volume_mixing_ratio(:, index_O2),        &
  !                       errmsg, errcode )
  !   if (errcode /= 0) return

  !   call set_O3_values( O3_profile, fixed_species_density(:, index_air), &
  !                       fixed_species_density(:, index_O3),              &
  !                       species_volume_mixing_ratio(:, index_O3),        & 
  !                       above_column_density(:, 1), errmsg, errcode )! TODO (jiwon) what is this '1'? 
  !   if (errcode /= 0) return

  !   ! TODO (jiwon) deallocate height delta
  ! end subroutine set_gas_species_values

end module musica_ccpp_tuvx_gas_species_profile