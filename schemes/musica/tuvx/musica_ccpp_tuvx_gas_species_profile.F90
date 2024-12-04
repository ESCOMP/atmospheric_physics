module musica_ccpp_tuvx_gas_species_profile
  use ccpp_kinds, only: kind_phys

  implicit none

  private
  public :: create_gas_species, configure_gas_species, create_gas_species_profile_group, &
  deallocate_gas_species, deallocate_gas_species_profile_group, set_gas_species_values

  !> Conversion factor from km to cm
  real(kind_phys), parameter, public :: km_to_cm = 1.0e5
  !> Conversion factor from m3 to cm3
  real(kind_phys), parameter, public :: m_3_to_cm_3 = 1.0e6
  !> Default value of number of wavelength grid bins
  integer, parameter :: DEFAULT_NUM_GRID_SECTIONS = 0
  !> Number of wavelength grid bins
  integer, protected :: num_midpoint_heights_ = DEFAULT_NUM_GRID_SECTIONS
  integer, protected :: num_interfaces_ = DEFAULT_NUM_GRID_SECTIONS
  real(kind_phys), allocatable :: height_delta_(:)

  type, public :: gas_species_t
    character(len=50) :: name
    character(len=20) :: unit
    real(kind_phys)   :: molar_mass
    real(kind_phys)   :: scale_height ! km
    integer           :: index_constituent_props
  end type gas_species_t

  type, public :: profile_group_t
    type(profile_t), pointer :: ptr
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
    allocate(gas_species_group(number_of_gas_species))
  
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

    type(grid_t),             intent(in)  :: height_grid
    type(gas_species_t),      intent(in)  :: gas_species_group(:)
    type(profile_t), pointer, intent(out) :: profile_group(:)
    character(len=512),       intent(out) :: errmsg
    integer,                  intent(out) :: errcode

    ! local variables
    type(error_t) :: error
    integer       :: number_of_gas_species
    integer       :: i_species

    number_of_gas_species = size(gas_species_group)
    allocate( profile_group( number_of_gas_species ) )

    do i_species = 1, number_of_gas_species
      write(*,*) "gas_species_group(i_species)%name: ", gas_species_group(i_species)%name
      profile_group(i_species) = profile_t( gas_species_group(i_species)%name, &
                            gas_species_group(i_species)%unit, height_grid, error )
      if ( has_error_occurred( error, errmsg, errcode ) ) return
      ! profile_group(i_species) = ptr
      ! TODO(jiwon) deallocate the previously allocated one?
    end do

    num_midpoint_heights_ = height_grid%number_of_sections( error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    num_interfaces_ = num_midpoint_heights_ + 1

    ! Set height delta values
    height_delta_ = get_height_delta( height_grid, num_interfaces_, errmsg, errcode )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    ! do i_species = 1, number_of_gas_species
    !   deallocate( profile_group(i_species))
    ! end do

  end subroutine create_gas_species_profile_group

  subroutine deallocate_gas_species(gas_species_group)
    type(gas_species_t), allocatable, intent(inout) :: gas_species_group(:)

    if (allocated( gas_species_group )) then
      write(*,*) "deallocate( gas_species_group )"
      deallocate( gas_species_group )
    end if

  end subroutine

  ! subroutine deallocate_gas_species_profile_group(profile_group)
  !   use musica_tuvx_profile, only: profile_t

  !   type(profile_t), allocatable, intent(inout) :: profile_group(:)

  !   ! local variables
  !   integer :: i_species
  !   type(profile_t), pointer :: ptr

  !   if ( allocated(profile_group) ) then
  !     write(*,*) "allocated(profile_group)"
  !     ! Deallocate each element of the pointer array
  !     do i_species = 1, size(profile_group)

  !       ! ptr = profile_group(i_species)
  !       if ( allocated( profile_group(i_species)) ) then
  !         deallocate( profile_group(i_species) )
  !         write(*,*) "deallocate( profile_group(i_species) ) ", i_species
  !       ! if ( associated( ptr) ) then
  !       !   write(*,*) "deallocate( profile_group(i_species) ) ", i_species
  !       !   deallocate( ptr)
  !       end if
  !     end do

  !   ! Deallocate the entire profile_group array
  !   ! if (allocated(profile_group)) then
  !   ! nullify(profile_group)
  !   ! end if

  !   else
  !     write(*,*) "allocated(profile_group) is false"
  !     return
  !   end if

  ! end subroutine


  subroutine deallocate_gas_species_profile_group(profile_group)
    use musica_util,         only: error_t
    use musica_tuvx_profile, only: profile_t
  
    type(profile_t), pointer, intent(inout) :: profile_group(:)
    character(len=512) :: errmsg
    integer :: errcode
  
    ! type(error_t) :: error
    integer :: i_species
    logical :: deallocation_error
    type(profile_t), pointer :: ptr
  
    write(*,*) "Called deallocate_gas_species_profile_group(profile_group)"
    deallocation_error = .false.
  
    ! Check if the profile_group is allocated
    ! if (.not. allocated(profile_group)) then
    !   write(*,*) 'profile_group is not allocated.'
    !   errcode = 1
    !   return
    ! end if
  
    ! Deallocate each element of profile_group (array of pointers)
    ! do i_species = 1, size(profile_group)
    !   if (allocated( profile_group(i_species)) ) then
    !     ! Deallocate each individual pointer
    !     deallocate(profile_group(i_species), stat=errcode)
    !     if (errcode /= 0) then
    !       ! errmsg = 'Error deallocating profile_group(' // trim(adjustl(i_species)) // ')'
    !       deallocation_error = .true.
    !       exit
    !     end if
    !   end if
    ! end do
    deallocate(profile_group, stat=errcode)
    !do i_species = 1, size(profile_group)
      !ptr => profile_group(i_species)
    !  deallocate(profile_group(i_species), stat=errcode)
    !end do 
    ! write(*,*) "size:L profile_group ", size(profile_group)
    ! ! Finally, deallocate the array of pointers itself
    ! if (.not. deallocation_error) then
    !   deallocate(profile_group, stat=errcode)
    !   if (errcode /= 0) then
    !     write(*,*) 'Error deallocating the profile_group array.'
    !   else
    !     write(*,*) 'profile_group successfully deallocated.'
    !   end if
    !end if
  
  end subroutine deallocate_gas_species_profile_group

  function get_height_delta(height_grid, num_interfaces, errmsg, errcode) result( height_delta )
    use musica_ccpp_util,    only: has_error_occurred
    use musica_tuvx_grid,    only: grid_t
    use musica_util,         only: error_t

    type(grid_t),     intent(in)  :: height_grid
    integer,          intent(in)  :: num_interfaces
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errcode

    ! local variables
    type(error_t)   :: error
    real(kind_phys) :: interfaces(num_interfaces)
    real(kind_phys) :: height_delta(num_interfaces-1)
    integer         :: i_elem

    call height_grid%get_edges(interfaces, error)
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    do i_elem = 1, size(height_delta)
      height_delta(i_elem) = interfaces(i_elem+1) - interfaces(i_elem)
    end do

  end function get_height_delta

  !> Sets TUV-x surface albedo values
  subroutine set_gas_species_values(profile, gas_species, constituent, errmsg, errcode)
    use ccpp_kinds,          only: kind_phys
    use musica_ccpp_util,    only: has_error_occurred
    use musica_tuvx_profile, only: profile_t
    use musica_util,         only: error_t

    type(profile_t),  intent(inout) :: profile
    type(gas_species_t), intent(in)  :: gas_species
    real(kind_phys),  intent(in)    :: constituent(:) ! mol m-3
    ! real(kind_phys),  intent(in)    :: height_delta(:)   ! km
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode

    ! local variables
    type(error_t)              :: error
    integer                    :: num_constituents
    real(kind_phys)            :: constituent_mol_per_cm_3(size(constituent)) ! mol cm-3
    real(kind_phys)            :: interfaces(size(constituent) + 2)
    real(kind_phys)            :: densities(size(constituent) + 1)

    num_constituents = size(constituent)
    constituent_mol_per_cm_3(:) = constituent(:) / m_3_to_cm_3

    ! TODO(jiwon) if ( is_parameterized ) then

    interfaces(1) = constituent_mol_per_cm_3(num_constituents)
    interfaces(2:num_constituents+1) = constituent_mol_per_cm_3(num_constituents:1:-1)
    interfaces(num_constituents+2) = constituent_mol_per_cm_3(1)

    densities(1:num_constituents+1) = &
        height_delta_(1:num_constituents+1) * km_to_cm &
        * sqrt(interfaces(1:num_constituents+1)) &
        * sqrt(interfaces(2:num_constituents+2))

    call profile%set_edge_values( interfaces, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    call profile%set_layer_densities( densities, error)
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    call profile%calculate_exo_layer_density( gas_species%scale_height, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end subroutine set_gas_species_values

end module musica_ccpp_tuvx_gas_species_profile