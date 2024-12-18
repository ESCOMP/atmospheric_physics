module musica_ccpp_species
  ! This module owns musica species

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: cleanup_musica_species, register_musica_species, initialize_musica_species_indices, &
            initialize_molar_mass_array, extract_subset_constituents, update_constituents, &
            check_initialization

  integer, parameter, public :: MUSICA_INT_UNASSIGNED = -99999

  !> Definition of the gas species type
  type, public :: musica_species_t
    character(len=:), allocatable :: name
    character(len=:), allocatable :: unit
    real(kind_phys)               :: molar_mass ! kg mol-1
    integer                       :: index_musica_species = MUSICA_INT_UNASSIGNED
    integer                       :: index_constituent_props = MUSICA_INT_UNASSIGNED
    logical                       :: profiled = .false. ! optional
    real(kind_phys)               :: scale_height = 0.0_kind_phys ! km, optional
  end type musica_species_t

  interface musica_species_t
    procedure species_t_constructor
  end interface musica_species_t

  type, public :: musica_species_ptr_t
    type(musica_species_t), pointer :: species
  end type musica_species_ptr_t

  ! Species are ordered to match the sequence of the MICM state array
  type(musica_species_t), allocatable, protected, public :: micm_species_set(:) ! index should match with the MICM state array
  type(musica_species_t), allocatable, protected, public :: tuvx_species_set(:)
  integer,                allocatable, protected, public :: micm_indices_constituent_props(:)
  integer,                allocatable, protected, public :: tuvx_indices_constituent_props(:)
  real(kind_phys),        allocatable, protected, public :: micm_molar_mass_array(:) ! kg mol-1
  integer,                             protected, public :: number_of_micm_species = MUSICA_INT_UNASSIGNED
  integer,                             protected, public :: number_of_tuvx_species = MUSICA_INT_UNASSIGNED

contains

  !> Constructor for musica_species_t object
  function species_t_constructor(name, unit, molar_mass, scale_height, &
      index_musica_species, index_constituent_props) result( this )
  
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: unit
    real(kind_phys),  intent(in) :: molar_mass   ! kg mol-1
    real(kind_phys),  intent(in) :: scale_height ! km
    integer,          intent(in) :: index_musica_species
    integer,          intent(in) :: index_constituent_props
    type(musica_species_t)       :: this

    this%name = name
    this%unit = unit
    this%molar_mass = molar_mass
    this%scale_height = scale_height
    this%index_musica_species = index_musica_species
    this%index_constituent_props = index_constituent_props

  end function species_t_constructor

  subroutine cleanup_musica_species()

    if (allocated( micm_species_set )) deallocate( micm_species_set )
    if (allocated( micm_indices_constituent_props )) deallocate( micm_indices_constituent_props )
    if (allocated( tuvx_species_set )) deallocate( tuvx_species_set )
    if (allocated( tuvx_indices_constituent_props )) deallocate( tuvx_indices_constituent_props )
    if (allocated( micm_molar_mass_array )) deallocate( micm_molar_mass_array )

  end subroutine  cleanup_musica_species

  subroutine register_musica_species(micm_species, tuvx_species)
    type(musica_species_t), intent(in)  :: micm_species(:)
    type(musica_species_t), intent(in)  :: tuvx_species(:)

    number_of_micm_species = size(micm_species)
    allocate( micm_species_set( number_of_micm_species ) )
    micm_species_set = micm_species

    number_of_tuvx_species = size(tuvx_species)
    allocate( tuvx_species_set( number_of_tuvx_species ) )
    tuvx_species_set = tuvx_species

  end subroutine register_musica_species

  !> Retrieve the species indices from the constituents array and store them
  subroutine find_musica_species_indices(constituent_props, musica_species_set, &
      indices_constituent_props, errmsg, errcode)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ccpp_const_utils,          only: ccpp_const_get_idx

    type(ccpp_constituent_prop_ptr_t), intent(in)    :: constituent_props(:)
    type(musica_species_t),            intent(inout) :: musica_species_set(:)
    integer,                           intent(inout) :: indices_constituent_props(:)
    character(len=512),                intent(out) :: errmsg
    integer,                           intent(out) :: errcode

    ! local variables
    integer :: i_elem, index_species

    do i_elem = 1, size(musica_species_set)
      call ccpp_const_get_idx(constituent_props, musica_species_set(i_elem)%name, &
          musica_species_set(i_elem)%index_constituent_props, errmsg, errcode)
      if (errcode /= 0) return

      index_species = musica_species_set(i_elem)%index_constituent_props
      if (index_species == MUSICA_INT_UNASSIGNED) then
        errmsg = "[MUSICA Error] Unable to find index for " // musica_species_set(i_elem)%name
        errcode = 1
        return
      end if
      indices_constituent_props(i_elem) = index_species
    end do

  end subroutine find_musica_species_indices

  !> Initialize arrays to store the species indices of the CCPP constituents
  subroutine initialize_musica_species_indices(constituent_props, errmsg, errcode)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    type(ccpp_constituent_prop_ptr_t), intent(in)  :: constituent_props(:)
    character(len=512),                intent(out) :: errmsg
    integer,                           intent(out) :: errcode

    if (.not. allocated( micm_species_set ) .or. &
        .not. allocated( tuvx_species_set )) then
      errmsg = "[MUSICA Error] The MUSICA species set(s) are not allocated."
      errcode = 1
      return
    end if

    allocate( micm_indices_constituent_props( size(micm_species_set) ) )
    call find_musica_species_indices(constituent_props, micm_species_set, &
                          micm_indices_constituent_props, errmsg, errcode)
    if (errcode /= 0) return

    allocate( tuvx_indices_constituent_props( size(tuvx_species_set) ) )
    call find_musica_species_indices(constituent_props, tuvx_species_set, &
                          tuvx_indices_constituent_props, errmsg, errcode)
    if (errcode /= 0) return

  end subroutine initialize_musica_species_indices

  subroutine initialize_molar_mass_array(constituent_props, errmsg, errcode)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    type(ccpp_constituent_prop_ptr_t), intent(in)  :: constituent_props(:)
    character(len=512),                intent(out) :: errmsg
    integer,                           intent(out) :: errcode

    ! local variables
    integer :: i_elem, index_species

    if (.not. allocated( micm_species_set )) then
      errmsg = "[MUSICA Error] The MICM species set is not allocated."
      errcode = 1
      return
    end if

    allocate( micm_molar_mass_array( size(micm_species_set) ) )
    do i_elem = 1, size(micm_species_set)
      call constituent_props( micm_species_set(i_elem)%index_constituent_props ) &
                            %molar_mass( micm_molar_mass_array(i_elem), errcode, errmsg )
      if (errcode /= 0) then
        errmsg = "[MUSICA Error] Unable to get molar mass for " // micm_species_set(i_elem)%name
        return
      end if
    end do

    ! Ask if this has been implemented
    ! TODO(jiwon) Check molar mass is non zero as it becomes a denominator for unit converison
    ! this code will be deleted when the framework does the check
    do i_elem = 1, size(micm_molar_mass_array)
      if (micm_molar_mass_array(i_elem) <= 0) then
        errcode = 1
        errmsg = "[MUSICA Error] Molar mass must be greater than zero for " &
                  // micm_species_set(i_elem)%name
        return
      end if
    end do

  end subroutine initialize_molar_mass_array

  !> Extract sub-constituents array using the indices from constituents array
  subroutine extract_subset_constituents(constituents, subset_constituents, errmsg, errcode)

    real(kind_phys),    intent(in)    :: constituents(:,:,:)        ! kg kg-1
    real(kind_phys),    intent(inout) :: subset_constituents(:,:,:) ! kg kg-1
    character(len=512), intent(out)   :: errmsg
    integer,            intent(out)   :: errcode
    
    ! local variables
    integer :: i_elem

    if ( size(subset_constituents, dim=3) == number_of_micm_species ) then
      do i_elem = 1, number_of_micm_species
        subset_constituents(:,:,i_elem) = constituents(:,:,micm_indices_constituent_props(i_elem))
      end do
    else if ( size(subset_constituents, dim=3) == number_of_tuvx_species ) then
      do i_elem = 1, number_of_tuvx_species
        subset_constituents(:,:,i_elem) = constituents(:,:,tuvx_indices_constituent_props(i_elem))
      end do
    else
      errmsg = "[MUSICA Error] The given dimension for the constituents &
                doesn't match the size of any species array."
      errcode = 1
      return
    end if

  end subroutine extract_subset_constituents

  subroutine update_constituents(subset_constituents, constituents, errmsg, errcode)

    real(kind_phys),    intent(in)    :: subset_constituents(:,:,:) ! kg kg-1
    real(kind_phys),    intent(inout) :: constituents(:,:,:)        ! kg kg-1
    character(len=512), intent(out)   :: errmsg
    integer,            intent(out)   :: errcode
  
    ! local variables
    integer :: i_elem

    if ( size(subset_constituents, dim=3) == number_of_micm_species ) then
      do i_elem = 1, number_of_micm_species
        constituents(:,:,micm_indices_constituent_props(i_elem)) = subset_constituents(:,:,i_elem)
      end do
    else if ( size(subset_constituents, dim=3) == number_of_tuvx_species ) then
      do i_elem = 1, number_of_tuvx_species
        constituents(:,:,tuvx_indices_constituent_props(i_elem)) = subset_constituents(:,:,i_elem)
      end do
    else
      errmsg = "[MUSICA Error] The given dimension for the constituents &
                doesn't match the size of any species array."
      errcode = 1
      return
    end if

  end subroutine update_constituents

  subroutine check_initialization(errmsg, errcode)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    if (.not. allocated( micm_species_set )) then 
      errmsg = "[MUSICA Error] MICM species set has not been allocated."
      errcode = 1
      return
    end if
    if (.not. allocated( micm_indices_constituent_props )) then 
      errmsg = "[MUSICA Error] MICM species indices array has not been allocated."
      errcode = 1
      return
    end if
    if (.not. allocated( tuvx_species_set )) then
      errmsg = "[MUSICA Error] TUV-X species set has not been allocated."
      errcode = 1
      return
    end if

    if (.not. allocated( tuvx_indices_constituent_props )) then
      errmsg = "[MUSICA Error] TUV-X species indices array has not been allocated."
      errcode = 1
      return
    end if
    if (.not. allocated( micm_molar_mass_array )) then 
      errmsg = "[MUSICA Error] MICM molar mass array has not been allocated."
      errcode = 1
      return
    end if
    if (number_of_micm_species == MUSICA_INT_UNASSIGNED) then
      errmsg = "[MUSICA Error] The 'number_of_micm_species' variable has not been initialized."
      errcode = 1
      return
    end if
    if (number_of_tuvx_species == MUSICA_INT_UNASSIGNED) then
      errmsg = "[MUSICA Error] The 'number_of_tuvx_species' variable has not been initialized."
      errcode = 1
      return
    end if

  end subroutine check_initialization

end module musica_ccpp_species