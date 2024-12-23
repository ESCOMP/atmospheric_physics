program test_musica_species

  use musica_ccpp_species

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  call test_register_musica_species()
	call test_initialize_musica_species_indices_and_molar_mass()
  call test_extract_and_update_subset_constituents()

contains

  subroutine test_register_musica_species()
    use ccpp_kinds,                		 only: kind_phys
    use ccpp_constituent_prop_mod, 		 only: ccpp_constituent_properties_t
    use musica_ccpp_species,       		 only: musica_species_t
    use musica_ccpp_tuvx_load_species, only: configure_tuvx_species, &
                                             check_tuvx_species_initialization

    integer, parameter                               :: NUM_MICM_SPECIES = 6
    integer, parameter                               :: NUM_TUVX_CONSTITUENTS = 4
    type(musica_species_t)                           :: micm_species(NUM_MICM_SPECIES)
    type(musica_species_t),              allocatable :: tuvx_species(:)
    type(ccpp_constituent_properties_t)              :: micm_constituent_props(NUM_MICM_SPECIES)
    type(ccpp_constituent_properties_t), allocatable :: tuvx_constituent_props(:)
    character(len=512)                               :: errmsg
    integer                                          :: errcode
    real(kind_phys)                                  :: molar_mass_group(NUM_MICM_SPECIES) = &
      [0.1_kind_phys, 0.2_kind_phys, 0.3_kind_phys, 0.4_kind_phys, 0.5_kind_phys, 0.6_kind_phys]
    integer                                          :: i_species
    character(len=512)                               :: species_names(NUM_MICM_SPECIES)

    species_names(1) = 'N2'
    species_names(2) = 'O2'  ! shared species
    species_names(3) = 'FOO'
    species_names(4) = 'O1D'
    species_names(5) = 'BAZ'
    species_names(6) = 'O3'  ! shared species

    do i_species = 1, NUM_MICM_SPECIES
      call micm_constituent_props(i_species)%instantiate( &
      std_name = trim(species_names(i_species)), &
      long_name = trim(species_names(i_species)), &
      units = 'kg kg-1', &
      vertical_dim = 'vertical_layer_dimension', &
      default_value = 0.0_kind_phys, &
      min_value = 0.0_kind_phys, &
      molar_mass = molar_mass_group(i_species), &
      advected = .true., &
      errcode = errcode, &
      errmsg = errmsg)

      micm_species(i_species) = musica_species_t( &
        name = species_names(i_species), &
        unit = 'kg kg-1', &
        molar_mass = molar_mass_group(i_species), &
        index_musica_species = i_species )
    end do

    call configure_tuvx_species(micm_species, tuvx_species, tuvx_constituent_props, &
                                errmsg, errcode)
    ASSERT(errcode == 0)

    call register_musica_species(micm_species, tuvx_species)
    ASSERT(allocated(micm_species_set))
    ASSERT(allocated(tuvx_species_set))
    ASSERT(number_of_micm_species == NUM_MICM_SPECIES)
    ASSERT(number_of_tuvx_species == NUM_TUVX_CONSTITUENTS)

    call check_tuvx_species_initialization(errmsg, errcode)
    ASSERT(errcode == 0)

    do i_species = 1, size(micm_species)
      call micm_species(i_species)%deallocate()
    end do
    do i_species = 1, size(tuvx_species)
      call tuvx_species(i_species)%deallocate()
    end do
    deallocate(tuvx_species)
    deallocate(tuvx_constituent_props)
  
    call cleanup_musica_species()

  end subroutine test_register_musica_species

  subroutine test_initialize_musica_species_indices_and_molar_mass()
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t, ccpp_constituent_prop_ptr_t
    use musica_ccpp, 							 only: musica_ccpp_register
    use musica_ccpp_namelist,      only: filename_of_micm_configuration, &
                                         filename_of_tuvx_configuration, &
                                         filename_of_tuvx_micm_mapping_configuration

    type(ccpp_constituent_properties_t), allocatable, target :: constituent_props(:)
    type(ccpp_constituent_prop_ptr_t),   allocatable         :: constituent_props_ptr(:)
    type(ccpp_constituent_properties_t), pointer             :: const_prop
    character(len=512)               							           :: errmsg
    integer                          							           :: errcode
    integer                                                  :: i

    filename_of_micm_configuration = 'musica_configurations/chapman/micm/config.json'
    filename_of_tuvx_configuration = 'musica_configurations/chapman/tuvx/config.json'
    filename_of_tuvx_micm_mapping_configuration = 'musica_configurations/chapman/tuvx_micm_mapping.json'

    call musica_ccpp_register(constituent_props, errmsg, errcode)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
    ASSERT(allocated(constituent_props))

    allocate(constituent_props_ptr(size(constituent_props)))
    do i = 1, size(constituent_props)
      const_prop => constituent_props(i)
      call constituent_props_ptr(i)%set(const_prop, errcode, errmsg)
    end do

    call initialize_musica_species_indices(constituent_props_ptr, errmsg, errcode)
    ASSERT(errcode == 0)
    ASSERT(allocated(micm_indices_constituent_props))
    ASSERT(allocated(tuvx_indices_constituent_props))

    call initialize_molar_mass_array(constituent_props_ptr, errmsg, errcode)
    ASSERT(errcode == 0)
    ASSERT(allocated(micm_molar_mass_array))

    do i = 1, size(micm_species_set)
      ASSERT(micm_species_set(i)%index_musica_species == i)
      ASSERT(micm_species_set(i)%index_constituent_props == micm_indices_constituent_props(i))
      ASSERT(micm_species_set(i)%molar_mass == micm_molar_mass_array(i))
    end do

    do i = 1, size(tuvx_species_set)
      ASSERT(tuvx_species_set(i)%index_musica_species == i)
      ASSERT(tuvx_species_set(i)%index_constituent_props == tuvx_indices_constituent_props(i))
    end do

    call check_initialization(errmsg, errcode)
    ASSERT(errcode == 0)

    call cleanup_musica_species()

  end subroutine test_initialize_musica_species_indices_and_molar_mass


  subroutine test_extract_and_update_subset_constituents()
    use ccpp_kinds, only: kind_phys

    integer, parameter :: NUM_COLUMNS = 3
    integer, parameter :: NUM_LAYERS = 3
    integer, parameter :: NUM_SPECIES = 5
    integer, parameter :: NUM_SUBSET_SPECIES = 3
    integer            :: indices_array(NUM_SUBSET_SPECIES) = [2, 4, 5]
    real(kind_phys)    :: constituents(NUM_COLUMNS, NUM_LAYERS, NUM_SPECIES)
    real(kind_phys)    :: subset_constituents(NUM_COLUMNS, NUM_LAYERS, NUM_SUBSET_SPECIES)
    real(kind_phys)    :: expected_constituents(NUM_COLUMNS, NUM_LAYERS, NUM_SPECIES)
    real(kind_phys)    :: expected_subset_constituents(NUM_COLUMNS, NUM_LAYERS, NUM_SUBSET_SPECIES)
    character(len=512) :: errmsg
    integer            :: errcode
    integer            :: i, j, k

    constituents = reshape([1.0, 2.0, 3.0, 4.0, 5.0, &
                            6.0, 7.0, 8.0, 9.0, 10.0, &
                            11.0, 12.0, 13.0, 14.0, 15.0, &
                            16.0, 17.0, 18.0, 19.0, 20.0, &
                            21.0, 22.0, 23.0, 24.0, 25.0, &
                            26.0, 27.0, 28.0, 29.0, 30.0, &
                            31.0, 32.0, 33.0, 34.0, 35.0, &
                            36.0, 37.0, 38.0, 39.0, 40.0, &
                            41.0, 42.0, 43.0, 44.0, 45.0, &
                            46.0, 47.0, 48.0, 49.0, 50.0], shape(constituents))
    expected_constituents(:,:,:) = constituents(:,:,:)
    do k = 1, NUM_SUBSET_SPECIES
      do j = 1, NUM_LAYERS
        do i = 1, NUM_COLUMNS
          expected_subset_constituents(i,j,k) = constituents(i,j,indices_array(k))
        end do
      end do
    end do

    call extract_subset_constituents(indices_array, constituents, &
                                     subset_constituents, errmsg, errcode)
    ASSERT(errcode == 0)
    do k = 1, NUM_SUBSET_SPECIES
      do j = 1, NUM_LAYERS
        do i = 1, NUM_COLUMNS
          ASSERT(expected_subset_constituents(i,j,k) == subset_constituents(i,j,k))
        end do
      end do
    end do

    do k = 1, NUM_SUBSET_SPECIES
      do j = 1, NUM_LAYERS
        do i = 1, NUM_COLUMNS
          subset_constituents(i,j,k) = subset_constituents(i,j,k) + 100.0_kind_phys
          expected_constituents(i,j,indices_array(k)) = &
                      expected_constituents(i,j,indices_array(k)) + 100.0_kind_phys
        end do
      end do
    end do

    call update_constituents(indices_array, subset_constituents, &
                            constituents, errmsg, errcode)
    ASSERT(errcode == 0)
    do k = 1, NUM_SPECIES
      do j = 1, NUM_LAYERS
        do i = 1, NUM_COLUMNS
          ASSERT(expected_constituents(i,j,k) == constituents(i,j,k))
        end do
      end do
    end do

  end subroutine test_extract_and_update_subset_constituents

end program test_musica_species