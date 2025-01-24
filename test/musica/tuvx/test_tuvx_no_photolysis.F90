program test_tuvx_no_photolysis

  use musica_ccpp_tuvx_no_photolysis_rate

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  call test_NO_and_N2_absence()
  call test_only_NO_exist()
  call test_only_N2_exist()
  call test_NO_and_N2_exist()
  call test_initialize_NO_constituents_indices_and_molar_mass()
  call test_calculate_NO_photolysis_rate()

contains

  subroutine test_NO_and_N2_absence()
    use ccpp_kinds,          only: kind_phys
    use musica_ccpp_species, only: musica_species_t

    integer, parameter     :: NUM_MICM_SPECIES = 3
    type(musica_species_t) :: micm_species(NUM_MICM_SPECIES)
    real(kind_phys)        :: molar_mass_group(NUM_MICM_SPECIES) = &
                              [0.1_kind_phys, 0.2_kind_phys, 0.3_kind_phys]
    integer                :: i_species
    character(len=512)     :: species_names(NUM_MICM_SPECIES)

    species_names(1) = 'O2'
    species_names(2) = 'O1D'
    species_names(3) = 'O3'

    do i_species = 1, NUM_MICM_SPECIES
      micm_species(i_species) = musica_species_t( &
        name = species_names(i_species), &
        unit = 'kg kg-1', &
        molar_mass = molar_mass_group(i_species), &
        index_musica_species = i_species )
    end do

    call check_NO_exist( micm_species )
    ASSERT(.not. is_NO_photolysis_active)

  end subroutine test_NO_and_N2_absence

  subroutine test_only_NO_exist()
    use ccpp_kinds,          only: kind_phys
    use musica_ccpp_species, only: musica_species_t

    integer, parameter     :: NUM_MICM_SPECIES = 3
    type(musica_species_t) :: micm_species(NUM_MICM_SPECIES)
    real(kind_phys)        :: molar_mass_group(NUM_MICM_SPECIES) = &
                              [0.1_kind_phys, 0.2_kind_phys, 0.3_kind_phys]
    integer                :: i_species
    character(len=512)     :: species_names(NUM_MICM_SPECIES)

    species_names(1) = 'NO'
    species_names(2) = 'O1D'
    species_names(3) = 'O3'

    do i_species = 1, NUM_MICM_SPECIES
      micm_species(i_species) = musica_species_t( &
        name = species_names(i_species), &
        unit = 'kg kg-1', &
        molar_mass = molar_mass_group(i_species), &
        index_musica_species = i_species )
    end do

    call check_NO_exist( micm_species )
    ASSERT(.not. is_NO_photolysis_active)

  end subroutine test_only_NO_exist

  subroutine test_only_N2_exist()
    use ccpp_kinds,          only: kind_phys
    use musica_ccpp_species, only: musica_species_t

    integer, parameter     :: NUM_MICM_SPECIES = 3
    type(musica_species_t) :: micm_species(NUM_MICM_SPECIES)
    real(kind_phys)        :: molar_mass_group(NUM_MICM_SPECIES) = &
                              [0.1_kind_phys, 0.2_kind_phys, 0.3_kind_phys]
    integer                :: i_species
    character(len=512)     :: species_names(NUM_MICM_SPECIES)

    species_names(1) = 'N2'
    species_names(2) = 'O1D'
    species_names(3) = 'O3'

    do i_species = 1, NUM_MICM_SPECIES
      micm_species(i_species) = musica_species_t( &
        name = species_names(i_species), &
        unit = 'kg kg-1', &
        molar_mass = molar_mass_group(i_species), &
        index_musica_species = i_species )
    end do

    call check_NO_exist( micm_species )
    ASSERT(.not. is_NO_photolysis_active)

  end subroutine test_only_N2_exist

  subroutine test_NO_and_N2_exist()
    use ccpp_kinds,          only: kind_phys
    use musica_ccpp_species, only: musica_species_t

    integer, parameter     :: NUM_MICM_SPECIES = 3
    type(musica_species_t) :: micm_species(NUM_MICM_SPECIES)
    real(kind_phys)        :: molar_mass_group(NUM_MICM_SPECIES) = &
                              [0.0280134_kind_phys, 0.2_kind_phys, 0.0300100_kind_phys]
    integer                :: i_species
    character(len=512)     :: species_names(NUM_MICM_SPECIES)

    species_names(1) = 'N2'
    species_names(2) = 'O2'
    species_names(3) = 'NO'

    do i_species = 1, NUM_MICM_SPECIES
      micm_species(i_species) = musica_species_t( &
        name = species_names(i_species), &
        unit = 'kg kg-1', &
        molar_mass = molar_mass_group(i_species), &
        index_musica_species = i_species )
    end do

    call check_NO_exist( micm_species )
    ASSERT(is_NO_photolysis_active)
    ASSERT(MOLAR_MASS_N2 == 0.0280134_kind_phys)
    ASSERT(MOLAR_MASS_NO == 0.0300100_kind_phys)

  end subroutine test_NO_and_N2_exist

  subroutine test_initialize_NO_constituents_indices_and_molar_mass()
    use ccpp_kinds,                only: kind_phys
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t, ccpp_constituent_prop_ptr_t
    use musica_ccpp, 							 only: musica_ccpp_register, musica_ccpp_final
    use musica_ccpp_namelist,      only: filename_of_micm_configuration, &
                                         filename_of_tuvx_configuration, &
                                         filename_of_tuvx_micm_mapping_configuration

    integer,                             parameter           :: expected_index_N2 = 1
    integer,                             parameter           :: expected_index_NO = 2
    real(kind_phys),                     parameter           :: expected_MOLAR_MASS_N2 = 0.0280134_kind_phys
    real(kind_phys),                     parameter           :: expected_MOLAR_MASS_NO = 0.0300100_kind_phys
    type(ccpp_constituent_properties_t), allocatable, target :: constituent_props(:)
    type(ccpp_constituent_prop_ptr_t),   allocatable         :: constituent_props_ptr(:)
    type(ccpp_constituent_properties_t), pointer             :: const_prop
    character(len=512)               							           :: errmsg
    integer                          							           :: errcode
    integer                                                  :: i

    filename_of_micm_configuration = 'musica_configurations/NO_photolysis/micm/config.json'
    filename_of_tuvx_configuration = 'musica_configurations/NO_photolysis/tuvx/config.json'
    filename_of_tuvx_micm_mapping_configuration = 'musica_configurations/NO_photolysis/tuvx_micm_mapping.json'

    call musica_ccpp_register( constituent_props, errmsg, errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
    ASSERT(allocated(constituent_props))
    ASSERT(is_NO_photolysis_active)
    ASSERT(MOLAR_MASS_N2 == expected_MOLAR_MASS_N2)
    ASSERT(MOLAR_MASS_NO == expected_MOLAR_MASS_NO)

    allocate(constituent_props_ptr(size(constituent_props)))
    do i = 1, size(constituent_props)
      const_prop => constituent_props(i)
      call constituent_props_ptr(i)%set( const_prop, errcode, errmsg )
    end do

    call set_NO_index_constituent_props(constituent_props_ptr, errmsg, errcode)
    ASSERT(errcode == 0)
    ASSERT(index_N2 == expected_index_N2)
    ASSERT(index_NO == expected_index_NO)
    ASSERT(NO_photolysis_indices_constituent_props(1) == expected_index_N2)
    ASSERT(NO_photolysis_indices_constituent_props(2) == expected_index_NO)

    call check_NO_initialization(errmsg, errcode)
    ASSERT(errcode == 0)

    call musica_ccpp_final(errmsg, errcode)
    ASSERT(errcode == 0)

  end subroutine test_initialize_NO_constituents_indices_and_molar_mass
  
  subroutine test_calculate_NO_photolysis_rate()
    use ccpp_kinds,                          only: kind_phys
    use ccpp_constituent_prop_mod,           only: ccpp_constituent_prop_ptr_t
    use ccpp_constituent_prop_mod,           only: ccpp_constituent_properties_t
    use musica_ccpp_micm,                    only: micm
    use musica_ccpp, 							           only: musica_ccpp_register, musica_ccpp_init, musica_ccpp_final
    use musica_ccpp_namelist,                only: filename_of_micm_configuration, &
                                                   filename_of_tuvx_configuration, &
                                                   filename_of_tuvx_micm_mapping_configuration
    use musica_ccpp_tuvx_load_species,       only: index_O2, index_O3
    use musica_ccpp_species,                 only: tuvx_indices_constituent_props, extract_subset_constituents
    use musica_ccpp_tuvx_no_photolysis_rate, only: NO_photolysis_indices_constituent_props
    use test_musica_data,                    only: get_wavelength_edges

    implicit none

    integer, parameter                                                  :: NUM_SPECIES = 10
    integer, parameter                                                  :: NUM_TUVX_CONSTITUENTS = 1     ! cloud liquid water content
    integer, parameter                                                  :: NUM_TUVX_ONLY_GAS_SPECIES = 1 ! dry air
    integer, parameter                                                  :: NUM_TUVX_CONSTITUENTS_INCLUDE_SHARED_SPECIES = 4
    integer, parameter                                                  :: NUM_NO_CONSTITUENTS = 2
    ! This test requires that the number of grid cells = 4, which is the default
    ! vector dimension for MICM. This restriction will be removed once
    ! https://github.com/NCAR/musica/issues/217 is finished.
    integer, parameter                                                  :: NUM_COLUMNS = 2
    integer, parameter                                                  :: NUM_LAYERS = 2
    integer, parameter                                                  :: NUM_WAVELENGTH_BINS = 102
    integer                                                             :: errcode
    character(len=512)                                                  :: errmsg
    real(kind_phys), dimension(NUM_WAVELENGTH_BINS+1)                   :: photolysis_wavelength_grid_interfaces        ! m
    integer, parameter                                                  :: num_photolysis_wavelength_grid_sections = 8  ! (count)
    real(kind_phys), dimension(num_photolysis_wavelength_grid_sections) :: extraterrestrial_flux                        ! photons cm-2 s-1 nm-1
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS)                  :: dry_air_density                              ! kg m-3
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS, &
           NUM_SPECIES+NUM_TUVX_CONSTITUENTS+NUM_TUVX_ONLY_GAS_SPECIES) :: constituents                                 ! kg kg-1
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS, &
                          NUM_TUVX_CONSTITUENTS_INCLUDE_SHARED_SPECIES) :: constituents_tuvx_species                    ! kg kg-1
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS, &
                                           NUM_NO_CONSTITUENTS)         :: constituents_NO_photolysis                   ! kg kg-1
    real(kind_phys), dimension(NUM_COLUMNS)                             :: solar_zenith_angle                           ! radians
    real(kind_phys), dimension(NUM_LAYERS+1)                            :: height_at_interfaces                         ! km
    real(kind_phys), dimension(NUM_LAYERS)                              :: photolysis_rate_NO(NUM_LAYERS)
    type(ccpp_constituent_prop_ptr_t),   allocatable                    :: constituent_props_ptr(:)
    type(ccpp_constituent_properties_t), allocatable, target            :: constituent_props(:)
    type(ccpp_constituent_properties_t), pointer                        :: const_prop
    real(kind_phys)                                                     :: base_conc
    character(len=:), allocatable                                       :: micm_species_name
    integer                                                             :: i, j, k, i_col
    integer                                                             :: O2_index, Ar_index, CO2_index, H2O_index, N_index
    integer                                                             :: NO_index, O_index, O1D_index, O3_index, N2_index

    filename_of_micm_configuration = 'musica_configurations/NO_photolysis/micm/config.json'
    filename_of_tuvx_configuration = 'musica_configurations/NO_photolysis/tuvx/config.json'
    filename_of_tuvx_micm_mapping_configuration = 'musica_configurations/NO_photolysis/tuvx_micm_mapping.json'

    call get_wavelength_edges(photolysis_wavelength_grid_interfaces)
    dry_air_density(:,1) = (/ 3.5_kind_phys, 4.5_kind_phys /)
    dry_air_density(:,2) = (/ 5.5_kind_phys, 6.5_kind_phys /)
    extraterrestrial_flux(:) = &
      (/ 1.5e13_kind_phys, 1.5e13_kind_phys, 1.4e13_kind_phys, 1.4e13_kind_phys, &
        1.3e13_kind_phys, 1.2e13_kind_phys, 1.1e13_kind_phys, 1.0e13_kind_phys /)
    solar_zenith_angle(:) = (/ 0.0_kind_phys, 2.1_kind_phys /)
    height_at_interfaces(:) = (/ 3.01_kind_phys, 1.01_kind_phys, 0.01_kind_phys /)

    call musica_ccpp_register( constituent_props, errmsg, errcode )
    ASSERT(errcode == 0)

    allocate(constituent_props_ptr(size(constituent_props)))
    do i = 1, size(constituent_props)
      const_prop => constituent_props(i)
      call constituent_props_ptr(i)%set( const_prop, errcode, errmsg )
      ASSERT(errcode == 0)
    end do

    call musica_ccpp_init( NUM_COLUMNS, NUM_LAYERS, NUM_LAYERS+1, photolysis_wavelength_grid_interfaces, &
                          constituent_props_ptr, errmsg, errcode )
    ASSERT(errcode == 0)

    do i = 1, micm%species_ordering%size()
      micm_species_name = micm%species_ordering%name(i)
      if (micm_species_name == "O2") then
        O2_index = i
        base_conc = 0.21_kind_phys
      else if (micm_species_name == "Ar") then
        Ar_index = i
        base_conc = 1.0e-9_kind_phys
      else if (micm_species_name == "CO2") then
        CO2_index = i
        base_conc = 1.0e-9_kind_phys
      else if (micm_species_name == "H2O") then
        H2O_index = i
        base_conc = 1.0e-9_kind_phys
      else if (micm_species_name == "N") then
        N_index = i
        base_conc = 1.0e-9_kind_phys
      else if (micm_species_name == "NO") then
        NO_index = i
        base_conc = 1.0e-9_kind_phys
      else if (micm_species_name == "O") then
        O_index = i
        base_conc = 1.0e-9_kind_phys
      else if (micm_species_name == "O1D") then
        O1D_index = i
        base_conc = 1.0e-9_kind_phys
      else if (micm_species_name == "O3") then
        O3_index = i
        base_conc = 1.0e-4_kind_phys
      else if (micm_species_name == "N2") then
        N2_index = i
        base_conc = 0.79_kind_phys
      else
        write(*,*) "Unknown species: ", micm_species_name
        stop 3
      endif
      do j = 1, NUM_COLUMNS
        do k = 1, NUM_LAYERS
          constituents(j,k,i) = base_conc * (1.0 + 0.1 * (j-1) + 0.01 * (k-1))
        end do
      end do
    end do
    ! set initial cloud liquid water mixing ratio to ~1e-3 kg kg-1
    do j = 1, NUM_COLUMNS
      do k = 1, NUM_LAYERS
        constituents(j,k,NUM_SPECIES+NUM_TUVX_CONSTITUENTS) = 1.0e-3_kind_phys * (1.0 + 0.1 * (j-1) + 0.01 * (k-1))
      end do
    end do
    do j = 1, NUM_COLUMNS
      do k = 1, NUM_LAYERS
        constituents(j,k,NUM_SPECIES+NUM_TUVX_CONSTITUENTS+NUM_TUVX_ONLY_GAS_SPECIES) = &
          1.0e-3_kind_phys * (1.0 + 0.1 * (j-1) + 0.01 * (k-1))
      end do
    end do

    call extract_subset_constituents(tuvx_indices_constituent_props, constituents, &
        constituents_tuvx_species, errmsg, errcode)
    ASSERT(errcode == 0)

    call extract_subset_constituents(NO_photolysis_indices_constituent_props, constituents, &
        constituents_NO_photolysis, errmsg, errcode)
    ASSERT(errcode == 0)

    do i_col = 1, size(dry_air_density, dim=1)
      photolysis_rate_NO = calculate_NO_photolysis_rate( solar_zenith_angle(i_col), & 
                                       extraterrestrial_flux, height_at_interfaces, & 
                                       dry_air_density(i_col,:),                    &
                                       constituents_tuvx_species(i_col,:,index_O2), &
                                       constituents_tuvx_species(i_col,:,index_O3), &
                                       constituents_NO_photolysis(i_col,:,:) )
    end do

    call musica_ccpp_final(errmsg, errcode)
    ASSERT(errcode == 0)

  end subroutine test_calculate_NO_photolysis_rate

  ! subroutine test_calculate_NO_photolysis_rate()
  !   use ccpp_kinds, only: kind_phys
  !   implicit none

  !   integer, parameter :: number_of_vertical_layers = 2

  !   real(kind_phys) :: solar_zenith_angle
  !   real(kind_phys), dimension(4) :: extraterrestrial_flux 
  !   real(kind_phys), dimension(1,number_of_vertical_layers,4) :: constituents ! (column, layers, constituents)  kg kg-1
  !   real(kind_phys), dimension(1,number_of_vertical_layers,2) :: constituents_NO_photolysis
  !   real(kind_phys), dimension(3) :: height_at_interfaces ! (layers + 1) km
  !   real(kind_phys), dimension(1,number_of_vertical_layers) :: dry_air_density ! (column, layers) kg m-3
  !   integer :: N2_index, O2_index, O3_index, NO_index
  !   real(kind_phys) :: molar_mass_N2, molar_mass_O2, molar_mass_O3, molar_mass_NO ! kg mol-1
  !   real(kind_phys), dimension(number_of_vertical_layers) :: jNO
  !   ! Some of this initialization data corresponds to values listed in the Minschwaner and Siskind (1993) paper
  !   ! see the musica_ccpp_tuvx_no_photolysis_rate.F90 file for the full citation

  !   ! Initialize test data
  !   solar_zenith_angle = 60.0_kind_phys
  !   extraterrestrial_flux = (/ 1.5e13_kind_phys, 1.4e13_kind_phys, 1.3e13_kind_phys, 1.2e13_kind_phys /)
  !   height_at_interfaces = reshape([40.0_kind_phys, 30.0_kind_phys, 20.0_kind_phys], shape(height_at_interfaces))
  !   dry_air_density = reshape([1.2_kind_phys, 1.1_kind_phys], shape(dry_air_density))
  !   N2_index = 1
  !   O2_index = 2
  !   O3_index = 3
  !   NO_index = 4
  !   molar_mass_N2 = 0.0280134_kind_phys
  !   molar_mass_O2 = 0.0319988_kind_phys
  !   molar_mass_O3 = 0.0479982_kind_phys
  !   molar_mass_NO = 0.3000610_kind_phys

  !   constituents(1,:,N2_index) = 0.78_kind_phys
  !   constituents(1,:,O2_index) = 0.21_kind_phys
  !   constituents(1,:,O3_index) = 20.0e-9_kind_phys
  !   constituents(1,:,NO_index) = 0.3e-9_kind_phys
  !   constituents_NO_photolysis(1,:,1) = constituents(1,:,N2_index)
  !   constituents_NO_photolysis(1,:,2) = constituents(1,:,NO_index)

  !   ! Call the function to test
  !   ! jNO = calculate_NO_photolysis_rate(size(constituents, dim=2), solar_zenith_angle, extraterrestrial_flux, constituents(1,:,:), height_at_interfaces, &
  !   !                                    dry_air_density(1,:), N2_index, O2_index, O3_index, NO_index, molar_mass_N2, molar_mass_O2, molar_mass_O3, molar_mass_NO)

  !   jNO = calculate_NO_photolysis_rate(size(constituents, dim=2), &
  !           solar_zenith_angle, extraterrestrial_flux, height_at_interfaces, &
  !           dry_air_density(1,:), constituents(1,:,O2_index), &
  !           constituents(1,:,O3_index), constituents_NO_photolysis(1,:,:))

  !   ! Validate the results
  !   print *, jNO
  !   ASSERT(jNO(1) .ne. 0.0_kind_phys)
  !   ASSERT(jNO(2) .ne. 0.0_kind_phys)
  !   ASSERT(jNO(1) .lt. 1.0_kind_phys)
  !   ASSERT(jNO(2) .lt. 1.0_kind_phys)

  ! end subroutine test_calculate_NO_photolysis_rate

end program test_tuvx_no_photolysis