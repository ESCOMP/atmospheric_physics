! Copyright (C) 2024-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
program test_tuvx_load_species

  use musica_ccpp_tuvx_load_species

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  call test_configure_shared_gas_species_tuvx_micm()
  call test_configure_partial_shared_gas_species()
  call test_configure_no_shared_gas_species()

contains

  subroutine test_configure_shared_gas_species_tuvx_micm()
    ! There are three gas species required for TUV-x: dry air, O2, and O3.
    ! This test focuses on configuring MUSICA species and constituent properties
    ! when the MICM species include all of these. Cloud liquid water content
    ! is the only component specific to TUVX.
    use ccpp_kinds,                only: kind_phys
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use musica_ccpp_species,       only: musica_species_t, MUSICA_INT_UNASSIGNED
    use musica_util,               only: error_t

    integer, parameter                               :: NUM_MICM_SPECIES = 6
    integer, parameter                               :: NUM_SHARED_SPECIES_BETWEEN_MICM_TUVX = 3
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
    character(len=512)                               :: name, unit, species_name
    real(kind_phys)                                  :: molar_mass   ! kg mol-1
    real(kind_phys)                                  :: scale_height ! km
    integer                                          :: index_musica, index_constituent_props
    logical                                          :: is_advected, tmp_bool, has_profile

    species_names(1) = 'N2'
    species_names(2) = 'O2'  ! shared species
    species_names(3) = 'FOO'
    species_names(4) = 'O1D'
    species_names(5) = 'air' ! shared species
    species_names(6) = 'O3'  ! shared species

    do i_species = 1, NUM_MICM_SPECIES
      call micm_constituent_props(i_species)%instantiate( &
      std_name = trim(species_names(i_species)), &
      long_name = trim(species_names(i_species)), &
      diag_name = trim(species_names(i_species)), &
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

    call configure_tuvx_species( micm_species, tuvx_species, tuvx_constituent_props, &
                                 errmsg, errcode )
    ASSERT(errcode == 0)
    ASSERT(allocated(tuvx_constituent_props))
    ASSERT(size(tuvx_constituent_props) == NUM_TUVX_CONSTITUENTS - NUM_SHARED_SPECIES_BETWEEN_MICM_TUVX)
    do i_species = 1, size(tuvx_constituent_props)
      ASSERT(tuvx_constituent_props(i_species)%is_instantiated(errcode, errmsg))
      call tuvx_constituent_props(i_species)%standard_name(species_name, errcode, errmsg)
      ASSERT(errcode == 0)
      call tuvx_constituent_props(i_species)%molar_mass(molar_mass, errcode, errmsg)
      ASSERT(errcode == 0)
      call tuvx_constituent_props(i_species)%is_advected(is_advected, errcode, errmsg)
      ASSERT(errcode == 0)
      tmp_bool = (trim(species_name) == 'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water' .and. &
                  molar_mass == 0.018_kind_phys .and. is_advected)
      ASSERT(tmp_bool)
    end do

    ASSERT(allocated(tuvx_species))
    ASSERT(size(tuvx_species) == NUM_TUVX_CONSTITUENTS)
    do i_species = 1, size(tuvx_species)
      name = tuvx_species(i_species)%name
      unit = tuvx_species(i_species)%unit
      molar_mass = tuvx_species(i_species)%molar_mass
      scale_height = tuvx_species(i_species)%scale_height
      index_musica = tuvx_species(i_species)%index_musica_species
      index_constituent_props = tuvx_species(i_species)%index_constituent_props
      has_profile = tuvx_species(i_species)%profiled
      tmp_bool = (trim(name) == 'air' .and. trim(unit) == 'molecule cm-3' .and. molar_mass == 0.0289644_kind_phys .and. &
                scale_height == 8.01_kind_phys .and. index_musica == 2 .and. index_constituent_props == MUSICA_INT_UNASSIGNED &
                .and. has_profile) .or.  &
                (trim(name) == 'O2' .and. trim(unit) == 'molecule cm-3' .and. molar_mass == 0.0319988_kind_phys  .and. &
                scale_height == 7.0_kind_phys .and. index_musica == 3 .and. index_constituent_props == MUSICA_INT_UNASSIGNED &
                .and. has_profile) .or.  &
                (trim(name) == 'O3' .and. trim(unit) == 'molecule cm-3' .and. molar_mass == 0.0479982_kind_phys .and. &
                scale_height == 7.0_kind_phys .and. index_musica == 4 .and. index_constituent_props == MUSICA_INT_UNASSIGNED &
                .and. has_profile) .or.  &
                (trim(name) == 'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water' .and. &
                trim(unit) == 'kg kg-1' .and. molar_mass == 0.018_kind_phys .and. &
                scale_height == 0.0_kind_phys .and. index_musica == 1 .and. index_constituent_props == MUSICA_INT_UNASSIGNED &
                .and. .not. has_profile)
      ASSERT(tmp_bool)
    end do

    call check_tuvx_species_initialization( errmsg, errcode )
    ASSERT(errcode == 0)
    ASSERT(index_cloud_liquid_water_content == 1)
    ASSERT(index_dry_air == 2)
    ASSERT(index_O2 == 3)
    ASSERT(index_O3 == 4)

    do i_species = 1, size(tuvx_species)
      call tuvx_species(i_species)%deallocate()
    end do
    deallocate(tuvx_species)
    deallocate(tuvx_constituent_props)

  end subroutine test_configure_shared_gas_species_tuvx_micm

  subroutine test_configure_partial_shared_gas_species()
    ! This test case applies when some gas species are registered in MICM.
    ! It checks which species are already registered and which are not,
    ! adding only the new species to the constituent properties.
    use ccpp_kinds,                only: kind_phys
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use musica_ccpp_species,       only: musica_species_t, MUSICA_INT_UNASSIGNED
    use musica_util,               only: error_t

    integer, parameter                               :: NUM_MICM_SPECIES = 6
    integer, parameter                               :: NUM_SHARED_SPECIES_BETWEEN_MICM_TUVX = 2
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
    character(len=512)                               :: name, unit, species_name
    real(kind_phys)                                  :: molar_mass   ! kg mol-1
    real(kind_phys)                                  :: scale_height ! km
    integer                                          :: index_musica, index_constituent_props
    logical                                          :: is_advected, tmp_bool, has_profile

    species_names(1) = 'N2'
    species_names(2) = 'O2' ! shared species
    species_names(3) = 'FOO'
    species_names(4) = 'O1D'
    species_names(5) = 'BAZ'
    species_names(6) = 'O3' ! shared species

    do i_species = 1, NUM_MICM_SPECIES
      call micm_constituent_props(i_species)%instantiate( &
      std_name = trim(species_names(i_species)), &
      long_name = trim(species_names(i_species)), &
      diag_name = trim(species_names(i_species)), &
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

    call configure_tuvx_species( micm_species, tuvx_species, tuvx_constituent_props, &
                                 errmsg, errcode )
    ASSERT(errcode == 0)
    ASSERT(allocated(tuvx_constituent_props))
    ASSERT(size(tuvx_constituent_props) == NUM_TUVX_CONSTITUENTS - NUM_SHARED_SPECIES_BETWEEN_MICM_TUVX)
    do i_species = 1, size(tuvx_constituent_props)
      ASSERT(tuvx_constituent_props(i_species)%is_instantiated(errcode, errmsg))
      call tuvx_constituent_props(i_species)%standard_name(species_name, errcode, errmsg)
      ASSERT(errcode == 0)
      call tuvx_constituent_props(i_species)%molar_mass(molar_mass, errcode, errmsg)
      ASSERT(errcode == 0)
      call tuvx_constituent_props(i_species)%is_advected(is_advected, errcode, errmsg)
      ASSERT(errcode == 0)
      tmp_bool = (trim(species_name) == 'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water' .and. &
                  molar_mass == 0.018_kind_phys .and. is_advected) .or. &
                 (trim(species_name) == 'air' .and. molar_mass == 0.0289644_kind_phys .and. .not. is_advected)
      ASSERT(tmp_bool)
    end do

    ASSERT(allocated(tuvx_species))
    ASSERT(size(tuvx_species) == NUM_TUVX_CONSTITUENTS)
    do i_species = 1, size(tuvx_species)
      name = tuvx_species(i_species)%name
      unit = tuvx_species(i_species)%unit
      molar_mass = tuvx_species(i_species)%molar_mass
      scale_height = tuvx_species(i_species)%scale_height
      index_musica = tuvx_species(i_species)%index_musica_species
      index_constituent_props = tuvx_species(i_species)%index_constituent_props
      has_profile = tuvx_species(i_species)%profiled
      tmp_bool = (trim(name) == 'air' .and. trim(unit) == 'molecule cm-3' .and. molar_mass == 0.0289644_kind_phys .and. &
                scale_height == 8.01_kind_phys .and. index_musica == 2 .and. index_constituent_props == MUSICA_INT_UNASSIGNED &
                .and. has_profile) .or.  &
                (trim(name) == 'O2' .and. trim(unit) == 'molecule cm-3' .and. molar_mass == 0.0319988_kind_phys  .and. &
                scale_height == 7.0_kind_phys .and. index_musica == 3 .and. index_constituent_props == MUSICA_INT_UNASSIGNED &
                .and. has_profile) .or.  &
                (trim(name) == 'O3' .and. trim(unit) == 'molecule cm-3' .and. molar_mass == 0.0479982_kind_phys .and. &
                scale_height == 7.0_kind_phys .and. index_musica == 4 .and. index_constituent_props == MUSICA_INT_UNASSIGNED &
                .and. has_profile) .or.  &
                (trim(name) == 'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water' .and. &
                trim(unit) == 'kg kg-1' .and. molar_mass == 0.018_kind_phys .and. &
                scale_height == 0.0_kind_phys .and. index_musica == 1 .and. index_constituent_props == MUSICA_INT_UNASSIGNED &
                .and. .not. has_profile)
      ASSERT(tmp_bool)
    end do

    call check_tuvx_species_initialization( errmsg, errcode )
    ASSERT(errcode == 0)
    ASSERT(index_cloud_liquid_water_content == 1)
    ASSERT(index_dry_air == 2)
    ASSERT(index_O2 == 3)
    ASSERT(index_O3 == 4)

    do i_species = 1, size(tuvx_species)
      call tuvx_species(i_species)%deallocate()
    end do
    deallocate(tuvx_species)
    deallocate(tuvx_constituent_props)

  end subroutine test_configure_partial_shared_gas_species

  subroutine test_configure_no_shared_gas_species()
    ! This test case applies when there are no shared species between MICM and TUV-x.
    ! All configured components are added to the constituent properties.
    use ccpp_kinds,                only: kind_phys
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use musica_ccpp_species,       only: musica_species_t, MUSICA_INT_UNASSIGNED
    use musica_util,               only: error_t

    integer, parameter                               :: NUM_MICM_SPECIES = 6
    integer, parameter                               :: NUM_SHARED_SPECIES_BETWEEN_MICM_TUVX = 0
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
    character(len=512)                               :: name, unit, species_name
    real(kind_phys)                                  :: molar_mass   ! kg mol-1
    real(kind_phys)                                  :: scale_height ! km
    integer                                          :: index_musica, index_constituent_props
    logical                                          :: is_advected, tmp_bool, has_profile

    species_names(1) = 'N2'
    species_names(2) = 'BAR'
    species_names(3) = 'FOO'
    species_names(4) = 'O1D'
    species_names(5) = 'BAZ'
    species_names(6) = 'BOB'

    do i_species = 1, NUM_MICM_SPECIES
      call micm_constituent_props(i_species)%instantiate( &
      std_name = trim(species_names(i_species)), &
      long_name = trim(species_names(i_species)), &
      diag_name = trim(species_names(i_species)), &
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

    call configure_tuvx_species( micm_species, tuvx_species, tuvx_constituent_props, &
                                 errmsg, errcode )
    ASSERT(errcode == 0)
    ASSERT(allocated(tuvx_constituent_props))
    ASSERT(size(tuvx_constituent_props) == NUM_TUVX_CONSTITUENTS - NUM_SHARED_SPECIES_BETWEEN_MICM_TUVX)
    do i_species = 1, size(tuvx_constituent_props)
      ASSERT(tuvx_constituent_props(i_species)%is_instantiated(errcode, errmsg))
      call tuvx_constituent_props(i_species)%standard_name(species_name, errcode, errmsg)
      ASSERT(errcode == 0)
      call tuvx_constituent_props(i_species)%molar_mass(molar_mass, errcode, errmsg)
      ASSERT(errcode == 0)
      call tuvx_constituent_props(i_species)%is_advected(is_advected, errcode, errmsg)
      ASSERT(errcode == 0)
      tmp_bool = (trim(species_name) == 'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water' .and. &
                  molar_mass == 0.018_kind_phys .and. is_advected) .or. &
                 (trim(species_name) == 'air' .and. molar_mass == 0.0289644_kind_phys .and. .not. is_advected) .or. &
                 (trim(species_name) == 'O2' .and. molar_mass == 0.0319988_kind_phys .and. .not. is_advected) .or. &
                 (trim(species_name) == 'O3' .and. molar_mass == 0.0479982_kind_phys .and. .not. is_advected)
      ASSERT(tmp_bool)
    end do

    ASSERT(allocated(tuvx_species))
    ASSERT(size(tuvx_species) == NUM_TUVX_CONSTITUENTS)
    do i_species = 1, size(tuvx_species)
      name = tuvx_species(i_species)%name
      unit = tuvx_species(i_species)%unit
      molar_mass = tuvx_species(i_species)%molar_mass
      scale_height = tuvx_species(i_species)%scale_height
      index_musica = tuvx_species(i_species)%index_musica_species
      index_constituent_props = tuvx_species(i_species)%index_constituent_props
      has_profile = tuvx_species(i_species)%profiled
      tmp_bool = (trim(name) == 'air' .and. trim(unit) == 'molecule cm-3' .and. molar_mass == 0.0289644_kind_phys .and. &
                scale_height == 8.01_kind_phys .and. index_musica == 2 .and. index_constituent_props == MUSICA_INT_UNASSIGNED &
                .and. has_profile) .or.  &
                (trim(name) == 'O2' .and. trim(unit) == 'molecule cm-3' .and. molar_mass == 0.0319988_kind_phys  .and. &
                scale_height == 7.0_kind_phys .and. index_musica == 3 .and. index_constituent_props == MUSICA_INT_UNASSIGNED &
                .and. has_profile) .or.  &
                (trim(name) == 'O3' .and. trim(unit) == 'molecule cm-3' .and. molar_mass == 0.0479982_kind_phys .and. &
                scale_height == 7.0_kind_phys .and. index_musica == 4 .and. index_constituent_props == MUSICA_INT_UNASSIGNED &
                .and. has_profile) .or.  &
                (trim(name) == 'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water' .and. &
                trim(unit) == 'kg kg-1' .and. molar_mass == 0.018_kind_phys .and. &
                scale_height == 0.0_kind_phys .and. index_musica == 1 .and. index_constituent_props == MUSICA_INT_UNASSIGNED &
                .and. .not. has_profile)
      ASSERT(tmp_bool)
    end do

    call check_tuvx_species_initialization( errmsg, errcode )
    ASSERT(errcode == 0)
    ASSERT(index_cloud_liquid_water_content == 1)
    ASSERT(index_dry_air == 2)
    ASSERT(index_O2 == 3)
    ASSERT(index_O3 == 4)

    do i_species = 1, size(tuvx_species)
      call tuvx_species(i_species)%deallocate()
    end do
    deallocate(tuvx_species)
    deallocate(tuvx_constituent_props)

  end subroutine test_configure_no_shared_gas_species

end program test_tuvx_load_species
