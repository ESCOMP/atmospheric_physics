module musica_ccpp_tuvx_load_species
  use ccpp_kinds,          only: kind_phys
  use musica_ccpp_species, only: MUSICA_INT_UNASSIGNED

  implicit none
  private

  public :: configure_tuvx_species, check_tuvx_species_initialization
  
  integer, protected, public :: index_cloud_liquid_water_content = MUSICA_INT_UNASSIGNED
  integer, protected, public :: index_dry_air = MUSICA_INT_UNASSIGNED
  integer, protected, public :: index_O2 = MUSICA_INT_UNASSIGNED
  integer, protected, public :: index_O3 = MUSICA_INT_UNASSIGNED

  ! Constants
  ! Cloud liquid water
  character(len=*), parameter, public :: CLOUD_LIQUID_WATER_CONTENT_LABEL = &
    'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water'
  character(len=*), parameter, public :: CLOUD_LIQUID_WATER_CONTENT_LONG_NAME = &
    'cloud water mass mixing ratio with respect to moist air plus all airborne condensates'
  character(len=*), parameter, public :: CLOUD_LIQUID_WATER_CONTENT_UNITS = 'kg kg-1'
  real(kind_phys),  parameter, public :: CLOUD_LIQUID_WATER_CONTENT_MOLAR_MASS = 0.018_kind_phys ! kg mol-1
  ! Gas species - dry air, O2, O3
  character(len=*), parameter, public :: DRY_AIR_LABEL = 'air'
  character(len=*), parameter, public :: O2_LABEL = 'O2'
  character(len=*), parameter, public :: O3_LABEL = 'O3'
  character(len=*), parameter, public :: TUVX_GAS_SPECIES_UNITS = 'molecule cm-3'
  real(kind_phys),  parameter, public :: SCALE_HEIGHT_DRY_AIR = 8.01_kind_phys ! km
  real(kind_phys),  parameter, public :: SCALE_HEIGHT_O2 = 7.0_kind_phys       ! km
  real(kind_phys),  parameter, public :: SCALE_HEIGHT_O3 = 7.0_kind_phys       ! km
  !> Molar mass value of dry air is obtained from 'CAM-SIMA/src/utils/std_atm_profile.F90'
  real(kind_phys),  parameter, public :: MOLAR_MASS_DRY_AIR = 0.0289644_kind_phys ! kg mol-1
  real(kind_phys),  parameter, public :: MOLAR_MASS_O2 = 0.0319988_kind_phys      ! kg mol-1
  real(kind_phys),  parameter, public :: MOLAR_MASS_O3 = 0.0479982_kind_phys      ! kg mol-1

contains

  !> Configures the TUV-x species and their constituent properties.
  ! If the MICM configuration includes any TUV-x gas species, constituent properties
  ! are not created; otherwise, new constituent properties are generated for each species.
  subroutine configure_tuvx_species(micm_species, tuvx_species, constituent_props, &
                                    errmsg, errcode)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use musica_ccpp_species,       only: musica_species_t
    use musica_util,               only: error_t

    type(musica_species_t),                           intent(inout) :: micm_species(:)
    type(musica_species_t),              allocatable, intent(out)   :: tuvx_species(:)
    type(ccpp_constituent_properties_t), allocatable, intent(out)   :: constituent_props(:)
    character(len=512),                               intent(out)   :: errmsg
    integer,                                          intent(out)   :: errcode

    ! local variables
    integer                             :: num_new_species = 4
    integer                             :: num_micm_species
    ! temp_constituents_props is used to store TUVx-specific constituents and gas species
    ! that are not registered by MICM. Its fixed array size represents the maximum number
    ! of possible constituents.
    type(ccpp_constituent_properties_t) :: temp_constituent_props(4)
    logical                             :: is_dry_air_registered = .false.
    logical                             :: is_O2_registered = .false.
    logical                             :: is_O3_registered = .false.
    integer                             :: i_new, i_species, i_tuvx_species

    num_micm_species = size(micm_species)
    is_dry_air_registered = .false.
    is_O2_registered = .false.
    is_O3_registered = .false.

    ! Register cloud liquid water content needed for cloud optics calculations
    i_new = 1
    call temp_constituent_props(i_new)%instantiate( &
      std_name = CLOUD_LIQUID_WATER_CONTENT_LABEL, &
      long_name = CLOUD_LIQUID_WATER_CONTENT_LONG_NAME, &
      units = CLOUD_LIQUID_WATER_CONTENT_UNITS, &
      vertical_dim = "vertical_layer_dimension", &
      default_value = 0.0_kind_phys, &
      min_value = 0.0_kind_phys, &
      molar_mass = CLOUD_LIQUID_WATER_CONTENT_MOLAR_MASS, &
      advected = .true., &
      errcode = errcode, &
      errmsg = errmsg )
    if (errcode /= 0) return

    ! Iterate through the MICM species to check if any TUV-x gas
    ! species are included; if present, updates the scale height and profiled status.
    do i_species = 1, num_micm_species
      if (is_dry_air_registered .and. is_O2_registered .and. is_O3_registered) exit

      if ( micm_species(i_species)%name == DRY_AIR_LABEL ) then
        is_dry_air_registered = .true.
        micm_species(i_species)%profiled = .true.
        micm_species(i_species)%scale_height = SCALE_HEIGHT_DRY_AIR
      else if ( micm_species(i_species)%name == O2_LABEL ) then
        is_O2_registered = .true.
        micm_species(i_species)%profiled = .true.
        micm_species(i_species)%scale_height = SCALE_HEIGHT_O2
      else if ( micm_species(i_species)%name == O3_LABEL ) then
        is_O3_registered = .true.
        micm_species(i_species)%profiled = .true.
        micm_species(i_species)%scale_height = SCALE_HEIGHT_O3
      end if
    end do

    if (.not. is_dry_air_registered) then
      i_new = i_new + 1
      call temp_constituent_props(i_new)%instantiate( &
        std_name = DRY_AIR_LABEL, &
        long_name = DRY_AIR_LABEL, &
        units = 'kg kg-1', &
        vertical_dim = "vertical_layer_dimension", &
        default_value = 0.0_kind_phys, &
        min_value = 0.0_kind_phys, &
        molar_mass = MOLAR_MASS_DRY_AIR, &
        advected = .false., &
        errcode = errcode, &
        errmsg = errmsg )
      if (errcode /= 0) return
    end if

    if (.not. is_O2_registered) then
      i_new = i_new + 1
      call temp_constituent_props(i_new)%instantiate( &
        std_name = O2_LABEL, &
        long_name = O2_LABEL, &
        units = 'kg kg-1', &
        vertical_dim = "vertical_layer_dimension", &
        default_value = 0.0_kind_phys, &
        min_value = 0.0_kind_phys, &
        molar_mass = MOLAR_MASS_O2, &
        advected = .false., &
        errcode = errcode, &
        errmsg = errmsg )
      if (errcode /= 0) return
    end if

    if (.not. is_O3_registered) then
      i_new = i_new + 1
      call temp_constituent_props(i_new)%instantiate( &
        std_name = O3_LABEL, &
        long_name = O3_LABEL, &
        units = 'kg kg-1', &
        vertical_dim = "vertical_layer_dimension", &
        default_value = 0.0_kind_phys, &
        min_value = 0.0_kind_phys, &
        molar_mass = MOLAR_MASS_O3, &
        advected = .false., &
        errcode = errcode, &
        errmsg = errmsg )
      if (errcode /= 0) return
    end if

    allocate( constituent_props(i_new) )
    constituent_props(:) = temp_constituent_props(1:i_new)

    allocate( tuvx_species(num_new_species) )
    i_tuvx_species = 1
    index_cloud_liquid_water_content = i_tuvx_species
    tuvx_species(i_tuvx_species) = musica_species_t( &
      name = CLOUD_LIQUID_WATER_CONTENT_LABEL, &
      unit = CLOUD_LIQUID_WATER_CONTENT_UNITS, &
      molar_mass = CLOUD_LIQUID_WATER_CONTENT_MOLAR_MASS, &
      index_musica_species = i_tuvx_species )

    i_tuvx_species = i_tuvx_species + 1
    index_dry_air = i_tuvx_species
    tuvx_species(i_tuvx_species) = musica_species_t( &
      name = DRY_AIR_LABEL, &
      unit = TUVX_GAS_SPECIES_UNITS, & ! TUV-x profile unit, different from molar mass unit
      molar_mass = MOLAR_MASS_DRY_AIR, & ! kg mol-1
      index_musica_species = i_tuvx_species, &
      profiled = .true., &
      scale_height = SCALE_HEIGHT_DRY_AIR )

    i_tuvx_species = i_tuvx_species + 1
    index_O2 = i_tuvx_species
    tuvx_species(i_tuvx_species) = musica_species_t( &
      name = O2_LABEL, &
      unit = TUVX_GAS_SPECIES_UNITS, & ! TUV-x profile unit, different from molar mass unit
      molar_mass = MOLAR_MASS_O2, & ! kg mol-1
      index_musica_species = i_tuvx_species, &
      profiled = .true., &
      scale_height = SCALE_HEIGHT_O2 )

    i_tuvx_species = i_tuvx_species + 1
    index_O3 = i_tuvx_species
    tuvx_species(i_tuvx_species) = musica_species_t( &
      name = O3_LABEL, &
      unit = TUVX_GAS_SPECIES_UNITS, & ! TUV-x profile unit, different from molar mass unit
      molar_mass = MOLAR_MASS_O3, & ! kg mol-1
      index_musica_species = i_tuvx_species, &
      profiled = .true., &
      scale_height = SCALE_HEIGHT_O3 )

  end subroutine configure_tuvx_species

  !> Ensures that the indices of all TUV-x species are initialized.
  ! This function is typically called during the initialization phase,
  ! so that the indices can be used during the run phase without the need
  ! for additional checks.
  subroutine check_tuvx_species_initialization(errmsg, errcode)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    errmsg = ''
    errcode = 0

    if ((index_cloud_liquid_water_content == MUSICA_INT_UNASSIGNED) .or. &
        (index_dry_air == MUSICA_INT_UNASSIGNED) .or. &
        (index_O2 == MUSICA_INT_UNASSIGNED) .or. &
        (index_O3== MUSICA_INT_UNASSIGNED)) then
      errmsg = "[MUSICA Error] TUV-x species index (or indices) has not been initialized."
      errcode = 1
      return
    end if

  end subroutine check_tuvx_species_initialization

end module musica_ccpp_tuvx_load_species