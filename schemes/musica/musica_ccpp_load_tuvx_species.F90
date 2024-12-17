module muscia_ccpp_load_tuvx_species
  use musica_ccpp_species

  implicit none
  private

  public :: configure_tuvx_species
contains
call configure_tuvx_species(constituent_props, musica_species, tuvx_specific_species, &
errmsg, errcode)
  ! Add constituent props and then create musica_species
  ! This is another reason in favor of moving all mechanism parsing of open atmos; 
  ! you could choose to do this and it would be valid. For micm, we always want 
  ! a full mechanism becauase that's what we need
  subroutine configure_tuvx_species(constituent_props, musica_species, tuvx_specific_species, &
                                    errmsg, errcode)
    use ccpp_kinds,                only: kind_phys
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use musica_ccpp_species,       only: musica_species_t
    use musica_util,               only: error_t

    type(ccpp_constituent_properties_t), allocatable, intent(out)   :: constituent_props(:)
    type(musica_species_t),              allocatable, intent(inout) :: musica_species(:)
    type(musica_species_t),              allocatable, intent(inout) :: tuvx_specific_species(:)
    character(len=512),                               intent(out)   :: errmsg
    integer,                                          intent(out)   :: errcode

    ! local variables
    integer                             :: num_new_species = 4
    integer                             :: num_registered_species = size(musica_species)
    type(ccpp_constituent_properties_t) :: temp_constituent_props(num_new_species)
    type(musica_species_t)              :: temp_musica_species(num_new_species)
    type(musica_species_t)              :: copy_musica_species(num_registered_species)
    logical                             :: is_O2_registered = .false.
    logical                             :: is_O3_registered = .false.
    logical                             :: is_dry_air_registered = .false.
    integer                             :: i_new, i_registered

    character(len=*), parameter :: CLOUD_LIQUID_WATER_CONTENT_LABEL = &
        'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water'
    character(len=*), parameter :: CLOUD_LIQUID_WATER_CONTENT_LONG_NAME = &
        'cloud water mass mixing ratio with respect to moist air plus all airborne condensates'
    character(len=*), parameter :: CLOUD_LIQUID_WATER_CONTENT_UNITS = 'kg kg-1'
    real(kind_phys),  parameter :: CLOUD_LIQUID_WATER_CONTENT_MOLAR_MASS = 0.018_kind_phys ! kg mol-1
    real(kind_phys),  parameter :: SCALE_HEIGHT_DRY_AIR = 8.01_kind_phys ! km
    real(kind_phys),  parameter :: SCALE_HEIGHT_O2 = 7.0_kind_phys       ! km
    real(kind_phys),  parameter :: SCALE_HEIGHT_O3 = 7.0_kind_phys       ! km
    !> Molar mass value of dry air is obtained from 'CAM-SIMA/src/utils/std_atm_profile.F90'
    ! TODO(jiwon) - how to make this an input argument?
    real(kind_phys),  parameter :: MOLAR_MASS_DRY_AIR = 0.0289644_kind_phys ! kg mol-1
    real(kind_phys),  parameter :: MOLAR_MASS_O2 = 0.0319988_kind_phys      ! kg mol-1
    real(kind_phys),  parameter :: MOLAR_MASS_O3 = 0.0479982_kind_phys      ! kg mol-1

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

    temp_musica_species(i_new) = musica_species_t( &
      name = CLOUD_LIQUID_WATER_CONTENT_LABEL, &
      unit = CLOUD_LIQUID_WATER_CONTENT_UNITS, &
      molar_mass = CLOUD_LIQUID_WATER_CONTENT_MOLAR_MASS, &
      index_musica_species = num_registered_species + i_new )

    ! Add gas species - dry air, O2, O3 - to be profiled
    ! iterate through all the registered species to 
    ! check if the species is already registered and if so
    ! update scale_height
    do i_registered = 1, num_registered_species
      if (is_dry_air_registered .and. is_O2_registered .and. is_O3_registered) exit

      if (musica_species(i)%name == "dry_air") then
        is_dry_air_registered = .true.
        musica_species(i_registered)%profiled = .true.
        musica_species(i_registered)%scale_height = SCALE_HEIGHT_DRY_AIR
      else if ( musica_species(i_registered)%name == "O2" ) then
        is_O2_registered = .true.
        musica_species(i_registered)%profiled = .true.
        musica_species(i_registered)%scale_height = SCALE_HEIGHT_O2
      else if (musica_species(i)%name == "O3") then 
        is_O3_registered = .true.
        musica_species(i_registered)%profiled = .true.
        musica_species(i_registered)%scale_height = SCALE_HEIGHT_O3
      end if
    end do

    if (.not. is_dry_air_registered) then
      i_new = i_new + 1

      call constituent_props(i_new)%instantiate( &
        std_name = 'dry_air', &
        long_name = 'dry_air', &
        units = 'kg kg-1', &
        vertical_dim = "vertical_layer_dimension", &
        default_value = 0.0_kind_phys, &
        min_value = 0.0_kind_phys, &
        molar_mass = MOLAR_MASS_DRY_AIR, &
        advected = .false., &
        errcode = errcode, &
        errmsg = errmsg )
      if (errcode /= 0) return

      temp_musica_species(i_new) = musica_species_t( &
        name = 'dry_air', &
        unit = "molecule cm-3", & ! TUVX profile unit, which can be different from molar mass unit
        molar_mass = MOLAR_MASS_DRY_AIR, & ! kg mol-1
        index_musica_species = num_registered_species + i_new, &
        profiled = .true., &
        scale_height = SCALE_HEIGHT_DRY_AIR )
    end if

    if (.not. is_O2_registered) then
      i_new = i_new + 1

      call constituent_props(i_new)%instantiate( &
        std_name = 'O2', &
        long_name = 'O2', &
        units = 'kg kg-1', &
        vertical_dim = "vertical_layer_dimension", &
        default_value = 0.0_kind_phys, &
        min_value = 0.0_kind_phys, &
        molar_mass = MOLAR_MASS_O2, &
        advected = .false., &
        errcode = errcode, &
        errmsg = errmsg )
      if (errcode /= 0) return

      temp_musica_species(i_new) = musica_species_t( &
        name = 'O2', &
        unit = "molecule cm-3", & ! TUVX profile unit, which can be different from molar mass unit
        molar_mass = MOLAR_MASS_DRY_O2, & ! kg mol-1
        index_musica_species = num_registered_species + i_new, &
        profiled = .true., &
        scale_height = SCALE_HEIGHT_O2 )
    end if

    if (.not. is_O3_registered) then
      i_new = i_new + 1

      call constituent_props(i_new)%instantiate( &
        std_name = 'O3', &
        long_name = 'O3', &
        units = 'kg kg-1', &
        vertical_dim = "vertical_layer_dimension", &
        default_value = 0.0_kind_phys, &
        min_value = 0.0_kind_phys, &
        molar_mass = MOLAR_MASS_O3, &
        advected = .false., &
        errcode = errcode, &
        errmsg = errmsg )
      if (errcode /= 0) return

      temp_musica_species(i_new) = musica_species_t( &
        name = 'O3', &
        unit = "molecule cm-3", & ! TUVX profile unit, which can be different from molar mass unit
        molar_mass = MOLAR_MASS_DRY_O3, & ! kg mol-1
        index_musica_species = num_registered_species + i_new, &
        profiled = .true., &
        scale_height = SCALE_HEIGHT_O3 )
    end if

    allocate( constituent_props( size(i_new) ) )
    constituent_props(:) = temp_musica_species(1:i_new)

    if (i_new > 0 ) then
      copy_musica_species = musica_species
      deallocate( musica_species )
      allocate( musica_species( num_registered_species + i_new )) 
      musica_species = [ copy_musica_species, temp_musica_species ]
    end if

  end subroutine configure_tuvx_species

end module muscia_ccpp_load_tuvx_species
