! Copyright (C) 2024-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

!> Photolysis rate constant calculation using TUV-x for MUSICA chemistry
module musica_ccpp_photolysis

  implicit none
  private

  public :: musica_ccpp_photolysis_run

contains

  !> \section arg_table_musica_ccpp_photolysis_run Argument Table
  !! \htmlinclude musica_ccpp_photolysis_run.html
  !!
  !! The standard name for the variable 'surface_temperature' is
  !! 'blackbody_temperature_at_surface' as this is the standard name
  !! used for 'cam_in%ts,' which represents the same quantity.
  subroutine musica_ccpp_photolysis_run(temperature, dry_air_density, constituents,          &
                                        geopotential_height_wrt_surface_at_midpoint,         &
                                        geopotential_height_wrt_surface_at_interface,        &
                                        surface_geopotential, surface_temperature,           &
                                        surface_albedo, extraterrestrial_flux,               &
                                        standard_gravitational_acceleration,                 &
                                        cloud_area_fraction, air_pressure_thickness,         &
                                        solar_zenith_angle, earth_sun_distance,              &
                                        rate_parameters, errmsg, errcode)
    use ccpp_kinds,          only: kind_phys
    use musica_ccpp_tuvx,    only: tuvx_run
    use musica_ccpp_util,    only: do_tuvx
    use musica_ccpp_species, only: number_of_tuvx_species, tuvx_indices_constituent_props, &
                                   extract_subset_constituents

    real(kind_phys), target, intent(in)  :: temperature(:,:)                                  ! K (column, layer)
    real(kind_phys), target, intent(in)  :: dry_air_density(:,:)                              ! kg m-3 (column, layer)
    real(kind_phys), target, intent(in)  :: constituents(:,:,:)                               ! kg kg-1 (column, layer, constituent)
    real(kind_phys),         intent(in)  :: geopotential_height_wrt_surface_at_midpoint(:,:)  ! m (column, layer)
    real(kind_phys),         intent(in)  :: geopotential_height_wrt_surface_at_interface(:,:) ! m (column, interface)
    real(kind_phys),         intent(in)  :: surface_geopotential(:)                           ! m2 s-2 (column)
    real(kind_phys),         intent(in)  :: surface_temperature(:)                            ! K (column)
    real(kind_phys),         intent(in)  :: surface_albedo(:)                                 ! fraction (column)
    real(kind_phys),         intent(in)  :: extraterrestrial_flux(:)                          ! photons cm-2 s-1 nm-1 (wavelength interface)
    real(kind_phys),         intent(in)  :: standard_gravitational_acceleration               ! m s-2
    real(kind_phys),         intent(in)  :: cloud_area_fraction(:,:)                          ! fraction (column, layer)
    real(kind_phys),         intent(in)  :: air_pressure_thickness(:,:)                       ! Pa (column, layer)
    real(kind_phys),         intent(in)  :: solar_zenith_angle(:)                             ! radians (column)
    real(kind_phys),         intent(in)  :: earth_sun_distance                                ! AU
    real(kind_phys),         intent(out) :: rate_parameters(:,:,:)                            ! various units (column, layer, parameter)
    character(len=512),      intent(out) :: errmsg
    integer,                 intent(out) :: errcode

    ! local variables
    real(kind_phys), dimension(size(constituents, dim=1), &
                               size(constituents, dim=2), &
                               number_of_tuvx_species) :: constituents_tuvx_species ! kg kg-1

    errmsg = ''
    errcode = 0

    ! Zero all rate parameters at initialization.
    ! TUV-x only writes the mapped photolysis slots, so unmapped slots
    ! (USER./EMIS./LOSS./SURF. reactions)
    ! would otherwise be uninitialized memory.
    !
    ! Schemes ordered between this scheme and
    ! musica_ccpp_chemistry may fill non-photolysis slots.
    rate_parameters(:,:,:) = 0.0_kind_phys

    if (.not. do_tuvx) return

    call extract_subset_constituents(tuvx_indices_constituent_props, constituents, &
                                     constituents_tuvx_species, errmsg, errcode)
    if (errcode /= 0) return

    call tuvx_run(temperature, dry_air_density,                 &
                  constituents_tuvx_species,                    &
                  geopotential_height_wrt_surface_at_midpoint,  &
                  geopotential_height_wrt_surface_at_interface, &
                  surface_geopotential, surface_temperature,    &
                  surface_albedo,                               &
                  extraterrestrial_flux,                        &
                  standard_gravitational_acceleration,          &
                  cloud_area_fraction,                          &
                  air_pressure_thickness,                       &
                  solar_zenith_angle,                           &
                  earth_sun_distance,                           &
                  rate_parameters,                              &
                  errmsg, errcode)
    if (errcode /= 0) return

  end subroutine musica_ccpp_photolysis_run

end module musica_ccpp_photolysis
