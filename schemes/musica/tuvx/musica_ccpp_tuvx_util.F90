module musica_ccpp_tuvx_util
  implicit none

  private
  public :: tuvx_deallocate

contains

  !> This is a helper subroutine created to deallocate objects associated with TUV-x
  subroutine tuvx_deallocate(grids, profiles, radiators, tuvx, height_grid, &
                             wavelength_grid, temperature_profile,          &
                             surface_albedo_profile, extraterrestrial_flux_profile)
    use musica_tuvx, only: tuvx_t, grid_map_t, profile_map_t, radiator_map_t, &
                           grid_t, profile_t

    type(grid_map_t),     pointer :: grids
    type(profile_map_t),  pointer :: profiles
    type(radiator_map_t), pointer :: radiators
    type(tuvx_t),         pointer :: tuvx
    type(grid_t),         pointer :: height_grid
    type(grid_t),         pointer :: wavelength_grid
    type(profile_t),      pointer :: temperature_profile
    type(profile_t),      pointer :: surface_albedo_profile
    type(profile_t),      pointer :: extraterrestrial_flux_profile

    if (associated( grids )) deallocate( grids )
    if (associated( profiles )) deallocate( profiles )
    if (associated( radiators )) deallocate( radiators )

    if (associated( tuvx )) then
      deallocate( tuvx )
      tuvx => null()
    end if

    if (associated( height_grid )) then
      deallocate( height_grid )
      height_grid => null()
    end if

    if (associated( wavelength_grid )) then
      deallocate( wavelength_grid )
      wavelength_grid => null()
    end if

    if (associated( temperature_profile )) then
      deallocate( temperature_profile )
      temperature_profile => null()
    end if

    if (associated( surface_albedo_profile )) then
      deallocate( surface_albedo_profile )
      surface_albedo_profile => null()
    end if

    if (associated( extraterrestrial_flux_profile )) then
      deallocate( extraterrestrial_flux_profile )
      extraterrestrial_flux_profile => null()
    end if

  end subroutine tuvx_deallocate

end module musica_ccpp_tuvx_util