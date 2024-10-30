module read_tuvx_data

  implicit none

  private
  public :: read_extraterrestrial_flux_from_file

  character(len=50), dimension(4), allocatable :: filepath_of_extraterrestrial_flux

contains

  subroutine read_extraterrestrial_flux_from_file(num_wavelength_grid_sections, & 
      wavelength_grid_interfaces, extraterrestrial_flux)
    use ccpp_kinds,          only: kind_phys

    real(kind_phys), intent(out) :: num_wavelength_grid_sections  ! (count)
    real(kind_phys), intent(out) :: wavelength_grid_interfaces(:) ! nm
    real(kind_phys), intent(out) :: extraterrestrial_flux(:)      ! photons cm-2 s-1 nm-1

    allocate(character(len=50), dimension(4))
    filepath_of_extraterrestrial_flux(1) = 'data/profiles/solar/susim_hi.flx'
    filepath_of_extraterrestrial_flux(2) = 'data/profiles/solar/atlas3_1994_317_a.dat'
    filepath_of_extraterrestrial_flux(3) = 'data/profiles/solar/sao2010.solref.converted'
    filepath_of_extraterrestrial_flux(4) = 'data/profiles/solar/neckel.flx'
      
    !TODO(jiwon) temporary
    num_wavelength_grid_sections = 4
    wavelength_grid_interfaces(:) = [2, 4]
    extraterrestrial_flux(:) = [2, 4]

  end subroutine read_extraterrestrial_flux_from_file

end module read_tuvx_data