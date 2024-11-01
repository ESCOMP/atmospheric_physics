module read_tuvx_data

  implicit none

  private
  public :: read_extraterrestrial_flux

  character(len=50), dimension(4), public :: filepath_of_extraterrestrial_flux

contains

  !> Reads a data file to retrieve the extraterrestrial radiation flux values. 
  ! This function is a temporary implementation and will be updated or
  ! replaced in future versions of the code.
  subroutine read_extraterrestrial_flux(num_wavelength_grid_sections, & 
      wavelength_grid_interfaces, extraterrestrial_flux)
    use ccpp_kinds,          only: kind_phys

    integer,                      intent(out) :: num_wavelength_grid_sections  ! (count)
    real(kind_phys), allocatable, intent(out) :: wavelength_grid_interfaces(:) ! nm
    real(kind_phys), allocatable, intent(out) :: extraterrestrial_flux(:)      ! photons cm-2 s-1 nm-1

    filepath_of_extraterrestrial_flux(1) = 'data/profiles/solar/susim_hi.flx'
    filepath_of_extraterrestrial_flux(2) = 'data/profiles/solar/atlas3_1994_317_a.dat'
    filepath_of_extraterrestrial_flux(3) = 'data/profiles/solar/sao2010.solref.converted'
    filepath_of_extraterrestrial_flux(4) = 'data/profiles/solar/neckel.flx'
      
    num_wavelength_grid_sections = 8

    allocate(wavelength_grid_interfaces(num_wavelength_grid_sections + 1))
    allocate(extraterrestrial_flux(num_wavelength_grid_sections))

    wavelength_grid_interfaces(:) = &
      [200.0_kind_phys, 210.0_kind_phys, 220.0_kind_phys, 230.0_kind_phys, &
       240.0_kind_phys, 250.0_kind_phys, 260.0_kind_phys, 270.0_kind_phys, 280.0_kind_phys]

    extraterrestrial_flux(:) = &
      [1.5e13_kind_phys, 1.5e13_kind_phys, 1.4e13_kind_phys, 1.4e13_kind_phys, &
       1.3e13_kind_phys, 1.2e13_kind_phys, 1.1e13_kind_phys, 1.0e13_kind_phys]

  end subroutine read_extraterrestrial_flux

end module read_tuvx_data