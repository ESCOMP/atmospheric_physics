module musica_ccpp_tuvx_no_photolysis_rate

  !> @file musica_ccpp_tuvx_no_photolysis_rate.F90
  !> @brief Module for calculating the NO photolysis rate
  ! This module calculates the NO photolysis rate using the same method in CAM
  ! https://github.com/ESCOMP/CAM/blob/ab476f9b7345cbefdc4cf67ff17f0fe85d8c7387/src/chemistry/mozart/mo_jshort.F90#L1796-L1870
  ! Much of the code is also based off of a PR which incorporated TUVx into CAM
  ! https://github.com/ESCOMP/CAM/blob/c12d1e46e0fdc1dccb0a651a6c9fefd6bb80b2ba/src/chemistry/mozart/mo_tuvx.F90#L1849-L1943
  ! The actual method is described in this paper:
  ! Minschwaner, K., Siskind, D.E., 1993. A new calculation of nitric oxide photolysis in the stratosphere, mesosphere, and lower thermosphere.
  ! Journal of Geophysical Research: Atmospheres 98, 20401–20412. https://doi.org/10.1029/93JD02007
  ! Acronyms:
  ! SRB: Schumann–Runge bands, a group of electronic transitions in molecular oxygen that absorb solar radiation
  ! NO: Nitric oxide

  implicit none

  private
  public :: calculate_NO_photolysis_rate

contains

  !> Calculates the NO photolysis rate
  function calculate_NO_photolysis_rate(solar_zenith_angle, extraterrestrial_flux, constituents, height_at_interfaces, &
    dry_air_density ) &
      result(jNO)
    use ccpp_kinds,          only: kind_phys
    ! inputs
    real(kind_phys),    intent(in)    :: solar_zenith_angle       ! degrees
    real(kind_phys),         intent(in)    :: extraterrestrial_flux(:) ! photons cm-2 s-1 nm-1
    real(kind_phys), target, intent(inout) :: constituents(:,:,:)      ! various (column, layer, constituent)
    real(kind_phys),    intent(in)    :: height_at_interfaces(:)  ! m
    real(kind_phys),    intent(in)         :: dry_air_density(:,:)     ! kg m-3 (column, layer)

    ! local variables
    real(kind_phys) :: n2_dens(size(constituents, dim=2)+1), o2_dens(size(constituents, dim=2)+1)
    real(kind_phys) :: o3_dens(size(constituents, dim=2)+1), no_dens(size(constituents, dim=2)+1)
    ! species slant column densities (molecule cm-2)
    real(kind_phys) :: o2_slant(size(constituents, dim=2)+1), o3_slant(size(constituents, dim=2)+1)
    real(kind_phys) :: no_slant(size(constituents, dim=2)+1)
    ! working photo rate array
    real(kind_phys) :: work_jno(size(constituents, dim=2)+1)
    ! parameters needed to calculate slant column densities
    ! (see sphers routine description for details)
    integer       :: nid(size(constituents, dim=2)+1)
    real(kind_phys) :: dsdh(0:size(constituents, dim=2)+1,size(constituents, dim=2)+1)
    ! layer thickness (cm)
    real(kind_phys) :: delz(size(constituents, dim=2)+1)
    ! conversion from km to cm
    real(kind_phys), parameter :: km2cm = 1.0e5_r8
    ! final photolysis rate
    real(kind_phys) :: jNO

    jNO = 0.0

  end function calculate_NO_photolysis_rate


end module musica_ccpp_tuvx_no_photolysis_rate