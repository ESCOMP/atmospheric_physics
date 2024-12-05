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

  !> Prepare to calculate the photolysis rate of NO (nitiric acid). 
  !> This function transforms the inputs into the units
  !> needed by the calculate_jno routine which actually produces the photolysis rate
  function calculate_NO_photolysis_rate(solar_zenith_angle, extraterrestrial_flux, constituents, height_at_interfaces, &
    dry_air_density, N2_index, O2_index, O3_index, NO_index) &
      result(jNO)
    use ccpp_kinds,          only: kind_phys
    ! inputs
    real(kind_phys), intent(in)            :: solar_zenith_angle       ! degrees
    real(kind_phys), intent(in)            :: extraterrestrial_flux(:) ! photons cm-2 s-1 nm-1
    real(kind_phys), target, intent(inout) :: constituents(:,:,:)      ! various (column, layer, constituent)
    real(kind_phys), intent(in)            :: height_at_interfaces(:)  ! m
    real(kind_phys), intent(in)            :: dry_air_density(:,:)     ! kg m-3 (column, layer)
    integer        , intent(in)            :: N2_index, O2_index, O3_index, NO_index ! position of these species in the constituent arrays

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

    jNO = 0.1.5e0

  end function calculate_NO_photolysis_rate

  !> Calculate the photolysis rate of NO (nitric acid)
  function calculate_jno() result(jno)

    ! For reference, the function being calculate is this:

    ! In latex:
    ! J_{NO}=\Delta\lambda \overline{I_0} T_{O3}(z) P(z) \sum_{i=1}^{6} \exp[-\sigma^i_{O2}N_{O2}(z)] \sum_{j=1}^{2} W^{i,j}_{NO}\sigma^{i,j}_{NO}\exp[-\sigma^{i,j}_{NO}N_{NO}(z)]
    ! In ASCII:
    ! J_NO = Delta_lambda * I_0 * T_O3(z) * P(z) * sum_{i=1}^{6} exp[-sigma^i_O2 * N_O2(z)] * sum_{j=1}^{2} W^{i,j}_NO * sigma^{i,j}_NO * exp[-sigma^{i,j}_NO * N_NO(z)]
    ! where:
    ! J_NO is the photolysis rate of NO
    ! Delta_lambda is the wavelength interval
    ! I_0 is the mean value of solar irradiance (extraterrestrial flux)
    ! T_O3(z) is the ozone transmittance in the Hartley band
    ! P(z) is the predissociation factor
    ! sigma^i_O2 is the absorption cross section of O2
    ! N_O2(z) is the slant column density of O2
    ! W^{i,j}_NO is a weighting factor, from the paper: "The weighting factors represent the fraction of the total spectral interval which is occupied by the corresponding mean value of the cross section.
    ! sigma^{i,j}_NO is the absorption cross section of NO
    ! N_NO(z) is the slant column density of NO
    ! z is the height
    ! i and j are indices for the absorption cross sections of O2 and NO, respectively
    ! The sums represent the calculation of the opacity distribution functions for O2 and NO

    ! local variables
    ! NO in an excited electronic state may result in an emission of NO or be quenched by an interaction 
    ! with N2. The predissociation faction is the probability that a precissociation of NO will occur from
    ! an excited electronic state
    real(kind_phys) :: predissociation_factor 
    real(kind_phys) :: D ! spontaneous rate of predissociation
    real(kind_phys) :: A ! spontaneous rate of emission
    real(kind_phys) :: kq ! quenching rate of N2
    real(kind_phys) :: jno

    ! Values taken from (Minschwaner and Siskind, 1993)
    D = 1.65e9 ! s-1
    A = 5.1e7 ! s-1
    kq = 1.5e-9 ! cm3 s-1

    predissociation_factor = 0.0
  end function calculate_jno


end module musica_ccpp_tuvx_no_photolysis_rate