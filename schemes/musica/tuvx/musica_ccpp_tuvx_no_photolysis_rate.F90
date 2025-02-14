module musica_ccpp_tuvx_no_photolysis_rate

  use ccpp_kinds,          only: dk => kind_phys
  use musica_ccpp_species, only: MUSICA_INT_UNASSIGNED

  implicit none

  !> @file musica_ccpp_tuvx_no_photolysis_rate.F90
  !> @brief Module for calculating the NO photolysis rate
  ! This module calculates the NO photolysis rate using the same method in CAM
  ! https://github.com/ESCOMP/CAM/blob/ab476f9b7345cbefdc4cf67ff17f0fe85d8c7387/src/chemistry/mozart/mo_jshort.F90#L1796-L1870
  ! Much of the code is also based off of a PR which incorporated TUVx into CAM
  ! https://github.com/ESCOMP/CAM/blob/c12d1e46e0fdc1dccb0a651a6c9fefd6bb80b2ba/src/chemistry/mozart/mo_tuvx.F90#L1849-L1943
  !
  ! The actual method is described in this paper:
  ! Minschwaner, K., Siskind, D.E., 1993. A new calculation of nitric oxide photolysis in the stratosphere, mesosphere, and lower thermosphere.
  ! Journal of Geophysical Research: Atmospheres 98, 20401–20412. https://doi.org/10.1029/93JD02007
  !
  ! Acronyms:
  ! SRB: Schumann–Runge bands, a group of electronic transitions in molecular oxygen that absorb solar radiation
  ! NO: Nitric oxide
  ! MS93: A reference to the Minschwaner and Siskind paper

  private
  public :: check_NO_exist, set_NO_index_constituent_props, check_NO_initialization, &
            calculate_NO_photolysis_rate, set_NO_photolysis_index

  character(len=*), parameter, public :: N2_label = 'N2'
  character(len=*), parameter, public :: NO_label = 'NO'
  integer,          parameter, public :: num_NO_photolysis_constituents = 2
  logical,          protected, public :: is_NO_photolysis_active = .false.
  integer,          protected, public :: index_NO_photolysis_rate = MUSICA_INT_UNASSIGNED
  integer,          protected, public :: index_N2 = MUSICA_INT_UNASSIGNED
  integer,          protected, public :: index_NO = MUSICA_INT_UNASSIGNED
  integer,          protected, public :: constituent_props_index_N2 = MUSICA_INT_UNASSIGNED
  integer,          protected, public :: constituent_props_index_NO = MUSICA_INT_UNASSIGNED
  integer,          protected, public :: NO_photolysis_indices_constituent_props(2) = &
                                        [MUSICA_INT_UNASSIGNED, MUSICA_INT_UNASSIGNED]
  real(dk),         protected, public :: MOLAR_MASS_N2 = MUSICA_INT_UNASSIGNED ! kg mol-1
  real(dk),         protected, public :: MOLAR_MASS_NO = MUSICA_INT_UNASSIGNED ! kg mol-1

contains

  !> Checks for the presence of NO and N2, and sets the flag to indicate
  !! whether the NO photolysis calculation is enabled or disabled
  subroutine check_NO_exist(micm_species, includes_no_photolysis)
    use musica_ccpp_species, only: musica_species_t

    type(musica_species_t), intent(in) :: micm_species(:)
    logical,                intent(in) :: includes_no_photolysis

    ! local variables
    logical :: is_N2_registered = .false.
    logical :: is_NO_registered = .false.
    integer :: num_micm_species
    integer :: i_species

    is_N2_registered = .false.
    is_NO_registered = .false.
    num_micm_species = size(micm_species)

    do i_species = 1, num_micm_species
      if ( is_N2_registered .and. is_NO_registered ) then
        is_NO_photolysis_active = .true.
        exit
      end if

      if ( micm_species(i_species)%name == N2_label ) then
        is_N2_registered = .true.
        MOLAR_MASS_N2 = micm_species(i_species)%molar_mass
      else if ( micm_species(i_species)%name == NO_label ) then
        is_NO_registered = .true.
        MOLAR_MASS_NO = micm_species(i_species)%molar_mass
      end if
    end do

    if (is_N2_registered .and. is_NO_registered .and. includes_no_photolysis) then
      is_NO_photolysis_active = .true.
    end if

  end subroutine check_NO_exist

  subroutine set_NO_photolysis_index(no_index)
    integer, intent(in) :: no_index

    index_NO_photolysis_rate = no_index

  end subroutine set_NO_photolysis_index

  !> Gets the indices of NO constituents from the CCPP constituent properties,
  !! and stores them in the indices array along with their relative positions
  !! in the NO constituents array
  subroutine set_NO_index_constituent_props(constituent_props, errmsg, errcode)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ccpp_const_utils,          only: ccpp_const_get_idx

    type(ccpp_constituent_prop_ptr_t), intent(in)  :: constituent_props(:)
    character(len=512),                intent(out) :: errmsg
    integer,                           intent(out) :: errcode

    call ccpp_const_get_idx( constituent_props, N2_label, &
              constituent_props_index_N2, errmsg, errcode )
    if (errcode /= 0) return
    index_N2 = 1

    call ccpp_const_get_idx( constituent_props, NO_label, &
              constituent_props_index_NO, errmsg, errcode )
    if (errcode /= 0) return
    index_NO = index_N2 + 1

    NO_photolysis_indices_constituent_props = &
      [constituent_props_index_N2, constituent_props_index_NO]

  end subroutine set_NO_index_constituent_props

  !> Checks if the variables associated with NO photolysis constituents
  !! have been initialized
  subroutine check_NO_initialization(errmsg, errcode)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    errmsg = ''
    errcode = 0

    if ((index_N2 == MUSICA_INT_UNASSIGNED) .or. &
        (index_NO == MUSICA_INT_UNASSIGNED) .or. &
        (constituent_props_index_N2 == MUSICA_INT_UNASSIGNED) .or. &
        (constituent_props_index_NO == MUSICA_INT_UNASSIGNED) .or. &
        (NO_photolysis_indices_constituent_props(1) == MUSICA_INT_UNASSIGNED) .or. &
        (NO_photolysis_indices_constituent_props(2) == MUSICA_INT_UNASSIGNED) .or. &
        (MOLAR_MASS_N2 == MUSICA_INT_UNASSIGNED) .or. &
        (MOLAR_MASS_NO == MUSICA_INT_UNASSIGNED)) then
      errmsg = "[MUSICA Error] The index (or indices) or molar mass values &
                for NO, N2 have not been initialized."
      errcode = 1
      return
    end if

  end subroutine check_NO_initialization

  !> Converts mixing ratio to molecule cm-3
  subroutine convert_mixing_ratio_to_molecule_cm3(mixing_ratio, &
      dry_air_density, molar_mass, molecule_cm3)
    use musica_ccpp_tuvx_gas_species, only: m_3_to_cm_3

    real(dk), intent(in)  :: mixing_ratio(:)    ! kg kg-1
    real(dk), intent(in)  :: dry_air_density(:) ! kg m-3
    real(dk), intent(in)  :: molar_mass         ! kg mol-1
    real(dk), intent(out) :: molecule_cm3(:)    ! molecule cm-3

    ! TODO: get from CAM-SIMA
    real(dk) :: avogadro = 6.022e23_dk

    molecule_cm3 = mixing_ratio * dry_air_density / molar_mass / m_3_to_cm_3 * avogadro

  end subroutine convert_mixing_ratio_to_molecule_cm3

  !> Prepares to calculate the photolysis rate of NO (nitiric acid).
  !! This function transforms the inputs into the units
  !! needed by the calculate_jno routine which actually produces the photolysis rate
  function calculate_NO_photolysis_rate(solar_zenith_angle, extraterrestrial_flux, &
      height_at_interfaces, dry_air_density, constituents_O2, constituents_O3,     &
      constituents_NO_photolysis) result(jNO)
    use musica_ccpp_tuvx_load_species, only: MOLAR_MASS_O2, MOLAR_MASS_O3
    use musica_ccpp_tuvx_gas_species,  only: km_to_cm

    real(dk), intent(in) :: solar_zenith_angle              ! degrees
    real(dk), intent(in) :: extraterrestrial_flux(:)        ! photons cm-2 s-1 nm-1 (layer)
    real(dk), intent(in) :: height_at_interfaces(:)         ! km (layer)
    real(dk), intent(in) :: dry_air_density(:)              ! kg m-3 (layer)
    real(dk), intent(in) :: constituents_O2(:)              ! kg kg-1 (layer)
    real(dk), intent(in) :: constituents_O3(:)              ! kg kg-1 (layer)
    real(dk), intent(in) :: constituents_NO_photolysis(:,:) ! kg kg-1 (layer, constituent)

    ! local variables
    integer  :: number_of_vertical_layers
    real(dk) :: N2_densities(size(dry_air_density)+1) ! mol cm-3
    real(dk) :: O2_densities(size(dry_air_density)+1) ! mol cm-3
    real(dk) :: O3_densities(size(dry_air_density)+1) ! mol cm-3
    real(dk) :: NO_densities(size(dry_air_density)+1) ! mol cm-3
    real(dk) :: O2_slant_colunm_densities(size(dry_air_density)+1) ! mol cm-2
    real(dk) :: O3_slant_colunm_densities(size(dry_air_density)+1) ! mol cm-2
    real(dk) :: NO_slant_colunm_densities(size(dry_air_density)+1) ! mol cm-2
    real(dk) :: current_jNO(size(dry_air_density)+1)

    ! parameters needed to calculate slant column densities
    ! (see sphers routine description for details)
    integer  :: number_of_crossed_layers(size(dry_air_density)+1)
    real(dk) :: slant_path(0:size(dry_air_density)+1, size(dry_air_density)+1)
    real(dk) :: delta_z(size(dry_air_density)+1) ! layer thickness (cm)
    real(dk) :: jNO(size(dry_air_density)) ! final photolysis rate

    number_of_vertical_layers = size(dry_air_density)
    ! TODO: what are these constants? scale heights?
    ! TODO: the values at index 1 appear to be for values above the model top in CAM, but how does that affect cam sima?
    call convert_mixing_ratio_to_molecule_cm3(constituents_NO_photolysis(:,index_N2), &
          dry_air_density, MOLAR_MASS_N2, N2_densities(2:))
    N2_densities(1) = N2_densities(2) * 0.9_dk

    call convert_mixing_ratio_to_molecule_cm3(constituents_O2, &
          dry_air_density, MOLAR_MASS_O2, O2_densities(2:))
    O2_densities(1) = O2_densities(2) * 7.0_dk / ( height_at_interfaces(1) - height_at_interfaces(2) )

    call convert_mixing_ratio_to_molecule_cm3(constituents_O3, &
          dry_air_density, MOLAR_MASS_O3, O3_densities(2:))
    O3_densities(1) = O3_densities(2) * 7.0_dk / ( height_at_interfaces(1) - height_at_interfaces(2) )

    call convert_mixing_ratio_to_molecule_cm3(constituents_NO_photolysis(:,index_NO), &
          dry_air_density, MOLAR_MASS_NO, NO_densities(2:))
    NO_densities(1) = NO_densities(2) * 0.9_dk

    ! calculate slant column densities
    call calculate_slant_path( number_of_vertical_layers, height_at_interfaces, &
          solar_zenith_angle, slant_path, number_of_crossed_layers )
    delta_z(1:number_of_vertical_layers) = km_to_cm *              &
                (height_at_interfaces(1:number_of_vertical_layers) &
                - height_at_interfaces(2:number_of_vertical_layers+1) )
    call calculate_slant_column_density( number_of_vertical_layers+1, delta_z,&
          slant_path, number_of_crossed_layers, O2_densities, O2_slant_colunm_densities )
    call calculate_slant_column_density( number_of_vertical_layers+1, delta_z, &
          slant_path, number_of_crossed_layers, O3_densities, O3_slant_colunm_densities )
    call calculate_slant_column_density( number_of_vertical_layers+1, delta_z, &
          slant_path, number_of_crossed_layers, NO_densities, NO_slant_colunm_densities )

    jNO = calculate_jno(number_of_vertical_layers, extraterrestrial_flux, N2_densities, &
                        O2_slant_colunm_densities, O3_slant_colunm_densities,           &
                        NO_slant_colunm_densities, current_jNO)

  end function calculate_NO_photolysis_rate

  !> Calculate the photolysis rate of NO (nitric acid)
  function calculate_jno(num_vertical_layers, extraterrestrial_flux, n2_dens, o2_slant, o3_slant, no_slant, work_jno) &
    result(jno)

    ! inputs
    integer, intent(in)     :: num_vertical_layers
    real(dk), intent(in)    :: extraterrestrial_flux(:) ! photons cm-2 s-1 nm-1
    real(dk), intent(in)    :: n2_dens(:)               ! molecule cm-3
    real(dk), intent(in)    :: o2_slant(:)              ! molecule cm-2
    real(dk), intent(in)    :: o3_slant(:)              ! molecule cm-2
    real(dk), intent(in)    :: no_slant(:)              ! molecule cm-2
    real(dk), intent(inout) :: work_jno(:)              ! various

    ! local variables
    ! NO in an excited electronic state may result in an emission of NO or be quenched by an interaction 
    ! with N2. The predissociation faction is the probability that a precissociation of NO will occur from
    ! an excited electronic state
    real(dk) :: o3_transmission_factor(size(o3_slant), 4)  ! Define the correct size
    real(dk) :: jno(num_vertical_layers)
    real(dk) :: jno100, jno90, jno50
    integer  :: wavelength_bins ! wavelength bins for MS, 93
    integer  :: idx 
    integer :: lev
    !------------------------------------------------------------------------------
    !   ... O3 SRB Cross Sections from WMO 1985, interpolated onto MS, 1993 grid
    !------------------------------------------------------------------------------
    ! TODO: What are the units? cm2  molecule-1 ?
    real(dk), save :: o3_cross_section(4) = (/ 7.3307600e-19_dk, 6.9660105E-19_dk, 5.9257699E-19_dk, 4.8372219E-19_dk /)

    !------------------------------------------------------------------------------
    !   ... O2 SRB Cross Sections for the six ODF regions, MS, 1993
    !------------------------------------------------------------------------------
    ! TODO: What are 250, 290, and 2100? 
    ! TODO: What are the units? cm2  molecule-1 ?
    real(dk), save :: cross_section250(6)  = (/ 1.117e-23_dk, 2.447e-23_dk, 7.188e-23_dk, 3.042e-22_dk, 1.748e-21_dk, 1.112e-20_dk /)
    real(dk), save :: cross_section290(6)  = (/ 1.350e-22_dk, 2.991e-22_dk, 7.334e-22_dk, 3.074e-21_dk, 1.689e-20_dk, 1.658e-19_dk /)
    real(dk), save :: cross_section2100(6) = (/ 2.968e-22_dk, 5.831e-22_dk, 2.053e-21_dk, 8.192e-21_dk, 4.802e-20_dk, 2.655e-19_dk /)

    !------------------------------------------------------------------------------
    !   ... O2 SRB Cross Sections for the six ODF regions, MS, 1993
    !------------------------------------------------------------------------------
    ! TODO: What are 250, 290, and 2100? 
    ! TODO: What are the units? cm2  molecule-1 ?
    real(dk), save :: o2_250_cross_section(6)  = (/ 1.117e-23_dk, 2.447e-23_dk, 7.188e-23_dk, 3.042e-22_dk, 1.748e-21_dk, 1.112e-20_dk /)
    real(dk), save :: o2_290_cross_section(6)  = (/ 1.350e-22_dk, 2.991e-22_dk, 7.334e-22_dk, 3.074e-21_dk, 1.689e-20_dk, 1.658e-19_dk /)
    real(dk), save :: o2_2100_cross_section(6) = (/ 2.968e-22_dk, 5.831e-22_dk, 2.053e-21_dk, 8.192e-21_dk, 4.802e-20_dk, 2.655e-19_dk /)

    !------------------------------------------------------------------------------
    !   ... delta wavelength of the MS, 1993 grid
    !------------------------------------------------------------------------------
    ! in nm
    ! TODO: What is the grid? The paper only has 3 bands listed in the table, but this has 4...
    real(dk), save :: delta_wavelength(4) = (/ 1.50_dk, 1.50_dk, 5.6_dk, 2.3_dk /)
    

    ! TODO: Are these correct names
    integer, parameter :: number_of_levels = 6
    integer, parameter :: number_of_intervals = 2

    real(dk), dimension(24) :: no_50_weighting_cross_sections, no_90_weighting_cross_sections, no_1000_weighting_cross_sections

    ! TODO: What are 50, 90, 100?
    real(dk), dimension(number_of_levels, number_of_intervals) :: no_50_weighting_factor
    real(dk), dimension(number_of_levels, number_of_intervals) :: no_90_weighting_factor
    real(dk), dimension(number_of_levels, number_of_intervals) :: no_100_weighting_factor
    real(dk), dimension(number_of_levels, number_of_intervals) :: no_50_cross_section
    real(dk), dimension(number_of_levels, number_of_intervals) :: no_90_cross_section
    real(dk), dimension(number_of_levels, number_of_intervals) :: no_100_cross_section

    !------------------------------------------------------------------------------
    !   ... 6 sub-intervals for O2 5-0 at 265K,
    !    2 sub-sub-intervals for NO 0-0 at 250K
    !------------------------------------------------------------------------------
    no_50_weighting_cross_sections(:) = (/    0._dk,       0._dk,       0._dk,       0._dk, &
                  5.12e-02_dk, 5.68e-03_dk, 1.32e-18_dk, 4.41e-17_dk, &
                  1.36e-01_dk, 1.52e-02_dk, 6.35e-19_dk, 4.45e-17_dk, &
                  1.65e-01_dk, 1.83e-02_dk, 7.09e-19_dk, 4.50e-17_dk, &
                  1.41e-01_dk, 1.57e-02_dk, 2.18e-19_dk, 2.94e-17_dk, &
                  4.50e-02_dk, 5.00e-03_dk, 4.67e-19_dk, 4.35e-17_dk /)
    
    !------------------------------------------------------------------------------
    !   ... sub-intervals for o2 9-0 band,
    !    2 sub-sub-intervals for no 1-0 at 250 k
    !------------------------------------------------------------------------------
    no_90_weighting_cross_sections(:) = (/        0._dk,       0._dk,       0._dk,       0._dk, &
                      0._dk,       0._dk,       0._dk,       0._dk, &
                1.93e-03_dk, 2.14e-04_dk, 3.05e-21_dk, 3.20e-21_dk, &
                9.73e-02_dk, 1.08e-02_dk, 5.76e-19_dk, 5.71e-17_dk, &
                9.75e-02_dk, 1.08e-02_dk, 2.29e-18_dk, 9.09e-17_dk, &
                3.48e-02_dk, 3.86e-03_dk, 2.21e-18_dk, 6.00e-17_dk /)
    
    !------------------------------------------------------------------------------
    ! ... sub-intervals for o2 10-0 band,
    !    2 sub-sub-intervals for no 1-0 at 250 k
    !------------------------------------------------------------------------------
    no_1000_weighting_cross_sections(:) = (/  4.50e-02_dk, 5.00e-03_dk, 1.80e-18_dk, 1.40e-16_dk, &
                1.80e-01_dk, 2.00e-02_dk, 1.50e-18_dk, 1.52e-16_dk, &
                2.25e-01_dk, 2.50e-02_dk, 5.01e-19_dk, 7.00e-17_dk, &
                2.25e-01_dk, 2.50e-02_dk, 7.20e-20_dk, 2.83e-17_dk, &
                1.80e-01_dk, 2.00e-02_dk, 6.72e-20_dk, 2.73e-17_dk, &
                4.50e-02_dk, 5.00e-03_dk, 1.49e-21_dk, 6.57e-18_dk /)
    
    no_50_weighting_factor (1:6,1) = no_50_weighting_cross_sections(1:24:4)
    no_50_weighting_factor (1:6,2) = no_50_weighting_cross_sections(2:24:4)
    no_50_cross_section (1:6,1) = no_50_weighting_cross_sections(3:24:4)
    no_50_cross_section (1:6,2) = no_50_weighting_cross_sections(4:24:4)
    no_90_weighting_factor (1:6,1) = no_90_weighting_cross_sections(1:24:4)
    no_90_weighting_factor (1:6,2) = no_90_weighting_cross_sections(2:24:4)
    no_90_cross_section (1:6,1) = no_90_weighting_cross_sections(3:24:4)
    no_90_cross_section (1:6,2) = no_90_weighting_cross_sections(4:24:4)
    no_100_weighting_factor(1:6,1) = no_1000_weighting_cross_sections(1:24:4)
    no_100_weighting_factor(1:6,2) = no_1000_weighting_cross_sections(2:24:4)
    no_100_cross_section(1:6,1) = no_1000_weighting_cross_sections(3:24:4)
    no_100_cross_section(1:6,2) = no_1000_weighting_cross_sections(4:24:4)


    jno = 0.0_dk
    jno100 = 0.0_dk
    jno90 = 0.0_dk
    jno50 = 0.0_dk
    wavelength_bins = 4

    ! ... derive O3 transmission for the three O2 SRB
    ! ... idx = 1,2, and 4 are used below for jno
    ! ... TODO: why is 3 not used?
    !------------------------------------------------------------------------------
    do idx = 1, wavelength_bins
      o3_transmission_factor(:,idx) = exp( -o3_cross_section(idx) * o3_slant(:) )
    end do

    !------------------------------------------------------------------------------
    !   ... Call calculate_photolysis_rate_for_interval Function to derive SR Band JNO contributions
    !         Called in order of wavelength interval (shortest first)
    ! TODO: what are 90, 100, and 50?
    !------------------------------------------------------------------------------
    do lev = 1, num_vertical_layers
      jno100   = calculate_photolysis_rate_for_interval( 1, o2_2100_cross_section, no_100_weighting_factor, no_100_cross_section )
      jno90    = calculate_photolysis_rate_for_interval( 2, o2_290_cross_section,  no_90_weighting_factor,  no_90_cross_section )
      jno50    = calculate_photolysis_rate_for_interval( 4, o2_250_cross_section,  no_50_weighting_factor,  no_50_cross_section )
      jno(lev) = jno50 + jno90 + jno100
    end do

    contains

    function calculate_photolysis_rate_for_interval( wavelength_interval, o2_cross_section, no_weighting_factor, no_cross_section ) result(result)
      !------------------------------------------------------------------------------
      !   ... uses xsec at center of g subinterval for o2
      !           uses mean values for no
      !------------------------------------------------------------------------------
      ! For reference, the function being calculate is this:
      ! In latex:
      ! J_{NO}=\Delta\lambda \overline{I_0} T_{O3}(z) P(z) \sum_{i=1}^{6} \exp[-\sigma^i_{O2}N_{O2}(z)] \sum_{j=1}^{2} W^{i,j}_{NO}\sigma^{i,j}_{NO}\exp[-\sigma^{i,j}_{NO}N_{NO}(z)]
      ! In ASCII:
      ! J_NO = Delta_lambda * I_0 * T_O3(z) * P(z) * sum_{i=1}^{6} exp[-sigma^i_O2 * N_O2(z)] * sum_{j=1}^{2} W^{i,j}_NO * sigma^{i,j}_NO * exp[-sigma^{i,j}_NO * N_NO(z)]
      ! where:
      ! J_NO is the photolysis rate of NO
      ! Delta_lambda is the wavelength interval
      ! I_0 is the mean value of solar irradiance (extraterrestrial flux)
      ! T_O3(z) is the ozone transmission factor in the Hartley band and is calculated as T_O3 = exp[-sigma_O3 N_O3]
      !  where sigma_O3 is the absorption cross section of O3 and N_O3 is the slant column density of O3
      ! P(z) is the predissociation factor
      ! sigma^i_O2 is the absorption cross section of O2
      ! N_O2(z) is the slant column density of O2
      ! W^{i,j}_NO is a weighting factor, from (Minschwaner and Siskind, 1993): "The weighting factors represent the fraction of the total spectral interval which is occupied by the corresponding mean value of the cross section."
      ! sigma^{i,j}_NO is the absorption cross section of NO
      ! N_NO(z) is the slant column density of NO
      ! z is the height
      ! i and j are indices for the absorption cross sections of O2 and NO, respectively
      ! The sums represent the calculation of the opacity distribution functions for O2 and NO
    
      !----------------------------------------------------------------
      !... Dummy arguments
      !----------------------------------------------------------------
      integer, intent(in)     :: wavelength_interval
      real(dk),    intent(in) :: o2_cross_section(6)
      real(dk),    intent(in) :: no_cross_section(number_of_levels,number_of_intervals)
      real(dk),    intent(in) :: no_weighting_factor(number_of_levels,number_of_intervals)
      
      !----------------------------------------------------------------
      !... Function declarations
      !----------------------------------------------------------------
      real(dk) :: result
      
      !----------------------------------------------------------------
      !... Local variables
      !----------------------------------------------------------------
      integer  ::  j, i
      real(dk) :: no_optical_depth
      real(dk) :: no_transmission
      real(dk) :: o2_transmission
      real(dk) :: o2_optical_depth
      real(dk) :: rate
      real(dk) :: jno1
      real(dk) :: D ! spontaneous rate of predissociation
      real(dk) :: A ! spontaneous rate of emission
      real(dk) :: kq ! quenching rate of N2
      real(dk) :: predissociation_factor
      
      !----------------------------------------------------------------
      !... derive the photolysis frequency for no within a given
      !         srb (i.e., 5-0, 9-0, 10-0)
      !----------------------------------------------------------------
      rate = 0._dk
      do i = 1,6
        o2_optical_depth = o2_slant(lev) * o2_cross_section(i) * 1e-6
        ! TODO: What is the value of 50? units? where does it come from?
        if( o2_optical_depth < 50._dk ) then
          o2_transmission = exp( -o2_optical_depth )
        else
          o2_transmission = 0._dk
         end if
        jno1 = 0._dk
        do j = 1,2
          no_optical_depth = no_slant(lev)*no_cross_section(i,j)
          ! TODO: What is the value of 50? units? where does it come from?
          if( no_optical_depth < 50._dk ) then
            no_transmission = exp( -no_optical_depth )
          else
            no_transmission = 0._dk
          end if
          jno1 = jno1 + no_cross_section(i,j) * no_weighting_factor(i,j) * no_transmission
        end do
        rate = rate + jno1 * o2_transmission
      end do

      result = delta_wavelength(wavelength_interval) * extraterrestrial_flux(wavelength_interval) &
                * o3_transmission_factor(lev, wavelength_interval) * rate

      ! Values taken from (Minschwaner and Siskind, 1993)
      D = 1.65e9_dk ! s-1
      A = 5.1e7_dk ! s-1
      kq = 1.5e-9_dk ! cm3 s-1

      !----------------------------------------------------------------
      !... correct for the predissociation of the deltq 1-0
      !         transition in the srb (5-0)
      !----------------------------------------------------------------
      if( wavelength_interval == 4 ) then
        predissociation_factor = D/(A + D + (kq * n2_dens(lev)))
        result = predissociation_factor * result
      end if
      
    end function calculate_photolysis_rate_for_interval

  end function calculate_jno

  subroutine calculate_slant_path( nlev, altitude, zenith_angle_degrees, slant_path, number_of_crossed_layers )
    !=============================================================================!
    !   Subroutine calculate_slant_path                                           !
    !=============================================================================!
    !   PURPOSE:                                                                  !
    !   Calculate slant path over vertical depth ds/dh in spherical geometry.     !
    !   Calculation is based on:  A.Dahlback, and K.Stamnes, A new spheric model  !
    !   for computing the radiation field available for photolysis and heating    !
    !   at twilight, Planet.Space Sci., v39, n5, pp. 671-683, 1991 (Appendix B)   !
    !=============================================================================!
    !   PARAMETERS:                                                               !
    !   nlev    - INTEGER, number of specified altitude levels in the working (I) !
    !             grid                                                            !
    !   altitude       - REAL, specified altitude working grid (km)           (I) !
    !   zenith_angle_degrees - REAL, solar zenith angle (degrees)             (I) !
    !   slant_path    - REAL, slant path of direct beam through each          (O) !
    !   layer crossed when travelling from the top of the atmosphere to layer i;  !
    !             slant_path(i,j), i = 0..nlev-1, j = 1..nlev-1                             !
    !   number_of_crossed_layers     - INTEGER, number of layers crossed by   (O) !
    !             the direct beam when travelling from the top of the atmosphere  !
    !             to layer i;                                                     !
    !             number_of_crossed_layers(i), i = 0..nlev-1                                             !
    !=============================================================================!
    !   EDIT HISTORY:                                                             !
    !   Original: Taken By Doug Kinnison from Sasha Madronich, TUV Code, V4.1a,   !
    !             on 1/1/02                                                       !
    !=============================================================================!
    ! taken from CAM: https://github.com/ESCOMP/CAM/blob/ab476f9b7345cbefdc4cf67ff17f0fe85d8c7387/src/chemistry/mozart/mo_jshort.F90#L1126-L1266
    
    !------------------------------------------------------------------------------
    !       ... Dummy arguments
    !------------------------------------------------------------------------------

    use musica_ccpp_util, only: DEGREE_TO_RADIAN, EARTH_RADIUS_KM

    integer,  intent(in)   :: nlev                             ! number model vertical levels
    integer,  intent(out)  :: number_of_crossed_layers(0:nlev) ! see above
    real(dk), intent (in)  :: zenith_angle_degrees             ! zenith_angle
    real(dk), intent (in)  :: altitude(nlev)                   ! geometric altitude (km)
    real(dk), intent (out) :: slant_path(0:nlev,nlev)          ! see above
    
    
    !------------------------------------------------------------------------------
    !       ... Local variables
    !------------------------------------------------------------------------------
    real(dk) :: radius_at_current_height
    real(dk) :: zenith_angle_radians
    real(dk) :: rpsinz
    real(dk) :: const0
    real(dk) :: rj
    real(dk) :: rjp1
    real(dk) :: dsj
    real(dk) :: dhj
    real(dk) :: ga
    real(dk) :: gb
    real(dk) :: slant_path_sign
    real(dk) :: zd(0:nlev-1)

    integer :: i
    integer :: j
    integer :: k
    integer :: id
    integer :: nlayer
    
    !------------------------------------------------------------------------------
    !       ... set zenith angle in radians
    !------------------------------------------------------------------------------
    zenith_angle_radians = zenith_angle_degrees*DEGREE_TO_RADIAN
    const0 = sin( zenith_angle_radians )
    
    !------------------------------------------------------------------------------
    !       ... set number of layers:
    !------------------------------------------------------------------------------
    nlayer = nlev - 1

    !------------------------------------------------------------------------------
    !       ... include the elevation above sea level to the radius of the earth:
    !------------------------------------------------------------------------------
    radius_at_current_height = EARTH_RADIUS_KM + altitude(nlev)

    !------------------------------------------------------------------------------
    !       ... inverse coordinate of z
    !------------------------------------------------------------------------------
    do k = 0,nlayer
      zd(k) = altitude(k+1) - altitude(nlev)
    end do
    
    !------------------------------------------------------------------------------
    !       ... initialize dsdh(i,j), nid(i)
    !------------------------------------------------------------------------------
    number_of_crossed_layers(:) = 0
    do j = 1,nlev
      slant_path(:,j) = 0._dk
    end do
    
    !------------------------------------------------------------------------------
    !       ... calculate ds/dh of every layer
    !------------------------------------------------------------------------------
    do i = 0,nlayer
      rpsinz = (radius_at_current_height + zd(i)) * const0
      if( zenith_angle_degrees <= 90._dk .or. rpsinz >= radius_at_current_height ) then
        !------------------------------------------------------------------------------
        ! Find index of layer in which the screening height lies
        !------------------------------------------------------------------------------
        id = i
        if( zenith_angle_degrees > 90._dk ) then
          do j = 1,nlayer
            if( rpsinz < (zd(j-1) + radius_at_current_height) .and.  rpsinz >= (zd(j) + radius_at_current_height) ) then
            id = j
              exit
            end if
          end do
        end if
        do j = 1,id
          slant_path_sign = 1._dk
          if( j == id .and. id == i .and. zenith_angle_degrees > 90._dk ) then
            slant_path_sign = -1._dk
          end if
          rj   = radius_at_current_height + zd(j-1)
          rjp1 = radius_at_current_height + zd(j)
          dhj  = zd(j-1) - zd(j)
          ga   = max( rj*rj - rpsinz*rpsinz,0._dk )
          gb   = max( rjp1*rjp1 - rpsinz*rpsinz,0._dk )
          if( id > i .and. j == id ) then
            dsj = sqrt( ga )
          else
            dsj = sqrt( ga ) - slant_path_sign*sqrt( gb )
          end if
          slant_path(i,j) = dsj / dhj
        end do
        number_of_crossed_layers(i) = id
      else
        number_of_crossed_layers(i) = -1
      end if
    end do
  end subroutine calculate_slant_path

  subroutine calculate_slant_column_density( nlev, delta_z, slant_path, number_of_crossed_layers, species_concentration, slant_column )
    !=============================================================================!
    !   PURPOSE:                                                                  !
    !   Derive Column                                                             !
    !=============================================================================!
    !   PARAMETERS:                                                               !
    !   NLEV   - INTEGER, number of specified altitude levels in the working  (I) !
    !            grid                                                             !
    !   delta_z   - REAL, specified altitude working grid (km)                (I) !
    !   slant_path- REAL, slant path of direct beam through each layer crossed (O)!
    !             when travelling from the top of the atmosphere to layer i;      !
    !             DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                             !
    !   number_of_crossed_layers    - INTEGER, number of layers crossed by the    !
    !             direct beam when travelling from the top of the atmosphere to   !
    !             layer i;                                                        !
    !             NID(i), i = 0..NZ-1                                             !
    !            specified altitude at each specified wavelength                  !
    !   species_concentration - REAL, absorber concentration, molecules cm-3      !
    !   slant_column   - REAL, absorber Slant Column, molecules cm-2              !
    !=============================================================================!
    !   EDIT HISTORY:                                                             !
    !   09/01  Read in profile from an input file, DEK                            !
    !   01/02  Taken from Sasha Madronich's TUV code                              !
    !=============================================================================!
    ! taken from CAM: https://github.com/ESCOMP/CAM/blob/ab476f9b7345cbefdc4cf67ff17f0fe85d8c7387/src/chemistry/mozart/mo_jshort.F90#L1268C7-L1369C31
    
    
    !------------------------------------------------------------------------------
    !       ... Dummy arguments
    !------------------------------------------------------------------------------
    integer,  intent(in)    :: nlev
    integer,  intent(in)    :: number_of_crossed_layers(0:nlev)       ! see above
    real(dk), intent(in)    :: delta_z(nlev)       ! layer thickness (cm)
    real(dk), intent(in)    :: slant_path(0:nlev,nlev) ! see above
    real(dk), intent(in)    :: species_concentration(nlev)      ! absorber concentration (molec. cm-3)
    real(dk), intent(out)   :: slant_column(nlev)     ! absorber Slant Column (molec. cm-2)
    
    !------------------------------------------------------------------------------
    !       ... Local variables
    !------------------------------------------------------------------------------
    real(dk), parameter :: largest = 1.e+36_dk

    real(dk) :: sum
    real(dk) :: hscale
    real(dk) :: numer, denom
    real(dk) :: cz(nlev)

    integer :: id
    integer :: j
    integer :: k
    
    !------------------------------------------------------------------------------
    !     ... compute column increments (logarithmic integrals)
    !------------------------------------------------------------------------------
    do k = 1,nlev-1
      if( species_concentration(k) /= 0._dk .and. species_concentration(k+1) /= 0._dk ) then
        cz(nlev-k) = (species_concentration(k) - species_concentration(k+1))/log( species_concentration(k)/species_concentration(k+1) ) * delta_z(k)
      else
        cz(nlev-k) = .5_dk*(species_concentration(k) + species_concentration(k+1)) * delta_z(k)
      end if
    end do
    
    !------------------------------------------------------------------------------
    !     ... Include exponential tail integral from infinity to model top
    !         specify scale height near top of data.For WACCM-X model, scale
    !         height needs to be increased for higher model top
    !------------------------------------------------------------------------------
    ! TODO: what should I do about nlev and pver here?
    ! if (nlev==pver) then
      ! if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
        ! hscale     = 20.e5_dk
      ! else
      !   hscale     = 10.e5_dk
      ! endif
      ! cz(nlev-1) = cz(nlev-1) + hscale * absden(1)
    ! endif
    
    !------------------------------------------------------------------------------
    !       ...  Calculate vertical and slant column from each level:
    !            work downward
    !------------------------------------------------------------------------------
    do id = 0,nlev-1
      sum = 0._dk
      if( number_of_crossed_layers(id) >= 0 ) then
        !------------------------------------------------------------------------------
        !       ...  Single pass layers:
        !------------------------------------------------------------------------------
        do j = 1, min(number_of_crossed_layers(id), id)
            sum = sum + cz(nlev-j)*slant_path(id,j)
        end do
        !------------------------------------------------------------------------------
        !       ...  Double pass layers:
        !------------------------------------------------------------------------------
        do j = min(number_of_crossed_layers(id),id)+1, number_of_crossed_layers(id)
            sum = sum + 2._dk*cz(nlev-j)*slant_path(id,j)
        end do
      else
        sum = largest
      end if
      slant_column(nlev-id) = sum
    end do
    slant_column(nlev) = .95_dk*slant_column(nlev-1)
  end subroutine calculate_slant_column_density

end module musica_ccpp_tuvx_no_photolysis_rate