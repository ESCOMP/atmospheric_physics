module musica_ccpp_tuvx_no_photolysis_rate

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

  use ccpp_kinds, only: kind_phys 
  dk => kind_phys

  implicit none

  private
  public :: calculate_NO_photolysis_rate, convert_mixing_ratio_to_molecule_cm3

  !------------------------------------------------------------------------------
  !   	... O3 SRB Cross Sections from WMO 1985, interpolated onto MS, 1993 grid
  !------------------------------------------------------------------------------
  ! in cm2
  real(kind_phys), save :: o3_cross_section(4) = (/ 7.3307600e-19_dk, 6.9660105E-19_dk, 5.9257699E-19_dk, 4.8372219E-19_dk /)

  !------------------------------------------------------------------------------
  !   	... delta wavelength of the MS, 1993 grid
  !------------------------------------------------------------------------------
  ! in nm
  ! TODO: What is the grid? The paper only has 3 bands listed in the table, but this has 4...
  real(kind_phys), save :: delta_wavelength(4) = (/ 1.50_dk, 1.50_dk, 5.6_dk, 2.3_dk /)
  
  !------------------------------------------------------------------------------
  !   	... O2 SRB Cross Sections for the six ODF regions, MS, 1993
  !------------------------------------------------------------------------------
  ! TODO: What are 250, 290, and 2100? 
  ! TODO: What are the units? cm2  molecule-1 ?
  real(kind_phys), save :: cross_section250(6)  = (/ 1.117e-23_dk, 2.447e-23_dk, 7.188e-23_dk, 3.042e-22_dk, 1.748e-21_dk, 1.112e-20_dk /)
  real(kind_phys), save :: cross_section290(6)  = (/ 1.350e-22_dk, 2.991e-22_dk, 7.334e-22_dk, 3.074e-21_dk, 1.689e-20_dk, 1.658e-19_dk /)
  real(kind_phys), save :: cross_section2100(6) = (/ 2.968e-22_dk, 5.831e-22_dk, 2.053e-21_dk, 8.192e-21_dk, 4.802e-20_dk, 2.655e-19_dk /)

  real(kind_phys), dimension(24) :: _a, _b, _c
  real(kind_phys) :: wtno50(6,2)
  real(kind_phys) :: wtno90(6,2)
  real(kind_phys) :: wtno100(6,2)

  !------------------------------------------------------------------------------
  !   	... 6 sub-intervals for O2 5-0 at 265K,
  !	    2 sub-sub-intervals for NO 0-0 at 250K
  !------------------------------------------------------------------------------
  _a(:) = (/    0._dk,       0._dk,       0._dk,       0._dk, &
                5.12e-02_dk, 5.68e-03_dk, 1.32e-18_dk, 4.41e-17_dk, &
                1.36e-01_dk, 1.52e-02_dk, 6.35e-19_dk, 4.45e-17_dk, &
                1.65e-01_dk, 1.83e-02_dk, 7.09e-19_dk, 4.50e-17_dk, &
                1.41e-01_dk, 1.57e-02_dk, 2.18e-19_dk, 2.94e-17_dk, &
                4.50e-02_dk, 5.00e-03_dk, 4.67e-19_dk, 4.35e-17_dk /)
  
  !------------------------------------------------------------------------------
  !   	... sub-intervals for o2 9-0 band,
  !	    2 sub-sub-intervals for no 1-0 at 250 k
  !------------------------------------------------------------------------------
  _b(:) = (/        0._dk,       0._dk,       0._dk,       0._dk, &
                    0._dk,       0._dk,       0._dk,       0._dk, &
              1.93e-03_dk, 2.14e-04_dk, 3.05e-21_dk, 3.20e-21_dk, &
              9.73e-02_dk, 1.08e-02_dk, 5.76e-19_dk, 5.71e-17_dk, &
              9.75e-02_dk, 1.08e-02_dk, 2.29e-18_dk, 9.09e-17_dk, &
              3.48e-02_dk, 3.86e-03_dk, 2.21e-18_dk, 6.00e-17_dk /)
  
  !------------------------------------------------------------------------------
  ! 	... sub-intervals for o2 10-0 band,
  !	    2 sub-sub-intervals for no 1-0 at 250 k
  !------------------------------------------------------------------------------
  _c(:) = (/  4.50e-02_dk, 5.00e-03_dk, 1.80e-18_dk, 1.40e-16_dk, &
              1.80e-01_dk, 2.00e-02_dk, 1.50e-18_dk, 1.52e-16_dk, &
              2.25e-01_dk, 2.50e-02_dk, 5.01e-19_dk, 7.00e-17_dk, &
              2.25e-01_dk, 2.50e-02_dk, 7.20e-20_dk, 2.83e-17_dk, &
              1.80e-01_dk, 2.00e-02_dk, 6.72e-20_dk, 2.73e-17_dk, &
              4.50e-02_dk, 5.00e-03_dk, 1.49e-21_dk, 6.57e-18_dk /)
  
  wtno50 (1:6,1) = _a(1:24:4)
  wtno50 (1:6,2) = _a(2:24:4)
  csno50 (1:6,1) = _a(3:24:4)
  csno50 (1:6,2) = _a(4:24:4)
  wtno90 (1:6,1) = _b(1:24:4)
  wtno90 (1:6,2) = _b(2:24:4)
  csno90 (1:6,1) = _b(3:24:4)
  csno90 (1:6,2) = _b(4:24:4)
  wtno100(1:6,1) = _c(1:24:4)
  wtno100(1:6,2) = _c(2:24:4)
  csno100(1:6,1) = _c(3:24:4)
  csno100(1:6,2) = _c(4:24:4)

contains

  !> Convert mixing ratio to molecule cm-3
  subroutine convert_mixing_ratio_to_molecule_cm3(mixing_ratio, dry_air_density, molar_mass, molecule_cm3)
    real(kind_phys), intent(in)  :: mixing_ratio(:)       ! kg kg-1
    real(kind_phys), intent(in)  :: dry_air_density(:)    ! kg m-3
    real(kind_phys), intent(in)  :: molar_mass            ! kg mol-1
    real(kind_phys), intent(out) :: molecule_cm3(:)       ! molecule cm-3

    integer :: i
    real(kind_phys), parameter :: avogadro_number = 6.02214076e23_dk ! mol-1

    do i = 1, size(mixing_ratio)
      molecule_cm3(i) = mixing_ratio(i) * dry_air_density(i) / molar_mass * avogadro_number * 1.0e-6
    end do
  end subroutine convert_mixing_ratio_to_molecule_cm3

  !> Prepare to calculate the photolysis rate of NO (nitiric acid). 
  !> This function transforms the inputs into the units
  !> needed by the calculate_jno routine which actually produces the photolysis rate
  function calculate_NO_photolysis_rate(solar_zenith_angle, extraterrestrial_flux, constituents, height_at_interfaces, &
    dry_air_density, N2_index, O2_index, O3_index, NO_index, molar_mass_N2, molar_mass_O2, molar_mass_O3, molar_mass_NO) &
      result(jNO)
    use ccpp_kinds,          only: kind_phys
    ! inputs
    real(kind_phys), intent(in)            :: solar_zenith_angle       ! degrees
    real(kind_phys), intent(in)            :: extraterrestrial_flux(:) ! photons cm-2 s-1 nm-1
    real(kind_phys), target, intent(inout) :: constituents(:,:,:)      ! various (column, layer, constituent)
    real(kind_phys), intent(in)            :: height_at_interfaces(:)  ! km
    real(kind_phys), intent(in)            :: dry_air_density(:,:)     ! kg m-3 (column, layer)
    integer        , intent(in)            :: N2_index, O2_index, O3_index, NO_index ! position of these species in the constituent arrays
    real(kind_phys), intent(in)            :: molar_mass_N2, molar_mass_O2, molar_mass_O3, molar_mass_NO

    ! local variables
    ! species column densities (molecule cm-3)
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
    integer       :: num_vertical_layers
    real(kind_phys) :: dsdh(0:size(constituents, dim=2)+1,size(constituents, dim=2)+1)
    ! layer thickness (cm)
    real(kind_phys) :: delz(size(constituents, dim=2)+1)
    ! conversion from km to cm
    real(kind_phys), parameter :: km2cm = 1.0e5_dk
    ! final photolysis rate
    real(kind_phys) :: jNO

    ! number of vertical levels
    num_vertical_layers = size(constituents, dim=2)

    ! TODO: what are these constants? scale heights?
    ! TODO: the values at index 1 appear to be for values above the model top in CAM, but how does that affect cam sima?
    call convert_mixing_ratio_to_molecule_cm3(constituents(:,:,N2_index), dry_air_density, molar_mass_N2, n2_dens(2:))
    n2_dens(1) = n2_dens(2) * 0.9_dk

    call convert_mixing_ratio_to_molecule_cm3(constituents(:,:,O2_index), dry_air_density, molar_mass_O2, o2_dens(2:))
    o2_dens(1) = o2_dens(2) * 7.0_dk / ( height_at_interfaces(1) - height_at_interfaces(2) )

    call convert_mixing_ratio_to_molecule_cm3(constituents(:,:,O3_index), dry_air_density, molar_mass_O3, o3_dens(2:))
    o3_dens(1) = o3_dens(2) * 7.0_dk / ( height_at_interfaces(1) - height_at_interfaces(2) )

    call convert_mixing_ratio_to_molecule_cm3(constituents(:,:,NO_index), dry_air_density, molar_mass_NO, no_dens(2:))
    no_dens(1) = no_dens(2) * 0.9_dk

    ! ================================
    ! calculate slant column densities
    ! ================================
    call sphers( num_vertical_layers+1, height_int, solar_zenith_angle, dsdh, nid )
    delz(1:num_vertical_layers) = km2cm * ( height_int(1:num_vertical_layers) - height_int(2:num_vertical_layers+1) )
    call slant_col( num_vertical_layers+1, delz, dsdh, nid, o2_dens, o2_slant )
    call slant_col( num_vertical_layers+1, delz, dsdh, nid, o3_dens, o3_slant )
    call slant_col( num_vertical_layers+1, delz, dsdh, nid, no_dens, no_slant )


    jNO = calculate_jno(num_vertical_layers, extraterrestrial_flux, n2_dens, o2_slant, o3_slant, no_slant, work_jno)

  end function calculate_NO_photolysis_rate

  !> Calculate the photolysis rate of NO (nitric acid)
  function calculate_jno(num_vertical_layers, extraterrestrial_flux, n2_dens, o2_slant, o3_slant, no_slant, work_jno) 
    result(jno)

    use ccpp_kinds, only: kind_phys
    ! inputs
    integer, intent(in)            :: num_vertical_layers
    real(kind_phys), intent(in)    :: extraterrestrial_flux(:) ! photons cm-2 s-1 nm-1
    real(kind_phys), intent(in)    :: n2_dens(:)              ! molecule cm-3
    real(kind_phys), intent(in)    :: o2_slant(:)             ! molecule cm-2
    real(kind_phys), intent(in)    :: o3_slant(:)             ! molecule cm-2
    real(kind_phys), intent(in)    :: no_slant(:)             ! molecule cm-2
    real(kind_phys), intent(inout) :: work_jno(:)             ! various

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

    ! local variables
    ! NO in an excited electronic state may result in an emission of NO or be quenched by an interaction 
    ! with N2. The predissociation faction is the probability that a precissociation of NO will occur from
    ! an excited electronic state
    real(kind_phys) :: o3_transmission_factor
    real(kind_phys) :: jno
    integer         :: wavelength_bins ! wavelength bins for MS, 93
    integer         :: idx 

    jno = 0.0_dk
    wavelength_bins = 4

    ! ... derive O3 transmission for the three O2 SRB
    ! ... idx = 1,2, and 4 are used below for jno
    ! ... TODO: why is 3 not used?
    !------------------------------------------------------------------------------
    do idx = 1,wavelength_bins
      o3_transmission_factor(:,idx) = exp( -o3_cross_section(idx)*o3_slant(:) )
    end do


    function pjno( w, cso2, wtno, csno )
      !------------------------------------------------------------------------------
      !   	... uses xsec at center of g subinterval for o2
      !           uses mean values for no
      !------------------------------------------------------------------------------
      
      !------------------------------------------------------------------------------
      !	... parameters
      !------------------------------------------------------------------------------
      integer, parameter :: ngint = 6
      integer, parameter :: nno = 2
    
      !----------------------------------------------------------------
      !	... Dummy arguments
      !----------------------------------------------------------------
      integer, intent(in)     :: w
      real(dk),    intent(in) :: cso2(ngint)
      real(dk),    intent(in) :: csno(ngint,nno)
      real(dk),    intent(in) :: wtno(ngint,nno)
      
      !----------------------------------------------------------------
      !	... Function declarations
      !----------------------------------------------------------------
      real(dk) :: pjno
      
      !----------------------------------------------------------------
      !	... Local variables
      !----------------------------------------------------------------
      integer  ::  jj, i, k
      real(dk) :: tauno
      real(dk) :: transno
      real(dk) :: transo2
      real(dk) :: tauo2
      real(dk) :: j
      real(dk) :: jno1
      real(kind_phys) :: D ! spontaneous rate of predissociation
      real(kind_phys) :: A ! spontaneous rate of emission
      real(kind_phys) :: kq ! quenching rate of N2
      
      !----------------------------------------------------------------
      !	... derive the photolysis frequency for no within a given
      !         srb (i.e., 5-0, 9-0, 10-0)
      !----------------------------------------------------------------
      j = 0._dk
      do k = 1,ngint
        tauo2 = o2scol(lev) * cso2(k)
        if( tauo2 < 50._dk ) then
          transo2 = exp( -tauo2 )
        else
          transo2 = 0._dk
         end if
        jno1 = 0._dk
        do jj = 1,nno
          tauno = noscol(lev)*csno(k,jj)
          if( tauno < 50._dk ) then
            transno = exp( -tauno )
          else
            transno = 0._dk
          end if
          jno1 = jno1 + csno(k,jj) * wtno(k,jj) * transno
        end do
        j = j + jno1*transo2
      end do
      
      pjno = wlintv_ms93(w)*etfphot_ms93(w)*tauo3(lev,w)*jno

      ! Values taken from (Minschwaner and Siskind, 1993)
      D = 1.65e9_dk ! s-1
      A = 5.1e7_dk ! s-1
      kq = 1.5e-9_dk ! cm3 s-1

      !----------------------------------------------------------------
      !	... correct for the predissociation of the deltq 1-0
      !         transition in the srb (5-0)
      !----------------------------------------------------------------
      if( w == 4 ) then
        pjno = D/(A + D + (kq*n2cc(nlev-lev+1)))*pjno
      end if
      
    end function pjno


  end function calculate_jno

  subroutine sphers( nlev, z, zenith_angle, dsdh, nid )
    !=============================================================================!
    !   Subroutine sphers                                                         !
    !=============================================================================!
    !   PURPOSE:                                                                  !
    !   Calculate slant path over vertical depth ds/dh in spherical geometry.     !
    !   Calculation is based on:  A.Dahlback, and K.Stamnes, A new spheric model  !
    !   for computing the radiation field available for photolysis and heating    !
    !   at twilight, Planet.Space Sci., v39, n5, pp. 671-683, 1991 (Appendix B)   !
    !=============================================================================!
    !   PARAMETERS:                                                               !
    !   NZ      - INTEGER, number of specified altitude levels in the working (I) !
    !             grid                                                            !
    !   Z       - REAL, specified altitude working grid (km)                  (I) !
    !   ZEN     - REAL, solar zenith angle (degrees)                          (I) !
    !   DSDH    - REAL, slant path of direct beam through each layer crossed  (O) !
    !             when travelling from the top of the atmosphere to layer i;      !
    !             DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                             !
    !   NID     - INTEGER, number of layers crossed by the direct beam when   (O) !
    !             travelling from the top of the atmosphere to layer i;           !
    !             NID(i), i = 0..NZ-1                                             !
    !=============================================================================!
    !   EDIT HISTORY:                                                             !
    !   Original: Taken By Doug Kinnison from Sasha Madronich, TUV Code, V4.1a,   !
    !             on 1/1/02                                                       !
    !=============================================================================!
    ! taken from CAM: https://github.com/ESCOMP/CAM/blob/ab476f9b7345cbefdc4cf67ff17f0fe85d8c7387/src/chemistry/mozart/mo_jshort.F90#L1126-L1266
    
    !------------------------------------------------------------------------------
    !       ... Dummy arguments
    !------------------------------------------------------------------------------
    integer,  intent(in)   :: nlev              ! number model vertical levels
    integer,  intent(out)  :: nid(0:nlev)       ! see above
    real(dk), intent (in)  :: zenith_angle		  ! zenith_angle
    real(dk), intent (in)  :: z(nlev)		        ! geometric altitude (km)
    real(dk), intent (out) :: dsdh(0:nlev,nlev) ! see above
    
    
    !------------------------------------------------------------------------------
    !       ... Local variables
    !------------------------------------------------------------------------------
    real(dk) :: radius
    real(dk) :: re
    real(dk) :: zenrad
    real(dk) :: rpsinz
    real(dk) :: const0
    real(dk) :: rj
    real(dk) :: rjp1
    real(dk) :: dsj
    real(dk) :: dhj
    real(dk) :: ga
    real(dk) :: gb
    real(dk) :: sm
    real(dk) :: zd(0:nlev-1)

    integer :: i
    integer :: j
    integer :: k
    integer :: id
    integer :: nlayer


    ! TODO: Get this from CAM-SIMA
    radius = 6.37100e9_dk ! radius earth (km)
    
    !------------------------------------------------------------------------------
    !       ... set zenith angle in radians
    !------------------------------------------------------------------------------
    zenrad = zenith_angle*d2r
    const0 = sin( zenrad )
    
    !------------------------------------------------------------------------------
    !       ... set number of layers:
    !------------------------------------------------------------------------------
    nlayer = nlev - 1

    !------------------------------------------------------------------------------
    !       ... include the elevation above sea level to the radius of the earth:
    !------------------------------------------------------------------------------
    re = radius + z(nlev)

    !------------------------------------------------------------------------------
    !       ... inverse coordinate of z
    !------------------------------------------------------------------------------
    do k = 0,nlayer
      zd(k) = z(k+1) - z(nlev)
    end do
    
    !------------------------------------------------------------------------------
    !       ... initialize dsdh(i,j), nid(i)
    !------------------------------------------------------------------------------
    nid(:) = 0
    do j = 1,nlev
      dsdh(:,j) = 0._dk
    end do
    
    !------------------------------------------------------------------------------
    !       ... calculate ds/dh of every layer
    !------------------------------------------------------------------------------
    do i = 0,nlayer
      rpsinz = (re + zd(i)) * const0
      if( zenith_angle <= 90._dk .or. rpsinz >= re ) then
    !------------------------------------------------------------------------------
    ! Find index of layer in which the screening height lies
    !------------------------------------------------------------------------------
          id = i
          if( zenith_angle > 90._dk ) then
            do j = 1,nlayer
                if( rpsinz < (zd(j-1) + re) .and.  rpsinz >= (zd(j) + re) ) then
      id = j
      exit
    end if
            end do
          end if

          do j = 1,id
            sm = 1._dk
            if( j == id .and. id == i .and. zenith_angle > 90._dk ) then
              sm = -1._dk
            end if
            rj   = re + zd(j-1)
            rjp1 = re + zd(j)
            dhj  = zd(j-1) - zd(j)
            ga   = max( rj*rj - rpsinz*rpsinz,0._dk )
            gb   = max( rjp1*rjp1 - rpsinz*rpsinz,0._dk )
            if( id > i .and. j == id ) then
              dsj = sqrt( ga )
            else
              dsj = sqrt( ga ) - sm*sqrt( gb )
            end if
            dsdh(i,j) = dsj / dhj
          end do
          nid(i) = id
      else
          nid(i) = -1
      end if
    end do
  end subroutine sphers

  subroutine slant_col( nlev, delz, dsdh, nid, absden, scol )
    !=============================================================================!
    !   PURPOSE:                                                                  !
    !   Derive Column
    !=============================================================================!
    !   PARAMETERS:                                                               !
    !   NLEV   - INTEGER, number of specified altitude levels in the working  (I) !
    !            grid                                                             !
    !   DELZ   - REAL, specified altitude working grid (km)                   (I) !
    !   DSDH   - REAL, slant path of direct beam through each layer crossed  (O)  !
    !             when travelling from the top of the atmosphere to layer i;      !
    !             DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                             !
    !   NID    - INTEGER, number of layers crossed by the direct beam when   (O)  !
    !             travelling from the top of the atmosphere to layer i;           !
    !             NID(i), i = 0..NZ-1                                             !
    !            specified altitude at each specified wavelength                  !
    !   absden - REAL, absorber concentration, molecules cm-3                     !
    !   SCOL   - REAL, absorber Slant Column, molecules cm-2                      !
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
    integer,  intent(in)    :: nid(0:nlev)       ! see above
    real(dk), intent(in)    :: delz(nlev)	       ! layer thickness (cm)
    real(dk), intent(in)    :: dsdh(0:nlev,nlev) ! see above
    real(dk), intent(in)    :: absden(nlev)      ! absorber concentration (molec. cm-3)
    real(dk), intent(out)   :: scol(nlev)		     ! absorber Slant Column (molec. cm-2)
    
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
      if( absden(k) /= 0._dk .and. absden(k+1) /= 0._dk ) then
        cz(nlev-k) = (absden(k) - absden(k+1))/log( absden(k)/absden(k+1) ) * delz(k)
      else
        cz(nlev-k) = .5_dk*(absden(k) + absden(k+1)) * delz(k)
      end if
    end do
    
    !------------------------------------------------------------------------------
    !     ... Include exponential tail integral from infinity to model top
    !         specify scale height near top of data.For WACCM-X model, scale
    !         height needs to be increased for higher model top
    !------------------------------------------------------------------------------
    if (nlev==pver) then
      ! if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
        hscale     = 20.e5_dk
      ! else
      !   hscale     = 10.e5_dk
      ! endif
      cz(nlev-1) = cz(nlev-1) + hscale * absden(1)
    endif
    
    !------------------------------------------------------------------------------
    !       ...  Calculate vertical and slant column from each level:
    !            work downward
    !------------------------------------------------------------------------------
    do id = 0,nlev-1
      sum = 0._dk
      if( nid(id) >= 0 ) then
        !------------------------------------------------------------------------------
        !       ...  Single pass layers:
        !------------------------------------------------------------------------------
        do j = 1, min(nid(id), id)
            sum = sum + cz(nlev-j)*dsdh(id,j)
        end do
        !------------------------------------------------------------------------------
        !       ...  Double pass layers:
        !------------------------------------------------------------------------------
        do j = min(nid(id),id)+1, nid(id)
            sum = sum + 2._dk*cz(nlev-j)*dsdh(id,j)
        end do
      else
        sum = largest
      end if
      scol(nlev-id) = sum
    end do
    scol(nlev) = .95_dk*scol(nlev-1)
  end subroutine slant_col

end module musica_ccpp_tuvx_no_photolysis_rate