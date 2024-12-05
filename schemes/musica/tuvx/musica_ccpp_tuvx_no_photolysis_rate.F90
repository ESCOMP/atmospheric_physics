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
  public :: calculate_NO_photolysis_rate, convert_mixing_ratio_to_molecule_cm3

contains

  !> Convert mixing ratio to molecule cm-3
  subroutine convert_mixing_ratio_to_molecule_cm3(mixing_ratio, dry_air_density, molar_mass, molecule_cm3)
    use ccpp_kinds, only: kind_phys
    real(kind_phys), intent(in)  :: mixing_ratio(:)       ! kg kg-1
    real(kind_phys), intent(in)  :: dry_air_density(:)    ! kg m-3
    real(kind_phys), intent(in)  :: molar_mass            ! kg mol-1
    real(kind_phys), intent(out) :: molecule_cm3(:)       ! molecule cm-3

    integer :: i
    real(kind_phys), parameter :: avogadro_number = 6.02214076e23 ! mol-1

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
    integer       :: pver
    real(kind_phys) :: dsdh(0:size(constituents, dim=2)+1,size(constituents, dim=2)+1)
    ! layer thickness (cm)
    real(kind_phys) :: delz(size(constituents, dim=2)+1)
    ! conversion from km to cm
    real(kind_phys), parameter :: km2cm = 1.0e5_r8
    ! final photolysis rate
    real(kind_phys) :: jNO

    ! number of vertical levels
    pver = size(constituents, dim=2)

    ! what are these constants? scale heights?
    call convert_mixing_ratio_to_molecule_cm3(constituents(:,:,N2_index), dry_air_density, molar_mass_N2, n2_dens(2:))
    n2_dens(1) = n2_dens(2) * 0.9_r8
    call convert_mixing_ratio_to_molecule_cm3(constituents(:,:,O2_index), dry_air_density, molar_mass_O2, o2_dens(2:))
    o2_dens(1) = o2_dens(2) * 7.0_r8 / ( height_at_interfaces(1) - height_at_interfaces(2) )
    call convert_mixing_ratio_to_molecule_cm3(constituents(:,:,O3_index), dry_air_density, molar_mass_O3, o3_dens(2:))
    o3_dens(1) = o3_dens(2) * 7.0_r8 / ( height_at_interfaces(1) - height_at_interfaces(2) )
    call convert_mixing_ratio_to_molecule_cm3(constituents(:,:,NO_index), dry_air_density, molar_mass_NO, no_dens(2:))
    no_dens(1) = no_dens(2) * 0.9_r8

    ! ================================
    ! calculate slant column densities
    ! ================================
    call sphers( pver+1, height_int, solar_zenith_angle, dsdh, nid )
    delz(1:pver) = km2cm * ( height_int(1:pver) - height_int(2:pver+1) )
    call slant_col( pver+1, delz, dsdh, nid, o2_dens, o2_slant )
    call slant_col( pver+1, delz, dsdh, nid, o3_dens, o3_slant )
    call slant_col( pver+1, delz, dsdh, nid, no_dens, no_slant )


    jNO = calculate_jno()

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
    real(kind_phys) :: predissociation_factor 
    real(kind_phys) :: D ! spontaneous rate of predissociation
    real(kind_phys) :: A ! spontaneous rate of emission
    real(kind_phys) :: kq ! quenching rate of N2
    real(kind_phys) :: jno

    jno = 0.0_r8

    ! Values taken from (Minschwaner and Siskind, 1993)
    D = 1.65e9_r8 ! s-1
    A = 5.1e7_r8 ! s-1
    kq = 1.5e-9_r8 ! cm3 s-1

    predissociation_factor = 0.0_r8

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
    real(r8), intent (in)  :: zenith_angle		  ! zenith_angle
    real(r8), intent (in)  :: z(nlev)		        ! geometric altitude (km)
    real(r8), intent (out) :: dsdh(0:nlev,nlev) ! see above
    
    
    !------------------------------------------------------------------------------
    !       ... Local variables
    !------------------------------------------------------------------------------
    real(r8) :: radius
    real(r8) :: re
    real(r8) :: zenrad
    real(r8) :: rpsinz
    real(r8) :: const0
    real(r8) :: rj
    real(r8) :: rjp1
    real(r8) :: dsj
    real(r8) :: dhj
    real(r8) :: ga
    real(r8) :: gb
    real(r8) :: sm
    real(r8) :: zd(0:nlev-1)

    integer :: i
    integer :: j
    integer :: k
    integer :: id
    integer :: nlayer


    radius = 6.37100e9_r8 ! radius earth (km)
    
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
      dsdh(:,j) = 0._r8
    end do
    
    !------------------------------------------------------------------------------
    !       ... calculate ds/dh of every layer
    !------------------------------------------------------------------------------
    do i = 0,nlayer
      rpsinz = (re + zd(i)) * const0
      if( zenith_angle <= 90._r8 .or. rpsinz >= re ) then
    !------------------------------------------------------------------------------
    ! Find index of layer in which the screening height lies
    !------------------------------------------------------------------------------
          id = i
          if( zenith_angle > 90._r8 ) then
            do j = 1,nlayer
                if( rpsinz < (zd(j-1) + re) .and.  rpsinz >= (zd(j) + re) ) then
      id = j
      exit
    end if
            end do
          end if

          do j = 1,id
            sm = 1._r8
            if( j == id .and. id == i .and. zenith_angle > 90._r8 ) then
              sm = -1._r8
            end if
            rj   = re + zd(j-1)
            rjp1 = re + zd(j)
            dhj  = zd(j-1) - zd(j)
            ga   = max( rj*rj - rpsinz*rpsinz,0._r8 )
            gb   = max( rjp1*rjp1 - rpsinz*rpsinz,0._r8 )
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
    real(r8), intent(in)    :: delz(nlev)	       ! layer thickness (cm)
    real(r8), intent(in)    :: dsdh(0:nlev,nlev) ! see above
    real(r8), intent(in)    :: absden(nlev)      ! absorber concentration (molec. cm-3)
    real(r8), intent(out)   :: scol(nlev)		     ! absorber Slant Column (molec. cm-2)
    
    !------------------------------------------------------------------------------
    !       ... Local variables
    !------------------------------------------------------------------------------
    real(r8), parameter :: largest = 1.e+36_r8

    real(r8) :: sum
    real(r8) :: hscale
    real(r8) :: numer, denom
    real(r8) :: cz(nlev)

    integer :: id
    integer :: j
    integer :: k
    
    !------------------------------------------------------------------------------
    !     ... compute column increments (logarithmic integrals)
    !------------------------------------------------------------------------------
    do k = 1,nlev-1
      if( absden(k) /= 0._r8 .and. absden(k+1) /= 0._r8 ) then
        cz(nlev-k) = (absden(k) - absden(k+1))/log( absden(k)/absden(k+1) ) * delz(k)
      else
        cz(nlev-k) = .5_r8*(absden(k) + absden(k+1)) * delz(k)
      end if
    end do
    
    !------------------------------------------------------------------------------
    !     ... Include exponential tail integral from infinity to model top
    !         specify scale height near top of data.For WACCM-X model, scale
    !         height needs to be increased for higher model top
    !------------------------------------------------------------------------------
    if (nlev==pver) then
      ! if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
        hscale     = 20.e5_r8
      ! else
      !   hscale     = 10.e5_r8
      ! endif
      cz(nlev-1) = cz(nlev-1) + hscale * absden(1)
    endif
    
    !------------------------------------------------------------------------------
    !       ...  Calculate vertical and slant column from each level:
    !            work downward
    !------------------------------------------------------------------------------
    do id = 0,nlev-1
      sum = 0._r8
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
            sum = sum + 2._r8*cz(nlev-j)*dsdh(id,j)
        end do
      else
        sum = largest
      end if
      scol(nlev-id) = sum
    end do
    scol(nlev) = .95_r8*scol(nlev-1)
  end subroutine slant_col

end module musica_ccpp_tuvx_no_photolysis_rate