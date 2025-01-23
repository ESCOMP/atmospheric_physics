! Copyright (C) 2002-2025 National Science Foundation-National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Compute tendencies from sedimentation of cloud liquid and ice particles
! Original authors: Byron Boville, September 2002 from code by P.J. Rasch
! CCPP-ized: Haipeng Lin, January 2025
module cloud_particle_sedimentation

  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  ! public CCPP-compliant subroutines
  public :: cloud_particle_sedimentation_init
  public :: cloud_particle_sedimentation_run

  ! tuning parameters for cloud liquid and ice particle sedimentation
  real(kind_phys), parameter :: vland = 1.5_kind_phys        ! liquid fall velocity over land [cm s-1]
  real(kind_phys), parameter :: vocean = 2.8_kind_phys       ! liquid fall velocity over ocean [cm s-1]
  real(kind_phys), parameter :: mxsedfac = 0.99_kind_phys    ! maximum sedimentation flux factor

  ! tuning parameters for Stokes terminal velocity
  logical,         parameter :: stokes = .true.              ! use Stokes velocity instead of McFarquhar and Heymsfield
  real(kind_phys)            :: cldsed_ice_stokes_fac        ! factor applied to ice fall vel
                                                             ! from Stokes terminal velocity
  real(kind_phys), parameter :: eta  = 1.7e-5_kind_phys      ! viscosity of air [kg m s-1]
  real(kind_phys), parameter :: r40  = 40._kind_phys         !  40 micron radius
  real(kind_phys), parameter :: r400 = 400._kind_phys        ! 400 micron radius
  real(kind_phys), parameter :: v400 = 1.00_kind_phys        ! fall velocity of 400 micron sphere [m s-1]
  real(kind_phys)            :: v40                          ! Stokes fall velocity of 40 micron sphere (m/s)
  ! v40 = (2._kind_phys/9._kind_phys) * rhoh2o * gravit/eta * r40**2 * 1.e-12_kind_phys
  real(kind_phys)            :: vslope                       ! linear slope for large particles [m s-1 micron-1]
  ! vslope = (v400 - v40)/(r400 -r40)

  ! tuning parameters for McFarquhar and Heymsfield terminal velocity
  real(kind_phys), parameter :: vice_small = 1._kind_phys    ! ice fall velocity for small concentration [cm s-1]
  real(kind_phys)            :: lbound                       ! lower bound for iciwc

  ! misc
  real(kind_phys), parameter :: um_to_m = 1.e-12_kind_phys

contains

!> \section arg_table_cloud_particle_sedimentation_init Argument Table
!! \htmlinclude arg_table_cloud_particle_sedimentation_init.html
  subroutine cloud_particle_sedimentation_init(&
    amIRoot, iulog, &
    cldsed_ice_stokes_fac_in, &
    rhoh2o, gravit, &
    errmsg, errflg)

    ! Input arguments
    logical,            intent(in)    :: amIRoot
    integer,            intent(in)    :: iulog          ! log output unit
    real(kind_phys),    intent(in)    :: cldsed_ice_stokes_fac_in
    real(kind_phys),    intent(in)    :: gravit
    real(kind_phys),    intent(in)    :: rhoh2o

    ! Output arguments
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    ! Local variables
    real(kind_phys)     :: ac, bc, cc                   ! constants for McFarquhar and Heymsfield method

    errmsg = ''
    errflg = 0
    cldsed_ice_stokes_fac = cldsed_ice_stokes_fac_in

    ! linear ramp variables for fall velocity
    v40 = (2._kind_phys/9._kind_phys)*rhoh2o*gravit/eta*r40**2*1.e-12_kind_phys
    vslope = (v400 - v40)/(r400 - r40)

    ! McFarquhar and Heymsfield lower bound for iciwc
    cc = 128.64_kind_phys
    bc = 53.242_kind_phys
    ac = 5.4795_kind_phys
    lbound = (-bc + sqrt(bc*bc - 4*ac*cc))/(2*ac)
    lbound = 10._kind_phys**lbound

    if(amIRoot) then
      write(iulog,*) 'cloud_particle_sedimentation_init: cldsed_ice_stokes_fac = ', cldsed_ice_stokes_fac
    endif

  end subroutine cloud_particle_sedimentation_init

  ! Compute gravitational sedimentation velocities for cloud liquid water
  ! and ice, based on Lawrence and Crutzen (1998).
  ! https://doi.org/10.3402/tellusb.v50i3.16129
  !
  ! and apply cloud particle gravitational sedimentation to condensate.
  !
  ! Original authors: mgl, March 1998 (MATCH-MPIC version 2.0)
  !                   adapted by P.J. Rasch
  !                   B.A. Boville, September 2002
  !                   P.J. Rasch, May 2003
!> \section arg_table_cloud_particle_sedimentation_run Argument Table
!! \htmlinclude arg_table_cloud_particle_sedimentation_run.html
  subroutine cloud_particle_sedimentation_run( &
    ncol, pver, pverp, dtime, &
    tmelt, gravit, latvap, latice, rair, rhoh2o, &
    icritc, & ! from prognostic_cloud_water namelist
    pint, pmid, pdel, t, cloud, &
    icefrac, landfrac, ocnfrac, &
    cldliq, cldice, snowh, landm, &
    ! below output for sedimentation velocity
    pvliq, pvice, &
    ! below output for sedimentation to condensate
    liqtend, icetend, wvtend, htend, sfliq, sfice, &
    errmsg, errflg)
    ! To-be-CCPPized dependencies
    use cloud_optical_properties,     only: reltab, reitab

    ! Input arguments
    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    integer,            intent(in)    :: pverp
    real(kind_phys),    intent(in)    :: dtime
    real(kind_phys),    intent(in)    :: tmelt
    real(kind_phys),    intent(in)    :: gravit
    real(kind_phys),    intent(in)    :: latvap
    real(kind_phys),    intent(in)    :: latice
    real(kind_phys),    intent(in)    :: rair
    real(kind_phys),    intent(in)    :: rhoh2o
    real(kind_phys),    intent(in)    :: icritc         ! tunable_parameter_for_autoconversion_of_cold_ice_for_rk_microphysics [kg kg-1]
    real(kind_phys),    intent(in)    :: pint(:,:)      ! air_pressure_at_interface [Pa]
    real(kind_phys),    intent(in)    :: pmid(:,:)      ! air_pressure [Pa]
    real(kind_phys),    intent(in)    :: pdel(:,:)      ! air_pressure_thickness [Pa]
    real(kind_phys),    intent(in)    :: t(:,:)         ! air_temperature [K]
    real(kind_phys),    intent(in)    :: cloud(:,:)     ! cloud_area_fraction [fraction]
    real(kind_phys),    intent(in)    :: icefrac(:)     ! sea_ice_area_fraction [fraction]
    real(kind_phys),    intent(in)    :: landfrac(:)    ! land_area_fraction [fraction]
    real(kind_phys),    intent(in)    :: ocnfrac(:)     ! ocean_area_fraction [fraction]
    real(kind_phys),    intent(in)    :: cldliq(:,:)    ! adv: cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
    real(kind_phys),    intent(in)    :: cldice(:,:)    ! adv: cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
    real(kind_phys),    intent(in)    :: snowh(:)       ! lwe_surface_snow_depth_over_land [m]
    real(kind_phys),    intent(in)    :: landm(:)       ! smoothed_land_area_fraction [fraction]

    ! Output arguments
    ! note: pvel is at the interfaces (loss from cell is based on pvel(k+1))
    real(kind_phys),    intent(out)   :: pvliq(:,:)     ! vertical velocity of cloud liquid drops [Pa s-1] pverp
    real(kind_phys),    intent(out)   :: pvice(:,:)     ! vertical velocity of cloud ice particles [Pa s-1] pverp
    real(kind_phys),    intent(out)   :: liqtend(:,:)   ! liquid condensate tendency -- to apply cldliq tendency
    real(kind_phys),    intent(out)   :: icetend(:,:)   ! ice condensate tendency -- to apply cldice tendency
    real(kind_phys),    intent(out)   :: wvtend(:,:)    ! water vapor tendency -- to apply wv tendency
    real(kind_phys),    intent(out)   :: htend(:,:)     ! heating rate [J kg-1 s-1] -- to apply s tendency
    real(kind_phys),    intent(out)   :: sfliq(:)       ! surface flux of liquid (rain) [kg m-2 s-1]
    real(kind_phys),    intent(out)   :: sfice(:)       ! lwe_cloud_ice_sedimentation_rate_at_surface_due_to_microphysics [m s-1]
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    ! Local variables
    real(kind_phys)     :: rho(ncol, pver)              ! air density [kg m-3]
    real(kind_phys)     :: vfall                        ! settling velocity of cloud particles [m s-1]
    real(kind_phys)     :: icice                        ! in cloud ice water content [kg kg-1]
    real(kind_phys)     :: iciwc                        ! in cloud ice water content [g m-3]
    real(kind_phys)     :: icefac
    real(kind_phys)     :: logiwc
    real(kind_phys)     :: rei(ncol, pver)              ! effective radius of ice particles [um]
    real(kind_phys)     :: rel(ncol, pver)              ! effective radius of liquid particles [um]
    real(kind_phys)     :: fxliq(ncol, pverp)           ! fluxes at interface, liquid (positive is down) [Pa]
    real(kind_phys)     :: fxice(ncol, pverp)           ! fluxes at interface, ice (positive is down) [Pa]
    real(kind_phys)     :: cldab(ncol)                  ! cloud in layer above [fraction]
    real(kind_phys)     :: evapliq                      ! evaporation of cloud liquid into environment [kg kg-1 s-1]
    real(kind_phys)     :: evapice                      ! evaporation of cloud ice into environment [kg kg-1 s-1]
    real(kind_phys)     :: cldovrl                      ! cloud overlap factor
    integer             :: i, k

    errmsg = ''
    errflg = 0

    !--------------------------------------------------
    ! COMPUTE GRAVITATIONAL SEDIMENTATION VELOCITIES
    !--------------------------------------------------
    ! NEED TO BE CAREFUL - VELOCITIES SHOULD BE AT THE *LOWER* INTERFACE
    ! (THAT IS, FOR K+1), FLUXES ARE ALSO AT THE LOWER INTERFACE (K+1),
    ! BUT MIXING RATIOS ARE AT THE MIDPOINTS (K)...
    !
    ! NOTE THAT PVEL IS ON INTERFACES, WHEREAS VFALL IS FOR THE CELL
    ! AVERAGES (I.E., MIDPOINTS); ASSUME THE FALL VELOCITY APPLICABLE TO THE
    ! LOWER INTERFACE (K+1) IS THE SAME AS THAT APPLICABLE FOR THE CELL (V(K))

    ! LIQUID:
    ! The fall velocities assume that droplets have a gamma distribution
    ! with effective radii for land and ocean as assessed by Han et al.;
    ! see Lawrence and Crutzen (1998) for a derivation.
    rho(:ncol,:) = pmid(:ncol,:)/(rair*t(:ncol,:))
    pvliq(:ncol,:) = 0._kind_phys

    ! get effective radius of liquid drop (rel)
    call reltab(ncol=ncol, pver=pver, tmelt=tmelt, t=t(:ncol,:), landfrac=landfrac(:ncol), landm=landm(:ncol), icefrac=icefrac(:ncol), snowh=snowh(:ncol), rel=rel(:ncol,:))

    do k = 1, pver
      do i = 1, ncol
        if (cloud(i, k) > 0._kind_phys .and. cldliq(i, k) > 0._kind_phys) then
          if (rel(i, k) < 40._kind_phys) then
            vfall = 2._kind_phys/9._kind_phys*rhoh2o*gravit*rel(i, k)**2/eta*um_to_m  ! um^2 -> m^2
          else
            vfall = v40 + vslope*(rel(i, k) - r40)    ! linear above 40 microns
          end if
          ! convert the fall speed to pressure units, but do not apply the traditional
          ! negative convention for pvel.
          pvliq(i, k + 1) = vfall*rho(i, k)*gravit      ! m s-1 to Pa s-1
        end if
      end do
    end do

    ! ICE:
    ! The fall velocities are based on data from McFarquhar and Heymsfield
    ! or on Stokes terminal velocity for spheres and the effective radius.
    pvice(:ncol,:) = 0._kind_phys

    if(stokes) then
      ! stokes terminal velocity < 40 microns

      ! get effective radius (rei)
      call reitab(ncol=ncol, pver=pver, t=t(:ncol,:), re=rei(:ncol,:))

      do k = 1, pver
        do i = 1, ncol
          if (cloud(i, k) > 0._kind_phys .and. cldice(i, k) > 0._kind_phys) then
            if (rei(i, k) < 40._kind_phys) then
              vfall = 2._kind_phys/9._kind_phys*rhoh2o*gravit*rei(i, k)**2/eta*um_to_m  ! microns^2 -> m^2
              vfall = vfall*cldsed_ice_stokes_fac
            else
              vfall = v40 + vslope*(rei(i, k) - r40)      ! linear above 40 microns
            end if

            ! convert the fall speed to pressure units, but do not apply the traditional
            ! negative convention for pvel.
            pvice(i, k + 1) = vfall*rho(i, k)*gravit      ! m s-1 to Pa s-1
          end if
        end do
      end do

    else
      ! McFarquhar and Heymsfield > icritc
      do k = 1, pver
        do i = 1, ncol
          if (cloud(i, k) > 0._kind_phys .and. cldice(i, k) > 0._kind_phys) then

            ! compute the in-cloud ice concentration (kg/kg)
            icice = cldice(i, k)/cloud(i, k)

            ! compute the ice water content in g/m3
            iciwc = icice*rho(i, k)*1.e3_kind_phys

            ! set the fall velocity (cm/s) to depend on the ice water content in g/m3,
            if (iciwc > lbound) then ! need this because of log10
              logiwc = log10(iciwc)
              !       Median -
              vfall = 128.64_kind_phys + 53.242_kind_phys*logiwc + 5.4795_kind_phys*logiwc**2
              !       Average -
              ! vfall = 122.63 + 44.111*logiwc + 4.2144*logiwc**2
            else
              vfall = 0._kind_phys
            end if

            ! set ice velocity to 1 cm/s if ice mixing ratio < icritc, ramp to value
            ! calculated above at 2*icritc
            if (icice <= icritc) then
              vfall = vice_small
            else if (icice < 2*icritc) then
              icefac = (icice - icritc)/icritc
              vfall = vice_small*(1._kind_phys - icefac) + vfall*icefac
            end if

            ! bound the terminal velocity of ice particles at high concentration
            vfall = min(100.0_kind_phys, vfall)

            ! convert the fall speed to pressure units, but do not apply the traditional
            ! negative convention for pvel.
            pvice(i, k + 1) = vfall &
                              *0.01_kind_phys &      ! cm to meters
                              *rho(i, k)*gravit      ! m s-1 to Pa s-1
          end if
        end do
      end do
    endif

    !--------------------------------------------------
    ! APPLY CLOUD PARTICLE GRAVITATIONAL SEDIMENTATION TO CONDENSATE
    !--------------------------------------------------
    ! Initialize variables
    fxliq(:ncol, :)   = 0._kind_phys ! flux at interfaces (liquid)
    fxice(:ncol, :)   = 0._kind_phys ! flux at interfaces (ice)
    liqtend(:ncol, :) = 0._kind_phys ! condensate tend (liquid)
    icetend(:ncol, :) = 0._kind_phys ! condensate tend (ice)
    wvtend(:ncol, :)  = 0._kind_phys ! environmental moistening
    htend(:ncol, :)   = 0._kind_phys ! evaporative cooling
    sfliq(:ncol)      = 0._kind_phys ! condensate sedimentation flux out bot of column (liquid)
    sfice(:ncol)      = 0._kind_phys ! condensate sedimentation flux out bot of column (ice)

    ! Get fluxes at interior points
    call getflx(ncol, pver, pverp, pint, cldliq, pvliq, dtime, fxliq, errmsg, errflg)
    if(errflg /= 0) then
      return
    endif
    call getflx(ncol, pver, pverp, pint, cldice, pvice, dtime, fxice, errmsg, errflg)
    if(errflg /= 0) then
      return
    endif

    ! Calculate fluxes at boundaries
    do i = 1, ncol
      fxliq(i, 1) = 0._kind_phys
      fxice(i, 1) = 0._kind_phys
      ! surface flux by upstream scheme
      fxliq(i, pverp) = cldliq(i, pver)*pvliq(i, pverp)*dtime
      fxice(i, pverp) = cldice(i, pver)*pvice(i, pverp)*dtime
    end do

    ! filter out any negative fluxes from the getflx routine
    ! (typical fluxes are of order > 1e-3 when clouds are present)
    do k = 2, pver
      fxliq(:ncol, k) = max(0._kind_phys, fxliq(:ncol, k))
      fxice(:ncol, k) = max(0._kind_phys, fxice(:ncol, k))
    end do

    ! Limit the flux out of the bottom of each cell to the water content in each phase.
    ! Apply mxsedfac to prevent generating very small negative cloud water/ice
    ! NOTE, REMOVED CLOUD FACTOR FROM AVAILABLE WATER. ALL CLOUD WATER IS IN CLOUDS.
    ! ***Should we include the flux in the top, to allow for thin surface layers?
    ! ***Requires simple treatment of cloud overlap, already included below.
    do k = 1, pver
      do i = 1, ncol
        fxliq(i, k + 1) = min(fxliq(i, k + 1), mxsedfac*cldliq(i, k)*pdel(i, k))
        fxice(i, k + 1) = min(fxice(i, k + 1), mxsedfac*cldice(i, k)*pdel(i, k))
        ! fxliq(i,k+1) = min( fxliq(i,k+1), cldliq(i,k) * pdel(i,k) + fxliq(i,k))
        ! fxice(i,k+1) = min( fxice(i,k+1), cldice(i,k) * pdel(i,k) + fxice(i,k))
        ! fxliq(i,k+1) = min( fxliq(i,k+1), cloud(i,k) * cldliq(i,k) * pdel(i,k) )
        ! fxice(i,k+1) = min( fxice(i,k+1), cloud(i,k) * cldice(i,k) * pdel(i,k) )
      end do
    end do

    ! Now calculate the tendencies assuming that condensate evaporates when
    ! it falls into environment, and does not when it falls into cloud.
    ! All flux out of the layer comes from the cloudy part.
    ! Assume maximum overlap for stratiform clouds
    !  if cloud above < cloud, all water falls into cloud below
    !  if cloud above > cloud, water split between cloud and environment
    do k = 1, pver
      cldab(:ncol) = 0._kind_phys
      do i = 1, ncol
        ! cloud overlap cloud factor
        cldovrl = min(cloud(i, k)/(cldab(i) + .0001_kind_phys), 1._kind_phys)
        cldab(i) = cloud(i, k)
        ! evaporation into environment cause moistening and cooling
        evapliq = fxliq(i, k)*(1._kind_phys - cldovrl)/(dtime*pdel(i, k))  ! into env [kg kg-1 s-1]
        evapice = fxice(i, k)*(1._kind_phys - cldovrl)/(dtime*pdel(i, k))  ! into env [kg kg-1 s-1]
        wvtend(i, k) = evapliq + evapice                                   ! evaporation into environment [kg kg-1 s-1]
        htend(i, k) = -latvap*evapliq - (latvap + latice)*evapice          ! evaporation [W kg-1]
        ! net flux into cloud changes cloud liquid/ice (all flux is out of cloud)
        liqtend(i, k) = (fxliq(i, k)*cldovrl - fxliq(i, k + 1))/(dtime*pdel(i, k))
        icetend(i, k) = (fxice(i, k)*cldovrl - fxice(i, k + 1))/(dtime*pdel(i, k))
      end do
    end do

    ! convert flux out the bottom to mass units Pa -> kg/m2/s
    sfliq(:ncol) = fxliq(:ncol, pverp)/(dtime*gravit)
    sfice(:ncol) = fxice(:ncol, pverp)/(dtime*gravit)

    ! Convert lwe_cloud_ice_sedimentation_rate_at_surface_due_to_microphysics from kg m-2 s-1 to precip units m s-1
    sfice(:ncol) = sfice(:ncol)/1000._kind_phys

  end subroutine cloud_particle_sedimentation_run

  ! Compute fluxes at interior points
  subroutine getflx(ncol, pver, pverp, &
    xw, phi, vel, deltat, flux, &
    errmsg, errflg)
    !.....xw1.......xw2.......xw3.......xw4.......xw5.......xw6
    !....psiw1.....psiw2.....psiw3.....psiw4.....psiw5.....psiw6
    !....velw1.....velw2.....velw3.....velw4.....velw5.....velw6
    !.........phi1......phi2.......phi3.....phi4.......phi5.......

    integer,         intent(in)  :: ncol
    integer,         intent(in)  :: pver
    integer,         intent(in)  :: pverp
    real(kind_phys), intent(in)  :: xw(ncol, pverp)
    real(kind_phys), intent(in)  :: phi(ncol, pver)
    real(kind_phys), intent(in)  :: vel(ncol, pverp)
    real(kind_phys), intent(in)  :: deltat

    real(kind_phys),    intent(out)   :: flux(ncol, pverp)
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    integer         :: i, k
    real(kind_phys) :: psi(ncol, pverp)
    real(kind_phys) :: fdot(ncol, pverp)
    real(kind_phys) :: xx(ncol)
    real(kind_phys) :: fxdot(ncol)
    real(kind_phys) :: fxdd(ncol)
    real(kind_phys) :: psistar(ncol)
    real(kind_phys) :: xxk(ncol, pver)

    do i = 1, ncol
      ! integral of phi
      psi(i, 1) = 0._kind_phys
      ! fluxes at boundaries
      flux(i, 1) = 0._kind_phys
      flux(i, pverp) = 0._kind_phys
    end do

    ! integral function
    do k = 2, pverp
      do i = 1, ncol
        psi(i, k) = phi(i, k - 1)*(xw(i, k) - xw(i, k - 1)) + psi(i, k - 1)
      end do
    end do

    ! Calculate the derivatives for the interpolating polynomial
    call cfdotmc(ncol, pver, pverp, xw, psi, fdot)

    !  NEW WAY
    !    calculate fluxes at interior pts
    do k = 2, pver
      do i = 1, ncol
        xxk(i, k) = xw(i, k) - vel(i, k)*deltat
      end do
    end do
    do k = 2, pver
      call cfint2(ncol, pverp, xw, psi, fdot, xxk(1, k), fxdot, fxdd, psistar, errmsg, errflg)
      if(errflg /= 0) then
        return
      endif
      do i = 1, ncol
        flux(i, k) = (psi(i, k) - psistar(i))
      end do
    end do
  end subroutine getflx

  ! Vertical level interpolation
  subroutine cfint2(ncol, pverp, &
    x, f, fdot, xin, fxdot, fxdd, psistar, &
    errmsg, errflg)

    integer,         intent(in)  :: ncol
    integer,         intent(in)  :: pverp
    real(kind_phys), intent(in)  :: x(ncol, pverp)
    real(kind_phys), intent(in)  :: f(ncol, pverp)
    real(kind_phys), intent(in)  :: fdot(ncol, pverp)
    real(kind_phys), intent(in)  :: xin(ncol)

    real(kind_phys), intent(out) :: fxdot(ncol)
    real(kind_phys), intent(out) :: fxdd(ncol)
    real(kind_phys), intent(out) :: psistar(ncol)
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    integer :: i, k
    integer :: intz(ncol)
    real(kind_phys) :: dx
    real(kind_phys) :: s
    real(kind_phys) :: c2
    real(kind_phys) :: c3
    real(kind_phys) :: xx
    real(kind_phys) :: xinf
    real(kind_phys) :: psi1, psi2, psi3, psim
    real(kind_phys) :: cfint
    real(kind_phys) :: cfnew
    real(kind_phys) :: xins(ncol)
    logical :: found_error

    do i = 1, ncol
      xins(i) = medan(x(i, 1), xin(i), x(i, pverp))
      intz(i) = 0
    end do

    ! first find the interval
    do k = 1, pverp - 1
      do i = 1, ncol
        if ((xins(i) - x(i, k))*(x(i, k + 1) - xins(i)) .ge. 0) then
          intz(i) = k
        end if
      end do
    end do

    found_error = .false.
    do i = 1, ncol
      if (intz(i) .eq. 0._kind_phys) found_error = .true.
    end do
    if (found_error) then
      do i = 1, ncol
        if (intz(i) .eq. 0._kind_phys) then
          write(errmsg,'(a,i0)') 'cloud_particle_sedimentation/cfint2: interval was not found for column ', i
          errflg = 1
        end if
      end do
    end if

    ! now interpolate
    do i = 1, ncol
      k = intz(i)
      dx = (x(i, k + 1) - x(i, k))
      s = (f(i, k + 1) - f(i, k))/dx
      c2 = (3*s - 2*fdot(i, k) - fdot(i, k + 1))/dx
      c3 = (fdot(i, k) + fdot(i, k + 1) - 2*s)/dx**2
      xx = (xins(i) - x(i, k))
      fxdot(i) = (3*c3*xx + 2*c2)*xx + fdot(i, k)
      fxdd(i) = 6*c3*xx + 2*c2
      cfint = ((c3*xx + c2)*xx + fdot(i, k))*xx + f(i, k)

      ! limit the interpolant
      psi1 = f(i, k) + (f(i, k + 1) - f(i, k))*xx/dx
      if (k .eq. 1) then
        psi2 = f(i, 1)
      else
        psi2 = f(i, k) + (f(i, k) - f(i, k - 1))*xx/(x(i, k) - x(i, k - 1))
      end if
      if (k + 1 .eq. pverp) then
        psi3 = f(i, pverp)
      else
        psi3 = f(i, k + 1) - (f(i, k + 2) - f(i, k + 1))*(dx - xx)/(x(i, k + 2) - x(i, k + 1))
      end if
      psim = medan(psi1, psi2, psi3)
      cfnew = medan(cfint, psi1, psim)
      !if (abs(cfnew - cfint)/(abs(cfnew) + abs(cfint) + 1.e-36_kind_phys) .gt. .03_kind_phys) then
      !    CHANGE THIS BACK LATER!!!
      !    $      .gt..1) then
      !    UNCOMMENT THIS LATER!!!
      !        write(iulog,*) ' cfint2 limiting important ', cfint, cfnew
      !end if
      psistar(i) = cfnew
    end do
  end subroutine cfint2

  ! Calculate the derivative for the interpolating polynomial, multi-column version
  subroutine cfdotmc(ncol, pver, pverp, x, f, fdot)
    ! assumed variable distribution
    !    x1.......x2.......x3.......x4.......x5.......x6    1,pverp points
    !    f1.......f2.......f3.......f4.......f5.......f6    1,pverp points
    !    ...sh1.......sh2......sh3......sh4......sh5....    1,pver points
    !    .........d2.......d3.......d4.......d5.........    2,pver points
    !    .........s2.......s3.......s4.......s5.........    2,pver points
    !    .............dh2......dh3......dh4.............    2,pver-1 points
    !    .............eh2......eh3......eh4.............    2,pver-1 points
    !    ..................e3.......e4..................    3,pver-1 points
    !    .................ppl3......ppl4................    3,pver-1 points
    !    .................ppr3......ppr4................    3,pver-1 points
    !    .................t3........t4..................    3,pver-1 points
    !    ................fdot3.....fdot4................    3,pver-1 points
    integer,         intent(in)  :: ncol
    integer,         intent(in)  :: pver
    integer,         intent(in)  :: pverp
    real(kind_phys), intent(in)  :: x(ncol, pverp)
    real(kind_phys), intent(in)  :: f(ncol, pverp)

    real(kind_phys), intent(out) :: fdot(ncol, pverp)       ! derivative at nodes

    integer :: i, k
    real(kind_phys) :: a, b, c                 ! work var
    real(kind_phys) :: s(ncol, pverp)         ! first divided differences at nodes
    real(kind_phys) :: sh(ncol, pverp)        ! first divided differences between nodes
    real(kind_phys) :: d(ncol, pverp)         ! second divided differences at nodes
    real(kind_phys) :: dh(ncol, pverp)        ! second divided differences between nodes
    real(kind_phys) :: e(ncol, pverp)         ! third divided differences at nodes
    real(kind_phys) :: eh(ncol, pverp)        ! third divided differences between nodes
    real(kind_phys) :: pp                      ! p prime
    real(kind_phys) :: ppl(ncol, pverp)       ! p prime on left
    real(kind_phys) :: ppr(ncol, pverp)       ! p prime on right
    real(kind_phys) :: qpl
    real(kind_phys) :: qpr
    real(kind_phys) :: ttt
    real(kind_phys) :: t
    real(kind_phys) :: tmin
    real(kind_phys) :: tmax
    real(kind_phys) :: delxh(ncol, pverp)

    do k = 1, pver
      ! first divided differences between nodes
      do i = 1, ncol
        delxh(i, k) = (x(i, k + 1) - x(i, k))
        sh(i, k) = (f(i, k + 1) - f(i, k))/delxh(i, k)
      end do

      ! first and second divided differences at nodes
      if (k .ge. 2) then
        do i = 1, ncol
          d(i, k) = (sh(i, k) - sh(i, k - 1))/(x(i, k + 1) - x(i, k - 1))
          s(i, k) = minmod(sh(i, k), sh(i, k - 1))
        end do
      end if
    end do

    ! second and third divided diffs between nodes
    do k = 2, pver - 1
      do i = 1, ncol
        eh(i, k) = (d(i, k + 1) - d(i, k))/(x(i, k + 2) - x(i, k - 1))
        dh(i, k) = minmod(d(i, k), d(i, k + 1))
      end do
    end do

    ! treat the boundaries
    do i = 1, ncol
      e(i, 2) = eh(i, 2)
      e(i, pver) = eh(i, pver - 1)
      ! outside level
      fdot(i, 1) = sh(i, 1) - d(i, 2)*delxh(i, 1) &
                   - eh(i, 2)*delxh(i, 1)*(x(i, 1) - x(i, 3))
      fdot(i, 1) = minmod(fdot(i, 1), 3*sh(i, 1))
      fdot(i, pverp) = sh(i, pver) + d(i, pver)*delxh(i, pver) &
                       + eh(i, pver - 1)*delxh(i, pver)*(x(i, pverp) - x(i, pver - 1))
      fdot(i, pverp) = minmod(fdot(i, pverp), 3*sh(i, pver))
      ! one in from boundary
      fdot(i, 2) = sh(i, 1) + d(i, 2)*delxh(i, 1) - eh(i, 2)*delxh(i, 1)*delxh(i, 2)
      fdot(i, 2) = minmod(fdot(i, 2), 3*s(i, 2))
      fdot(i, pver) = sh(i, pver) - d(i, pver)*delxh(i, pver) &
                      - eh(i, pver - 1)*delxh(i, pver)*delxh(i, pver - 1)
      fdot(i, pver) = minmod(fdot(i, pver), 3*s(i, pver))
    end do

    do k = 3, pver - 1
      do i = 1, ncol
        e(i, k) = minmod(eh(i, k), eh(i, k - 1))
      end do
    end do

    do k = 3, pver - 1
      do i = 1, ncol
        ! p prime at k-0.5
        ppl(i, k) = sh(i, k - 1) + dh(i, k - 1)*delxh(i, k - 1)
        ! p prime at k+0.5
        ppr(i, k) = sh(i, k) - dh(i, k)*delxh(i, k)

        t = minmod(ppl(i, k), ppr(i, k))

        ! derivate from parabola thru f(i,k-1), f(i,k), and f(i,k+1)
        pp = sh(i, k - 1) + d(i, k)*delxh(i, k - 1)

        ! quartic estimate of fdot
        fdot(i, k) = pp &
                     - delxh(i, k - 1)*delxh(i, k) &
                     *(eh(i, k - 1)*(x(i, k + 2) - x(i, k)) &
                       + eh(i, k)*(x(i, k) - x(i, k - 2)) &
                       )/(x(i, k + 2) - x(i, k - 2))

        ! now limit it
        qpl = sh(i, k - 1) &
              + delxh(i, k - 1)*minmod(d(i, k - 1) + e(i, k - 1)*(x(i, k) - x(i, k - 2)), &
                                       d(i, k) - e(i, k)*delxh(i, k))
        qpr = sh(i, k) &
              + delxh(i, k)*minmod(d(i, k) + e(i, k)*delxh(i, k - 1), &
                                   d(i, k + 1) + e(i, k + 1)*(x(i, k) - x(i, k + 2)))

        fdot(i, k) = medan(fdot(i, k), qpl, qpr)

        ttt = minmod(qpl, qpr)
        tmin = min(0._kind_phys, 3*s(i, k), 1.5_kind_phys*t, ttt)
        tmax = max(0._kind_phys, 3*s(i, k), 1.5_kind_phys*t, ttt)

        fdot(i, k) = fdot(i, k) + minmod(tmin - fdot(i, k), tmax - fdot(i, k))
      end do
    end do
  end subroutine cfdotmc

  pure function minmod(a, b) result(res)
      real(kind_phys), intent(in) :: a, b
      real(kind_phys) :: res
      res = 0.5_kind_phys*(sign(1._kind_phys, a) + sign(1._kind_phys, b))*min(abs(a), abs(b))
  end function minmod

  pure function medan(a, b, c) result(res)
      real(kind_phys), intent(in) :: a, b, c
      real(kind_phys) :: res
      res = a + minmod(b - a, c - a)
  end function medan

end module cloud_particle_sedimentation
