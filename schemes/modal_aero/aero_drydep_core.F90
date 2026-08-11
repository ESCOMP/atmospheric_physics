!===============================================================================
! Aerosol dry deposition
!   Portable science routines split from modal_aero/aero_model.F90 and
!   aer_drydep_mod.F90: surface deposition velocities of particles
!   (Zhang et al. 2001) and the aerodynamic resistance / friction velocity
!   patch over ocean and sea ice. Host constants and the landuse fractions
!   are passed as arguments; array sizing is by ncol/pver runtime arguments.
!===============================================================================
module aero_drydep_core

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: modal_aero_depvel_part
  public :: calcram

contains

  !=============================================================================
  !=============================================================================
  subroutine modal_aero_depvel_part(ncol, t, pmid, ram1, fv, vlc_dry, vlc_trb, vlc_grv, &
                                    radius_part, density_part, sig_part, moment, &
                                    pver, top_lev, n_land_type, fraction_landuse, &
                                    pi, boltz, gravit, rair, aspherical)   ! dmleung added aspherical flag 20 Oct 2025

!    calculates surface deposition velocity of particles
!    L. Zhang, S. Gong, J. Padro, and L. Barrie
!    A size-seggregated particle dry deposition scheme for an atmospheric aerosol module
!    Atmospheric Environment, 35, 549-560, 2001.
!
!    Authors: X. Liu

    ! !ARGUMENTS:
    !
    !
    real(kind_phys), intent(in) :: t(:, :)       !atm temperature (K)
    real(kind_phys), intent(in) :: pmid(:, :)    !atm pressure (Pa)
    real(kind_phys), intent(in) :: fv(:)           !friction velocity (m/s)
    real(kind_phys), intent(in) :: ram1(:)         !aerodynamical resistance (s/m)
    real(kind_phys), intent(in) :: radius_part(:, :)    ! mean (volume/number) particle radius (m)
    real(kind_phys), intent(in) :: density_part(:, :)   ! density of particle material (kg/m3)
    real(kind_phys), intent(in) :: sig_part(:, :)       ! geometric standard deviation of particles
    integer, intent(in) :: moment ! moment of size distribution (0 for number, 2 for surface area, 3 for volume)
    integer, intent(in) :: ncol
    integer, intent(in) :: pver                ! number of vertical levels
    integer, intent(in) :: top_lev             ! top level for modal aerosols
    integer, intent(in) :: n_land_type         ! number of land use types
    real(kind_phys), intent(in) :: fraction_landuse(:, :)  ! land use fractions (ncol, n_land_type)
    real(kind_phys), intent(in) :: pi                  ! host model constants
    real(kind_phys), intent(in) :: boltz               ! Boltzmann constant (J/K)
    real(kind_phys), intent(in) :: gravit              ! gravitational acceleration (m/s2)
    real(kind_phys), intent(in) :: rair                ! gas constant for dry air (J/K/kg)

    real(kind_phys), intent(out) :: vlc_trb(:)       !Turbulent deposn velocity (m/s)
    real(kind_phys), intent(out) :: vlc_grv(:, :)       !grav deposn velocity (m/s)
    real(kind_phys), intent(out) :: vlc_dry(:, :)       !dry deposn velocity (m/s)
    logical, intent(in), optional :: aspherical   ! dmleung: asphericity is strong for coarse-mode interstitial
    ! aerosols only, mostly dust and seasalt. For coarse mode aerosols, asphericity reduces coarse-mode gravitational
    ! settling velocity by 20 % following Fig. 4 of Yue Huang et al. (2020).
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! Local Variables
    integer  :: m, i, k, ix                !indices
    real(kind_phys) :: rho     !atm density (kg/m**3)
    real(kind_phys) :: vsc_dyn_atm(ncol, pver)   ![kg m-1 s-1] Dynamic viscosity of air
    real(kind_phys) :: vsc_knm_atm(ncol, pver)   ![m2 s-1] Kinematic viscosity of atmosphere
    real(kind_phys) :: shm_nbr       ![frc] Schmidt number
    real(kind_phys) :: stk_nbr       ![frc] Stokes number
    real(kind_phys) :: mfp_atm(ncol, pver)       ![m] Mean free path of air
    real(kind_phys) :: dff_aer       ![m2 s-1] Brownian diffusivity of particle
    real(kind_phys) :: slp_crc(ncol, pver) ![frc] Slip correction factor
    real(kind_phys) :: rss_trb       ![s m-1] Resistance to turbulent deposition
    real(kind_phys) :: rss_lmn       ![s m-1] Quasi-laminar layer resistance
    real(kind_phys) :: brownian      ! collection efficiency for Browning diffusion
    real(kind_phys) :: impaction     ! collection efficiency for impaction
    real(kind_phys) :: interception  ! collection efficiency for interception
    real(kind_phys) :: stickfrac     ! fraction of particles sticking to surface
    real(kind_phys) :: radius_moment(ncol, pver) ! median radius (m) for moment
    real(kind_phys) :: lnsig         ! ln(sig_part)
    real(kind_phys) :: dispersion    ! accounts for influence of size dist dispersion on bulk settling velocity
    ! assuming radius_part is number mode radius * exp(1.5 ln(sigma))

    integer  :: lt
    real(kind_phys) :: lnd_frc
    real(kind_phys) :: wrk1, wrk2, wrk3

    ! constants

    real(kind_phys), parameter :: asphericaldust_drydep = 0.8_kind_phys ! dmleung added 20 Oct 2025: aspherical dust reduces
    ! gravitational settling velocity by 15-20 %. Yue Huang et al. (2020)
    ! Climate Models and Remote Sensing Retrievals Neglect Substantial Desert Dust Asphericity

    real(kind_phys) :: gamma(11)      ! exponent of schmidt number
!   data gamma/0.54d+00,  0.56d+00,  0.57d+00,  0.54d+00,  0.54d+00, &
!              0.56d+00,  0.54d+00,  0.54d+00,  0.54d+00,  0.56d+00, &
!              0.50d+00/
    data gamma/0.56e+00_kind_phys, 0.54e+00_kind_phys, 0.54e+00_kind_phys, 0.56e+00_kind_phys, 0.56e+00_kind_phys, &
      0.56e+00_kind_phys, 0.50e+00_kind_phys, 0.54e+00_kind_phys, 0.54e+00_kind_phys, 0.54e+00_kind_phys, &
      0.54e+00_kind_phys/
    save gamma

    real(kind_phys) :: alpha(11)      ! parameter for impaction
!   data alpha/50.00d+00,  0.95d+00,  0.80d+00,  1.20d+00,  1.30d+00, &
!               0.80d+00, 50.00d+00, 50.00d+00,  2.00d+00,  1.50d+00, &
!             100.00d+00/
    data alpha/1.50e+00_kind_phys, 1.20e+00_kind_phys, 1.20e+00_kind_phys, 0.80e+00_kind_phys, 1.00e+00_kind_phys, &
      0.80e+00_kind_phys, 100.00e+00_kind_phys, 50.00e+00_kind_phys, 2.00e+00_kind_phys, 1.20e+00_kind_phys, &
      50.00e+00_kind_phys/
    save alpha

    real(kind_phys) :: radius_collector(11) ! radius (m) of surface collectors
!   data radius_collector/-1.00d+00,  5.10d-03,  3.50d-03,  3.20d-03, 10.00d-03, &
!                          5.00d-03, -1.00d+00, -1.00d+00, 10.00d-03, 10.00d-03, &
!                         -1.00d+00/
    data radius_collector/10.00e-03_kind_phys, 3.50e-03_kind_phys, 3.50e-03_kind_phys, 5.10e-03_kind_phys, 2.00e-03_kind_phys, &
      5.00e-03_kind_phys, -1.00e+00_kind_phys, -1.00e+00_kind_phys, 10.00e-03_kind_phys, 3.50e-03_kind_phys, &
      -1.00e+00_kind_phys/
    save radius_collector

    integer            :: iwet(11) ! flag for wet surface = 1, otherwise = -1
!   data iwet/1,   -1,   -1,   -1,   -1,  &
!            -1,   -1,   -1,    1,   -1,  &
!             1/
    data iwet/-1, -1, -1, -1, -1, &
      -1, 1, -1, 1, -1, &
      -1/
    save iwet

    vlc_trb = 0._kind_phys
    vlc_grv = 0._kind_phys
    vlc_dry = 0._kind_phys

    !------------------------------------------------------------------------
    do k = top_lev, pver ! radius_part is not defined above top_lev
      do i = 1, ncol

        lnsig = log(sig_part(i, k))
! use a maximum radius of 50 microns when calculating deposition velocity
        radius_moment(i, k) = min(50.0e-6_kind_phys, radius_part(i, k))* &
                              exp((float(moment) - 1.5_kind_phys)*lnsig*lnsig)
        dispersion = exp(2._kind_phys*lnsig*lnsig)

        rho = pmid(i, k)/rair/t(i, k)

        ! Quasi-laminar layer resistance: call rss_lmn_get
        ! Size-independent thermokinetic properties
        vsc_dyn_atm(i, k) = 1.72e-5_kind_phys*((t(i, k)/273.0_kind_phys)**1.5_kind_phys)*393.0_kind_phys/ &
                            (t(i, k) + 120.0_kind_phys)      ![kg m-1 s-1] RoY94 p. 102
        mfp_atm(i, k) = 2.0_kind_phys*vsc_dyn_atm(i, k)/ &   ![m] SeP97 p. 455
                        (pmid(i, k)*sqrt(8.0_kind_phys/(pi*rair*t(i, k))))
        vsc_knm_atm(i, k) = vsc_dyn_atm(i, k)/rho ![m2 s-1] Kinematic viscosity of air

        slp_crc(i, k) = 1.0_kind_phys + mfp_atm(i, k)* &
                        (1.257_kind_phys + 0.4_kind_phys*exp(-1.1_kind_phys*radius_moment(i, k)/(mfp_atm(i, k))))/ &
                        radius_moment(i, k)   ![frc] Slip correction factor SeP97 p. 464
        vlc_grv(i, k) = (4.0_kind_phys/18.0_kind_phys)*radius_moment(i, k)*radius_moment(i, k)*density_part(i, k)* &
                        gravit*slp_crc(i, k)/vsc_dyn_atm(i, k) ![m s-1] Stokes' settling velocity SeP97 p. 466
        vlc_grv(i, k) = vlc_grv(i, k)*dispersion

        ! dmleung edited 20 Oct 2025 based on Longlei Li's edits ++
        ! asphericity reduces gravitational settling velocity of coarse-mode aerosols by 20 %.
        ! scale flag is only true for coarse mode (m == n_coarse_dust).
        if (present(aspherical)) then
          if (aspherical) then
            vlc_grv(i, k) = vlc_grv(i, k)*asphericaldust_drydep
          end if
        end if
        ! dmleung --

        vlc_dry(i, k) = vlc_grv(i, k)
      end do
    end do
    k = pver  ! only look at bottom level for next part
    do i = 1, ncol
      dff_aer = boltz*t(i, k)*slp_crc(i, k)/ &    ![m2 s-1]
                (6.0_kind_phys*pi*vsc_dyn_atm(i, k)*radius_moment(i, k)) !SeP97 p.474
      shm_nbr = vsc_knm_atm(i, k)/dff_aer                        ![frc] SeP97 p.972

      wrk2 = 0._kind_phys
      wrk3 = 0._kind_phys
      do lt = 1, n_land_type
        lnd_frc = fraction_landuse(i, lt)
        if (lnd_frc /= 0._kind_phys) then
          brownian = shm_nbr**(-gamma(lt))
          if (radius_collector(lt) > 0.0_kind_phys) then
!       vegetated surface
            stk_nbr = vlc_grv(i, k)*fv(i)/(gravit*radius_collector(lt))
            interception = 2.0_kind_phys*(radius_moment(i, k)/radius_collector(lt))**2.0_kind_phys
          else
!       non-vegetated surface
            stk_nbr = vlc_grv(i, k)*fv(i)*fv(i)/(gravit*vsc_knm_atm(i, k))  ![frc] SeP97 p.965
            interception = 0.0_kind_phys
          end if
          impaction = (stk_nbr/(alpha(lt) + stk_nbr))**2.0_kind_phys

          if (iwet(lt) > 0) then
            stickfrac = 1.0_kind_phys
          else
            stickfrac = exp(-sqrt(stk_nbr))
            if (stickfrac < 1.0e-10_kind_phys) stickfrac = 1.0e-10_kind_phys
          end if
          rss_lmn = 1.0_kind_phys/(3.0_kind_phys*fv(i)*stickfrac*(brownian + interception + impaction))
          rss_trb = ram1(i) + rss_lmn + ram1(i)*rss_lmn*vlc_grv(i, k)

          wrk1 = 1.0_kind_phys/rss_trb
          wrk2 = wrk2 + lnd_frc*(wrk1)
          wrk3 = wrk3 + lnd_frc*(wrk1 + vlc_grv(i, k))
        end if
      end do  ! n_land_type
      vlc_trb(i) = wrk2
      vlc_dry(i, k) = wrk3
    end do !ncol

  end subroutine modal_aero_depvel_part

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subroutine Calcram
!
! !INTERFACE:
!

  subroutine calcram(ncol, landfrac, icefrac, ocnfrac, obklen, &
                     ustar, ram1in, ram1, t, pmid, &
                     pdel, fvin, fv, rair, gravit)
    !
    ! !DESCRIPTION:
    !
    ! Calc aerodynamic resistance over oceans and sea ice (comes in from land model)
    ! from Seinfeld and Pandis, p.963.
    !
    ! Author: Natalie Mahowald
    !
    integer, intent(in) :: ncol
    real(kind_phys), intent(in) :: ram1in(:)         !aerodynamical resistance (s/m)
    real(kind_phys), intent(in) :: fvin(:)                 ! sfc frc vel from land
    real(kind_phys), intent(out) :: ram1(:)         !aerodynamical resistance (s/m)
    real(kind_phys), intent(out) :: fv(:)                 ! sfc frc vel from land
    real(kind_phys), intent(in) :: obklen(:)                 ! obklen
    real(kind_phys), intent(in) :: ustar(:)                  ! sfc fric vel
    real(kind_phys), intent(in) :: landfrac(:)               ! land fraction
    real(kind_phys), intent(in) :: icefrac(:)                ! ice fraction
    real(kind_phys), intent(in) :: ocnfrac(:)                ! ocean fraction
    real(kind_phys), intent(in) :: t(:)       !atm temperature (K)
    real(kind_phys), intent(in) :: pmid(:)    !atm pressure (Pa)
    real(kind_phys), intent(in) :: pdel(:)    !atm pressure (Pa)
    real(kind_phys), intent(in) :: rair       ! gas constant for dry air (J/K/kg)
    real(kind_phys), intent(in) :: gravit     ! gravitational acceleration (m/s2)
    real(kind_phys), parameter :: zzocen = 0.0001_kind_phys   ! Ocean aerodynamic roughness length
    real(kind_phys), parameter :: zzsice = 0.0400_kind_phys   ! Sea ice aerodynamic roughness length
    real(kind_phys), parameter :: xkar = 0.4_kind_phys      ! Von Karman constant

    ! local variables
    real(kind_phys) :: z, psi, psi0, nu, nu0, temp, ram
    integer :: i
    !    write(iulog,*) rair,zzsice,zzocen,gravit,xkar

    do i = 1, ncol
      z = pdel(i)*rair*t(i)/pmid(i)/gravit/2.0_kind_phys   !use half the layer height like Ganzefeld and Lelieveld, 1995
      if (obklen(i) == 0) then
        psi = 0._kind_phys
        psi0 = 0._kind_phys
      else
        psi = min(max(z/obklen(i), -1.0_kind_phys), 1.0_kind_phys)
        psi0 = min(max(zzocen/obklen(i), -1.0_kind_phys), 1.0_kind_phys)
      end if
      temp = z/zzocen
      if (icefrac(i) > 0.5_kind_phys) then
        if (obklen(i) > 0) then
          psi0 = min(max(zzsice/obklen(i), -1.0_kind_phys), 1.0_kind_phys)
        else
          psi0 = 0.0_kind_phys
        end if
        temp = z/zzsice
      end if
      if (psi > 0._kind_phys) then
        ram = 1/xkar/ustar(i)*(log(temp) + 4.7_kind_phys*(psi - psi0))
      else
        nu = (1.00_kind_phys - 15.000_kind_phys*psi)**(.25_kind_phys)
        nu0 = (1.000_kind_phys - 15.000_kind_phys*psi0)**(.25_kind_phys)
        if (ustar(i) /= 0._kind_phys) then
          ram = 1/xkar/ustar(i)*(log(temp) &
                                 + log(((nu0**2 + 1.00_kind_phys)*(nu0 + 1.0_kind_phys)**2)/((nu**2 + 1.0_kind_phys)*(nu + 1.00_kind_phys)**2)) &
                                 + 2.0_kind_phys*(atan(nu) - atan(nu0)))
        else
          ram = 0._kind_phys
        end if
      end if
      if (landfrac(i) < 0.000000001_kind_phys) then
        fv(i) = ustar(i)
        ram1(i) = ram
      else
        fv(i) = fvin(i)
        ram1(i) = ram1in(i)
      end if
      !          write(iulog,*) i,pdel(i),t(i),pmid(i),gravit,obklen(i),psi,psi0,icefrac(i),nu,nu0,ram,ustar(i),&
      !             log(((nu0**2+1.00)*(nu0+1.0)**2)/((nu**2+1.0)*(nu+1.00)**2)),2.0*(atan(nu)-atan(nu0))

    end do

    ! fvitt -- fv == 0 causes a floating point exception in
    ! dry dep of sea salts and dust
    where (fv(:ncol) == 0._kind_phys)
      fv(:ncol) = 1.e-12_kind_phys
    end where

  end subroutine calcram

end module aero_drydep_core
