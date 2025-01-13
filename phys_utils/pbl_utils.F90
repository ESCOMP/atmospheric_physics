module pbl_utils

    use ccpp_kinds, only: kind_phys

    implicit none
    private

    public :: calc_rrho
    public :: calc_ustar
    public :: calc_kinematic_heat_flux
    public :: calc_kinematic_water_vapor_flux
    public :: calc_kinematic_buoyancy_flux
    public :: calc_obukhov_length
    public :: calc_virtual_temperature
    public :: austausch_atm
    public :: austausch_atm_free

    real(r8), parameter :: ustar_min = 0.01_r8
    real(r8), parameter :: zkmin     = 0.01_r8 ! Minimum kneutral*f(ri). CCM1 2.f.14

contains

    pure elemental function calc_rrho(rair, t, pmid) result(rrho)
        ! air density reciprocal

        real(r8), intent(in) :: t     ! surface temperature
        real(r8), intent(in) :: pmid  ! midpoint pressure (bottom level)
        real(r8), intent(in) :: rair  ! gas constant for dry air

        real(r8) :: rrho              ! 1./bottom level density

        rrho = rair * t / pmid
    end function calc_rrho

    pure elemental function calc_ustar(taux, tauy, rrho) result(ustar)
        ! https://glossary.ametsoc.org/wiki/Friction_velocity
        ! NOTE: taux,tauy come form the expansion of the Reynolds stress
        !
        ! Also found in:
        ! Stull, Roland B. An Introduction to Boundary Layer Meteorology. Springer Kluwer Academic Publishers, 1988. Print.
        ! DOI: https://doi.org/10.1007/978-94-009-3027-8
        ! Equation 2.10b, page 181

        real(r8), intent(in) :: taux ! surface u stress [N/m2]
        real(r8), intent(in) :: tauy ! surface v stress [N/m2]
        real(r8), intent(in) :: rrho ! 1./bottom level density

        real(r8) :: ustar            ! surface friction velocity [m/s]

        ustar = max( sqrt( sqrt(taux**2 + tauy**2)*rrho ), ustar_min )
    end function calc_ustar

    pure elemental function calc_kinematic_heat_flux(shflx, rrho, cpair) result(khfs)
        real(r8), intent(in)  :: shflx  ! surface heat flux (W/m2)
        real(r8), intent(in)  :: rrho   ! 1./bottom level density [ m3/kg ]
        real(r8), intent(in)  :: cpair  ! specific heat of dry air

        real(r8) :: khfs                ! sfc kinematic heat flux [mK/s]

        khfs = shflx*rrho/cpair
    end function calc_kinematic_heat_flux

    pure elemental function calc_kinematic_water_vapor_flux(qflx, rrho) result(kqfs)
        real(r8), intent(in)  :: qflx ! water vapor flux (kg/m2/s)
        real(r8), intent(in)  :: rrho ! 1./bottom level density [ m3/kg ]

        real(r8) :: kqfs              ! sfc kinematic water vapor flux [m/s]

        kqfs = qflx*rrho
    end function calc_kinematic_water_vapor_flux

    pure elemental function calc_kinematic_buoyancy_flux(khfs, zvir, ths, kqfs) result(kbfs)
        real(r8), intent(in) :: khfs  ! sfc kinematic heat flux [mK/s]
        real(r8), intent(in) :: zvir  ! rh2o/rair - 1
        real(r8), intent(in) :: ths   ! potential temperature at surface [K]
        real(r8), intent(in) :: kqfs  ! sfc kinematic water vapor flux [m/s]

        real(r8) :: kbfs              ! sfc kinematic buoyancy flux [m^2/s^3]

        kbfs = khfs + zvir*ths*kqfs
    end function calc_kinematic_buoyancy_flux

    pure elemental function calc_obukhov_length(thvs, ustar, g, vk, kbfs) result(obukhov_length)
        ! Stull, Roland B. An Introduction to Boundary Layer Meteorology. Springer Kluwer Academic Publishers, 1988. Print.
        ! DOI: https://doi.org/10.1007/978-94-009-3027-8
        ! Equation 5.7c, page 181

        real(r8), intent(in)  :: thvs              ! virtual potential temperature at surface
        real(r8), intent(in)  :: ustar             ! Surface friction velocity [ m/s ]
        real(r8), intent(in)  :: g                 ! acceleration of gravity
        real(r8), intent(in)  :: vk                ! Von Karman's constant
        real(r8), intent(in)  :: kbfs              ! sfc kinematic buoyancy flux [m*K/s]

        real(r8)              :: obukhov_length(n) ! Obukhov length

        ! Added sign(...) term to prevent division by 0 and using the fact that
        ! `kbfs = \overline{(w' \theta_v')}_s`
        obukhov_length = -thvs * ustar**3 / (g*vk*(kbfs + sign(1.e-10_r8,kbfs)))
    end function calc_obukhov_length

    pure elemental function calc_virtual_temperature(t, q, zvir) result(virtem)
        ! Williamson, D., Kiehl, J., Ramanathan, V., Dickinson, R., & Hack, J. (1987).
        ! Description of the NCAR Community Climate Model (CCM1).
        ! University Corporation for Atmospheric Research. https://doi.org/10.5065/D6TB14WH (Original work published 1987)
        ! Equation 2.a.7

        real(r8), intent(in) :: t, q    ! temperature and specific humidity
        real(r8), intent(in) :: zvir    ! rh2o/rair - 1

        real(r8)             :: virtual_temperature

        virtual_temperature = t * (1.0_r8 + zvir*q)
    end function calc_virtual_temperature

    pure function austausch_atm(pcols, ncol, pver, ntop, nbot, ml2, ri, s2) result(kvf)
        !---------------------------------------------------------------------- !
        !                                                                       !
        ! Purpose: Computes exchange coefficients for free turbulent flows.     !
        !                                                                       !
        ! Method:                                                               !
        !                                                                       !
        ! The free atmosphere diffusivities are based on standard mixing length !
        ! forms for the neutral diffusivity multiplied by functns of Richardson !
        ! number. K = l^2 * |dV/dz| * f(Ri). The same functions are used for    !
        ! momentum, potential temperature, and constitutents.                   !
        !                                                                       !
        ! The stable Richardson num function (Ri>0) is taken from Holtslag and  !
        ! Beljaars (1989), ECMWF proceedings. f = 1 / (1 + 10*Ri*(1 + 8*Ri))    !
        ! The unstable Richardson number function (Ri<0) is taken from  CCM1.   !
        ! f = sqrt(1 - 18*Ri)                                                   !
        !                                                                       !
        ! Author: B. Stevens (rewrite, August 2000)                             !
        !                                                                       !
        !---------------------------------------------------------------------- !

        integer,  intent(in)  :: pcols              ! Atmospheric columns dimension size
        integer,  intent(in)  :: ncol               ! Number of atmospheric columns
        integer,  intent(in)  :: pver               ! Number of atmospheric layers
        integer,  intent(in)  :: ntop               ! Top layer for calculation
        integer,  intent(in)  :: nbot               ! Bottom layer for calculation
        real(r8), intent(in)  :: ml2(pver+1)        ! Mixing lengths squared
        real(r8), intent(in)  :: s2(pcols,pver)     ! Shear squared
        real(r8), intent(in)  :: ri(pcols,pver)     ! Richardson number

        real(r8)              :: kvf(pcols,pver+1)  ! Eddy diffusivity for heat and tracers

        real(r8)              :: fofri              ! f(ri)
        real(r8)              :: kvn                ! Neutral Kv
        integer               :: i                  ! Longitude index
        integer               :: k                  ! Vertical index

        kvf(:ncol,:) = 0.0_r8

        ! Compute the atmosphere free vertical diffusion coefficients: kvh = kvq = kvm.
        do k = ntop, nbot - 1
            do i = 1, ncol
                if( ri(i,k) < 0.0_r8 ) then
                    fofri = unstable_gradient_richardson_stability_parameter(ri(i,k))
                else
                    fofri = stable_gradient_richardson_stability_parameter(ri(i,k))
                end if
                kvn = neutral_exchange_coefficient(ml2(k), s2(i,k))
                kvf(i,k+1) = max( zkmin, kvn * fofri )
            end do
        end do
        
    end function austausch_atm
        
    pure function austausch_atm_free(pcols, ncol, pver, ntop, nbot, ml2, ri, s2) result(kvf)
        !---------------------------------------------------------------------- !
        !                                                                       !
        ! same as austausch_atm but only mixing for Ri<0                        !
        ! i.e. no background mixing and mixing for Ri>0                         !
        !                                                                       !
        !---------------------------------------------------------------------- !
        
        integer,  intent(in)  :: pcols              ! Atmospheric columns dimension size
        integer,  intent(in)  :: ncol               ! Number of atmospheric columns
        integer,  intent(in)  :: pver               ! Number of atmospheric layers
        integer,  intent(in)  :: ntop               ! Top layer for calculation
        integer,  intent(in)  :: nbot               ! Bottom layer for calculation
        real(r8), intent(in)  :: ml2(pver+1)        ! Mixing lengths squared
        real(r8), intent(in)  :: s2(pcols,pver)     ! Shear squared
        real(r8), intent(in)  :: ri(pcols,pver)     ! Richardson no

        real(r8)              :: kvf(pcols,pver+1)  ! Eddy diffusivity for heat and tracers

        real(r8)              :: fofri              ! f(ri)
        real(r8)              :: kvn                ! Neutral Kv
        integer               :: i                  ! Longitude index
        integer               :: k                  ! Vertical index

        kvf(:ncol,:) = 0.0_r8

        ! Compute the free atmosphere vertical diffusion coefficients: kvh = kvq = kvm.
        do k = ntop, nbot - 1
            do i = 1, ncol
                if( ri(i,k) < 0.0_r8 ) then
                    fofri = unstable_gradient_richardson_stability_parameter(ri(i,k))
                    kvn = neutral_exchange_coefficient(ml2(k), s2(i,k))
                    kvf(i,k+1) = kvn * fofri
                end if
            end do
        end do
    end function austausch_atm_free

    pure elemental function unstable_gradient_richardson_stability_parameter(richardson_number) result(modifier)
        ! Williamson, D., Kiehl, J., Ramanathan, V., Dickinson, R., & Hack, J. (1987).
        ! Description of the NCAR Community Climate Model (CCM1).
        ! University Corporation for Atmospheric Research. https://doi.org/10.5065/D6TB14WH (Original work published 1987)
        ! Equation 2.f.13

        real(r8), intent(in)  :: richardson_number

        real(r8)              :: modifier

        modifier = sqrt( 1._r8 - 18._r8 * richardson_number )
    end function unstable_gradient_richardson_stability_parameter

    pure elemental function stable_gradient_richardson_stability_parameter(richardson_number) result(modifier)
        ! ECMWF 74888 Surface flux parameterization schemes: developments and experiences at KNMI (eq. 20)
        ! Originally used published equation from in CCM1 2.f.12

        real(r8), intent(in)  :: richardson_number

        real(r8)              :: modifier

        modifier = 1.0_r8 / ( 1.0_r8 + 10.0_r8 * richardson_number * ( 1.0_r8 + 8.0_r8 * richardson_number ) )
    end function stable_gradient_richardson_stability_parameter

    pure elemental function neutral_exchange_coefficient(mixing_length, shear_squared) result(neutral_k)
        ! Williamson, D., Kiehl, J., Ramanathan, V., Dickinson, R., & Hack, J. (1987).
        ! Description of the NCAR Community Climate Model (CCM1).
        ! University Corporation for Atmospheric Research. https://doi.org/10.5065/D6TB14WH (Original work published 1987)
        ! Equation 2.f.15

        real(r8), intent(in) :: mixing_length
        real(r8), intent(in) :: shear_squared

        real(r8)             :: neutral_k

        neutral_k = mixing_length * sqrt(shear_squared)
    end function neutral_exchange_coefficient
end module pbl_utils