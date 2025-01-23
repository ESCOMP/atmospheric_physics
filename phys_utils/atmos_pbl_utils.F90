module atmos_phys_pbl_utils

    use ccpp_kinds, only: kind_phys

    implicit none
    private

    public :: calc_rrho
    public :: calc_friction_velocity
    public :: calc_kinematic_heat_flux
    public :: calc_kinematic_water_vapor_flux
    public :: calc_kinematic_buoyancy_flux
    public :: calc_obukhov_length
    public :: calc_virtual_temperature
    public :: austausch_atm
    public :: austausch_atm_free

    real(kind_phys), parameter :: ustar_min = 0.01_kind_phys
    real(kind_phys), parameter :: zkmin     = 0.01_kind_phys ! Minimum kneutral*f(ri). CCM1 2.f.14

contains

    pure elemental function calc_rrho(rair, surface_temperature, pmid) result(rrho)
        ! air density reciprocal
        real(kind_phys), intent(in) :: rair  ! gas constant for dry air
        real(kind_phys), intent(in) :: surface_temperature
        real(kind_phys), intent(in) :: pmid  ! midpoint pressure (bottom level)

        real(kind_phys)             :: rrho  ! 1./bottom level density

        rrho = rair * surface_temperature / pmid
    end function calc_rrho

    pure elemental function calc_friction_velocity(taux, tauy, rrho) result(friction_velocity)
        ! https://glossary.ametsoc.org/wiki/Friction_velocity
        ! NOTE: taux,tauy come form the expansion of the Reynolds stress
        !
        ! Also found in:
        ! Stull, Roland B. An Introduction to Boundary Layer Meteorology. Springer Kluwer Academic Publishers, 1988. Print.
        ! DOI: https://doi.org/10.1007/978-94-009-3027-8
        ! Equation 2.10b, page 67

        real(kind_phys), intent(in) :: taux              ! surface u stress [N/m2]
        real(kind_phys), intent(in) :: tauy              ! surface v stress [N/m2]
        real(kind_phys), intent(in) :: rrho              ! 1./bottom level density

        real(kind_phys)             :: friction_velocity ! surface friction velocity [m/s]

        friction_velocity = max( sqrt( sqrt(taux**2 + tauy**2)*rrho ), ustar_min )
    end function calc_friction_velocity

    pure elemental function calc_kinematic_heat_flux(shflx, rrho, cpair) result(khfs)
        real(kind_phys), intent(in)  :: shflx  ! surface heat flux (W/m2)
        real(kind_phys), intent(in)  :: rrho   ! 1./bottom level density [ m3/kg ]
        real(kind_phys), intent(in)  :: cpair  ! specific heat of dry air

        real(kind_phys)              :: khfs   ! sfc kinematic heat flux [mK/s]

        khfs = shflx*rrho/cpair
    end function calc_kinematic_heat_flux

    pure elemental function calc_kinematic_water_vapor_flux(qflx, rrho) result(kqfs)
        real(kind_phys), intent(in)  :: qflx ! water vapor flux (kg/m2/s)
        real(kind_phys), intent(in)  :: rrho ! 1./bottom level density [ m3/kg ]

        real(kind_phys)              :: kqfs ! sfc kinematic water vapor flux [m/s]

        kqfs = qflx*rrho
    end function calc_kinematic_water_vapor_flux

    pure elemental function calc_kinematic_buoyancy_flux(khfs, zvir, ths, kqfs) result(kbfs)
        real(kind_phys), intent(in) :: khfs  ! sfc kinematic heat flux [mK/s]
        real(kind_phys), intent(in) :: zvir  ! rh2o/rair - 1
        real(kind_phys), intent(in) :: ths   ! potential temperature at surface [K]
        real(kind_phys), intent(in) :: kqfs  ! sfc kinematic water vapor flux [m/s]

        real(kind_phys)             :: kbfs  ! sfc kinematic buoyancy flux [mK/s] (`kbfs = \overline{(w' \theta'_v)}_s`)

        kbfs = khfs + zvir*ths*kqfs
    end function calc_kinematic_buoyancy_flux

    pure elemental function calc_obukhov_length(thvs, ustar, g, vk, kbfs) result(obukhov_length)
        ! Stull, Roland B. An Introduction to Boundary Layer Meteorology. Springer Kluwer Academic Publishers, 1988. Print.
        ! DOI: https://doi.org/10.1007/978-94-009-3027-8
        ! Equation 5.7c, page 181

        real(kind_phys), intent(in)  :: thvs              ! virtual potential temperature at surface
        real(kind_phys), intent(in)  :: ustar             ! Surface friction velocity [ m/s ]
        real(kind_phys), intent(in)  :: g                 ! acceleration of gravity
        real(kind_phys), intent(in)  :: vk                ! Von Karman's constant
        real(kind_phys), intent(in)  :: kbfs              ! sfc kinematic buoyancy flux [m*K/s]

        real(kind_phys)              :: obukhov_length    ! Obukhov length

        ! Added sign(...) term to prevent division by 0 and using the fact that
        ! `kbfs = \overline{(w' \theta_v')}_s`
        obukhov_length = -thvs * ustar**3 /                          &
                         (g*vk*(kbfs + sign(1.e-10_kind_phys,kbfs)))
    end function calc_obukhov_length

    pure elemental function calc_virtual_temperature(temperature, specific_humidity, zvir) result(virtual_temperature)
        ! Williamson, D., Kiehl, J., Ramanathan, V., Dickinson, R., & Hack, J. (1987).
        ! Description of the NCAR Community Climate Model (CCM1).
        ! University Corporation for Atmospheric Research. https://doi.org/10.5065/D6TB14WH (Original work published 1987)
        ! Equation 2.a.7

        real(kind_phys), intent(in) :: temperature
        real(kind_phys), intent(in) :: specific_humidity
        real(kind_phys), intent(in) :: zvir    ! rh2o/rair - 1

        real(kind_phys)             :: virtual_temperature

        virtual_temperature = temperature * (1.0_kind_phys + zvir*specific_humidity)
    end function calc_virtual_temperature

    pure elemental function austausch_atm(mixing_length_squared, &
                                          richardson_number,     &
                                          shear_squared)         &
                                          result(kvf)
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

        real(kind_phys), intent(in)  :: mixing_length_squared
        real(kind_phys), intent(in)  :: richardson_number
        real(kind_phys), intent(in)  :: shear_squared

        real(kind_phys)              :: kvf    ! Eddy diffusivity for heat and tracers

        real(kind_phys)              :: fofri  ! f(ri)
        real(kind_phys)              :: kvn    ! Neutral Kv

        if( richardson_number < 0.0_kind_phys ) then
            fofri = unstable_gradient_richardson_stability_parameter(richardson_number)
        else
            fofri = stable_gradient_richardson_stability_parameter(richardson_number)
        end if
        kvn = neutral_exchange_coefficient(mixing_length_squared, shear_squared)
        kvf = max( zkmin, kvn * fofri )
    end function austausch_atm

    pure elemental function austausch_atm_free(mixing_length_squared, &
                                               richardson_number,     &
                                               shear_squared)         &
                                               result(kvf)
        !---------------------------------------------------------------------- !
        !                                                                       !
        ! same as austausch_atm but only mixing for Ri<0                        !
        ! i.e. no background mixing and mixing for Ri>0                         !
        !                                                                       !
        !---------------------------------------------------------------------- !
        real(kind_phys), intent(in)  :: mixing_length_squared
        real(kind_phys), intent(in)  :: richardson_number
        real(kind_phys), intent(in)  :: shear_squared

        real(kind_phys)              :: kvf

        real(kind_phys)              :: fofri ! f(ri)
        real(kind_phys)              :: kvn   ! Neutral Kv

        kvf = 0.0_kind_phys
        if( richardson_number < 0.0_kind_phys ) then
            fofri = unstable_gradient_richardson_stability_parameter(richardson_number)
            kvn = neutral_exchange_coefficient(mixing_length_squared, shear_squared)
            kvf = kvn * fofri
        end if
    end function austausch_atm_free

    pure elemental function unstable_gradient_richardson_stability_parameter(richardson_number) result(modifier)
        ! Williamson, D., Kiehl, J., Ramanathan, V., Dickinson, R., & Hack, J. (1987).
        ! Description of the NCAR Community Climate Model (CCM1).
        ! University Corporation for Atmospheric Research. https://doi.org/10.5065/D6TB14WH (Original work published 1987)
        ! Equation 2.f.13

        real(kind_phys), intent(in)  :: richardson_number

        real(kind_phys)              :: modifier

        modifier = sqrt( 1._kind_phys - 18._kind_phys * richardson_number )
    end function unstable_gradient_richardson_stability_parameter

    pure elemental function stable_gradient_richardson_stability_parameter(richardson_number) result(modifier)
        ! Holtslag, A. A. M., and Beljaars A. C. M. , 1989: Surface flux parameterization schemes: Developments and experiences at KNMI.
        ! ECMWF Workshop on Parameterization of Fluxes and Land Surface, Reading, United Kingdom, ECMWF, 121–147.
        ! equation 20, page 140
        ! Originally used published equation from CCM1, 2.f.12, page 11

        real(kind_phys), intent(in)  :: richardson_number

        real(kind_phys)              :: modifier

        modifier = 1.0_kind_phys /                                                                                              &
                 ( 1.0_kind_phys + 10.0_kind_phys * richardson_number * ( 1.0_kind_phys + 8.0_kind_phys * richardson_number ) )
    end function stable_gradient_richardson_stability_parameter

    pure elemental function neutral_exchange_coefficient(mixing_length_squared, shear_squared) result(neutral_k)
        ! Williamson, D., Kiehl, J., Ramanathan, V., Dickinson, R., & Hack, J. (1987).
        ! Description of the NCAR Community Climate Model (CCM1).
        ! University Corporation for Atmospheric Research. https://doi.org/10.5065/D6TB14WH (Original work published 1987)
        ! Equation 2.f.15, page 12
        ! NOTE: shear_squared vriable currently (01/2025) computed in hb_diff.F90 (trbintd) matches referenced equation.

        real(kind_phys), intent(in) :: mixing_length_squared
        real(kind_phys), intent(in) :: shear_squared

        real(kind_phys)             :: neutral_k

        neutral_k = mixing_length_squared * sqrt(shear_squared)
    end function neutral_exchange_coefficient
end module atmos_phys_pbl_utils
