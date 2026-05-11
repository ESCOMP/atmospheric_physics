module atmos_phys_pbl_utils
    ! Planetary boundary layer related functions used for vertical diffusion schemes.

    use ccpp_kinds, only: kind_phys

    implicit none
    private

    public :: calc_ideal_gas_rrho
    public :: calc_friction_velocity
    public :: calc_kinematic_heat_flux
    public :: calc_kinematic_water_vapor_flux
    public :: calc_kinematic_buoyancy_flux
    public :: calc_obukhov_length
    public :: calc_virtual_temperature
    public :: calc_eddy_flux_coefficient
    public :: calc_free_atm_eddy_flux_coefficient
    public :: calc_brunt_vaisala_frequency_squared
    public :: calc_bulk_richardson_number
    public :: calc_vertical_shear_squared

    real(kind_phys), parameter :: minimum_friction_velocity     = 0.01_kind_phys ! Assuming minimum for coarse grids

contains

    pure elemental function calc_ideal_gas_rrho(rair, surface_temperature, pmid) result(rrho)
        ! air density reciprocal
        ! Taken from https://glossary.ametsoc.org/wiki/Equation_of_state
        ! where \alpha = rrho
        real(kind_phys), intent(in) :: rair  ! gas constant for dry air         [ J kg-1 K-1 ]
        real(kind_phys), intent(in) :: surface_temperature                    ! [ K          ]
        real(kind_phys), intent(in) :: pmid  ! midpoint pressure (bottom level) [ Pa         ]

        real(kind_phys)             :: rrho  ! 1./bottom level density          [ m3 kg-1    ]

        rrho = rair * surface_temperature / pmid
    end function calc_ideal_gas_rrho

    pure elemental function calc_friction_velocity(taux, tauy, rrho) result(friction_velocity)
        ! https://glossary.ametsoc.org/wiki/Friction_velocity
        ! NOTE: taux,tauy come form the expansion of the Reynolds stress
        !
        ! Also found in:
        ! Stull, Roland B. An Introduction to Boundary Layer Meteorology. Springer Kluwer Academic Publishers, 1988. Print.
        ! DOI: https://doi.org/10.1007/978-94-009-3027-8
        ! Equation 2.10b, page 67

        real(kind_phys), intent(in) :: taux              ! surface u stress          [ N m-2   ]
        real(kind_phys), intent(in) :: tauy              ! surface v stress          [ N m-2   ]
        real(kind_phys), intent(in) :: rrho              ! 1./bottom level density   [ m3 kg-1 ]

        real(kind_phys)             :: friction_velocity ! surface friction velocity [ m s-1   ]

        friction_velocity = max( sqrt( sqrt(taux**2 + tauy**2)*rrho ), minimum_friction_velocity )
    end function calc_friction_velocity

    pure elemental function calc_kinematic_heat_flux(shflx, rrho, cpair) result(khfs)
        real(kind_phys), intent(in)  :: shflx  ! surface heat flux            [ W m-2       ]
        real(kind_phys), intent(in)  :: rrho   ! 1./bottom level density      [ m3 kg-1     ]
        real(kind_phys), intent(in)  :: cpair  ! specific heat of dry air     [ J kg-1 K-1  ]

        real(kind_phys)              :: khfs   ! surface kinematic heat flux  [ m K s-1     ]

        khfs = shflx*rrho/cpair
    end function calc_kinematic_heat_flux

    pure elemental function calc_kinematic_water_vapor_flux(qflx, rrho) result(kqfs)
        real(kind_phys), intent(in)  :: qflx ! water vapor flux                   [ kg m-2 s-1    ]
        real(kind_phys), intent(in)  :: rrho ! 1./bottom level density            [ m3 kg-1       ]

        real(kind_phys)              :: kqfs ! surface kinematic water vapor flux [ kg kg-1 m s-1 ]

        kqfs = qflx*rrho
    end function calc_kinematic_water_vapor_flux

    pure elemental function calc_kinematic_buoyancy_flux(khfs, zvir, ths, kqfs) result(kbfs)
        real(kind_phys), intent(in) :: khfs  ! surface kinematic heat flux          [ m K s-1       ]
        real(kind_phys), intent(in) :: zvir  ! rh2o/rair - 1
        real(kind_phys), intent(in) :: ths   ! potential temperature at surface     [ K             ]
        real(kind_phys), intent(in) :: kqfs  ! surface kinematic water vapor flux   [ kg kg-1 m s-1 ]

        ! (`kbfs = \overline{(w' \theta'_v)}_s`)
        real(kind_phys)             :: kbfs  ! surface kinematic buoyancy flux      [ m K s-1 ]

        kbfs = khfs + zvir*ths*kqfs
    end function calc_kinematic_buoyancy_flux

    pure elemental function calc_obukhov_length(thvs, ustar, g, karman, kbfs) result(obukhov_length)
        ! Stull, Roland B. An Introduction to Boundary Layer Meteorology. Springer Kluwer Academic Publishers, 1988. Print.
        ! DOI: https://doi.org/10.1007/978-94-009-3027-8
        ! Equation 5.7c, page 181
        ! \frac{-\theta*u_*^3}{g*k*\overline{(w' \theta_v')}_s} = frac{-\theta*u_*^3}{g*k*kbfs}

        real(kind_phys), intent(in)  :: thvs              ! virtual potential temperature at surface [ K       ]
        real(kind_phys), intent(in)  :: ustar             ! Surface friction velocity                [ m s-1   ]
        real(kind_phys), intent(in)  :: g                 ! acceleration of gravity                  [ m s-2   ]
        real(kind_phys), intent(in)  :: karman            ! Von Karman's constant (unitless)
        real(kind_phys), intent(in)  :: kbfs              ! surface kinematic buoyancy flux          [ m K s-1 ]

        real(kind_phys)              :: obukhov_length    ! Obukhov length                           [ m       ]

        ! Added sign(...) term to prevent division by 0 and using the fact that
        ! `kbfs = \overline{(w' \theta_v')}_s`
        obukhov_length = -thvs * ustar**3 /                             &
                         (g*karman*(kbfs + sign(1.e-10_kind_phys,kbfs)))
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

    pure elemental function calc_eddy_flux_coefficient(mixing_length_squared, &
                                                       richardson_number,     &
                                                       shear_squared)         &
                                                       result(kvf)
        ! Computes exchange coefficients for turbulent flows.
        !
        ! The stable case (Richardson number, Ri>0) is taken from Holtslag and
        ! Beljaars (1989), ECMWF proceedings.
        ! The unstable case (Richardson number, Ri<0) is taken from  CCM1.

        real(kind_phys), intent(in)  :: mixing_length_squared     ! [ m2     ]
        real(kind_phys), intent(in)  :: richardson_number         ! [ 1      ]
        real(kind_phys), intent(in)  :: shear_squared             ! [ s-2    ]

        real(kind_phys)              :: kvf    ! Eddy diffusivity ! [ m2 s-1 ]

        real(kind_phys)              :: fofri  ! f(ri)
        real(kind_phys)              :: kvn    ! Neutral Kv

        if( richardson_number < 0.0_kind_phys ) then
            fofri = unstable_gradient_richardson_stability_parameter(richardson_number)
        else
            fofri = stable_gradient_richardson_stability_parameter(richardson_number)
        end if
        kvn = neutral_exchange_coefficient(mixing_length_squared, shear_squared)
        kvf = kvn * fofri
    end function calc_eddy_flux_coefficient

    pure elemental function calc_free_atm_eddy_flux_coefficient(mixing_length_squared, &
                                                                richardson_number,     &
                                                                shear_squared)         &
                                                                result(kvf)
        ! same as austausch_atm but only mixing for Ri<0
        ! i.e. no background mixing and mixing for Ri>0

        real(kind_phys), intent(in)  :: mixing_length_squared ! [ m2     ]
        real(kind_phys), intent(in)  :: richardson_number     ! [ 1      ]
        real(kind_phys), intent(in)  :: shear_squared         ! [ s-2    ]

        real(kind_phys)              :: kvf                   ! [ m2 s-1 ]

        real(kind_phys)              :: fofri ! f(ri)
        real(kind_phys)              :: kvn   ! Neutral Kv

        kvf = 0.0_kind_phys
        if( richardson_number < 0.0_kind_phys ) then
            fofri = unstable_gradient_richardson_stability_parameter(richardson_number)
            kvn = neutral_exchange_coefficient(mixing_length_squared, shear_squared)
            kvf = kvn * fofri
        end if
    end function calc_free_atm_eddy_flux_coefficient

    pure elemental function unstable_gradient_richardson_stability_parameter(richardson_number) result(modifier)
        ! Williamson, D., Kiehl, J., Ramanathan, V., Dickinson, R., & Hack, J. (1987).
        ! Description of the NCAR Community Climate Model (CCM1).
        ! University Corporation for Atmospheric Research. https://doi.org/10.5065/D6TB14WH (Original work published 1987)
        ! Equation 2.f.13
        ! \sqrt{ 1-18*Ri }

        real(kind_phys), intent(in)  :: richardson_number

        real(kind_phys)              :: modifier

        modifier = sqrt( 1._kind_phys - 18._kind_phys * richardson_number )
    end function unstable_gradient_richardson_stability_parameter

    pure elemental function stable_gradient_richardson_stability_parameter(richardson_number) result(modifier)
        ! Holtslag, A. A. M., and Beljaars A. C. M. , 1989: Surface flux parameterization schemes: Developments and experiences at KNMI.
        ! ECMWF Workshop on Parameterization of Fluxes and Land Surface, Reading, United Kingdom, ECMWF, 121–147.
        ! equation 20, page 140
        ! Originally used published equation from CCM1, 2.f.12, page 11
        ! \frac{1}{1+10*Ri*(1+8*Ri)}

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
        ! NOTE: shear_squared variable currently (01/2025) computed in hb_diff.F90 (s2 in trbintd(...)) matches referenced equation.

        real(kind_phys), intent(in) :: mixing_length_squared ! [ m2     ]
        real(kind_phys), intent(in) :: shear_squared         ! [ s-2    ]

        real(kind_phys)             :: neutral_k

        neutral_k = mixing_length_squared * sqrt(shear_squared)
    end function neutral_exchange_coefficient

    pure elemental function calc_brunt_vaisala_frequency_squared(thv1, thv2, z1, z2, g) result(n2)
        ! Vallis, G. K. (2006).
        ! Atmospheric and Oceanic Fluid Dynamics.  Cambridge University Press.
        ! https://doi.org/10.1017/CBO9780511790447
        ! Page 94, equation 2.225
        !
        ! Note:  The brunt-vaisala frequency is the same as the "buoyancy frequency" in the atmosphere.
        real(kind_phys), intent(in) :: thv1 ! (virtual) potential temperature of top layer [K]
        real(kind_phys), intent(in) :: thv2 ! (virtual) potential temperature of bottom layer [K]
        real(kind_phys), intent(in) :: z1   ! height of top layer [m]
        real(kind_phys), intent(in) :: z2   ! height of bottom layer [m]
        real(kind_phys), intent(in) :: g    ! gravitational acceleration [m s-2]
        real(kind_phys) :: n2               ! Brunt-Vaisala frequency squared [s-2]

        n2 = g*2.0_kind_phys*(thv1-thv2)/((thv1+thv2)*(z1-z2))
    end function calc_brunt_vaisala_frequency_squared

    pure elemental function calc_bulk_richardson_number(thv1, thv2, z1, z2, s2, g) result(ri)
        ! https://glossary.ametsoc.org/wiki/bulk-richardson-number/
        !
        ! Also found in:
        ! Stull, Roland B. An Introduction to Boundary Layer Meteorology. Springer Kluwer Academic Publishers, 1988. Print.
        ! DOI: https://doi.org/10.1007/978-94-009-3027-8
        ! Equation 5.6.3, page 177
        real(kind_phys), intent(in) :: thv1  ! (virtual) potential temperature of top layer [K]
        real(kind_phys), intent(in) :: thv2  ! (virtual) potential temperature of bottom layer [K]
        real(kind_phys), intent(in) :: z1    ! height of top layer [m]
        real(kind_phys), intent(in) :: z2    ! height of bottom layer [m]
        real(kind_phys), intent(in) :: s2    ! vertical shear squared [s-2]
        real(kind_phys), intent(in) :: g     ! gravitational acceleration [m s-2]
        real(kind_phys)             :: ri    ! bulk Richardson number [1]
        real(kind_phys)             :: n2    ! Brunt-Vaisala frequency [s-2]

        n2 = calc_brunt_vaisala_frequency_squared(thv1, thv2, z1, z2, g)
        ri = n2/s2
    end function calc_bulk_richardson_number

    pure elemental function calc_vertical_shear_squared(u1, u2, v1, v2, z1, z2, minimum_velocity_shear_squared) result(s2)
        ! Hack, J., Boville, B., Briegleb, B., Kiehl, J., & Williamson, D. (1993).
        ! Description of the NCAR Community Climate Model (CCM2). University Corporation for Atmospheric Research.
        ! https://doi.org/10.5065/D6QZ27XV (Original work published 1993)
        ! Page 71, Equation 4.e.10
        real(kind_phys), intent(in) :: u1
        real(kind_phys), intent(in) :: u2
        real(kind_phys), intent(in) :: v1
        real(kind_phys), intent(in) :: v2
        real(kind_phys), intent(in) :: z1
        real(kind_phys), intent(in) :: z2
        real(kind_phys), intent(in) :: minimum_velocity_shear_squared
        real(kind_phys) :: s2

        s2 = max((u1-u2)**2 + (v1-v2)**2, minimum_velocity_shear_squared)/((z1-z2)**2)
    end function calc_vertical_shear_squared
end module atmos_phys_pbl_utils
