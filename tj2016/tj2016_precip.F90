module TJ2016_precip

    use ccpp_kinds, only: kind_phys
  
    implicit none
    private
    save

    public :: tj2016_precip_run

contains

    !=======================================================================
    ! Moist processes   
    !=======================================================================

    !> \section arg_table_tj2016_precip_run  Argument Table
    !! \htmlinclude tj2016_precip_run.html
    subroutine tj2016_precip_run(ncol, pver, gravit, rairv,                                      &
        cpairv, latvap, rh2o, epsilo, rhoh2o, ps0, etamid, dtime,                                &
        pmid, pdel, T, qv, relhum, precl, tendency_of_air_enthalpy, scheme_name, errmsg, errflg)
        !------------------------------------------------
        !   Input / output parameters
        !------------------------------------------------

        integer,  intent(in)    :: ncol                    ! number of columns
        integer,  intent(in)    :: pver                    ! number of vertical levels
        
        real(kind_phys), intent(in)    :: gravit      ! g: gravitational acceleration (m/s2)
        real(kind_phys), intent(in)    :: rairv(:,:)  ! Rd: dry air gas constant (J/K/kg)
        real(kind_phys), intent(in)    :: cpairv(:,:) ! cp: specific heat of dry air (J/K/kg)
        real(kind_phys), intent(in)    :: latvap      ! L: latent heat of vaporization (J/kg)
        real(kind_phys), intent(in)    :: rh2o        ! Rv: water vapor gas constant (J/K/kg)
        real(kind_phys), intent(in)    :: epsilo      ! ratio of h2o to dry air molecular weights
        real(kind_phys), intent(in)    :: rhoh2o      ! density of liquid water (kg/m3)
        real(kind_phys), intent(in)    :: ps0         ! Base state surface pressure (Pa)
        real(kind_phys), intent(in)    :: etamid(:)   ! hybrid coordinate - midpoints

        real(kind_phys), intent(in)    :: dtime       ! time step (s)
        real(kind_phys), intent(in)    :: pmid(:,:)   ! mid-point pressure (Pa)
        real(kind_phys), intent(in)    :: pdel(:,:)   ! layer thickness (Pa)

        real(kind_phys), intent(in)    :: T(:,:)      ! air temperature (K)
        real(kind_phys), intent(inout) :: qv(:,:)     ! specific humidity Q (kg/kg)

        real(kind_phys), intent(out)   :: relhum(:,:)                   ! relative humidity
        real(kind_phys), intent(out)   :: precl(:)                      ! large-scale precipitation rate (m/s)
        real(kind_phys), intent(out)   :: tendency_of_air_enthalpy(:,:) ! change in air enthalpy (J/kg/s)
        character(len=512), intent(out):: scheme_name
        character(len=512), intent(out):: errmsg
        integer,            intent(out):: errflg

        !------------------------------------------------
        !   Local variables
        !------------------------------------------------

        ! Simple physics specific constants and variables

        real(kind_phys), parameter :: T0=273.16_kind_phys    ! control temperature (K) for calculation of qsat
        real(kind_phys), parameter :: e0=610.78_kind_phys    ! saturation vapor pressure (Pa) at T0 for calculation of qsat

        ! Variables for condensation, precipitation, and temperature
        real(kind_phys) :: qsat                       ! saturation value for Q (kg/kg)
        real(kind_phys) :: tmp, tmp_t, tmp_q
        real(kind_phys) :: localT(ncol, pver)                ! air temperature (K) _before_ run
        ! Loop variables
        integer  :: i, k

        !==========================================================================
        ! Set intial total, convective, and large scale precipitation rates to zero
        !==========================================================================
        scheme_name = "tj2016_precip_run"
        errmsg = ' '
        errflg = 0

        precl   = 0.0_kind_phys
        tendency_of_air_enthalpy    = 0.0_kind_phys

        !=========================================================================
        ! Placeholder location for an optional deep convection parameterization (not included here)
        !=========================================================================
        ! An example could be the simplified Betts-Miller (SBM) convection
        ! parameterization described in Frierson (JAS, 2007).
        ! The parameterization is expected to update 
        ! and the temporary state variables
        ! T and qv. T and qv will then be updated again with the 
        ! large-scale condensation process below.

        !=========================================================================
        ! Large-Scale Condensation and Precipitation without cloud stage 
        !=========================================================================
        localT = T
        do k = 1, pver
            do i = 1, ncol
                qsat = epsilo*e0/pmid(i,k)*exp(-latvap/rh2o*((1._kind_phys/T(i,k))-1._kind_phys/T0)) ! saturation value for Q
                if (qv(i,k) > qsat) then
                ! if > 100% relative humidity rain falls out
                    tmp         = 1._kind_phys/dtime*(qv(i,k)-qsat)/(1._kind_phys+(latvap/cpairv(i,k))*(epsilo*latvap*qsat/(rairv(i,k)*T(i,k)**2))) ! condensation rate
                    tmp_t       = latvap/cpairv(i,k)*tmp       ! dT/dt tendency from large-scale condensation
                    tmp_q       = -tmp                   ! dqv/dt tendency from large-scale condensation
                    precl(i)    = precl(i) + tmp*pdel(i,k)/(gravit*rhoh2o) ! large-scale precipitation rate (m/s)
                    localT(i,k) = localT(i,k)   + tmp_t*dtime ! update T (temperature)
                    qv(i,k)     = qv(i,k)  + tmp_q*dtime ! update qv (specific humidity)
                    ! recompute qsat with updated T
                    qsat = epsilo*e0/pmid(i,k)*exp(-latvap/rh2o*((1._kind_phys/localT(i,k))-1._kind_phys/T0)) ! saturation value for Q
                end if

                relhum(i,k) = qv(i,k) / qsat * 100._kind_phys ! in percent
                tendency_of_air_enthalpy(i,k) = (localT(i,k) - T(i,k)) / dtime * cpairv(i,k)
            end do
        end do

    end subroutine tj2016_precip_run
end module TJ2016_precip
