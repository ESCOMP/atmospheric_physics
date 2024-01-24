! subroutine test_micm_ccpp_api()
program test_micm_ccpp_api
    !-----------------------------------------------------------------------
    !
    ! Purpose: Test MICM CCPP API
    !          ......
    !
    !-----------------------------------------------------------------------
    use iso_c_binding
    use micm
    use ccpp_kinds,                only: kind_phys
    use ccpp_constituent_prop_mod
    ! use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    implicit none

    real(kind_phys)                                     :: time_step         ! s
    real(kind_phys),                   dimension(2,1)   :: temperature       ! K
    real(kind_phys),                   dimension(2,1)   :: pressure          ! Pa
    real(kind_phys),                   dimension(2,1)   :: dry_air_density   ! kg m-3 
    type(ccpp_constituent_prop_ptr_t), dimension(5)     :: constituent_props_ptr
    type(ccpp_constituent_properties_t), dimension(5)   :: constituent_props
    real(kind_phys),                   dimension(5)     :: molar_mass_arr    ! kg mol-1
    real(kind_phys),                   dimension(2,1,5) :: constituents      ! kg kg-1
    integer                                             :: iulog
    integer                                             :: errcode
    character(len=512)                                  :: errmsg

    time_step = 500d0

    temperature(1,1) = 10._kind_phys
    temperature(2,1) = 10._kind_phys

    pressure(1,1) = 100._kind_phys
    pressure(2,1) = 100._kind_phys

    dry_air_density(1,1) = 0.0627_kind_phys
    dry_air_density(2,1) = 0.0627_kind_phys

    molar_mass_arr = (/ 0.1_kind_phys, 0.1_kind_phys, 0.1_kind_phys, 0.1_kind_phys, 0.1_kind_phys /)

    constituents(1,1,1:5) = (/ 1._kind_phys, 2._kind_phys, 3._kind_phys, &
                            4._kind_phys, 5._kind_phys /)
    constituents(2,1,1:5) = (/ 10._kind_phys, 20._kind_phys, 30._kind_phys, &
                             40._kind_phys, 50._kind_phys /)

    iulog = 6

    !------------------------------------------------------------
    call constituent_props(1)%instantiate( &
        std_name="water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water",    &
        long_name="water vapor mixing ratio w.r.t moist air and condensed_water", &
        units="kg kg-1",                                                          &
        default_value=0._kind_phys,                                               &
        vertical_dim="vertical_layer_dimension",                                  &
        advected=.true.,                                                          &
        errcode=errcode, errmsg=errmsg)
    call constituent_props(2)%instantiate( &
        std_name="water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water",    &
        long_name="water vapor mixing ratio w.r.t moist air and condensed_water", &
        units="kg kg-1",                                                          &
        default_value=0._kind_phys,                                               &
        vertical_dim="vertical_layer_dimension",                                  &
        advected=.true.,                                                          &
        ! molar_mass=2000._kind_phys                                                     &
        errcode=errcode, errmsg=errmsg)
    call constituent_props(3)%instantiate( &
        std_name="water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water",    &
        long_name="water vapor mixing ratio w.r.t moist air and condensed_water", &
        units="kg kg-1",                                                          &
        default_value=0._kind_phys,                                               &
        vertical_dim="vertical_layer_dimension",                                  &
        advected=.true.,                                                          &
        errcode=errcode, errmsg=errmsg)
    call constituent_props(4)%instantiate( &
        std_name="water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water",    &
        long_name="water vapor mixing ratio w.r.t moist air and condensed_water", &
        units="kg kg-1",                                                          &
        default_value=0._kind_phys,                                               &
        vertical_dim="vertical_layer_dimension",                                  &
        advected=.true.,                                                          &
        errcode=errcode, errmsg=errmsg)
    call constituent_props(5)%instantiate( &
        std_name="water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water",    &
        long_name="water vapor mixing ratio w.r.t moist air and condensed_water", &
        units="kg kg-1",                                                          &
        default_value=0._kind_phys,                                               &
        vertical_dim="vertical_layer_dimension",                                  &
        advected=.true.,                                                          &
        errcode=errcode, errmsg=errmsg)

    constituent_props=>constituent_props_ptr(1)

    write(*,*) "  * [ATM] Creating MICM"
    call micm_init("chapman", iulog, errcode, errmsg)
    
    write(*,*) "  * [ATM] Initial time_step", time_step
    write(*,*) "  * [ATM] Initial temp", temperature
    write(*,*) "  * [ATM] Initial pressure", pressure
    write(*,*) "  * [ATM] Initial concentrations", constituents

    if (errcode == 0) then 
        write(*,*) "  * [ATM] Starting to solve"
        call micm_run(time_step, temperature, pressure, dry_air_density, constituent_props, &
            molar_mass_arr, constituents, iulog, errcode, errmsg)

        write(*,*) "  * [ATM] After solving, conentrations", constituents
    else
        write(*,*) "  * [ATM] Exiting due to the error in creating solver"
        stop 3
    endif

    write(*,*) "  * [ATM] Completed solving"
    call micm_final(iulog, errcode, errmsg)

! end subroutine
end program

! program run_test_micm_ccpp
! implicit none
!     call test_micm_ccpp_api()
! end program 