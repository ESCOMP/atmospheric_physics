subroutine test_micm_ccpp_api()
  use iso_c_binding
  use micm
  use ccpp_kinds,                only: kind_phys
  use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
  use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t

  implicit none

  real(kind_phys)                                     :: time_step             ! s
  real(kind_phys),                   dimension(2,1)   :: temperature           ! K
  real(kind_phys),                   dimension(2,1)   :: pressure              ! Pa
  real(kind_phys),                   dimension(2,1)   :: dry_air_density       ! kg m-3
  type(ccpp_constituent_prop_ptr_t), dimension(5)     :: constituent_props_ptr
  real(kind_phys),                   dimension(5)     :: molar_mass_arr        ! kg mol-1
  real(kind_phys),                   dimension(2,1,5) :: constituents          ! kg kg-1
  integer                                             :: iulog
  integer                                             :: errcode
  character(len=512)                                  :: errmsg

  ! local variables
  type(ccpp_constituent_properties_t), dimension(5), target   :: constituent_props
  type(ccpp_constituent_properties_t), pointer                :: const_prop
  integer                                                     :: i

  time_step = 60._kind_phys

  temperature(1,1) = 206.6374207_kind_phys
  temperature(2,1) = 206.6374207_kind_phys

  pressure(1,1) = 6152.049805_kind_phys
  pressure(2,1) = 6152.049805_kind_phys

  dry_air_density(1,1) = 3.580_kind_phys
  dry_air_density(2,1) = 3.580_kind_phys

  molar_mass_arr = (/ 200._kind_phys, 200._kind_phys, 200._kind_phys, 200._kind_phys, 200._kind_phys /)

  constituents(1,1,1:5) = (/ 0.75_kind_phys, 8.1e-6_kind_phys, 2.42e-17_kind_phys, &
                          1.15e-5_kind_phys, 6.61e-9_kind_phys /)
  constituents(2,1,1:5) = (/ 0.75_kind_phys, 8.1e-6_kind_phys, 2.42e-17_kind_phys, &
                          1.15e-5_kind_phys, 6.61e-9_kind_phys /)

  iulog = 6

  call constituent_props(1)%instantiate( &
    std_name="water_vapor_mixing_ratio",      &
    long_name="water_vapor_mixing_ratio",     &
    units="kg kg-1",                          &
    default_value=0._kind_phys,               &
    vertical_dim="vertical_layer_dimension",  &
    advected=.true.,                          &
    molar_mass=molar_mass_arr(1),             &
    errcode=errcode, errmsg=errmsg)
  call constituent_props(2)%instantiate( &
    std_name="water_vapor_mixing_ratio_a",    &
    long_name="water_vapor_mixing_ratio",     &
    units="kg kg-1",                          &
    default_value=0._kind_phys,               &
    vertical_dim="vertical_layer_dimension",  &
    advected=.true.,                          &
    molar_mass=molar_mass_arr(2),             &
    errcode=errcode, errmsg=errmsg)
  call constituent_props(3)%instantiate( &
    std_name="water_vapor_mixing_ratio_b",    &
    long_name="water_vapor_mixing_ratio",     &
    units="kg kg-1",                          &
    default_value=0._kind_phys,               &
    vertical_dim="vertical_layer_dimension",  &
    advected=.true.,                          &
    molar_mass=molar_mass_arr(3),             &
    errcode=errcode, errmsg=errmsg)
  call constituent_props(4)%instantiate( &
    std_name="water_vapor_mixing_ratio_c",    &
    long_name="water_vapor_mixing_ratio",     &
    units="kg kg-1",                          &
    default_value=0._kind_phys,               &
    vertical_dim="vertical_layer_dimension",  &
    advected=.true.,                          &
    molar_mass=molar_mass_arr(4),             &
    errcode=errcode, errmsg=errmsg)
  call constituent_props(5)%instantiate( &
    std_name="water_vapor_mixing_ratio_d",    &
    long_name="water_vapor_mixing_ratio",     &
    units="kg kg-1",                          &
    default_value=0._kind_phys,               &
    vertical_dim="vertical_layer_dimension",  &
    advected=.true.,                          &
    molar_mass=molar_mass_arr(5),             &
    errcode=errcode, errmsg=errmsg)

  do i = 1, 5
    const_prop => constituent_props(i)
    call constituent_props_ptr(i)%set(const_prop, errcode, errmsg)
  end do

  call micm_init("chapman", iulog, errcode, errmsg)

  write(*,*) "    -- Initial time_step", time_step
  write(*,*) "    -- Initial temp", temperature
  write(*,*) "    -- Initial pressure", pressure
  write(*,*) "    -- Initial concentrations", constituents

  if (errcode == 0) then
    call micm_run(time_step, temperature, pressure, dry_air_density, constituent_props_ptr, &
                constituents, iulog, errcode, errmsg)

    write(*,*) "    -- After solving, conentrations", constituents
  else
    write(*,*) "    -- Exiting due to the error in creating solver"
    stop 3
  endif

  write(*,*) "    -- Completed solving"
  call micm_final(iulog, errcode, errmsg)

end subroutine test_micm_ccpp_api

program run_test_micm_ccpp
implicit none
  call test_micm_ccpp_api()
end program