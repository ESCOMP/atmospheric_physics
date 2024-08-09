subroutine test_musica_ccpp_api()

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif

  use iso_c_binding
  use musica_ccpp
  use ccpp_kinds,                only: kind_phys
  use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
  use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t

  implicit none

  real(kind_phys)                                :: time_step                 ! s
  real(kind_phys),              dimension(2,1)   :: temperature               ! K
  real(kind_phys),              dimension(2,1)   :: pressure                  ! Pa
  real(kind_phys),              dimension(2,1)   :: dry_air_density           ! kg m-3
  type(ccpp_constituent_prop_ptr_t), allocatable :: constituent_props_ptr(:)
  real(kind_phys),              dimension(5)     :: molar_mass_arr            ! kg mol-1
  real(kind_phys),              dimension(2,1,5) :: constituents              ! kg kg-1
  real(kind_phys),              dimension(2,1)   :: height                    ! km
  real(kind_phys),              dimension(4,2,1) :: photolysis_rate_constants ! s-1
  integer                                        :: n_vertical_levels
  integer                                        :: errcode
  character(len=512)                             :: errmsg

  ! local variables
  type(ccpp_constituent_properties_t), allocatable, target :: constituent_props(:)
  type(ccpp_constituent_properties_t), pointer             :: const_prop
  integer                                                  :: i
  character(len=512)                                       :: species_name, units
  logical                                                  :: tmp_bool, is_advected
  real(kind_phys)                                          :: molar_mass

  n_vertical_levels = 2
  time_step = 60._kind_phys
  temperature(:,1) = (/ 206.6374207_kind_phys, 206.6374207_kind_phys /)
  pressure(:,1) = (/ 6152.049805_kind_phys, 6152.049805_kind_phys /)
  dry_air_density(:,1) = (/ 3.580_kind_phys, 3.580_kind_phys /)

  molar_mass_arr = (/ 200._kind_phys, 200._kind_phys, 200._kind_phys, 200._kind_phys, 200._kind_phys /)

  constituents(1,1,1:5) = (/ 0.75_kind_phys, 8.1e-6_kind_phys, 2.42e-17_kind_phys, &
                          1.15e-5_kind_phys, 6.61e-9_kind_phys /)
  constituents(2,1,1:5) = (/ 0.75_kind_phys, 8.1e-6_kind_phys, 2.42e-17_kind_phys, &
                          1.15e-5_kind_phys, 6.61e-9_kind_phys /)

  call musica_ccpp_register(constituent_props, errmsg, errcode)
  ASSERT(allocated(constituent_props))
  ASSERT(size(constituent_props) == 5)
  do i = 1, size(constituent_props)
    ASSERT(constituent_props(i)%is_instantiated(errcode, errmsg))
    ASSERT(errcode == 0)
    call constituent_props(i)%standard_name(species_name, errcode, errmsg)
    ASSERT(errcode == 0)
    call constituent_props(i)%molar_mass(molar_mass, errcode, errmsg)
    ASSERT(errcode == 0)
    call constituent_props(i)%is_advected(is_advected, errcode, errmsg)
    ASSERT(errcode == 0)
    tmp_bool = (trim(species_name) == "M" .and. molar_mass == 0.029_kind_phys .and. .not. is_advected) .or. &
               (trim(species_name) == "O2" .and. molar_mass == 0.032_kind_phys .and. .not. is_advected) .or. &
               (trim(species_name) == "O" .and. molar_mass == 0.016_kind_phys .and. .not. is_advected) .or. &
               (trim(species_name) == "O1D" .and. molar_mass == 0.016_kind_phys .and. .not. is_advected) .or. &
               (trim(species_name) == "O3" .and. molar_mass == 0.048_kind_phys .and. is_advected)
    ASSERT(tmp_bool)
    call constituent_props(i)%units(units, errcode, errmsg)
    ASSERT(errcode == 0)
    ASSERT(trim(units) == 'kg kg-1')
  end do
  if (errcode /= 0) then
    write(*,*) errcode, trim(errmsg)
    stop 3
  end if

  allocate(constituent_props_ptr(size(constituent_props)))
  do i = 1, size(constituent_props)
    const_prop => constituent_props(i)
    call constituent_props_ptr(i)%set(const_prop, errcode, errmsg)
  end do

  call musica_ccpp_init(n_vertical_levels, errmsg, errcode)

  if (errcode /= 0) then
    write(*,*) trim(errmsg)
    stop 3
  endif

  write(*,*) "    -- Initial time_step", time_step
  write(*,*) "    -- Initial temp", temperature
  write(*,*) "    -- Initial pressure", pressure
  write(*,*) "    -- Initial concentrations", constituents

  call musica_ccpp_run(time_step, temperature, pressure, dry_air_density, constituent_props_ptr, &
                       constituents, height, photolysis_rate_constants, errmsg, errcode)

  if (errcode /= 0) then
    write(*,*) trim(errmsg)
    stop 3
  endif

  write(*,*) "    -- After solving, concentrations", constituents

  call musica_ccpp_final(errmsg, errcode)

  if (errcode /= 0) then
    write(*,*) trim(errmsg)
    stop 3
  endif

end subroutine test_musica_ccpp_api

program run_test_musica_ccpp
implicit none
  call test_musica_ccpp_api()
end program
