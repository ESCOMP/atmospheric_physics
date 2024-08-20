subroutine test_musica_ccpp_api()

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif

  use iso_c_binding
  use musica_ccpp
  use musica_micm,               only: Rosenbrock, RosenbrockStandardOrder
  use ccpp_kinds,                only: kind_phys
  use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
  use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t

  implicit none
  integer                                        :: solver_type
  integer                                        :: num_column, num_vertical_layer, num_grid_cells
  real(kind_phys)                                :: time_step                 ! s
  real(kind_phys),              dimension(2,2)   :: temperature               ! K
  real(kind_phys),              dimension(2,2)   :: pressure                  ! Pa
  real(kind_phys),              dimension(2,2)   :: dry_air_density           ! kg m-3
  type(ccpp_constituent_prop_ptr_t), allocatable :: constituent_props_ptr(:)
  real(kind_phys),              dimension(5)     :: molar_mass_arr            ! kg mol-1
  real(kind_phys),              dimension(2,2,5) :: constituents              ! kg kg-1
  real(kind_phys),              dimension(2,2)   :: height                    ! km
  integer                                        :: errcode
  character(len=512)                             :: errmsg

  ! local variables
  type(ccpp_constituent_properties_t), allocatable, target :: constituent_props(:)
  type(ccpp_constituent_properties_t), pointer             :: const_prop
  character(len=512)                                       :: species_name, units
  logical                                                  :: tmp_bool, is_advected
  real(kind_phys)                                          :: molar_mass
  integer                                                  :: i

  solver_type = Rosenbrock
  num_column = 2
  num_vertical_layer = 2
  num_grid_cells = num_column * num_vertical_layer
  time_step = 60._kind_phys

  temperature(:,1) = (/ 100._kind_phys, 200._kind_phys /)
  temperature(:,2) = (/ 300._kind_phys, 400._kind_phys /)

  pressure(:,1) = (/ 6000.04_kind_phys, 7000.04_kind_phys /)
  pressure(:,2) = (/ 8000.04_kind_phys, 9000.04_kind_phys /)

  dry_air_density(:,1) = (/ 3.5_kind_phys, 4.5_kind_phys /)
  dry_air_density(:,2) = (/ 5.5_kind_phys, 6.5_kind_phys /)

  molar_mass_arr = (/ 200._kind_phys, 200._kind_phys, 200._kind_phys, 200._kind_phys, 200._kind_phys /)

  constituents(1,1,1:5) = (/ 0.1_kind_phys, 0.2_kind_phys, 0.3_kind_phys, 0.4_kind_phys, 0.5_kind_phys /)
  constituents(1,2,1:5) = (/ 0.41_kind_phys, 0.42_kind_phys, 0.43_kind_phys, 0.44_kind_phys, 0.45_kind_phys /)
  constituents(2,1,1:5) = (/ 0.21_kind_phys, 0.22_kind_phys, 0.23_kind_phys, 0.24_kind_phys, 0.25_kind_phys /)
  constituents(2,2,1:5) = (/ 0.31_kind_phys, 0.32_kind_phys, 0.33_kind_phys, 0.34_kind_phys, 0.35_kind_phys /)

  call musica_ccpp_register(constituent_props, solver_type, num_grid_cells, errmsg, errcode)
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

  call musica_ccpp_init(errmsg, errcode)

  if (errcode /= 0) then
    write(*,*) trim(errmsg)
    stop 3
  endif

  write(*,*) "[musica] Initial Time Step"
  write(*,fmt="(1x,f10.2)") time_step
  write(*,*) "[musica] Initial Temperature"
  write(*,fmt="(2(1x,f10.4))") temperature
  write(*,*) "[musica] Initial Pressure"
  write(*,fmt="(2(1x,f10.4))") pressure
  write(*,*) "[musica] Initial Concentrations"
  write(*,fmt="(4(3x,e13.6))") constituents

  call musica_ccpp_run(time_step, temperature, pressure, dry_air_density, constituent_props_ptr, &
                       constituents, height, errmsg, errcode)
  if (errcode /= 0) then
    write(*,*) trim(errmsg)
    stop 3
  endif

  write(*,*) "[musica] Solved Concentrations"
  write(*,fmt="(4(3x,e13.6))") constituents
  
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
