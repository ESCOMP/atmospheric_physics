module initialize_constituents

implicit none
private

public :: initialize_constituents_register
public :: initialize_constituents_run

contains

!> \section arg_table_initialize_constituents_register  Argument Table
!! \htmlinclude initialize_constituents_register.html
subroutine initialize_constituents_register(constituents, errmsg, errcode)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use cam_initfiles,             only: initial_file_get_id
    use pio,                       only: pio_inquire, file_desc_t, pio_inq_varname
    use ccpp_kinds,                only: kind_phys
    ! Dummy variables
    type(ccpp_constituent_properties_t), allocatable, intent(out) :: constituents(:)
    character(len=512), intent(out) :: errmsg
    integer, intent(out) :: errcode
    ! Local variables
    type(file_desc_t),     pointer       :: ncdata
    integer :: num_variables
    integer :: ierr
    integer :: var_index
    integer :: constituent_index
    integer :: known_const_index
    integer :: found_const_count
    logical :: known_constituent
    character(len=256) :: variable_name
    character(len=512) :: alloc_err_msg
    character(len=256), allocatable :: constituent_names(:)
    character(len=65), parameter :: water_species_std_names(6) = &
      (/'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water       ', &
        'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water', &
        'rain_mixing_ratio_wrt_moist_air_and_condensed_water              ', &
        'cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water         ', &
        'snow_mixing_ratio_wrt_moist_air_and_condensed_water              ', &
        'graupel_water_mixing_ratio_wrt_moist_air_and_condensed_water     '/)

    character(len=11), parameter :: const_file_names(6) = (/'cnst_Q     ', &
                                                 'cnst_CLDLIQ', &
                                                 'cnst_RAINQM', &
                                                 'cnst_CLDICE', &
                                                 'cnst_SNOWQM', &
                                                 'cnst_GRAUQM'/)
    character(len=11), parameter :: water_species_number_concentrations(5) = &
                                      (/'cnst_NUMLIQ', &
                                        'cnst_NUMRAI', &
                                        'cnst_NUMICE', &
                                        'cnst_NUMSNO', &
                                        'cnst_NUMGRA'/)

    errcode = 0
    errmsg = ''
    constituent_index = 0
    found_const_count = 0

    ncdata => initial_file_get_id()
    ! See how many variables are present on the file
    ierr = pio_inquire(ncdata, nVariables=num_variables)
    allocate(constituent_names(num_variables), stat=ierr, errmsg=alloc_err_msg)
    if (ierr /= 0) then
       errcode = 1
       write(errmsg,*) 'Failed to allocate "constituent_names": ', trim(alloc_err_msg)
       return
    end if

    ! Loop over all variables in the file and add each constituent to the
    !  dynamic constituent array
    do var_index = 1, num_variables
       ierr = pio_inq_varname(ncdata, var_index, variable_name)
       known_constituent = .false.
       if (index(variable_name, 'cnst_') > 0) then
          constituent_index = constituent_index + 1
          ! Replace with standard name if known, to avoid duplicates
          if (found_const_count < size(water_species_std_names)) then
             do known_const_index = 1, size(const_file_names)
                if (trim(const_file_names(known_const_index)) == trim(variable_name)) then
                   constituent_names(constituent_index) = water_species_std_names(known_const_index)
                   found_const_count = found_const_count + 1
                   known_constituent = .true.
                   exit
                end if
             end do
          end if
          if (.not. known_constituent) then
             constituent_names(constituent_index) = variable_name
          end if
       end if
    end do

    allocate(constituents(constituent_index), stat=ierr, errmsg=alloc_err_msg)
    if (ierr /= 0) then
       errcode = 1
       write(errmsg,*) 'Failed to allocate "constituents": ', trim(alloc_err_msg)
       return
    end if

    do var_index = 1, size(constituents)
       if (any(water_species_number_concentrations == trim(constituent_names(var_index)))) then
          call constituents(var_index)%instantiate(     &
             std_name = constituent_names(var_index),   &
             long_name = constituent_names(var_index),  &
             units = 'kg-1',                            &
             vertical_dim = 'vertical_layer_dimension', &
             min_value = 0.0_kind_phys,                 &
             advected = .true.,                         &
             water_species = .true.,                    &
             mixing_ratio_type = 'wet',                 &
             errcode = errcode,                         &
             errmsg = errmsg)
       else if (any(water_species_std_names == trim(constituent_names(var_index)))) then
          call constituents(var_index)%instantiate(     &
             std_name = constituent_names(var_index),   &
             long_name = constituent_names(var_index),  &
             units = 'kg kg-1',                         &
             vertical_dim = 'vertical_layer_dimension', &
             min_value = 0.0_kind_phys,                 &
             advected = .true.,                         &
             water_species = .true.,                    &
             mixing_ratio_type = 'wet',                 &
             errcode = errcode,                         &
             errmsg = errmsg)
       else
          call constituents(var_index)%instantiate(     &
             std_name = constituent_names(var_index),   &
             long_name = constituent_names(var_index),  &
             units = 'kg kg-1',                         &
             vertical_dim = 'vertical_layer_dimension', &
             min_value = 0.0_kind_phys,                 &
             advected = .true.,                         &
             errcode = errcode,                         &
             errmsg = errmsg)
       end if
       if (errcode /= 0) then
          exit
       end if
    end do

  end subroutine initialize_constituents_register

!> \section arg_table_initialize_constituents_run  Argument Table
!! \htmlinclude initialize_constituents_run.html
  subroutine initialize_constituents_run(errmsg, errcode)
     character(len=512), intent(out) :: errmsg
     integer,            intent(out) :: errcode

     errcode = 0
     errmsg = ''
  end subroutine initialize_constituents_run
end module initialize_constituents
