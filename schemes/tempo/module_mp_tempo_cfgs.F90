!> This module contains the default configuration for the TEMPO microphysics scheme in the form of Fortran derived types,
!> overriding the upstream counterpart.
module module_mp_tempo_cfgs
    implicit none

    private
    public :: ty_tempo_cfgs, ty_tempo_table_cfgs
    public :: set_tempo_table_names

    integer, parameter :: table_name_len = 100

    character(table_name_len), target :: ccn_table_name_target = ''
    character(table_name_len), target :: qrqg_table_name_target = ''
    character(table_name_len), target :: qrqs_table_name_target = ''
    character(table_name_len), target :: freezewater_table_name_target = ''

    !> \section arg_table_ty_tempo_cfgs Argument Table
    !! \htmlinclude ty_tempo_cfgs.html
    type :: ty_tempo_cfgs
        logical :: aerosolaware_flag = .true. !! flag to run aerosol-aware microphysics
        logical :: hailaware_flag = .true. !! flag to run hail-aware microphysics
        logical :: ml_for_bl_nc_flag = .false. !! flag to run machine-learning prediction for subgrid cloud number concentration
        logical :: ml_for_nc_flag = .false. !! flag to run machine-learning prediction for tempo cloud number concentration
        logical :: semi_sedi_flag = .false. !! flag for semi-lagrangian sedimentation
        logical :: refl10cm_from_melting_flag = .false. !! flag to calculate reflectivity for melting snow and graupel
        logical :: turn_off_micro_flag = .false. !! flag to turn off all microphysical processes
        logical :: cloud_condensation_flag = .true. !! flag to control cloud condensation
        logical :: verbose = .false. !! flag to turn on verbose print statements
        logical :: check_tables = .false. !! flag to turn on check for lookup tables
        character(32) :: single_moment_nc_opt = 'land' !! option for single moment cloud number concentration when land input is not present (options include 'land' that uses the parameter nt_c_l, 'ocean' that uses nt_c_o, or string value in m^-3, e.g. '10.e6')
        ! flags to turn on/off diagnostic output
        logical :: refl10cm_flag = .true. !! flag to output 10cm reflectivity
        logical :: re_cloud_flag = .true. !! flag to output cloud effective radius
        logical :: re_ice_flag = .true. !! flag to output ice effective radius
        logical :: re_snow_flag = .true. !! flag to output snow effective radius
        logical :: max_hail_diameter_flag = .true. !! flag to output maximum hail diameter
        logical :: rain_med_vol_diam_flag = .false. !! flag to output median volume diameter for rain
        logical :: graupel_med_vol_diam_flag = .false. !! flag to output median volume diameter for graupel
        logical :: cloud_number_mixing_ratio_flag = .false. !! flag to output cloud number mixing ratio
    contains
        procedure :: get_nc_val => resolve_nc_value
    end type

    type :: ty_tempo_table_cfgs
        !! ccn table name
        character(table_name_len), pointer :: ccn_table_name => ccn_table_name_target
        !! rain-graupel collection table name
        character(table_name_len), pointer :: qrqg_table_name => qrqg_table_name_target
        !! rain-snow collection table name
        character(table_name_len), pointer :: qrqs_table_name => qrqs_table_name_target
        !! freeze water collection table name
        character(table_name_len), pointer :: freezewater_table_name => freezewater_table_name_target
    end type
contains
    subroutine set_tempo_table_names( &
            ccn_table_name, qrqg_table_name, qrqs_table_name, freezewater_table_name, &
            errmsg, errflg)
        character(*), intent(in) :: ccn_table_name
        character(*), intent(in) :: qrqg_table_name
        character(*), intent(in) :: qrqs_table_name
        character(*), intent(in) :: freezewater_table_name
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        logical :: exist

        ccn_table_name_target = trim(adjustl(ccn_table_name))
        qrqg_table_name_target = trim(adjustl(qrqg_table_name))
        qrqs_table_name_target = trim(adjustl(qrqs_table_name))
        freezewater_table_name_target = trim(adjustl(freezewater_table_name))

        inquire(file=trim(adjustl(ccn_table_name_target)), iomsg=errmsg, iostat=errflg, exist=exist)

        if (errflg /= 0) then
            return
        end if

        if (.not. exist) then
            errmsg = 'set_tempo_table_names: CCN activation table does not exist at ' // &
                trim(adjustl(ccn_table_name_target))
            errflg = 1

            return
        end if

        inquire(file=trim(adjustl(qrqg_table_name_target)), iomsg=errmsg, iostat=errflg, exist=exist)

        if (errflg /= 0) then
            return
        end if

        if (.not. exist) then
            errmsg = 'set_tempo_table_names: Rain-graupel collection table does not exist at ' // &
                trim(adjustl(qrqg_table_name_target))
            errflg = 1

            return
        end if

        inquire(file=trim(adjustl(qrqs_table_name_target)), iomsg=errmsg, iostat=errflg, exist=exist)

        if (errflg /= 0) then
            return
        end if

        if (.not. exist) then
            errmsg = 'set_tempo_table_names: Rain-snow collection table does not exist at ' // &
                trim(adjustl(qrqs_table_name_target))
            errflg = 1

            return
        end if

        inquire(file=trim(adjustl(freezewater_table_name_target)), iomsg=errmsg, iostat=errflg, exist=exist)

        if (errflg /= 0) then
            return
        end if

        if (.not. exist) then
            errmsg = 'set_tempo_table_names: Freeze water collection table does not exist at ' // &
                trim(adjustl(freezewater_table_name_target))
            errflg = 1

            return
        end if

        errmsg = ''
        errflg = 0
    end subroutine set_tempo_table_names

    pure function resolve_nc_value(this, val_land, val_ocean)
        use ccpp_kinds, only: kind_phys

        class(ty_tempo_cfgs), intent(in) :: this
        real(kind_phys), intent(in) :: val_land, val_ocean
        real(kind_phys) :: resolve_nc_value

        integer :: ierr

        select case (trim(adjustl(this % single_moment_nc_opt)))
            case ('land')
                resolve_nc_value = val_land
            case ('ocean')
                resolve_nc_value = val_ocean
            case default
                read(this % single_moment_nc_opt, *, iostat=ierr) resolve_nc_value

                if (ierr /= 0) then
                    resolve_nc_value = val_land
                end if
        end select
    end function resolve_nc_value
end module module_mp_tempo_cfgs
