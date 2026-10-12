!> This module contains the default configuration for the TEMPO microphysics scheme in the form of Fortran derived types,
!> overriding the upstream counterpart.
module module_mp_tempo_cfgs
    implicit none

    private
    public :: ty_tempo_cfgs, ty_tempo_table_cfgs

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
        character(100) :: ccn_table_name = 'ccn_activate.bin' !! ccn table name
        character(100) :: qrqg_table_name = 'qr_acr_qg_data_tempo_v3' !! rain-graupel collection table name
        character(100) :: qrqs_table_name = 'qr_acr_qs_data_tempo_v3' !! rain-snow collection table name
        character(100) :: freezewater_table_name = 'freeze_water_data_tempo_v3' !! freeze water collection table name
    end type
contains
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
