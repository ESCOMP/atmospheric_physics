module rrtmgp_variables

 implicit none
 private

 public :: rrtmgp_variables_init
 public :: rrtmgp_variables_timestep_init
 public :: rrtmgp_variables_run

CONTAINS

!> \section arg_table_rrtmgp_variables_init Argument Table
!! \htmlinclude rrtmgp_variables_init.html
!!
   subroutine rrtmgp_variables_init(unset_real, active_call_array, tlev,       &
                   fluxlwup_Jac, rad_heat, use_tlev, snow_exists, grau_exists, &
                   errmsg, errflg)
      use ccpp_kinds, only: kind_phys
      real(kind_phys),     intent(in) :: unset_real           ! Definition of "unset"
      logical,            intent(out) :: active_call_array(:) ! Diagnostic subcycles
      real(kind_phys),    intent(out) :: tlev(:,:)            ! Air temperature at interfaces [K]
      real(kind_phys),    intent(out) :: fluxlwup_Jac(:,:)    ! Surface temperature flux Jacobian [W m-2 K-1]
      real(kind_phys),    intent(out) :: rad_heat(:,:)        ! Tendency of dry air enthalpy [J kg-1 s-1]
      logical,            intent(out) :: use_tlev             ! Flag to use temperature at interfaces in radiation calculation
      logical,            intent(out) :: snow_exists          ! Flag to include snow cloud area fraction
      logical,            intent(out) :: grau_exists          ! Flag to include graupel cloud area fraction
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Initialize error variables
      errflg = 0
      errmsg = ''

      ! Initialize the active call array
      active_call_array = .true.

      ! Set tlev & fluxlwup_Jac to unset values; not used by default in CAM-SIMA
      use_tlev = .false.
      tlev = unset_real
      fluxlwup_Jac = unset_real

      ! Initialize rad_heat
      rad_heat = unset_real

      ! REMOVECAM: The grau_exists and snow_exists flags should be set to .true. by
      ! schemes that introduce graupel and/or snow
      ! Set the snow and graupel flags to the values needed for the snapshot test
      snow_exists = .true.
      grau_exists = .false.

   end subroutine rrtmgp_variables_init

!> \section arg_table_rrtmgp_variables_timestep_init Argument Table
!! \htmlinclude rrtmgp_variables_timestep_init.html
!!
   subroutine rrtmgp_variables_timestep_init(ncol, nday, rrtmgp_phys_blksz_lw, &
                 rrtmgp_phys_blksz_sw, errmsg, errflg)
      integer,             intent(in) :: nday                 ! Daytime points dimension
      integer,             intent(in) :: ncol                 ! Total horizontal gridpoints
      integer,            intent(out) :: rrtmgp_phys_blksz_lw ! Number of LW columns to process at once
      integer,            intent(out) :: rrtmgp_phys_blksz_sw ! Number of SW columns to process at once
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errflg = 0
      errmsg = ''

      rrtmgp_phys_blksz_lw = ncol
      rrtmgp_phys_blksz_sw = nday

   end subroutine rrtmgp_variables_timestep_init

!> \section arg_table_rrtmgp_variables_run Argument Table
!! \htmlinclude rrtmgp_variables_run.html
!!
  subroutine rrtmgp_variables_run(graupel_in_rad, grau_exists, snow_exists, &
                  do_grau, do_snow, tiny_rad, errmsg, errflg)
     use ccpp_kinds,              only: kind_phys
     ! Inputs
     logical,                          intent(in) :: graupel_in_rad  ! Namelist control of whether to include graupel in radiation calculation
     logical,                          intent(in) :: grau_exists     ! Flag to include graupel cloud area fraction
     logical,                          intent(in) :: snow_exists     ! Flag to include snow cloud area fraction

     ! Outputs
     logical,                         intent(out) :: do_grau         ! Flag to use graupel in radiation calcuation
     logical,                         intent(out) :: do_snow         ! Flag to use snow in radiation calculation
     real(kind_phys),                 intent(out) :: tiny_rad        ! Definition of tiny for RRTMGP
     character(len=512),              intent(out) :: errmsg
     integer,                         intent(out) :: errflg

     ! Set error variables
     errflg = 0
     errmsg = ''

     ! Set definition of tiny for radiation
     tiny_rad = 1.e-80_kind_phys

     ! Initialize graupel flag
     do_grau = .false.

     ! Determine if we should include graupel in the radiation calculation
     if (graupel_in_rad .and. grau_exists) then
        do_grau = .true.
     end if

     ! Snow included if it exists
     do_snow = snow_exists

  end subroutine rrtmgp_variables_run
end module rrtmgp_variables
