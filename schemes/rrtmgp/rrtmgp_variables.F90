module rrtmgp_variables

 implicit none
 private

 public :: rrtmgp_variables_run

CONTAINS

!> \section arg_table_rrtmgp_variables_run Argument Table
!! \htmlinclude rrtmgp_variables_run.html
!!
  subroutine rrtmgp_variables_run(cldfsnow, cldfgrau, degrau, icgrauwp, unset_real, graupel_in_rad, do_snow, &
      do_grau, grau_associated, tiny_rad, errmsg, errflg)
     use ccpp_kinds,              only: kind_phys
     ! Inputs
     real(kind_phys), dimension(:,:),  intent(in) :: cldfsnow
     real(kind_phys), dimension(:,:),  intent(in) :: cldfgrau
     real(kind_phys), dimension(:,:),  intent(in) :: degrau
     real(kind_phys), dimension(:,:),  intent(in) :: icgrauwp
     real(kind_phys),                  intent(in) :: unset_real
     logical,                          intent(in) :: graupel_in_rad

     ! Outputs
     logical,                         intent(out) :: do_snow
     logical,                         intent(out) :: do_grau
     logical,                         intent(out) :: grau_associated
     real(kind_phys),                 intent(out) :: tiny_rad
     character(len=512),              intent(out) :: errmsg
     integer,                         intent(out) :: errflg

     ! Set error variables
     errflg = 0
     errmsg = ''

     ! Set definition of tiny for radiation
     tiny_rad = 1.e-80_kind_phys

     ! Initialize flags
     do_snow = .false.
     do_grau = .false.
     grau_associated = .false.

     ! Determine if the snow cloud fraction variable is set to something
     if (cldfsnow(1,1) /= unset_real) then
        do_snow = .true.
     end if

     ! Determine if the graupel cloud fraction variable is set to something
     if (cldfgrau(1,1) /= unset_real) then
        grau_associated = .true.
     end if

     ! Determine if we should include graupel in the radiation calculation
     if (graupel_in_rad .and. ((cldfgrau(1,1) /= unset_real) .and. (degrau(1,1) /= unset_real) .and. (icgrauwp(1,1) /= unset_real))) then
        do_grau = .true.
     end if


  end subroutine rrtmgp_variables_run
end module rrtmgp_variables
