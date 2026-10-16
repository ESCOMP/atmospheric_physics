! Helper scheme to create subcycles that run only on the first timestep.
! SDFs do not have a conditional construct, so this helper is used to set trip count
! of any subcycles with loop="number_of_first_timestep_only_subcycles" to 1 only
! on the first timestep of an initial run.
module set_first_timestep_only_subcycles

  implicit none
  private

  public :: set_first_timestep_only_subcycles_run

contains

!> \section arg_table_set_first_timestep_only_subcycles_run Argument Table
!! \htmlinclude set_first_timestep_only_subcycles_run.html
  subroutine set_first_timestep_only_subcycles_run(is_first_timestep, num_first_timestep_subcycles, &
                                                   errmsg, errflg)

    ! Input arguments
    logical,          intent(in)  :: is_first_timestep            ! flag for first timestep of an initial run [flag]

    ! Output arguments
    integer,          intent(out) :: num_first_timestep_subcycles ! trip count of first-timestep-only subcycles [count]
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    if (is_first_timestep) then
      num_first_timestep_subcycles = 1
    else
      num_first_timestep_subcycles = 0
    end if

  end subroutine set_first_timestep_only_subcycles_run

end module set_first_timestep_only_subcycles
