! Module to load namelist options for vertical diffusion solver scheme.
module vertical_diffusion_options

  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  ! CCPP-compliant public interfaces
  public :: vertical_diffusion_options_init

contains

!> \section arg_table_vertical_diffusion_options_init Argument Table
!! \htmlinclude arg_table_vertical_diffusion_options_init.html
  subroutine vertical_diffusion_options_init( &
    amIRoot, iulog, &
    do_iss_in, &
    am_correction_in, &
    errmsg, errflg)

    ! Input arguments
    logical,            intent(in)    :: amIRoot           ! are we on the MPI root task?
    integer,            intent(in)    :: iulog             ! log output unit
    logical,            intent(in)    :: do_iss_in         ! Use implicit turbulent surface stress computation
    logical,            intent(in)    :: am_correction_in  ! Do angular momentum conservation correction

    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    errmsg = ''
    errflg = 0

    if(amIRoot) then
      write(iulog,*) "vertical diffusion solver: do_iss ", do_iss_in
      write(iulog,*) "vertical diffusion solver: am_correction ", am_correction_in
    endif

  end subroutine vertical_diffusion_options_init

end module vertical_diffusion_options
