! Module to load namelist options for the Holstlag-Boville boundary layer scheme
module holtslag_boville_diff_options
  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  ! public CCPP-compliant subroutines
  public :: holtslag_boville_diff_options_init

contains

  subroutine holtslag_boville_diff_options_init( &
    amIRoot, iulog, &
    is_hbr_pbl_scheme, &
    errmsg, errflg)

    ! Input arguments
    logical,            intent(in)    :: amIRoot           ! are we on the MPI root task?
    integer,            intent(in)    :: iulog             ! log output unit
    logical,            intent(in)    :: is_hbr_pbl_scheme ! is HBR = true; is HB = false [flag]

    ! Output arguments
    character(len=512), intent(out)   :: errmsg            ! error message
    integer,            intent(out)   :: errflg            ! error flag

    errmsg = ''
    errflg = 0

    if(amIRoot) then
      if(is_hbr_pbl_scheme) then
        write(iulog,*) 'Holtslag-Boville PBL: initializing as HBR PBL scheme.'
      else
        write(iulog,*) 'Holtslag-Boville PBL: initializing as HB PBL scheme.'
      endif
    endif

  end subroutine holtslag_boville_diff_options_init

end module holtslag_boville_diff_options
