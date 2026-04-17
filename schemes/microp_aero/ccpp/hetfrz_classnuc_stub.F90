! Stub for heterogeneous freezing by classical nucleation theory (hetfrz_classnuc).
!
! This scheme exists to read the use_hetfrz_classnuc namelist flag and expose it
! via its CCPP standard name so downstream schemes (e.g. nucleate_ice_ccpp,
! eventually PUMAS) can see the flag. It performs no heterogeneous-freezing
! computation -- the real hetfrz_classnuc_ccpp is not yet ported to CAM-SIMA.
!
! Safe to use with use_hetfrz_classnuc=.true. for snapshot-testing flows where
! the flag must propagate downstream; a warning is printed at init time.
module hetfrz_classnuc_stub
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: hetfrz_classnuc_stub_init

contains

!> \section arg_table_hetfrz_classnuc_stub_init Argument Table
!! \htmlinclude hetfrz_classnuc_stub_init.html
  subroutine hetfrz_classnuc_stub_init( &
    amIRoot, iulog, &
    use_hetfrz_classnuc, &
    errmsg, errflg)

    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    logical,            intent(in)  :: use_hetfrz_classnuc

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    if (amIRoot) then
      if (use_hetfrz_classnuc) then
        write(iulog, '(A)') 'WARNING: hetfrz_classnuc_stub: use_hetfrz_classnuc=.true., but ' // &
                            'heterogeneous freezing by classical nucleation theory is not yet ' // &
                            'implemented in CAM-SIMA (stub scheme).'
        write(iulog, '(A)') 'WARNING: hetfrz_classnuc_stub: The flag is still propagated to ' // &
                            'downstream schemes (e.g. nucleate_ice_ccpp) via its CCPP standard ' // &
                            'name, but no hetfrz_classnuc computation will occur.'
      else
        write(iulog, '(A)') 'hetfrz_classnuc_stub: use_hetfrz_classnuc=.false. (stub scheme, ' // &
                            'no-op).'
      end if
    end if

  end subroutine hetfrz_classnuc_stub_init

end module hetfrz_classnuc_stub
