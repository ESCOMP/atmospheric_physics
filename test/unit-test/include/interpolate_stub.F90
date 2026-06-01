module interpolate_data
   use ccpp_kinds, only: kind_phys
   public :: lininterp_init
   public :: lininterp

   integer, public, parameter :: extrap_method_bndry = 1
   public :: interp_type
   type interp_type
     real(kind_phys), pointer :: wgts(:)
     real(kind_phys), pointer :: wgtn(:)
     integer, pointer  :: jjm(:)
     integer, pointer  :: jjp(:)
  end type interp_type

contains
   subroutine lininterp_init(yin, nin, yout, nout, extrap_method, interp_wgts)
    integer, intent(in) :: nin
    integer, intent(in) :: nout
    real(kind_phys), intent(in) :: yin(:)           ! input mesh
    real(kind_phys), intent(in) :: yout(:)         ! output mesh
    integer, intent(in) :: extrap_method       ! if 0 set values outside output grid to 0
                                               ! if 1 set to boundary value
                                               ! if 2 set to cyclic boundaries
    type (interp_type), intent(out) :: interp_wgts

    ! routine is a stub
  end subroutine lininterp_init

  subroutine lininterp(arrin, nin, arrout, nout, interp_wgts)
    integer, intent(in) :: nin, nout
    real(kind_phys), intent(in) :: arrin(nin)
    real(kind_phys), intent(out) :: arrout(nout)
    type (interp_type) :: interp_wgts

    ! routine is a stub
  end subroutine lininterp

end module interpolate_data
