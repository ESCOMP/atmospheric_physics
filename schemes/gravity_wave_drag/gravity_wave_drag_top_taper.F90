module gravity_wave_drag_top_taper
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: gravity_wave_drag_top_taper_init

contains

  ! Creates vertical ramp coefficients to taper gravity wave drag
  ! at the top levels of the model
!> \section arg_table_gravity_wave_drag_top_taper_init Argument Table
!! \htmlinclude gravity_wave_drag_top_taper_init.html
  subroutine gravity_wave_drag_top_taper_init( &
    pver, &
    amIRoot, iulog, &
    gw_top_taper, &
    nbot_gravity_wave_top_taper, &
    pref_edge, pref_mid, &
    vramp, &
    errmsg, errflg)

    ! Input arguments
    integer,          intent(in)  :: pver
    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog

    logical,          intent(in)  :: gw_top_taper                           ! Apply gravity wave top tapering [flag]
    integer,          intent(in)  :: nbot_gravity_wave_top_taper            ! Bottom level to taper gravity waves at top of model [index]
    real(kind_phys),  intent(in)  :: pref_edge(:)                           ! Reference pressure at layer interfaces [Pa]
    real(kind_phys),  intent(in)  :: pref_mid(:)                            ! Reference pressure at layer midpoints [Pa]

    ! Output arguments
    real(kind_phys),  intent(out) :: vramp(:)                               ! Gravity wave drag tapering coefficients [1]
    
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer :: k                                                            ! Vertical level index
    integer :: topndx                                                       ! Top level index for tapering [index]
    integer :: botndx                                                       ! Bottom level index for tapering [index]

    errmsg = ''
    errflg = 0

    ! Initialize vramp to unity (no tapering by default)
    vramp(:pver) = 1._kind_phys

    if (gw_top_taper) then
      topndx = 1
      botndx = nbot_gravity_wave_top_taper
      
      if (botndx > 1) then
        ! Compute tapering coefficients from bottom index to top
        do k = botndx, topndx, -1
          vramp(k) = vramp(k + 1) / (pref_edge(k + 1) / pref_edge(k))
        end do
        
        if (amIRoot) then
          write(iulog, '(A)') 'gravity_wave_drag_top_taper_init: GW taper coef (vramp):'
          do k = 1, pver
            write(iulog, "('k: ',I4,' taper coef, press(Pa): ',F12.8,E12.4)") k, vramp(k), pref_mid(k)
          end do
        end if
      end if
    end if
  end subroutine gravity_wave_drag_top_taper_init

end module gravity_wave_drag_top_taper
