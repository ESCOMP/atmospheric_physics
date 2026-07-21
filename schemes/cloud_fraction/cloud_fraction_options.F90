! Reads the cloud fraction (cldfrc_nl) namelist options without a run phase.
! This carries the namelist parameters so they can be read independently of
! whether compute_cloud_fraction runs in a given suite
! e.g., in CAM7 the parameters are shared by compute_cloud_fraction_two_moment
! but compute_cloud_fraction itself is not run).
module cloud_fraction_options
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: cloud_fraction_options_init

contains
!> \section arg_table_cloud_fraction_options_init Argument Table
!! \htmlinclude cloud_fraction_options_init.html
  subroutine cloud_fraction_options_init( &
    amIRoot, iulog, &
    inversion_cld_off, cldfrc_freeze_dry, cldfrc_ice, iceopt, &
    rhminl, rhminl_adj_land, rhminh, &
    premit, premib, icecrit, &
    errmsg, errflg)

    ! Input arguments
    logical,         intent(in)  :: amIRoot
    integer,         intent(in)  :: iulog
    logical,         intent(in)  :: inversion_cld_off ! turn off stratification-based cloud fraction (inversion_cld)
    logical,         intent(in)  :: cldfrc_freeze_dry ! flag for Vavrus correction
    logical,         intent(in)  :: cldfrc_ice        ! flag to compute ice cloud fraction
    integer,         intent(in)  :: iceopt            ! option for ice cloud closure
    real(kind_phys), intent(in)  :: rhminl            ! minimum rh for low stable clouds
    real(kind_phys), intent(in)  :: rhminl_adj_land   ! rhminl adjustment for snowfree land
    real(kind_phys), intent(in)  :: rhminh            ! minimum rh for high stable clouds
    real(kind_phys), intent(in)  :: premit            ! top pressure bound for mid level cloud
    real(kind_phys), intent(in)  :: premib            ! bottom pressure bound for mid level cloud
    real(kind_phys), intent(in)  :: icecrit           ! critical RH for ice clouds in Wilson and Ballard closure

    ! Output arguments
    character(len=512), intent(out) :: errmsg            ! error message
    integer,            intent(out) :: errflg            ! error flag

    errflg = 0
    errmsg = ''

    if (amIRoot) then
      write(iulog,*) 'tuning parameters cloud_fraction_options_init: inversion_cld_off ', inversion_cld_off, &
                     ' cldfrc_freeze_dry ', cldfrc_freeze_dry, ' cldfrc_ice ', cldfrc_ice
      write(iulog,*) 'tuning parameters cloud_fraction_options_init: rhminl ', rhminl, ' rhminl_adj_land ', rhminl_adj_land, &
                     ' rhminh ', rhminh, ' premit ', premit, ' premib ', premib
      write(iulog,*) 'tuning parameters cloud_fraction_options_init: iceopt ', iceopt, ' icecrit ', icecrit
    end if

  end subroutine cloud_fraction_options_init

end module cloud_fraction_options
