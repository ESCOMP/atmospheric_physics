! Reads ZM namelist options without run phase
module zm_conv_options
  use ccpp_kinds, only: kind_phys
  implicit none
  private

  public :: zm_conv_options_init

contains
!> \section arg_table_zm_conv_options_init Argument Table
!! \htmlinclude zm_conv_options_init.html
  subroutine zm_conv_options_init(masterproc, iulog, &
    zmconv_c0_lnd, zmconv_c0_ocn, zmconv_ke, zmconv_ke_lnd, &
    zmconv_momcu, zmconv_momcd, zmconv_num_cin, &
    no_deep_pbl_in, zmconv_tiedke_add, &
    zmconv_capelmt, zmconv_dmpdz, zmconv_parcel_pbl, zmconv_tau)

    integer, intent(in)                  :: zmconv_num_cin  ! Number of negative buoyancy regions that are allowed
                                                             ! before the convection top and CAPE calculations are completed.
    real(kind_phys),intent(in)           :: zmconv_c0_lnd
    real(kind_phys),intent(in)           :: zmconv_c0_ocn
    real(kind_phys),intent(in)           :: zmconv_ke
    real(kind_phys),intent(in)           :: zmconv_ke_lnd
    real(kind_phys),intent(in)           :: zmconv_momcu
    real(kind_phys),intent(in)           :: zmconv_momcd
    logical, intent(in)                  :: no_deep_pbl_in  ! no_deep_pbl = .true. eliminates ZM convection entirely within PBL
    real(kind_phys),intent(in)           :: zmconv_tiedke_add
    real(kind_phys),intent(in)           :: zmconv_capelmt
    real(kind_phys),intent(in)           :: zmconv_dmpdz
    logical, intent(in)                  :: zmconv_parcel_pbl ! Should the parcel properties include PBL mixing?
    real(kind_phys),intent(in)           :: zmconv_tau
    logical, intent(in)                  :: masterproc
    integer, intent(in)                  :: iulog

    if ( masterproc ) then
      write(iulog,*) 'tuning parameters zm_conv_options_init: tau',zmconv_tau
      write(iulog,*) 'tuning parameters zm_conv_options_init: c0_lnd',zmconv_c0_lnd, ', c0_ocn', zmconv_c0_ocn
      write(iulog,*) 'tuning parameters zm_conv_options_init: num_cin', zmconv_num_cin
      write(iulog,*) 'tuning parameters zm_conv_options_init: ke',zmconv_ke
      write(iulog,*) 'tuning parameters zm_conv_options_init: no_deep_pbl',no_deep_pbl_in
      write(iulog,*) 'tuning parameters zm_conv_options_init: zm_capelmt', zmconv_capelmt
      write(iulog,*) 'tuning parameters zm_conv_options_init: zm_dmpdz', zmconv_dmpdz
      write(iulog,*) 'tuning parameters zm_conv_options_init: zm_tiedke_add', zmconv_tiedke_add
      write(iulog,*) 'tuning parameters zm_conv_options_init: zm_parcel_pbl', zmconv_parcel_pbl
    endif

  end subroutine zm_conv_options_init

end module zm_conv_options
