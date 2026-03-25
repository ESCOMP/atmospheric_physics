! Rebins 4 sea salt bins
!   sslt1, sslt2, sslt3, sslt4
! into 2 bins for radiation:
!   SSLTA (accumulation mode) and SSLTC (coarse mode).
!
! Original author: Francis Vitt, 9 May 2008
! CCPPized: Haipeng Lin, Mar 2026
module rebin_seasalt
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  ! public CCPP-compliant subroutines
  public :: rebin_seasalt_register
  public :: rebin_seasalt_init
  public :: rebin_seasalt_run

  ! Fraction of total sea salt mass in coarse mode [1]
  real(kind_phys), parameter :: wgt_sscm = 6.0_kind_phys / 7.0_kind_phys

  ! Constituent indices (set during init)
  integer :: idx_sslt1, idx_sslt2, idx_sslt3, idx_sslt4
  integer :: idx_sslta, idx_ssltc

  ! Flag for whether sea salt constituents are available
  logical :: has_sslt = .false.

contains

  ! Register SSLTA and SSLTC as non-advected CCPP constituents.
!> \section arg_table_rebin_seasalt_register Argument Table
!! \htmlinclude rebin_seasalt_register.html
  subroutine rebin_seasalt_register( &
    seasalt_constituents, &
    errmsg, errflg)

    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t

    ! Output arguments:
    type(ccpp_constituent_properties_t), allocatable, intent(out) :: seasalt_constituents(:)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    character(len=*), parameter :: subname = 'rebin_seasalt_register'

    errmsg = ''
    errflg = 0

    allocate(seasalt_constituents(2), stat=errflg, errmsg=errmsg)
    if (errflg /= 0) then
      errmsg = subname // ': ' // trim(errmsg)
      return
    end if

    ! SSLTA: sea salt accumulation mode
    call seasalt_constituents(1)%instantiate( &
         std_name          = 'SSLTA', &
         diag_name         = 'SSLTA', &
         long_name         = 'sea salt accumulation mode mixing ratio', &
         units             = 'kg kg-1', &
         vertical_dim      = 'vertical_layer_dimension', &
         min_value         = 1.e-36_kind_phys, &
         advected          = .false., &
         water_species     = .false., &
         mixing_ratio_type = 'dry', &
         errcode           = errflg, &
         errmsg            = errmsg)
    if (errflg /= 0) return

    ! SSLTC: sea salt coarse mode
    call seasalt_constituents(2)%instantiate( &
         std_name          = 'SSLTC', &
         diag_name         = 'SSLTC', &
         long_name         = 'sea salt coarse mode mixing ratio', &
         units             = 'kg kg-1', &
         vertical_dim      = 'vertical_layer_dimension', &
         min_value         = 1.e-36_kind_phys, &
         advected          = .false., &
         water_species     = .false., &
         mixing_ratio_type = 'dry', &
         errcode           = errflg, &
         errmsg            = errmsg)
    if (errflg /= 0) return

  end subroutine rebin_seasalt_register

  ! Initialize constituent indices for sea salt bins and rebinned modes.
!> \section arg_table_rebin_seasalt_init Argument Table
!! \htmlinclude rebin_seasalt_init.html
  subroutine rebin_seasalt_init( &
    amIRoot, iulog, &
    pcnst, const_props, &
    errmsg, errflg)

    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ccpp_const_utils,          only: ccpp_const_get_idx

    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    integer,            intent(in)  :: pcnst
    type(ccpp_constituent_prop_ptr_t), &
                        intent(in)  :: const_props(:)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    character(len=*), parameter :: subname = 'rebin_seasalt_init'

    errmsg = ''
    errflg = 0

    ! Look up the 4 input sea salt bin constituent indices.
    ! These are registered by prescribed_aerosols (bulk aerosol specifier).
    call ccpp_const_get_idx(const_props, 'sslt1', idx_sslt1, errmsg, errflg)
    if (errflg /= 0) return
    call ccpp_const_get_idx(const_props, 'sslt2', idx_sslt2, errmsg, errflg)
    if (errflg /= 0) return
    call ccpp_const_get_idx(const_props, 'sslt3', idx_sslt3, errmsg, errflg)
    if (errflg /= 0) return
    call ccpp_const_get_idx(const_props, 'sslt4', idx_sslt4, errmsg, errflg)
    if (errflg /= 0) return

    if(idx_sslt1 < 0 .or. idx_sslt2 < 0 .or. &
       idx_sslt3 < 0 .or. idx_sslt4 < 0) then
      if(amIRoot) then
        write(iulog,*) subname//': did not find sslt1/2/3/4 and will not rebin seasalt to SSLTA/SSLTC'
      end if
      return
    end if

    ! Look up the output rebinned constituent indices (registered by this scheme)
    call ccpp_const_get_idx(const_props, &
         'SSLTA', idx_sslta, errmsg, errflg)
    if (errflg /= 0) return
    call ccpp_const_get_idx(const_props, &
         'SSLTC', idx_ssltc, errmsg, errflg)
    if (errflg /= 0) return

    has_sslt = .true.

    if (amIRoot) then
      write(iulog,*) subname//': initialized; coarse mode weight = ', wgt_sscm
    end if

  end subroutine rebin_seasalt_init

  ! Rebin 4 sea salt bins into accumulation and coarse modes.
!> \section arg_table_rebin_seasalt_run Argument Table
!! \htmlinclude rebin_seasalt_run.html
  subroutine rebin_seasalt_run( &
    ncol, pver, pcnst, &
    constituents, &
    errmsg, errflg)

    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    integer,            intent(in)    :: pcnst
    real(kind_phys),    intent(inout) :: constituents(:,:,:) ! constituent mixing ratios [kg kg-1]
    character(len=512), intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    ! Local variables
    integer         :: i, k
    real(kind_phys) :: sslt_sum  ! sum of 4 sea salt bins [kg kg-1]

    errmsg = ''
    errflg = 0

    if (.not. has_sslt) return

    do k = 1, pver
      do i = 1, ncol
        sslt_sum = constituents(i, k, idx_sslt1) &
                 + constituents(i, k, idx_sslt2) &
                 + constituents(i, k, idx_sslt3) &
                 + constituents(i, k, idx_sslt4)

        ! Accumulation mode: fraction (1 - wgt_sscm) = 1/7 of total
        constituents(i, k, idx_sslta) = (1.0_kind_phys - wgt_sscm) * sslt_sum

        ! Coarse mode: fraction wgt_sscm = 6/7 of total
        constituents(i, k, idx_ssltc) = wgt_sscm * sslt_sum
      end do
    end do

  end subroutine rebin_seasalt_run

end module rebin_seasalt
