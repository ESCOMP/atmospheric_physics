!------------------------------------------------------------------------------------
! Ion drag on the horizontal winds from a globally uniform ion drag tensor,
! used by WACCM configurations without ion chemistry.
!
! Ported from the ghg path of CAM's iondrag module (src/physics/waccm/iondrag.F90,
! ghg_init and iondrag_calc_ghg; B. Foster, Feb 2004).
!
! Subroutines:
!   ion_drag_ghg_init: initialize & print settings
!   ion_drag_ghg_run
!------------------------------------------------------------------------------------
module ion_drag_ghg

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  ! public CCPP-compliant subroutines
  public :: ion_drag_ghg_init
  public :: ion_drag_ghg_run

  integer, parameter :: plevtiod = 97

  real(kind_phys) alamxx(plevtiod)
  real(kind_phys) alamxy(plevtiod)

  real(kind_phys) pshtiod(plevtiod)           ! TIME pressure scale height

  real(kind_phys), allocatable :: alamxxi(:)  ! alamxx interpolated to waccm grid
  real(kind_phys), allocatable :: alamxyi(:)  ! alamxy interpoalted to waccm grid

  integer :: ntop_lev = 1
  integer :: nbot_lev = 0

  logical :: doiodrg = .false.

  !
  ! Data statement for ALAMXX
  !
  data alamxx /                                                     &
       0.13902E-17_kind_phys, 0.22222E-17_kind_phys, 0.34700E-17_kind_phys, 0.53680E-17_kind_phys, 0.83647E-17_kind_phys, &
       0.13035E-16_kind_phys, 0.20254E-16_kind_phys, 0.31415E-16_kind_phys, 0.48944E-16_kind_phys, 0.75871E-16_kind_phys, &
       0.11584E-15_kind_phys, 0.17389E-15_kind_phys, 0.25786E-15_kind_phys, 0.37994E-15_kind_phys, 0.58088E-15_kind_phys, &
       0.95179E-15_kind_phys, 0.19052E-14_kind_phys, 0.47869E-14_kind_phys, 0.14284E-13_kind_phys, 0.45584E-13_kind_phys, &
       0.14756E-12_kind_phys, 0.48154E-12_kind_phys, 0.14844E-11_kind_phys, 0.39209E-11_kind_phys, 0.83886E-11_kind_phys, &
       0.14213E-10_kind_phys, 0.20304E-10_kind_phys, 0.27449E-10_kind_phys, 0.39276E-10_kind_phys, 0.59044E-10_kind_phys, &
       0.83683E-10_kind_phys, 0.11377E-09_kind_phys, 0.14655E-09_kind_phys, 0.19059E-09_kind_phys, 0.28338E-09_kind_phys, &
       0.46326E-09_kind_phys, 0.73966E-09_kind_phys, 0.11785E-08_kind_phys, 0.18789E-08_kind_phys, 0.31037E-08_kind_phys, &
       0.53919E-08_kind_phys, 0.97251E-08_kind_phys, 0.17868E-07_kind_phys, 0.33041E-07_kind_phys, 0.61265E-07_kind_phys, &
       0.11406E-06_kind_phys, 0.20912E-06_kind_phys, 0.39426E-06_kind_phys, 0.76691E-06_kind_phys, 0.15113E-05_kind_phys, &
       0.29545E-05_kind_phys, 0.55644E-05_kind_phys, 0.97208E-05_kind_phys, 0.16733E-04_kind_phys, 0.28101E-04_kind_phys, &
       0.36946E-04_kind_phys, 0.44277E-04_kind_phys, 0.50982E-04_kind_phys, 0.57526E-04_kind_phys, 0.64190E-04_kind_phys, &
       0.71471E-04_kind_phys, 0.80311E-04_kind_phys, 0.96121E-04_kind_phys, 0.11356E-03_kind_phys, 0.14131E-03_kind_phys, &
       0.18695E-03_kind_phys, 0.26058E-03_kind_phys, 0.36900E-03_kind_phys, 0.50812E-03_kind_phys, 0.66171E-03_kind_phys, &
       0.80763E-03_kind_phys, 0.92583E-03_kind_phys, 0.10038E-02_kind_phys, 0.10382E-02_kind_phys, 0.10333E-02_kind_phys, &
       0.99732E-03_kind_phys, 0.93994E-03_kind_phys, 0.86984E-03_kind_phys, 0.79384E-03_kind_phys, 0.71691E-03_kind_phys, &
       0.64237E-03_kind_phys, 0.57224E-03_kind_phys, 0.50761E-03_kind_phys, 0.44894E-03_kind_phys, 0.39624E-03_kind_phys, &
       0.34929E-03_kind_phys, 0.30767E-03_kind_phys, 0.27089E-03_kind_phys, 0.23845E-03_kind_phys, 0.20985E-03_kind_phys, &
       0.18462E-03_kind_phys, 0.16233E-03_kind_phys, 0.14260E-03_kind_phys, 0.12510E-03_kind_phys, 0.10955E-03_kind_phys, &
       0.95699E-04_kind_phys, 0.83347E-04_kind_phys/

  !
  ! Data statement for ALAMXY
  !
  data alamxy /                                                     &
       0.74471E-24_kind_phys, 0.22662E-23_kind_phys, 0.69004E-23_kind_phys, 0.20345E-22_kind_phys, 0.58465E-22_kind_phys, &
       0.16542E-21_kind_phys, 0.46240E-21_kind_phys, 0.12795E-20_kind_phys, 0.35226E-20_kind_phys, 0.96664E-20_kind_phys, &
       0.26650E-19_kind_phys, 0.76791E-19_kind_phys, 0.25710E-18_kind_phys, 0.10897E-17_kind_phys, 0.56593E-17_kind_phys, &
       0.30990E-16_kind_phys, 0.16792E-15_kind_phys, 0.85438E-15_kind_phys, 0.40830E-14_kind_phys, 0.18350E-13_kind_phys, &
       0.79062E-13_kind_phys, 0.33578E-12_kind_phys, 0.13348E-11_kind_phys, 0.45311E-11_kind_phys, 0.12443E-10_kind_phys, &
       0.27052E-10_kind_phys, 0.49598E-10_kind_phys, 0.86072E-10_kind_phys, 0.15807E-09_kind_phys, 0.30480E-09_kind_phys, &
       0.55333E-09_kind_phys, 0.96125E-09_kind_phys, 0.15757E-08_kind_phys, 0.25896E-08_kind_phys, 0.48209E-08_kind_phys, &
       0.96504E-08_kind_phys, 0.18494E-07_kind_phys, 0.34296E-07_kind_phys, 0.61112E-07_kind_phys, 0.10738E-06_kind_phys, &
       0.18747E-06_kind_phys, 0.32054E-06_kind_phys, 0.52872E-06_kind_phys, 0.83634E-06_kind_phys, 0.12723E-05_kind_phys, &
       0.18748E-05_kind_phys, 0.26362E-05_kind_phys, 0.36986E-05_kind_phys, 0.52079E-05_kind_phys, 0.72579E-05_kind_phys, &
       0.98614E-05_kind_phys, 0.12775E-04_kind_phys, 0.15295E-04_kind_phys, 0.18072E-04_kind_phys, 0.20959E-04_kind_phys, &
       0.19208E-04_kind_phys, 0.16285E-04_kind_phys, 0.13628E-04_kind_phys, 0.11784E-04_kind_phys, 0.11085E-04_kind_phys, &
       0.11916E-04_kind_phys, 0.14771E-04_kind_phys, 0.20471E-04_kind_phys, 0.29426E-04_kind_phys, 0.42992E-04_kind_phys, &
       0.62609E-04_kind_phys, 0.90224E-04_kind_phys, 0.12870E-03_kind_phys, 0.18281E-03_kind_phys, 0.26029E-03_kind_phys, &
       0.37224E-03_kind_phys, 0.53254E-03_kind_phys, 0.75697E-03_kind_phys, 0.10623E-02_kind_phys, 0.14660E-02_kind_phys, &
       0.19856E-02_kind_phys, 0.26393E-02_kind_phys, 0.34473E-02_kind_phys, 0.44327E-02_kind_phys, 0.56254E-02_kind_phys, &
       0.70672E-02_kind_phys, 0.88174E-02_kind_phys, 0.10960E-01_kind_phys, 0.13613E-01_kind_phys, 0.16934E-01_kind_phys, &
       0.21137E-01_kind_phys, 0.26501E-01_kind_phys, 0.33388E-01_kind_phys, 0.42263E-01_kind_phys, 0.53716E-01_kind_phys, &
       0.68491E-01_kind_phys, 0.87521E-01_kind_phys, 0.11196E+00_kind_phys, 0.14320E+00_kind_phys, 0.18295E+00_kind_phys, &
       0.23321E+00_kind_phys, 0.29631E+00_kind_phys/

contains

!> \section arg_table_ion_drag_ghg_init Argument Table
!! \htmlinclude arg_table_ion_drag_ghg_init.html
  subroutine ion_drag_ghg_init(pver, pref_mid, masterproc, iulog, do_ion_drag, errmsg, errflg)

    use ccpp_scheme_utils, only: ccpp_constituent_index
    use interpolate_data,  only: lininterp

    !------------------Input arguments---------------------------------------

    integer,            intent(in)  :: pver
    real(kind_phys),    intent(in)  :: pref_mid(:)       ! model ref pressure at midpoint [Pa]
    logical,            intent(in)  :: masterproc
    integer,            intent(in)  :: iulog
    logical,            intent(out) :: do_ion_drag       ! ion drag is applied
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    !-----------------local workspace---------------------------------------
    integer k
    integer kinv
    integer m
    integer const_idx
    integer ierr
    integer n_ions_found

    real(kind_phys) rpsh                             ! ref pressure scale height
    real(kind_phys) pshccm(pver)                     ! CCM pressure scale height

    real(kind_phys), parameter :: preftgcm = 5.e-5_kind_phys       ! TIME GCM reference pressure (Pa)

    ! Ion species whose joint presence selects the ions path of CAM's iondrag module
    character(len=3), parameter :: ion_species(7) = (/ 'e  ', 'Op ', 'O2p', 'NOp', 'O  ', 'O2 ', 'N  ' /)

    !------------------------------------------------------------------------

    errmsg = ''
    errflg = 0
    do_ion_drag = .false.

    ! CAM uses the ions path, not this tensor, when all seven ion species are
    ! present. Short-lived species are not constituents and are not seen here.
    ! CAM reports a missing constituent as an error code (its unassigned index
    ! is positive), CAM-SIMA as a negative index with no error; both tests are
    ! needed.
    n_ions_found = 0
    do m = 1, size(ion_species)
       call ccpp_constituent_index(trim(ion_species(m)), const_idx, errcode=ierr)
       if (ierr == 0 .and. const_idx > 0) n_ions_found = n_ions_found + 1
    end do
    if (n_ions_found == size(ion_species)) then
       errflg = 1
       errmsg = 'ion_drag_ghg_init: all ion species are present; this configuration needs the ions ion drag path'
       return
    end if

    ! Activation gate of CAM's iondrag_init, which runs before ghg_init there.
    ! Duplicated here so that the standalone scheme activates as CAM does.
    ! nbot_lev is reset so that a repeated init does not start from a stale value.
    doiodrg = .false.
    nbot_lev = 0

    !-------------------------------------------------------------------------------
    ! find lower bnd for iondrag
    !-------------------------------------------------------------------------------
    if( pref_mid(1) < 0.1_kind_phys ) then
       do k = 1, pver
          if (pref_mid(k) < 50._kind_phys) nbot_lev = k
       end do
    end if

    if (nbot_lev > 0) then
       doiodrg = .true.
    endif

    if ( .not. doiodrg ) then
       if ( masterproc ) then
          write(iulog,*) ' '
          write(iulog,*) 'ion_drag_ghg_init: Does not have waccm level. Ion drag does not apply. '
          write(iulog,*) ' '
       end if
       return
    endif

    ! With the defualt values of nbot_lev and ntop_lev, ion drag calcualtion are NOT carried out
    nbot_lev=0
    ntop_lev=1

    do k = 1, pver
       rpsh=log(1e5_kind_phys/pref_mid(k))
       if (rpsh .gt. 14._kind_phys) nbot_lev  = k
    end do
    if (nbot_lev .gt. ntop_lev) doiodrg=.true.
    if (masterproc) then
       write(iulog,fmt='(a15)') 'From IONDRAGI:'
       write(iulog,fmt='(1a12,1i10)') 'NTOP_LEV  =',ntop_lev
       write(iulog,fmt='(1a12,1i10)') 'NBOT_LEV  =',nbot_lev
       write(iulog,*) 'IONDRAG flag is',doiodrg
    endif
    do_ion_drag = doiodrg
    if (.not.doiodrg) return

    if (allocated(alamxxi)) deallocate(alamxxi)
    if (allocated(alamxyi)) deallocate(alamxyi)
    allocate(alamxxi(pver), alamxyi(pver), stat=ierr, errmsg=errmsg)
    if (ierr /= 0) then
       errflg = ierr
       errmsg = 'ion_drag_ghg_init: allocate of alamxxi, alamxyi failed: '//trim(errmsg)
       return
    end if

    !     obtain TIME/GCM pressure scale height
    pshtiod(1)=-17._kind_phys
    do k=2,plevtiod
       pshtiod(k)=pshtiod(k-1)+0.25_kind_phys
    enddo

    !     map TIME-psh into CCM-psh
    pshtiod=pshtiod-log(preftgcm/1E5_kind_phys)

    !     CCM psh
    !     note that vertical indexing is inverted with respect to CCM standard
    do k=1,pver
       kinv=pver-k+1
       pshccm(kinv)=log(1e5_kind_phys/pref_mid(k))
    enddo

    !     vertical interpolation
    write(iulog,*) ' '
    write(iulog,*) 'iondragi: before lininterp for alamxx'
    write(iulog,*) '          nlatin,nlatout =',plevtiod,pver
    write(iulog,*) '          yin'
    write(iulog,'(1p,5g15.8)') pshtiod
    write(iulog,*) '          yout'
    write(iulog,'(1p,5g15.8)') pshccm
    write(iulog,*) ' '

    call lininterp (alamxx  ,pshtiod,plevtiod, alamxxi   ,pshccm,pver)

    call lininterp (alamxy  ,pshtiod,plevtiod, alamxyi   ,pshccm,pver)

    !     invert indeces back to CCM convention
    alamxxi(1:pver)=alamxxi(pver:1:-1)
    alamxyi(1:pver)=alamxyi(pver:1:-1)

  end subroutine ion_drag_ghg_init

  !=========================================================================

!> \section arg_table_ion_drag_ghg_run Argument Table
!! \htmlinclude arg_table_ion_drag_ghg_run.html
  subroutine ion_drag_ghg_run(ncol, clat, u, v, dudt, dvdt, errmsg, errflg)

    !
    !     This subroutine calculates ion drag using globally uniform
    !     ion drag tensor:
    !
    !                |alamxx       alamxy   |
    !                |                      |
    !         lambda=|                      |
    !                |                      |
    !                |alamyx       alamyy   |
    !
    !     alamxx and alamxy are provided is data statements
    !     alamyy is obtaine from alamxx:
    !
    !
    !       alamyy = alamxx (sin(DIP_ANGLE))**2
    !
    !     where
    !
    !       DIP_ANGLE = arctan(2.*tan(clat))
    !

    !--------------------Input arguments------------------------------------

    integer,            intent(in)  :: ncol
    real(kind_phys),    intent(in)  :: clat(:)       ! latitude [rad]
    real(kind_phys),    intent(in)  :: u(:,:)        ! eastward wind [m s-1]
    real(kind_phys),    intent(in)  :: v(:,:)        ! northward wind [m s-1]
    real(kind_phys),    intent(out) :: dudt(:,:)     ! tendency of eastward wind due to ion drag [m s-2]
    real(kind_phys),    intent(out) :: dvdt(:,:)     ! tendency of northward wind due to ion drag [m s-2]
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    !---------------------Local workspace-------------------------------------

    real(kind_phys) alamyyi                          ! ALAMYY
    real(kind_phys) dipan                            ! dip angle

    integer i
    integer k

    !-------------------------------------------------------------------------

    errmsg = ''
    errflg = 0

    if (.not.doiodrg) then
       dudt(:,:) = 0.0_kind_phys
       dvdt(:,:) = 0.0_kind_phys
       return
    end if

    !     calculate zonal wind drag
    dudt(:,:)=0.0_kind_phys
    do i=1,ncol
       do k=ntop_lev,nbot_lev
          dudt(i,k)=-alamxyi(k)*v(i,k)-alamxxi(k)*u(i,k)
       enddo
    enddo

    !     calculate meridional wind drag
    dvdt(:,:)=0.0_kind_phys
    do i=1,ncol
       dipan=atan(2._kind_phys*tan(clat(i)))
       do k=ntop_lev,nbot_lev
          alamyyi=alamxxi(k)*(sin(dipan))**2._kind_phys
          dvdt(i,k)=+alamxyi(k)*u(i,k)-alamyyi*v(i,k)
       enddo
    enddo

  end subroutine ion_drag_ghg_run

end module ion_drag_ghg
