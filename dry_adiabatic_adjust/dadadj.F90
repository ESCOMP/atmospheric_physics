module dadadj
  !=======================================================================
  ! GFDL style dry adiabatic adjustment
  !
  ! Method:
  ! if stratification is unstable, adjustment to the dry adiabatic lapse
  ! rate is forced subject to the condition that enthalpy is conserved.
  !=======================================================================

  use ccpp_kinds, only:  kind_phys

  implicit none
  private
  save

  public :: dadadj_init ! init routine
  public :: dadadj_run  ! main routine

  integer  :: nlvdry  ! number of layers from top of model to apply the adjustment
  integer  :: niter   ! number of iterations for convergence

CONTAINS

  !> \section arg_table_dadadj_init Argument Table
  !! \htmlinclude dadadj_init.html
  subroutine dadadj_init(dadadj_nlvdry, dadadj_niter, nz, errmsg, errflg)
    !------------------------------------------------
    !   Input / output parameters
    !------------------------------------------------
    integer,            intent(in)  :: dadadj_nlvdry
    integer,            intent(in)  :: dadadj_niter
    integer,            intent(in)  :: nz
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    if (dadadj_nlvdry >= nz .or. dadadj_nlvdry < 0) then
       errflg = 1
       write(errmsg,*) 'dadadj_init: dadadj_nlvdry=',dadadj_nlvdry,' but must be less than the number of vertical levels ',&
       '(',nz,'), and must be a positive integer.`
    end if

    nlvdry = dadadj_nlvdry
    niter = dadadj_niter

  end subroutine dadadj_init

  !> \section arg_table_dadadj_run Argument Table
  !! \htmlinclude dadadj_run.html
  subroutine dadadj_run( &
       ncol, dt, pmid, pint, pdel, state_t, state_q, cappa, t_tend, &
       q_tend, dadpdf, scheme_name, errmsg, errflg)

    !------------------------------------------------
    !   Input / output parameters
    !------------------------------------------------
    integer,  intent(in)        :: ncol  ! number of atmospheric columns
    real(kind_phys), intent(in) :: dt          ! physics timestep
    real(kind_phys), intent(in) :: pmid(:,:)   ! pressure at model levels
    real(kind_phys), intent(in) :: pint(:,:)   ! pressure at model interfaces
    real(kind_phys), intent(in) :: pdel(:,:)   ! vertical delta-p
    real(kind_phys), intent(in) :: cappa(:,:)  ! variable Kappa
    real(kind_phys), intent(in) :: state_t(:,:)   ! temperature (K)
    real(kind_phys), intent(in) :: state_q(:,:)   ! specific humidity
    real(kind_phys), intent(out), target :: t_tend(:,:)   ! temperature tendency
    real(kind_phys), intent(out), target :: q_tend(:,:)   ! specific humidity tendency
    real(kind_phys), intent(out) :: dadpdf(:,:)  ! PDF of where adjustments happened

    character(len=64),  intent(out) :: scheme_name
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    !------------------------------------------------
    !   Local variables
    !------------------------------------------------

    integer                      :: i,k      ! longitude, level indices
    integer                      :: jiter    ! iteration index
    real(kind_phys), allocatable :: c1dad(:) ! intermediate constant
    real(kind_phys), allocatable :: c2dad(:) ! intermediate constant
    real(kind_phys), allocatable :: c3dad(:) ! intermediate constant
    real(kind_phys), allocatable :: c4dad(:) ! intermediate constant
    real(kind_phys) :: gammad    ! dry adiabatic lapse rate (deg/Pa)
    real(kind_phys) :: zeps      ! convergence criterion (deg/Pa)
    real(kind_phys) :: rdenom    ! reciprocal of denominator of expression
    real(kind_phys) :: dtdp      ! delta-t/delta-p
    real(kind_phys) :: zepsdp    ! zeps*delta-p
    real(kind_phys) :: zgamma    ! intermediate constant
    real(kind_phys) :: qave      ! mean q between levels
    real(kind_phys) :: cappaint  ! Kappa at level intefaces
    real(kind_phys), pointer :: t(:,:)
    real(kind_phys), pointer :: q(:,:)

    logical :: ilconv          ! .TRUE. ==> convergence was attained
    logical :: dodad(ncol)     ! .TRUE. ==> do dry adjustment

    !-----------------------------------------------------------------------

    zeps = 2.0e-5_kind_phys           ! set convergence criteria
    errmsg = ''
    errflg = 0
    scheme_name = 'DADADJ'

    allocate(c1dad(nlvdry), stat=ierr)
    if (ierr /= 0) then
       errcode = ierr
       errmsg = trim(scheme_name)//': Allocate of c1dad(nlvdry) failed'
       return
    end if
    allocate(c2dad(nlvdry), stat=ierr)
    if (ierr /= 0) then
       errcode = ierr
       errmsg = trim(scheme_name)//': Allocate of c2dad(nlvdry) failed'
       return
    end if
    allocate(c3dad(nlvdry), stat=ierr)
    if (ierr /= 0) then
       errcode = ierr
       errmsg = trim(scheme_name)//': Allocate of c3dad(nlvdry) failed'
       return
    end if
    allocate(c4dad(nlvdry), stat=ierr)
    if (ierr /= 0) then
       errcode = ierr
       errmsg = trim(scheme_name)//': Allocate of c4dad(nlvdry) failed'
       return
    end if

    ! t_tend< and tend_dtdq used as workspace until needed to calculate tendencies
    t => t_tend
    q => q_tend

    t = state_t
    q = state_q

    ! Find gridpoints with unstable stratification

    do i = 1, ncol
       cappaint = 0.5_kind_phys*(cappa(i,2) + cappa(i,1))
       gammad = cappaint*0.5_kind_phys*(t(i,2) + t(i,1))/pint(i,2)
       dtdp = (t(i,2) - t(i,1))/(pmid(i,2) - pmid(i,1))
       dodad(i) = (dtdp + zeps) > gammad
    end do

    dadpdf(:ncol,:) = 0._kind_phys
    do k= 2, nlvdry
       do i = 1, ncol
         cappaint = 0.5_kind_phys*(cappa(i,k+1) + cappa(i,k))
         gammad = cappaint*0.5_kind_phys*(t(i,k+1) + t(i,k))/pint(i,k+1)
         dtdp = (t(i,k+1) - t(i,k))/(pmid(i,k+1) - pmid(i,k))
         dodad(i) = dodad(i) .or. (dtdp + zeps) > gammad
         if ((dtdp + zeps) > gammad) then
            dadpdf(i,k) = 1._kind_phys
         end if
      end do
   end do

   ! Make a dry adiabatic adjustment
   ! Note: nlvdry ****MUST**** be < pver

   COL: do i = 1, ncol

      if (dodad(i)) then

         zeps = 2.0e-5_kind_phys

         do k = 1, nlvdry
            cappaint = 0.5_kind_phys*(cappa(i,k+1) + cappa(i,k))
            c1dad(k) = cappaint*0.5_kind_phys*(pmid(i,k+1)-pmid(i,k))/pint(i,k+1)
            c2dad(k) = (1._kind_phys - c1dad(k))/(1._kind_phys + c1dad(k))
            rdenom = 1._kind_phys/(pdel(i,k)*c2dad(k) + pdel(i,k+1))
            c3dad(k) = rdenom*pdel(i,k)
            c4dad(k) = rdenom*pdel(i,k+1)
         end do

         ilconv = .false.

         DBLZEP: do while (.not. ilconv)

            do jiter = 1, niter
               ilconv = .true.

               do k = 1, nlvdry
                  zepsdp = zeps*(pmid(i,k+1) - pmid(i,k))
                  zgamma = c1dad(k)*(t(i,k) + t(i,k+1))

                  if ((t(i,k+1)-t(i,k)) >= (zgamma+zepsdp)) then
                     ilconv = .false.
                     t(i,k+1) = t(i,k)*c3dad(k) + t(i,k+1)*c4dad(k)
                     t(i,k) = c2dad(k)*t(i,k+1)
                     qave = (pdel(i,k+1)*q(i,k+1) + pdel(i,k)*q(i,k))/(pdel(i,k+1)+ pdel(i,k))
                     q(i,k+1) = qave
                     q(i,k) = qave
                  end if

               end do

               if (ilconv) cycle COL ! convergence => next longitude

            end do

            zeps = zeps + zeps
            if (zeps > 1.e-4_kind_phys) then
               errflg = i
               errmsg = trim(scheme_name)//': Convergence failure, zeps > 1.e-4'
               return                ! error return
            end if
         end do DBLZEP

      end if

   end do COL

   deallocate(c1dad, c2dad, c3dad, c4dad)

   t_tend = (t - state_t)/dt
   q_tend = (q - state_q)/dt

 end subroutine dadadj_run

end module dadadj
