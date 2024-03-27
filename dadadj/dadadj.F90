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
  subroutine dadadj_init(dadadj_nlvdry, dadadj_ninter, vertical_layer_dimension, errmsg, errflg)
    !------------------------------------------------
    !   Input / output parameters
    !------------------------------------------------
    integer,            intent(in)  :: dadadj_nlvdry
    integer,            intent(in)  :: dadadj_ninter
    integer,            intent(in)  :: vertical_layer_dimension
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    if (dadadj_nlvdry >= vertical_layer_dimension .or. dadadj_nlvdry < 0) then
       errflg = 1
       write(errmsg,*) 'dadadj_init: dadadj_nlvdry=',dadadj_nlvdry,' but must be less than the number of vertical levels '
    end if

    nlvdry = dadadj_nlvdry
    niter = dadadj_ninter

  end subroutine dadadj_init

  !> \section arg_table_dadadj_run Argument Table
  !! \htmlinclude dadadj_run.html
  subroutine dadadj_run( &
       horizontal_loop_extent, dt, pmid, pint, pdel, state_t, state_q, cappav, tend_t, &
       tend_q, dadpdf, scheme_name, errmsg, errflg)

    !------------------------------------------------
    !   Input / output parameters
    !------------------------------------------------
    integer,  intent(in)        :: horizontal_loop_extent  ! number of atmospheric columns
    real(kind_phys), intent(in) :: dt          ! physics timestep
    real(kind_phys), intent(in) :: pmid(:,:)   ! pressure at model levels
    real(kind_phys), intent(in) :: pint(:,:)   ! pressure at model interfaces
    real(kind_phys), intent(in) :: pdel(:,:)   ! vertical delta-p
    real(kind_phys), intent(in) :: cappav(:,:) ! variable Kappa
    real(kind_phys), intent(in) :: state_t(:,:)   ! temperature (K)
    real(kind_phys), intent(in) :: state_q(:,:)   ! specific humidity
    real(kind_phys), intent(out), TARGET :: tend_t(:,:)   ! temperature tendency
    real(kind_phys), intent(out), TARGET :: tend_q(:,:)   ! specific humidity tendency
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
    real(kind_phys) :: cappa     ! Kappa at level intefaces
    real(kind_phys), pointer :: t(:,:) 
    real(kind_phys), pointer :: q(:,:) 

    logical :: ilconv          ! .TRUE. ==> convergence was attained
    logical :: dodad(horizontal_loop_extent)     ! .TRUE. ==> do dry adjustment

    !-----------------------------------------------------------------------

    zeps = 2.0e-5_kind_phys           ! set convergence criteria
    errmsg = ''
    errflg = 0
    scheme_name = 'DADADJ'

    allocate(c1dad(nlvdry), c2dad(nlvdry), c3dad(nlvdry), c4dad(nlvdry))

    ! tend_t and tend_q used as workspace until needed to calculate tendencies
    t => tend_t
    q => tend_q
    
    t = state_t
    q = state_q

    ! Find gridpoints with unstable stratification

    do i = 1, horizontal_loop_extent
       cappa = 0.5_kind_phys*(cappav(i,2) + cappav(i,1))
       gammad = cappa*0.5_kind_phys*(t(i,2) + t(i,1))/pint(i,2)
       dtdp = (t(i,2) - t(i,1))/(pmid(i,2) - pmid(i,1))
       dodad(i) = (dtdp + zeps) .gt. gammad
    end do

    dadpdf(:horizontal_loop_extent,:) = 0._kind_phys
    do k= 2, nlvdry
       do i = 1, horizontal_loop_extent
         cappa = 0.5_kind_phys*(cappav(i,k+1) + cappav(i,k))
         gammad = cappa*0.5_kind_phys*(t(i,k+1) + t(i,k))/pint(i,k+1)
         dtdp = (t(i,k+1) - t(i,k))/(pmid(i,k+1) - pmid(i,k))
         dodad(i) = dodad(i) .or. (dtdp + zeps).gt.gammad
         if ((dtdp + zeps).gt.gammad) then
            dadpdf(i,k) = 1._kind_phys
         end if
      end do
   end do

   ! Make a dry adiabatic adjustment
   ! Note: nlvdry ****MUST**** be < pver

   i=1
   do while(errflg==0 .and. i <= horizontal_loop_extent)
      if (dodad(i)) then
         do k = 1, nlvdry
            c1dad(k) = cappa*0.5_kind_phys*(pmid(i,k+1)-pmid(i,k))/pint(i,k+1)
            c2dad(k) = (1._kind_phys - c1dad(k))/(1._kind_phys + c1dad(k))
            rdenom = 1._kind_phys/(pdel(i,k)*c2dad(k) + pdel(i,k+1))
            c3dad(k) = rdenom*pdel(i,k)
            c4dad(k) = rdenom*pdel(i,k+1)
         end do

         jiter = 1
         ilconv = .false.
         do while (.not. ilconv .and. zeps <= 1.e-4_kind_phys)

            if (jiter == 1) ilconv = .true.

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

            ! if reach niter double convergence criterion
            if (jiter == niter ) then
               zeps = zeps + zeps
               jiter = 1
            else
               jiter = jiter + 1
            end if

         end do

         if (zeps > 1.e-4_kind_phys) then
            errflg = i
            write(errmsg,*) 'dadadj_run: No convergence in column ',i
            exit
         end if

      end if
      i=i+1
   end do

   deallocate(c1dad, c2dad, c3dad, c4dad)

   tend_t = (t - state_t)/dt
   tend_q = (q - state_q)/dt

 end subroutine dadadj_run

end module dadadj
