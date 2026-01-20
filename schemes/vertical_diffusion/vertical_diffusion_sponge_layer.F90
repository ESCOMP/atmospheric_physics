! Vertical diffusion sponge layer scheme
! Adds artificial sponge layer vertical diffusion at the top of the atmosphere
!
! CCPP-ized: Haipeng Lin, June 2025
module vertical_diffusion_sponge_layer
  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  ! CCPP-compliant public interfaces
  public :: vertical_diffusion_sponge_layer_init
  public :: vertical_diffusion_sponge_layer_run
  public :: vertical_diffusion_sponge_layer_final

  ! Module variables for sponge layer parameters
  real(kind_phys), allocatable :: kvm_sponge(:)  ! sponge layer diffusion coefficients [m^2 s-1]

contains

!> \section arg_table_vertical_diffusion_sponge_layer_init Argument Table
!! \htmlinclude vertical_diffusion_sponge_layer_init.html
  subroutine vertical_diffusion_sponge_layer_init( &
    amIRoot, iulog, &
    ptop_ref, &
    diff_sponge_fac, &
    errmsg, errflg)

    ! Input arguments
    logical,            intent(in)  :: amIRoot           ! are we on the MPI root task?
    integer,            intent(in)  :: iulog             ! log output unit
    real(kind_phys),    intent(in)  :: ptop_ref          ! reference top pressure [Pa]
    real(kind_phys),    intent(in)  :: diff_sponge_fac

    ! Output arguments
    character(len=512), intent(out) :: errmsg            ! error message
    integer,            intent(out) :: errflg            ! error flag

    ! Local variables
    integer :: k

    errmsg = ''
    errflg = 0

    if (diff_sponge_fac > 0) then
       !
       ! Add sponge layer vertical diffusion based on model top pressure
       !
       ! This code follows the spectral-element dynamical core that has
       ! a hardcoded vertical profile for del2 sponge diffusion.
       !
       ! Note: kvm_sponge values listed below are also listed in CAM's namelist_definition.xml file,
       !       so if changed here then please also change them in namelist_definition.xml in CAM
       !
       if (ptop_ref < 1.e-3_kind_phys) then
          !
          ! WACCM7 (ptop=2.04E-4Pa) or higher top
          !
          allocate(kvm_sponge(6), stat=errflg, errmsg=errmsg)
          if (errflg /= 0) then
             return
          end if
          kvm_sponge(1) = 2.e5_kind_phys
          kvm_sponge(2) = 2.e5_kind_phys
          kvm_sponge(3) = 1.5e5_kind_phys
          kvm_sponge(4) = 1.0e5_kind_phys
          kvm_sponge(5) = 0.5e5_kind_phys
          kvm_sponge(6) = 0.1e5_kind_phys
       else
          !
          ! CAM7 MT (ptop 0.42Pa) and LT (205.48Pa)
          !
          allocate(kvm_sponge(4), stat=errflg, errmsg=errmsg)
          if (errflg /= 0) then
             return
          end if
          kvm_sponge(1) = 2.e4_kind_phys
          kvm_sponge(2) = 2.e4_kind_phys
          kvm_sponge(3) = 0.5e4_kind_phys
          kvm_sponge(4) = 0.1e4_kind_phys
       end if
    end if
    if (amIRoot) then
       write(iulog, *) 'Sponge layer vertical diffusion factor: ', diff_sponge_fac
       write(iulog, *) '(ptop_ref = ', ptop_ref, ' Pa)'
      if (allocated(kvm_sponge)) then
         write(iulog, *) 'Artificial sponge layer vertical diffusion added:'
        do k = 1, size(kvm_sponge(:))
          write(iulog, '(a44,i2,a17,e7.2,a8)') 'vertical diffusion coefficient at interface ', k, &
                                              ' is increased by ', kvm_sponge(k), ' m2 s-1'
        end do
      else
        write(iulog, *) 'No sponge layer vertical diffusion applied'
      end if
    end if

  end subroutine vertical_diffusion_sponge_layer_init

!> \section arg_table_vertical_diffusion_sponge_layer_run Argument Table
!! \htmlinclude vertical_diffusion_sponge_layer_run.html
  subroutine vertical_diffusion_sponge_layer_run( &
    ncol, pverp, &
    diff_sponge_fac, &
    kvm, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pverp
    real(kind_phys),    intent(in)    :: diff_sponge_fac

    ! Input/output arguments
    real(kind_phys),    intent(inout) :: kvm(:,:)         ! Eddy diffusivity for momentum [m^2 s-1], interfaces

    ! Output arguments
    character(len=512), intent(out)   :: errmsg           ! error message
    integer,            intent(out)   :: errflg           ! error flag

    ! Local variables
    integer :: k

    errmsg = ''
    errflg = 0

    ! Add sponge layer vertical diffusion
    if (allocated(kvm_sponge)) then
      do k = 1, size(kvm_sponge(:))
        kvm(:ncol, k) = kvm(:ncol, k) + diff_sponge_fac * kvm_sponge(k)
      end do
    end if

  end subroutine vertical_diffusion_sponge_layer_run

!> \section arg_table_vertical_diffusion_sponge_layer_final Argument Table
!! \htmlinclude vertical_diffusion_sponge_layer_final.html
  subroutine vertical_diffusion_sponge_layer_final( &
    errmsg, errflg)

    ! Output arguments
    character(len=512), intent(out) :: errmsg            ! error message
    integer,            intent(out) :: errflg            ! error flag

    errmsg = ''
    errflg = 0

    if (allocated(kvm_sponge)) then
      deallocate(kvm_sponge)
    end if

  end subroutine vertical_diffusion_sponge_layer_final

end module vertical_diffusion_sponge_layer
