module convect_shallow_sum_to_deep
  use ccpp_kinds, only: kind_phys
  implicit none
  private

  public :: convect_shallow_sum_to_deep_run

contains
  ! Combines shallow and deep convection quantities into total convective quantities.
!> \section arg_table_convect_shallow_sum_to_deep_run Argument Table
!! \htmlinclude arg_table_convect_shallow_sum_to_deep_run.html
  subroutine convect_shallow_sum_to_deep_run( &
       ncol,                                  &
       pver,                                  &
       pverp,                                 &
       ! All convection input+outputs (input is from deep)
       cmfmc_total,                           &
       qc_total,                              &
       rliq_total,                            &
       ! Deep convection inputs
       rprddp,                                &
       ! Shallow convection inputs
       cmfmc_sh,                              &
       qc_sh,                                 &
       rliq_sh,                               &
       rprdsh,                                & ! cmfdqr
       ! Total convection outputs
       rprdtot,                               &
       errmsg, errflg)

      ! Input arguments
      integer,           intent(in)    :: ncol                   ! Number of columns
      integer,           intent(in)    :: pver                   ! Number of model layers
      integer,           intent(in)    :: pverp                  ! pver + 1

      ! Deep convective inputs
      real(kind_phys),   intent(in)    :: rprddp(:,:)            ! Deep convection convective rainout Q tendency [kg kg-1 s-1]

      ! Shallow convective inputs
      real(kind_phys),   intent(in)    :: cmfmc_sh(:,:)          ! Shallow convection cloud mass flux [kg m-2 s-1]
      real(kind_phys),   intent(in)    :: qc_sh(:,:)             ! Shallow convection cloud water tendency [kg kg-1 s-1]
      real(kind_phys),   intent(in)    :: rliq_sh(:)             ! Shallow convection reserved liquid [kg m-2]
      real(kind_phys),   intent(in)    :: rprdsh(:,:)            ! Shallow convection convective rainout Q tendency [kg kg-1 s-1]

      ! Input/output (total) arguments
      ! In: deep only, Out: deep+shallow
      real(kind_phys),   intent(inout) :: cmfmc_total(:,:)     ! Total convective mass flux [kg m-2 s-1]
      real(kind_phys),   intent(inout) :: qc_total(:,:)        ! Total convection cloud water tendency [kg kg-1 s-1]
      real(kind_phys),   intent(inout) :: rliq_total(:)        ! Total convective reserved liquid [kg m-2]

      ! Output arguments
      real(kind_phys),   intent(out)   :: rprdtot(:,:)           ! Total convective rainout Q tendency [kg kg-1 s-1]
      character(len=512),intent(out)   :: errmsg
      integer,           intent(out)   :: errflg

      errmsg = ''
      errflg = 0

      ! Sum convective mass flux at interfaces
      cmfmc_total(:ncol,:pverp) = cmfmc_total(:ncol,:pverp) + cmfmc_sh(:ncol,:pverp)

      ! Sum detrainment of water at midpoints
      qc_total(:ncol,:pver) = qc_total(:ncol,:pver) + qc_sh(:ncol,:pver)

      ! Sum vertically-integrated reserved liquid
      rliq_total(:ncol) = rliq_total(:ncol) + rliq_sh(:ncol)

      ! Sum RPRDTOT = RPRDSH + RPRDDP
      ! CMFDQR is the shallow convective rainout only. The total is stored in RPRDTOT.
      rprdtot(:ncol,:pver) = rprdsh(:ncol,:pver) + rprddp(:ncol,:pver)

  end subroutine convect_shallow_sum_to_deep_run
end module convect_shallow_sum_to_deep
