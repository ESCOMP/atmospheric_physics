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
       ! Deep convection inputs
       cmfmc_deep,                            &
       qc_deep,                               &
       rliq_deep,                             &
       ! Shallow convection inputs
       cmfmc_sh,                              &
       qc_sh,                                 &
       rliq_sh,                               &
       ! Total convection outputs
       cmfmc_total,                           &
       qc_total,                              &
       rliq_total,                            &
       errmsg, errflg)

      ! Input arguments
      integer,           intent(in)  :: ncol                   ! Number of columns
      integer,           intent(in)  :: pver                   ! Number of model layers
      integer,           intent(in)  :: pverp                  ! pver + 1
      real(kind_phys),   intent(in)  :: cmfmc_deep(:,:)        ! Deep convection cloud mass flux [kg m-2 s-1]
      real(kind_phys),   intent(in)  :: qc_deep(:,:)           ! Deep convection cloud water tendency [kg kg-1 s-1]
      real(kind_phys),   intent(in)  :: rliq_deep(:)           ! Deep convection reserved liquid [kg m-2]
      real(kind_phys),   intent(in)  :: cmfmc_sh(:,:)          ! Shallow convection cloud mass flux [kg m-2 s-1]
      real(kind_phys),   intent(in)  :: qc_sh(:,:)             ! Shallow convection cloud water tendency [kg kg-1 s-1]
      real(kind_phys),   intent(in)  :: rliq_sh(:)             ! Shallow convection reserved liquid [kg m-2]

      ! Output arguments
      real(kind_phys),   intent(out) :: cmfmc_total(:,:)       ! Total convective mass flux [kg m-2 s-1]
      real(kind_phys),   intent(out) :: qc_total(:,:)          ! Total convective cloud water tendency [kg kg-1 s-1]
      real(kind_phys),   intent(out) :: rliq_total(:)          ! Total convective reserved liquid [kg m-2]
      character(len=*),  intent(out) :: errmsg
      integer,           intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! Sum mass flux at interfaces
      cmfmc_total(:ncol,:pverp) = cmfmc_deep(:ncol,:pverp) + cmfmc_sh(:ncol,:pverp)

      ! Sum cloud water tendency at midpoints
      qc_total(:ncol,:pver) = qc_deep(:ncol,:pver) + qc_sh(:ncol,:pver)

      ! Sum vertically-integrated reserved liquid
      rliq_total(:ncol) = rliq_deep(:ncol) + rliq_sh(:ncol)

  end subroutine convect_shallow_sum_to_deep_run
end module convect_shallow_sum_to_deep
