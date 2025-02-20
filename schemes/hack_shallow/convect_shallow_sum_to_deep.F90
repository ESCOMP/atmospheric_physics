! Copyright (C) 2024-2025 National Science Foundation-National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Note: There is an implicit dependency in the shallow convective code
! that deep convection runs *before* shallow convection.
! Efforts have been made to make this dependency clear in the form that
! deep convective quantities are merged with the shallow convective
! quantities (incl. determining the final cloud top/bottom, cloud mass flux, ...)
! in this module, and shallow convection modules only using the quantities
! with a _from_shallow_convection suffix.
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
       pver, pverp,                           &
       pmid,                                  &
       ! All convection input+outputs (input is from deep)
       rliq_total,                            &
       ! Deep convection inputs
       cmfmc_deep,                            & ! see notes.
       rprddp,                                &
       qc_deep,                               &
       cnt_deep,                              &
       cnb_deep,                              &
       ! Shallow convection inputs
       cmfmc_sh,                              &
       qc_sh,                                 &
       rliq_sh,                               &
       rprdsh,                                & ! cmfdqr
       cnt_sh,                                &
       cnb_sh,                                &
       ! Total convection outputs
       cmfmc_total,                           &
       rprdtot,                               &
       qc_total,                              &
       cnt, cnb, p_cnt, p_cnb,                &
       errmsg, errflg)

      ! Input arguments
      integer,           intent(in)    :: ncol                   ! # of columns
      integer,           intent(in)    :: pver                   ! # of model layers
      integer,           intent(in)    :: pverp                  ! # of layer interfaces (pver + 1)
      real(kind_phys),   intent(in)    :: pmid(:,:)              ! air_pressure [Pa]

      ! Deep convective inputs
      real(kind_phys),   intent(in)    :: cmfmc_deep(:,:)        ! Deep convection cloud mass flux [kg m-2 s-1]
      real(kind_phys),   intent(in)    :: rprddp(:,:)            ! Deep convection convective rainout Q tendency [kg kg-1 s-1]
      real(kind_phys),   intent(in)    :: qc_deep(:,:)           ! Deep convection cloud water tendency [kg kg-1 s-1]
      integer,           intent(in)    :: cnt_deep(:)            ! Deep convection cloud top index [index] (jctop / jt)
      integer,           intent(in)    :: cnb_deep(:)            ! Deep convection cloud base index [index] (jcbot / maxg)

      ! Shallow convective inputs
      real(kind_phys),   intent(in)    :: cmfmc_sh(:,:)          ! Shallow convection cloud mass flux [kg m-2 s-1]
      real(kind_phys),   intent(in)    :: qc_sh(:,:)             ! Shallow convection cloud water tendency [kg kg-1 s-1]
      real(kind_phys),   intent(in)    :: rliq_sh(:)             ! Shallow convection reserved liquid [m s-1]
      real(kind_phys),   intent(in)    :: rprdsh(:,:)            ! Shallow convection convective rainout Q tendency [kg kg-1 s-1]
      integer,           intent(in)    :: cnt_sh(:)              ! Shallow convection cloud top index [index]
      integer,           intent(in)    :: cnb_sh(:)              ! Shallow convection cloud base index [index]

      ! Input/output (total) arguments
      ! At this point, "total" is only deep because deep convection has written
      ! to these quantities but shallow convection has not.
      ! In: deep only, Out: deep+shallow
      real(kind_phys),   intent(inout) :: rliq_total(:)          ! Total convective reserved liquid [m s-1]

      ! Output arguments
      real(kind_phys),   intent(out)   :: cmfmc_total(:,:)       ! Total convective mass flux [kg m-2 s-1]
      real(kind_phys),   intent(out)   :: rprdtot(:,:)           ! Total convective rainout Q tendency [kg kg-1 s-1]
      real(kind_phys),   intent(out)   :: qc_total(:,:)          ! Total convection cloud water tendency [kg kg-1 s-1]
      integer,           intent(out)   :: cnt(:)                 ! Convective cloud top index [index]
      integer,           intent(out)   :: cnb(:)                 ! Convective cloud base index [index]
      real(kind_phys),   intent(out)   :: p_cnt(:)               ! Convective cloud top pressure [Pa]
      real(kind_phys),   intent(out)   :: p_cnb(:)               ! Convective cloud base pressure [Pa]
      character(len=512),intent(out)   :: errmsg
      integer,           intent(out)   :: errflg

      ! Local variables
      integer                          :: i

      errmsg = ''
      errflg = 0

      ! Sum convective mass flux at interfaces
      cmfmc_total(:ncol,:pverp) = cmfmc_deep(:ncol,:pverp) + cmfmc_sh(:ncol,:pverp)

      ! Sum detrainment of water at midpoints
      qc_total(:ncol,:pver) = qc_deep(:ncol,:pver) + qc_sh(:ncol,:pver)

      ! Sum vertically-integrated reserved liquid
      rliq_total(:ncol) = rliq_total(:ncol) + rliq_sh(:ncol)

      ! Sum RPRDTOT = RPRDSH + RPRDDP
      ! CMFDQR is the shallow convective rainout only. The total is stored in RPRDTOT.
      rprdtot(:ncol,:pver) = rprdsh(:ncol,:pver) + rprddp(:ncol,:pver)

      ! Merge the shallow and deep convective "cloud top/base" and compute pressures
      ! Note: Indices decrease with height
      cnt(:ncol) = cnt_deep(:ncol)
      cnb(:ncol) = cnb_deep(:ncol)
      do i = 1, ncol
        ! if shallow cloud top is higher then use shallow cloud top
        if(cnt_sh(i) < cnt(i)) cnt(i) = cnt_sh(i)

        ! if shallow cloud base is lower then use shallow cloud base
        if(cnb_sh(i) > cnb(i)) cnb(i) = cnb_sh(i)

        ! if cloud base is at model top then set it to cloud top
        if(cnb(i) == 1) cnb(i) = cnt(i)

        ! set pressures
        p_cnt(i) = pmid(i, cnt(i))
        p_cnb(i) = pmid(i, cnb(i))
      enddo


  end subroutine convect_shallow_sum_to_deep_run
end module convect_shallow_sum_to_deep
