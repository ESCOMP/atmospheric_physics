! Copyright (C) 2024-2026 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

!> Gas-phase chemistry solve using MICM for MUSICA chemistry
module musica_ccpp_chemistry

  implicit none
  private

  public :: musica_ccpp_chemistry_run

contains

  !> \section arg_table_musica_ccpp_chemistry_run Argument Table
  !! \htmlinclude musica_ccpp_chemistry_run.html
  subroutine musica_ccpp_chemistry_run(time_step, temperature, pressure, dry_air_density, &
                                       rate_parameters, constituents, log_output_unit,    &
                                       errmsg, errcode)
    use ccpp_kinds,       only: kind_phys
    use musica_ccpp_micm, only: micm_run

    real(kind_phys),         intent(in)    :: time_step            ! s
    real(kind_phys), target, intent(in)    :: temperature(:,:)     ! K (column, layer)
    real(kind_phys), target, intent(in)    :: pressure(:,:)        ! Pa (column, layer)
    real(kind_phys), target, intent(in)    :: dry_air_density(:,:) ! kg m-3 (column, layer)
    real(kind_phys), target, intent(in)    :: rate_parameters(:,:,:) ! various units (column, layer, parameter)
    real(kind_phys), target, intent(inout) :: constituents(:,:,:)  ! kg kg-1 (column, layer, constituent)
    integer,                 intent(in)    :: log_output_unit      ! file unit number for logging
    character(len=*),        intent(out)   :: errmsg
    integer,                 intent(out)   :: errcode

    errcode = 0
    errmsg = ''

    call micm_run(time_step, temperature, pressure, dry_air_density, rate_parameters, &
                  constituents, log_output_unit, errmsg, errcode)
    if (errcode /= 0) return

  end subroutine musica_ccpp_chemistry_run

end module musica_ccpp_chemistry
