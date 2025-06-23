! Copyright (C) 2024-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

! Stub for auto-generated MUSICA namelist module
module musica_ccpp_namelist

  implicit none

  private

  character(len=250), public :: micm_solver_type = 'Rosenbrock'
  character(len=250), public :: filename_of_micm_configuration = 'musica_configurations/chapman/micm/config.json'
  character(len=250), public :: filename_of_tuvx_configuration = 'musica_configurations/chapman/tuvx/config.json'
  character(len=250), public :: filename_of_tuvx_micm_mapping_configuration = 'musica_configurations/chapman/tuvx_micm_mapping.json'

end module musica_ccpp_namelist
