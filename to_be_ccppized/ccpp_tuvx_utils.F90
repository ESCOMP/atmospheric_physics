! Copyright (C) 2024-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
module ccpp_tuvx_utils

  implicit none

  private
  public :: rebin, read_extraterrestrial_flux

contains

  !> Regrids normalized flux data to match a specified wavelength grid
  ! This function is copied from CAM/src/chemistry/utils/mo_util.F90
  subroutine rebin( source_dimension, target_dimension, source_coordinates, &
    target_coordinates, source, target )
  use ccpp_kinds, only: kind_phys

  integer,         intent(in)  :: source_dimension
  integer,         intent(in)  :: target_dimension
  real(kind_phys), intent(in)  :: source_coordinates(source_dimension+1)
  real(kind_phys), intent(in)  :: target_coordinates(target_dimension+1)
  real(kind_phys), intent(in)  :: source(source_dimension)
  real(kind_phys), intent(out) :: target(target_dimension)

  ! local variables
  integer         :: i, si, si1, sil, siu
  real(kind_phys) :: y, sl, su, tl, tu

  do i = 1, target_dimension
    tl = target_coordinates(i)
    if( tl < source_coordinates( source_dimension + 1) ) then
      do sil = 1, source_dimension + 1
        if( tl <= source_coordinates( sil ) ) then
          exit
        end if
      end do
      tu = target_coordinates( i + 1 )
      do siu = 1, source_dimension + 1
        if( tu <= source_coordinates( siu ) ) then
          exit
        end if
      end do
      y   = 0._kind_phys
      sil = max( sil, 2 )
      siu = min( siu, source_dimension + 1 )
      do si = sil, siu
        si1 = si - 1
        sl  = max( tl, source_coordinates( si1 ) )
        su  = min( tu, source_coordinates( si ) )
        y   = y + ( su - sl ) * source( si1 )
      end do
      target(i) = y / (target_coordinates( i + 1 ) - target_coordinates( i ) )
    else
      target(i) = 0._kind_phys
    end if
  end do

  end subroutine rebin

  !> Reads a data file to retrieve the extraterrestrial radiation flux values. 
  ! This function is a temporary implementation and will be replaced in
  ! future versions of the code.
  subroutine read_extraterrestrial_flux()
    character(len=50), dimension(4) :: filepath_of_extraterrestrial_flux

    filepath_of_extraterrestrial_flux(1) = 'musica_configurations/chapman/tuvx/data/profiles/solar/susim_hi.flx'
    filepath_of_extraterrestrial_flux(2) = 'musica_configurations/chapman/tuvx/data/profiles/solar/atlas3_1994_317_a.dat'
    filepath_of_extraterrestrial_flux(3) = 'musica_configurations/chapman/tuvx/data/profiles/solar/sao2010.solref.converted'
    filepath_of_extraterrestrial_flux(4) = 'musica_configurations/chapman/tuvx/data/profiles/solar/neckel.flx'

  end subroutine read_extraterrestrial_flux

end module ccpp_tuvx_utils