! Copyright (C) 2026 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
module musica_ccpp_micm_surface

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: surface_reaction_parameters_from_mode
  public :: find_surface_reactions, find_surface_reactions_in_ordering
  public :: mode_from_label

  !> MICM surface reaction carries two rate parameters per grid cell:
  !! '<label>.effective radius [m]' and '<label>.particle number concentration [# m-3]'
  character(len=*), parameter, public :: RADIUS_SUFFIX = '.effective radius [m]'
  character(len=*), parameter, public :: NUMBER_SUFFIX = '.particle number concentration [# m-3]'
  integer,          parameter, public :: SURFACE_LABEL_LENGTH = 256

contains

  !> Convert one aerosol mode's surface area density and wet number-median
  !! diameter into the two rate parameters of a MICM surface reaction
  !!
  !! MICM computes a surface reaction rate as
  !!   k = 4 pi r^2 N / ( r / D_g + 4 / (v gamma) )
  !! and CAM computes the per-mode heterogeneous rate (mo_usrrxt hetrxtrate) as
  !!   k_m = A_m / ( (d_m / 2) / D_g + 4 / (v gamma) )
  !! so setting
  !!   r = d_m / 2
  !!   N = A_m / (4 pi r^2)
  !! reproduces CAM's per-mode rate. N is a surface-equivalent number, not the
  !! mode's particle number concentration: it carries the lognormal surface
  !! enhancement and the mass fraction of the mapped species types already
  !! included in A_m.
  !!
  !! Both parameters are zero at and above the tropopause layer, where CAM does
  !! not compute tropospheric surface area density, and wherever the mode has
  !! no surface or no diameter.
  pure subroutine surface_reaction_parameters_from_mode(surface_area_density, &
      diameter, tropopause_layer, radius, number_concentration)

    use musica_ccpp_util, only: PI

    real(kind_phys), intent(in)  :: surface_area_density(:,:) ! cm2 cm-3 (column, layer)
    real(kind_phys), intent(in)  :: diameter(:,:)             ! cm (column, layer)
    integer,         intent(in)  :: tropopause_layer(:)       ! tropopause vertical layer index (column)
    real(kind_phys), intent(out) :: radius(:,:)               ! m (column, layer)
    real(kind_phys), intent(out) :: number_concentration(:,:) ! # m-3 (column, layer)

    real(kind_phys), parameter :: CM_TO_M = 1.0e-2_kind_phys
    real(kind_phys), parameter :: CM2_CM3_TO_M2_M3 = 1.0e2_kind_phys
    integer :: i_column, i_layer

    radius(:,:) = 0.0_kind_phys
    number_concentration(:,:) = 0.0_kind_phys

    do i_layer = 1, size(surface_area_density, 2)
      do i_column = 1, size(surface_area_density, 1)
        if (i_layer <= tropopause_layer(i_column)) cycle
        if (surface_area_density(i_column, i_layer) <= 0.0_kind_phys) cycle
        if (diameter(i_column, i_layer) <= 0.0_kind_phys) cycle
        radius(i_column, i_layer) = 0.5_kind_phys * diameter(i_column, i_layer) * CM_TO_M
        number_concentration(i_column, i_layer) = &
            surface_area_density(i_column, i_layer) * CM2_CM3_TO_M2_M3 &
            / (4.0_kind_phys * PI * radius(i_column, i_layer)**2)
      end do
    end do

  end subroutine surface_reaction_parameters_from_mode

  !> Find the MICM surface reactions among the rate parameters and the slots of
  !! their effective radius and particle number concentration
  !!
  !! names(i) is the label of the rate parameter stored at slot indices(i).
  !! A reaction is identified by its effective radius parameter; a reaction
  !! without the matching particle number concentration parameter is an error.
  subroutine find_surface_reactions(names, indices, labels, radius_index, &
      number_index, errmsg, errcode)

    character(len=*),                                  intent(in)  :: names(:)
    integer,                                           intent(in)  :: indices(:)
    character(len=SURFACE_LABEL_LENGTH), allocatable,  intent(out) :: labels(:)
    integer,                             allocatable,  intent(out) :: radius_index(:)
    integer,                             allocatable,  intent(out) :: number_index(:)
    character(len=*),                                  intent(out) :: errmsg
    integer,                                           intent(out) :: errcode

    integer :: i_param, i_partner, i_reaction, num_reactions, label_length

    errmsg = ''
    errcode = 0

    num_reactions = 0
    do i_param = 1, size(names)
      if (ends_with(trim(names(i_param)), RADIUS_SUFFIX)) num_reactions = num_reactions + 1
    end do
    allocate(labels(num_reactions), radius_index(num_reactions), number_index(num_reactions))

    i_reaction = 0
    do i_param = 1, size(names)
      if (.not. ends_with(trim(names(i_param)), RADIUS_SUFFIX)) cycle
      i_reaction = i_reaction + 1
      label_length = len_trim(names(i_param)) - len(RADIUS_SUFFIX)
      if (label_length > SURFACE_LABEL_LENGTH) then
        errcode = 1
        errmsg = '[MUSICA Error] MICM surface reaction label too long: '//trim(names(i_param))
        return
      end if
      labels(i_reaction) = names(i_param)(1:label_length)
      radius_index(i_reaction) = indices(i_param)

      number_index(i_reaction) = -1
      do i_partner = 1, size(names)
        if (trim(names(i_partner)) == trim(labels(i_reaction))//NUMBER_SUFFIX) then
          number_index(i_reaction) = indices(i_partner)
          exit
        end if
      end do
      if (number_index(i_reaction) < 0) then
        errcode = 1
        errmsg = '[MUSICA Error] MICM surface reaction "'//trim(labels(i_reaction))// &
                 '" has no particle number concentration rate parameter'
        return
      end if
    end do

  end subroutine find_surface_reactions

  !> find_surface_reactions over a MICM rate parameter ordering
  subroutine find_surface_reactions_in_ordering(ordering, labels, radius_index, &
      number_index, errmsg, errcode)

    use musica_util, only: mappings_t

    type(mappings_t),                                  intent(in)  :: ordering
    character(len=SURFACE_LABEL_LENGTH), allocatable,  intent(out) :: labels(:)
    integer,                             allocatable,  intent(out) :: radius_index(:)
    integer,                             allocatable,  intent(out) :: number_index(:)
    character(len=*),                                  intent(out) :: errmsg
    integer,                                           intent(out) :: errcode

    ! long enough for the longest label plus the longer suffix
    character(len=SURFACE_LABEL_LENGTH + len(NUMBER_SUFFIX)), allocatable :: names(:)
    integer, allocatable :: indices(:)
    character(len=:),     allocatable :: name
    integer :: i_param

    allocate(names(ordering%size()), indices(ordering%size()))
    do i_param = 1, size(names)
      name = ordering%name(i_param)
      if (len(name) > len(names)) then
        errcode = 1
        errmsg = '[MUSICA Error] MICM rate parameter label too long: '//name
        return
      end if
      names(i_param) = name
      indices(i_param) = ordering%index(i_param)
    end do

    call find_surface_reactions(names, indices, labels, radius_index, number_index, &
                                errmsg, errcode)

  end subroutine find_surface_reactions_in_ordering

  !> Aerosol mode a surface reaction label maps to, by the '_a<mm>' suffix
  !! convention (CAM's per-mode naming): 'SURF.N2O5_hydrolysis_a1' -> 1.
  !! Returns 0 when the label has no such suffix.
  pure integer function mode_from_label(label) result(mode)

    character(len=*), intent(in) :: label

    integer :: i_char, i_digit, label_length

    mode = 0
    label_length = len_trim(label)
    i_char = label_length
    do while (i_char > 0)
      if (label(i_char:i_char) < '0' .or. label(i_char:i_char) > '9') exit
      i_char = i_char - 1
    end do
    ! need one to three digits, preceded by '_a'
    if (i_char == label_length .or. label_length - i_char > 3 .or. i_char < 2) return
    if (label(i_char-1:i_char) /= '_a') return
    do i_digit = i_char + 1, label_length
      mode = 10 * mode + (ichar(label(i_digit:i_digit)) - ichar('0'))
    end do

  end function mode_from_label

  pure logical function ends_with(str, suffix)
    character(len=*), intent(in) :: str, suffix
    ends_with = .false.
    if (len(str) >= len(suffix)) then
      ends_with = str(len(str) - len(suffix) + 1:) == suffix
    end if
  end function ends_with

end module musica_ccpp_micm_surface
