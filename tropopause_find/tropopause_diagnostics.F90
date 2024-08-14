module tropopause_output
  ! ... output tropopause diagnostics
  ! this will be moved to cam_diagnostics when History is available. (hplin, 8/14/24)
  use ccpp_kinds,           only : kind_phys

  implicit none
  private

  ! CCPP-compliant subroutines
  public :: tropopause_diagnostics_init
  public :: tropopause_diagnostics_run

contains
  ! Initialize the output history fields.
  subroutine tropopause_diagnostics_init(errmsg, errflg)

    !use cam_history,         only: history_add_field
    !use cam_history_support, only: horiz_only

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ' '
    errflg = 0

    ! Define the output fields.
    ! call addfld('TROP_P',          horiz_only,  'A', 'Pa',          'Tropopause Pressure',              flag_xyfill=.True.)
    ! call addfld('TROP_T',          horiz_only,  'A', 'K',           'Tropopause Temperature',           flag_xyfill=.True.)
    ! call addfld('TROP_Z',          horiz_only,  'A', 'm',           'Tropopause Height',                flag_xyfill=.True.)
    ! call addfld('TROP_DZ',         (/ 'lev' /), 'A', 'm',           'Relative Tropopause Height')
    ! call addfld('TROP_PD',         (/ 'lev' /), 'A', 'probability', 'Tropopause Probabilty')
    ! call addfld('TROP_FD',         horiz_only,  'A', 'probability', 'Tropopause Found')

    ! call addfld('TROPP_P',         horiz_only,  'A', 'Pa',          'Tropopause Pressure (primary)',    flag_xyfill=.True.)
    ! call addfld('TROPP_T',         horiz_only,  'A', 'K',           'Tropopause Temperature (primary)', flag_xyfill=.True.)
    ! call addfld('TROPP_Z',         horiz_only,  'A', 'm',           'Tropopause Height (primary)',      flag_xyfill=.True.)
    ! call addfld('TROPP_DZ',        (/ 'lev' /), 'A', 'm',           'Relative Tropopause Height (primary)')
    ! call addfld('TROPP_PD',        (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (primary)')
    ! call addfld('TROPP_FD',        horiz_only,  'A', 'probability', 'Tropopause Found (primary)')

    ! call addfld('TROPF_P',         horiz_only,  'A',  'Pa',         'Tropopause Pressure (cold point)',    flag_xyfill=.True.)
    ! call addfld('TROPF_T',         horiz_only,  'A',  'K',          'Tropopause Temperature (cold point)', flag_xyfill=.True.)
    ! call addfld('TROPF_Z',         horiz_only,  'A',  'm',          'Tropopause Height (cold point)',      flag_xyfill=.True.)
    ! call addfld('TROPF_DZ',        (/ 'lev' /),  'A', 'm',          'Relative Tropopause Height (cold point)', flag_xyfill=.True.)
    ! call addfld('TROPF_PD',        (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (cold point)')
    ! call addfld('TROPF_FD',        horiz_only,  'A', 'probability', 'Tropopause Found (cold point)')

    ! call addfld( 'hstobie_trop',   (/ 'lev' /), 'I',  'fraction of model time', 'Lowest level with stratospheric chemsitry')
    ! call addfld( 'hstobie_linoz',  (/ 'lev' /), 'I',  'fraction of model time', 'Lowest possible Linoz level')
    ! call addfld( 'hstobie_tropop', (/ 'lev' /), 'I', 'fraction of model time', &
    !      'Troposphere boundary calculated in chemistry' )

    ! ! If requested, be prepared to output results from all of the methods.
    ! if (output_all) then
    !   call addfld('TROPA_P',  horiz_only,  'A',  'Pa',          'Tropopause Pressure (analytic)',        flag_xyfill=.True.)
    !   call addfld('TROPA_T',  horiz_only,  'A',  'K',          'Tropopause Temperature (analytic)',      flag_xyfill=.True.)
    !   call addfld('TROPA_Z',  horiz_only,  'A',  'm',          'Tropopause Height (analytic)',           flag_xyfill=.True.)
    !   call addfld('TROPA_PD', (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (analytic)')
    !   call addfld('TROPA_FD', horiz_only,  'A', 'probability', 'Tropopause Found (analytic)')

    !   call addfld('TROPC_P',  horiz_only,  'A',  'Pa',         'Tropopause Pressure (climatology)',      flag_xyfill=.True.)
    !   call addfld('TROPC_T',  horiz_only,  'A',  'K',          'Tropopause Temperature (climatology)',   flag_xyfill=.True.)
    !   call addfld('TROPC_Z',  horiz_only,  'A',  'm',          'Tropopause Height (climatology)',        flag_xyfill=.True.)
    !   call addfld('TROPC_PD', (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (climatology)')
    !   call addfld('TROPC_FD', horiz_only,  'A', 'probability', 'Tropopause Found (climatology)')

    !   call addfld('TROPS_P',  horiz_only,  'A',  'Pa',         'Tropopause Pressure (stobie)',           flag_xyfill=.True.)
    !   call addfld('TROPS_T',  horiz_only,  'A',  'K',          'Tropopause Temperature (stobie)',        flag_xyfill=.True.)
    !   call addfld('TROPS_Z',  horiz_only,  'A',  'm',          'Tropopause Height (stobie)',             flag_xyfill=.True.)
    !   call addfld('TROPS_PD', (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (stobie)')
    !   call addfld('TROPS_FD', horiz_only,  'A', 'probability', 'Tropopause Found (stobie)')

    !   call addfld('TROPT_P',  horiz_only,  'A',  'Pa',         'Tropopause Pressure (twmo)',             flag_xyfill=.True.)
    !   call addfld('TROPT_T',  horiz_only,  'A',  'K',          'Tropopause Temperature (twmo)',          flag_xyfill=.True.)
    !   call addfld('TROPT_Z',  horiz_only,  'A',  'm',          'Tropopause Height (twmo)',               flag_xyfill=.True.)
    !   call addfld('TROPT_PD', (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (twmo)')
    !   call addfld('TROPT_FD', horiz_only,  'A', 'probability', 'Tropopause Found (twmo)')

    !   call addfld('TROPW_P',  horiz_only,  'A',  'Pa',         'Tropopause Pressure (WMO)',              flag_xyfill=.True.)
    !   call addfld('TROPW_T',  horiz_only,  'A',  'K',          'Tropopause Temperature (WMO)',           flag_xyfill=.True.)
    !   call addfld('TROPW_Z',  horiz_only,  'A',  'm',          'Tropopause Height (WMO)',                flag_xyfill=.True.)
    !   call addfld('TROPW_PD', (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (WMO)')
    !   call addfld('TROPW_FD', horiz_only,  'A', 'probability', 'Tropopause Found (WMO)')

    !   call addfld('TROPH_P',  horiz_only,  'A',  'Pa',         'Tropopause Pressure (Hybrid Stobie)',    flag_xyfill=.True.)
    !   call addfld('TROPH_T',  horiz_only,  'A',  'K',          'Tropopause Temperature (Hybrid Stobie)', flag_xyfill=.True.)
    !   call addfld('TROPH_Z',  horiz_only,  'A',  'm',          'Tropopause Height (Hybrid Stobie)',      flag_xyfill=.True.)
    !   call addfld('TROPH_PD', (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (Hybrid Stobie)')
    !   call addfld('TROPH_FD', horiz_only,  'A', 'probability', 'Tropopause Found (Hybrid Stobie)')
    ! end if
  end subroutine tropopause_diagnostics_init

  ! Output the tropopause pressure and temperature to the history files. Two sets
  ! of output will be generated, one for the default algorithm and another one
  ! using the default routine, but backed by a climatology when the default
  ! algorithm fails.
  !> \section arg_table_tropopause_diagnostics_run Argument Table
  !! \htmlinclude tropopause_diagnostics_run.html
  subroutine tropopause_diagnostics_run(ncol, pver, lat, pint, pmid, t, zi, zm, phis)
    !use cam_history,  only : outfld

    real(kind_phys), intent(in)         :: ncol          ! Number of atmospheric columns
    real(kind_phys), intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: lat(:,:)      ! Latitudes (radians)
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa), pverp
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m), pverp
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m), pver
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)

    ! Local Variables
    integer       :: i
    integer       :: alg
    integer       :: ncol                     ! number of cloumns in the chunk
    !integer       :: lchnk                    ! chunk identifier
    integer       :: tropLev(pcols)           ! tropopause level index
    real(kind_phys)      :: tropP(pcols)             ! tropopause pressure (Pa)
    real(kind_phys)      :: tropT(pcols)             ! tropopause temperature (K)
    real(kind_phys)      :: tropZ(pcols)             ! tropopause height (m)
    real(kind_phys)      :: tropFound(pcols)         ! tropopause found
    real(kind_phys)      :: tropDZ(pcols, pver)      ! relative tropopause height (m)
    real(kind_phys)      :: tropPdf(pcols, pver)     ! tropopause probability distribution

    ! Information about the chunk.
    !lchnk = pstate%lchnk

    ! Find the tropopause using the default algorithm backed by the climatology.
    call tropopause_find_run(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP=tropP, tropT=tropT, tropZ=tropZ)

    tropPdf(:,:) = 0._kind_phys
    tropFound(:) = 0._kind_phys
    tropDZ(:,:) = fillvalue
    do i = 1, ncol
      if (tropLev(i) /= NOTFOUND) then
        tropPdf(i, tropLev(i)) = 1._kind_phys
        tropFound(i) = 1._kind_phys
        tropDZ(i,:) = zm(i,:) - tropZ(i)
      end if
    end do

    ! call outfld('TROP_P',   tropP(:ncol),      ncol, lchnk)
    ! call outfld('TROP_T',   tropT(:ncol),      ncol, lchnk)
    ! call outfld('TROP_Z',   tropZ(:ncol),      ncol, lchnk)
    ! call outfld('TROP_DZ',  tropDZ(:ncol, :), ncol, lchnk)
    ! call outfld('TROP_PD',  tropPdf(:ncol, :), ncol, lchnk)
    ! call outfld('TROP_FD',  tropFound(:ncol),  ncol, lchnk)


    ! Find the tropopause using just the primary algorithm.
    call tropopause_find_run(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP=tropP, tropT=tropT, tropZ=tropZ, backup=TROP_ALG_NONE)

    tropPdf(:,:) = 0._kind_phys
    tropFound(:) = 0._kind_phys
    tropDZ(:,:) = fillvalue

    do i = 1, ncol
      if (tropLev(i) /= NOTFOUND) then
        tropPdf(i, tropLev(i)) = 1._kind_phys
        tropFound(i) = 1._kind_phys
        tropDZ(i,:) = zm(i,:) - tropZ(i)
      end if
    end do

    ! call outfld('TROPP_P',   tropP(:ncol),      ncol, lchnk)
    ! call outfld('TROPP_T',   tropT(:ncol),      ncol, lchnk)
    ! call outfld('TROPP_Z',   tropZ(:ncol),      ncol, lchnk)
    ! call outfld('TROPP_DZ',  tropDZ(:ncol, :), ncol, lchnk)
    ! call outfld('TROPP_PD',  tropPdf(:ncol, :), ncol, lchnk)
    ! call outfld('TROPP_FD',  tropFound(:ncol),  ncol, lchnk)


    ! Find the tropopause using just the cold point algorithm.
    call tropopause_find_run(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP=tropP, tropT=tropT, tropZ=tropZ, primary=TROP_ALG_CPP, backup=TROP_ALG_NONE)

    tropPdf(:,:) = 0._kind_phys
    tropFound(:) = 0._kind_phys
    tropDZ(:,:) = fillvalue

    do i = 1, ncol
      if (tropLev(i) /= NOTFOUND) then
        tropPdf(i, tropLev(i)) = 1._kind_phys
        tropFound(i) = 1._kind_phys
        tropDZ(i,:) = zm(i,:) - tropZ(i)
      end if
    end do

    ! call outfld('TROPF_P',   tropP(:ncol),      ncol, lchnk)
    ! call outfld('TROPF_T',   tropT(:ncol),      ncol, lchnk)
    ! call outfld('TROPF_Z',   tropZ(:ncol),      ncol, lchnk)
    ! call outfld('TROPF_DZ',  tropDZ(:ncol, :), ncol, lchnk)
    ! call outfld('TROPF_PD',  tropPdf(:ncol, :), ncol, lchnk)
    ! call outfld('TROPF_FD',  tropFound(:ncol),  ncol, lchnk)


    ! If requested, do all of the algorithms.
    if (output_all) then

      do alg = 2, TROP_NALG

        ! Find the tropopause using just the analytic algorithm.
        call tropopause_find_run(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP=tropP, tropT=tropT, tropZ=tropZ, primary=alg, backup=TROP_ALG_NONE)

        tropPdf(:,:) = 0._kind_phys
        tropFound(:) = 0._kind_phys

        do i = 1, ncol
          if (tropLev(i) /= NOTFOUND) then
            tropPdf(i, tropLev(i)) = 1._kind_phys
            tropFound(i) = 1._kind_phys
          end if
        end do

        ! call outfld('TROP' // TROP_LETTER(alg) // '_P',   tropP(:ncol),      ncol, lchnk)
        ! call outfld('TROP' // TROP_LETTER(alg) // '_T',   tropT(:ncol),      ncol, lchnk)
        ! call outfld('TROP' // TROP_LETTER(alg) // '_Z',   tropZ(:ncol),      ncol, lchnk)
        ! call outfld('TROP' // TROP_LETTER(alg) // '_PD',  tropPdf(:ncol, :), ncol, lchnk)
        ! call outfld('TROP' // TROP_LETTER(alg) // '_FD',  tropFound(:ncol),  ncol, lchnk)
      end do
    end if

    return
  end subroutine tropopause_output
end module tropopause_output