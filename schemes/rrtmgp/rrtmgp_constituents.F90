module rrtmgp_constituents

   implicit none
   private

   public :: rrtmgp_constituents_init
   public :: rrtmgp_constituents_run

   ! Constituent array index for each gas in the radiatively active gas list,
   ! resolved once during initialization. Gases left unresolved keep
   ! int_unassigned and are given zero concentrations in radiation.
   integer, allocatable :: rad_gas_indices(:)

contains

!> \section arg_table_rrtmgp_constituents_init Argument Table
!! \htmlinclude rrtmgp_constituents_init.html
!!
   subroutine rrtmgp_constituents_init(amIRoot, iulog, rad_climate, gaslist, const_props, errmsg, errflg)
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t, int_unassigned
      use ccpp_scheme_utils,         only: ccpp_constituent_index
      use radiation_utils,           only: parse_rad_climate_entry

      logical,            intent(in)  :: amIRoot         ! are we on the MPI root task?
      integer,            intent(in)  :: iulog           ! log output unit
      character(len=256), intent(in)  :: rad_climate(:)  ! (namelist) list of radiatively active gases and sources
      character(len=*),   intent(in)  :: gaslist(:)      ! Radiatively active gas list
      type(ccpp_constituent_prop_ptr_t), intent(in) :: const_props(:) ! Constituent properties
      character(len=*),   intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables
      character(len=8)   :: source
      character(len=32)  :: identifier
      character(len=32)  :: gas_name
      character(len=256) :: diag_name
      character(len=256) :: alloc_errmsg
      logical            :: entry_found
      logical            :: is_advected
      integer            :: gas_idx, entry_idx, const_idx, ierr

      errmsg = ''
      errflg = 0

      ! Initialization should be idempotent:
      if (allocated(rad_gas_indices)) then
         return
      end if

      allocate(rad_gas_indices(size(gaslist)), stat=ierr, errmsg=alloc_errmsg)
      if (ierr /= 0) then
         write(errmsg, *) 'rrtmgp_constituents_init: Unable to allocate rad_gas_indices - message: ', trim(alloc_errmsg)
         errflg = 1
         return
      end if
      rad_gas_indices = int_unassigned

      ! Map each radiatively active gas to the constituent that provides it.
      !
      ! H2O is always taken from the water vapor constituent (it needs no rad_climate entry);
      ! every other gas must have a rad_climate entry ("flag:identifier:gas_name"):
      ! 'A'/'N' entries are resolved by matching the identifier against the constituent
      ! diagnostic names (unresolved is fatal)
      ! 'Z' entries give the gas zero concentration in radiation.
      !
      ! A gas with no entry at all is an error.
      gas_loop: do gas_idx = 1, size(gaslist)

         if (trim(gaslist(gas_idx)) == 'H2O') then
            call ccpp_constituent_index('water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water', &
                 rad_gas_indices(gas_idx), errflg, errmsg)
            if (errflg /= 0) then
               return
            end if
            if (rad_gas_indices(gas_idx) == int_unassigned) then
               write(errmsg, *) 'rrtmgp_constituents_init: water vapor constituent', &
                    ' (water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water) not found for radiatively active gas "H2O"'
               errflg = 1
               return
            end if
            ! Tolerated for CAM namelist compatibility; the water vapor constituent is always used for H2O.
            if (amIRoot) then
               write(iulog, *) 'rrtmgp_constituents_init: H2O is taken from the water vapor constituent;', &
                    ' any rad_climate entry for H2O is ignored'
            end if
            cycle gas_loop
         end if

         entry_found = .false.
         entry_loop: do entry_idx = 1, size(rad_climate)
            if (len_trim(rad_climate(entry_idx)) == 0) then
               exit entry_loop
            end if

            call parse_rad_climate_entry(rad_climate(entry_idx), source, identifier, gas_name, errmsg, errflg)
            if (errflg /= 0) then
               return
            end if

            if (trim(gas_name) /= trim(gaslist(gas_idx))) then
               cycle entry_loop
            end if
            entry_found = .true.

            if (source == 'Z') then
               ! Zero concentration: leave the gas unresolved; the run phase
               ! fills unresolved gases with zero.
               if (amIRoot) then
                  write(iulog, *) 'rrtmgp_constituents_init: gas "', trim(gaslist(gas_idx)), &
                       '" has a "Z" rad_climate entry; it will have zero concentration in radiation'
               end if
            else if (source == 'A' .or. source == 'N') then
               ! Find the constituent whose diagnostic name matches the identifier.
               const_loop: do const_idx = 1, size(const_props)
                  call const_props(const_idx)%diagnostic_name(diag_name, errcode=errflg, errmsg=errmsg)
                  if (errflg /= 0) then
                     return
                  end if
                  if (trim(diag_name) == trim(identifier)) then
                     call const_props(const_idx)%const_index(rad_gas_indices(gas_idx), errflg, errmsg)
                     if (errflg /= 0) then
                        return
                     end if

                     ! The A/N flag does not affect behavior (everything in SIMA is a constituent)
                     ! but warn if the flag is inconsistent w/ the constituent that actually provides the gas.
                     call const_props(const_idx)%is_advected(is_advected, errflg, errmsg)
                     if (errflg /= 0) then
                        return
                     end if
                     if (amIRoot .and. (is_advected .neqv. (source == 'A'))) then
                        write(iulog, *) 'rrtmgp_constituents_init: WARNING: rad_climate flag "', trim(source), &
                             '" for gas "', trim(gaslist(gas_idx)), '" does not match constituent "', &
                             trim(identifier), '" (advected: ', is_advected, ')'
                     end if
                     exit const_loop
                  end if
               end do const_loop
               if (rad_gas_indices(gas_idx) == int_unassigned) then
                  write(errmsg, *) 'rrtmgp_constituents_init: no constituent with diagnostic name "', &
                       trim(identifier), '" found for radiatively active gas "', trim(gaslist(gas_idx)), &
                       '" (rad_climate entry "', trim(rad_climate(entry_idx)), '")'
                  errflg = 1
                  return
               end if
            else
               write(errmsg, *) 'rrtmgp_constituents_init: invalid gas source flag "', trim(source), &
                    '" in rad_climate entry "', trim(rad_climate(entry_idx)), '" (must be A, N, or Z)'
               errflg = 1
               return
            end if
            exit entry_loop
         end do entry_loop

         if (.not. entry_found) then
            write(errmsg, *) 'rrtmgp_constituents_init: radiatively active gas "', trim(gaslist(gas_idx)), &
                 '" has no rad_climate entry; all radiation gases must be specified', &
                 ' (use a "Z" entry to give a gas zero concentration)'
            errflg = 1
            return
         end if

      end do gas_loop

   end subroutine rrtmgp_constituents_init

!> \section arg_table_rrtmgp_constituents_run Argument Table
!! \htmlinclude rrtmgp_constituents_run.html
!!
   subroutine rrtmgp_constituents_run(const_array, rad_const_array, errmsg, errflg)
       use ccpp_constituent_prop_mod, only: int_unassigned
       use ccpp_kinds,                only: kind_phys
       real(kind_phys),           intent(in) :: const_array(:,:,:)     ! Constituents array
       real(kind_phys),          intent(out) :: rad_const_array(:,:,:) ! Radiatively active constituent mixing ratios
       integer,                  intent(out) :: errflg
       character(len=*),         intent(out) :: errmsg

       ! Local variables
       integer :: gas_idx

       errflg = 0
       errmsg = ''

       rad_const_array = 0._kind_phys

       do gas_idx = 1, size(rad_gas_indices)
          if (rad_gas_indices(gas_idx) /= int_unassigned) then
             rad_const_array(:,:,gas_idx) = const_array(:,:,rad_gas_indices(gas_idx))
          end if
       end do

   end subroutine rrtmgp_constituents_run

end module rrtmgp_constituents
