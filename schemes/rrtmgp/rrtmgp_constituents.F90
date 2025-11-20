module rrtmgp_constituents

   public :: rrtmgp_constituents_register
   public :: rrtmgp_constituents_run

contains

!> \section arg_table_rrtmgp_constituents_register Argument Table
!! \htmlinclude rrtmgp_constituents_register.html
!!
   subroutine rrtmgp_constituents_register(rad_climate, rrtmgp_dyn_consts, errmsg, errflg)
      use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
      use ccpp_kinds,                only: kind_phys
      type(ccpp_constituent_properties_t), allocatable, intent(out) :: rrtmgp_dyn_consts(:) ! Runtime constituent properties
      character(len=256), intent(in)  :: rad_climate(:)  ! (namelist) list of radiatively active gases and sources
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables
      character(len=1)  :: source
      character(len=32) :: long_name
      character(len=32) :: stdname
      character(len=256) :: tmpstr, alloc_errmsg
      integer            :: gas_idx, strlen, ipos, ierr, idx

      errmsg = ''
      errflg = 0

      ! Allocate the dynamic constituents array
      allocate(rrtmgp_dyn_consts(size(rad_climate)), stat=ierr, errmsg=alloc_errmsg)
      if (ierr /= 0) then
         write(errmsg, *) 'rrtmgp_constituents_register: Unable to allocate rrtmgp_dyn_consts - message: ', alloc_errmsg
         errflg = 1
         return
      end if

      ! Parse gases, long names, and sources from rad_climate
      parse_loop: do gas_idx = 1, size(rad_climate)
      if ( len_trim(rad_climate(gas_idx)) == 0 ) then
         exit parse_loop
      endif

      ! There are no fields in the input strings in which a blank character is allowed.
      ! To simplify the parsing go through the input strings and remove blanks.
      tmpstr = adjustl(rad_climate(gas_idx))
      do
         strlen = len_trim(tmpstr)
         ipos = index(tmpstr, ' ')
         if (ipos == 0 .or. ipos > strlen) exit
         tmpstr = tmpstr(:ipos-1) // tmpstr(ipos+1:strlen)
      end do

      ! Locate the ':' separating source from long name.
      idx = index(tmpstr, ':')
      source = tmpstr(:idx-1)
      tmpstr = tmpstr(idx+1:)

      ! locate the ':' separating long name from rad gas ("standard") name
      idx = scan(tmpstr, ':')

      long_name = tmpstr(:idx-1)
      stdname = tmpstr(idx+1:)

      ! Register the constituent based on the source
      if (source == 'A') then
          call rrtmgp_dyn_consts(gas_idx)%instantiate(     &
             std_name = stdname,   &
             long_name = long_name,  &
             units = 'kg kg-1',                            &
             vertical_dim = 'vertical_layer_dimension', &
             min_value = 0.0_kind_phys,                 &
             advected = .true.,                         &
             water_species = .false.,                    &
             mixing_ratio_type = 'dry',                 &
             errcode = errflg,                         &
             errmsg = errmsg)
      else if (source == 'N') then
          call rrtmgp_dyn_consts(gas_idx)%instantiate(     &
             std_name = stdname,   &
             long_name = long_name,  &
             units = 'kg kg-1',                            &
             vertical_dim = 'vertical_layer_dimension', &
             min_value = 0.0_kind_phys,                 &
             advected = .false.,                         &
             water_species = .false.,                    &
             mixing_ratio_type = 'dry',                 &
             errcode = errflg,                         &
             errmsg = errmsg)
      else if (source == 'Z') then
          call rrtmgp_dyn_consts(gas_idx)%instantiate(     &
             std_name = stdname,   &
             long_name = long_name,  &
             units = 'kg kg-1',                            &
             vertical_dim = 'vertical_layer_dimension', &
             min_value = 0.0_kind_phys,                 &
             default_value = 0.0_kind_phys,             &
             advected = .false.,                         &
             water_species = .false.,                    &
             mixing_ratio_type = 'dry',                 &
             errcode = errflg,                         &
             errmsg = errmsg)
      else
         write(errmsg,*) 'rrtmgp_constituent_register: invalid gas source "', source, '" for radiation', &
            ' constituent "', stdname, '"'
         errflg = 1
         return
      end if

      end do parse_loop

   end subroutine rrtmgp_constituents_register

!> \section arg_table_rrtmgp_constituents_run Argument Table
!! \htmlinclude rrtmgp_constituents_run.html
!!
   subroutine rrtmgp_constituents_run(gaslist, const_array, rad_const_array, errmsg, errflg)
       use ccpp_constituent_prop_mod, only: int_unassigned
       use ccpp_scheme_utils,         only: ccpp_constituent_index
       use ccpp_kinds,                only: kind_phys
       character(len=5),          intent(in) :: gaslist(:)             ! Radiatively active gas list
       real(kind_phys),           intent(in) :: const_array(:,:,:)     ! Constituents array
       real(kind_phys),          intent(out) :: rad_const_array(:,:,:) ! Radiatively active constituent mixing ratios
       integer,                  intent(out) :: errflg
       character(len=512),       intent(out) :: errmsg

       ! Local variables
       integer :: gas_idx
       integer :: const_idx

       errflg = 0
       errmsg = ''

       rad_const_array = 0._kind_phys

       do gas_idx = 1, size(gaslist)
          ! Find the index of the current gas in the constituents array
          if (trim(gaslist(gas_idx)) == 'H2O') then
             call ccpp_constituent_index('water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water', const_idx, errflg, errmsg)
          else
             call ccpp_constituent_index(trim(gaslist(gas_idx)), const_idx, errflg, errmsg)
          end if
          if (errflg /= 0) then
             return
          end if
          if (const_idx /= int_unassigned) then
             rad_const_array(:,:,gas_idx) = const_array(:,:,const_idx)
          end if
       end do

   end subroutine rrtmgp_constituents_run

end module rrtmgp_constituents
