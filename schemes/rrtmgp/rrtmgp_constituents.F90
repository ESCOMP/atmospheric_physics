module rrtmgp_constituents

   public :: rrtmgp_constituents_register

contains

!> \section arg_table_rrtmgp_constituents_register Argument Table
!! \htmlinclude rrtmgp_constituents_register.html
!!
   subroutine rrtmgp_constituents_register(rad_climate, rrtmgp_dyn_consts, errmsg, errcode)
      use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
      type(ccpp_constituent_properties_t), allocatable, intent(out) :: rrtmgp_dyn_consts(:)
      character(len=256), intent(in)  :: rad_climate(:)
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errcode

      ! Local variables
      character(len=1)  :: source
      character(len=32) :: long_name
      character(len=32) :: stdname
      character(len=256) :: tmpstr, alloc_errmsg
      integer            :: gas_idx, strlen, ipos, ierr

      errmsg = ''
      errcode = 0

      ! Allocate the dynamic constituents array
      allocate(rrtmgp_dyn_consts(size(rad_climate)), stat=ierr, errmsg=alloc_errmsg)
      if (ierr /= 0) then
         write(errmsg, *) 'rrtmgp_constituents_register: Unable to allocate rrtmgp_dyn_consts - message: ', alloc_errmsg
         errcode = 1
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
      source = tmpstr(:jdx-1)
      tmpstr = tmpstr(jdx+1:)

      ! locate the ':' separating long name from rad gas ("standard") name
      idx = scan(tmpstr, ':')

      long_name = tmpstr(:jdx-1)
      stdname = tmpstr(jdx+1:)

      ! Register the constituent based on the source
      if (source == 'A') then
          call rrtmgp_dyn_consts(gas_idx)%instantiate(     &
             std_name = stdname,   &
             long_name = long_name,  &
             units = 'kg-1',                            &
             vertical_dim = 'vertical_layer_dimension', &
             min_value = 0.0_kind_phys,                 &
             advected = .true.,                         &
             water_species = .false.,                    &
             mixing_ratio_type = 'dry',                 &
             errcode = errcode,                         &
             errmsg = errmsg)
      else if (source == 'N') then
          call rrtmgp_dyn_consts(gas_idx)%instantiate(     &
             std_name = stdname,   &
             long_name = long_name,  &
             units = 'kg-1',                            &
             vertical_dim = 'vertical_layer_dimension', &
             min_value = 0.0_kind_phys,                 &
             advected = .false.,                         &
             water_species = .false.,                    &
             mixing_ratio_type = 'dry',                 &
             errcode = errcode,                         &
             errmsg = errmsg)
      else if (source == 'Z') then
          call rrtmgp_dyn_consts(gas_idx)%instantiate(     &
             std_name = stdname,   &
             long_name = long_name,  &
             units = 'kg-1',                            &
             vertical_dim = 'vertical_layer_dimension', &
             min_value = 0.0_kind_phys,                 &
             default_value = 0.0_kind_phys,             &
             advected = .false.,                         &
             water_species = .false.,                    &
             mixing_ratio_type = 'dry',                 &
             errcode = errcode,                         &
             errmsg = errmsg)
      else
         write(errmsg,*) 'rrtmgp_constituent_register: invalid gas source "', source, '" for radiation', &
            ' constituent "', stdname, '"'
         errcode = 1
         return
      end if

      end do parse_loop

   end subroutine rrtmgp_constituents_register

!> \section arg_table_rrtmgp_constituents_init Argument Table
!! \htmlinclude rrtmgp_constituents_init.html
!!
   subroutine rrtmgp_constituents_init(gaslist, errmsg, errcode)
      character(len=5),   intent(out) :: gaslist(:)
      integer,            intent(out) :: errcode
      character(len=512), intent(out) :: errmsg

      errcode = 0
      errmsg = ''

      gaslist =  (/'H2O  ','O3   ', 'O2   ', 'CO2  ', 'N2O  ', 'CH4  ', 'CFC11', 'CFC12'/)

   end subroutine rrtmgp_constituents_init

end module rrtmgp_constituents
