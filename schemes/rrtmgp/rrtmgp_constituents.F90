module rrtmgp_constituents

   public :: rrtmgp_constituents_register

contains

!> \section arg_table_rrtmgp_constituents_register Argument Table
!! \htmlinclude rrtmgp_constituents_register.html
!!
   subroutine rrtmgp_constituents_register(nradgas, rad_climate, rrtmgp_dyn_consts, errmsg, errcode)
      use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
      use ccpp_kinds,                only: kind_phys
      integer,            intent(in)  :: nradgas
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
      allocate(rrtmgp_dyn_consts(nradgas), stat=ierr, errmsg=alloc_errmsg)
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
!! \htmlinclude rrtmgp_constituents_int.html
!!
   subroutine rrtmgp_constituents_init(ndiag, ncol, unset_real, diag_cur, active_call_array, &
      rrtmgp_phys_blksz, tlev, fluxlwup_Jac, rad_heat, fsnt, fsns, is_first_restart_step,    &
      use_tlev, top_at_one, errmsg, errcode)
      use ccpp_kinds, only: kind_phys
      integer,             intent(in) :: ndiag
      integer,             intent(in) :: ncol
      real(kind_phys),     intent(in) :: unset_real
      integer,            intent(out) :: diag_cur
      logical,            intent(out) :: active_call_array(:)
      integer,            intent(out) :: rrtmgp_phys_blksz
      real(kind_phys),    intent(out) :: tlev(:,:)
      real(kind_phys),    intent(out) :: fluxlwup_Jac(:,:)
      real(kind_phys),    intent(out) :: rad_heat(:,:)
      real(kind_phys),    intent(out) :: fsnt(:)
      real(kind_phys),    intent(out) :: fsns(:)
      logical,            intent(out) :: is_first_restart_step
      logical,            intent(out) :: use_tlev
      logical,            intent(out) :: top_at_one
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errcode

      errcode = 1
      errmsg = ''

      active_call_array = .true.
      is_first_restart_step = .false.
      top_at_one = .true.

      diag_cur = 1
      rrtmgp_phys_blksz = ncol
      ! Set tlev & fluxlwup_Jac to unset values; not used by default in CAM-SIMA
      use_tlev = .false.
      tlev = unset_real
      fluxlwup_Jac = unset_real
      rad_heat = unset_real

      ! PEVERWHEE - remove when shortwave is done
      fsnt = 0.0_kind_phys
      fsns = 0.0_kind_phys

   end subroutine rrtmgp_constituents_init
!> \section arg_table_rrtmgp_constituents_run Argument Table
!! \htmlinclude rrtmgp_constituents_run.html
!!
   subroutine rrtmgp_constituents_run(gaslist, const_array, rad_const_array, errmsg, errcode)
       use ccpp_constituent_prop_mod, only: int_unassigned
       use ccpp_scheme_utils,         only: ccpp_constituent_index
       use ccpp_kinds,                only: kind_phys
       character(len=5),          intent(in) :: gaslist(:)
       real(kind_phys),           intent(in) :: const_array(:,:,:)
       real(kind_phys),          intent(out) :: rad_const_array(:,:,:)
       integer,                  intent(out) :: errcode
       character(len=512),       intent(out) :: errmsg

       ! Local variables
       integer :: gas_idx
       integer :: const_idx

       errcode = 0
       errmsg = ''

       rad_const_array = 0._kind_phys

       do gas_idx = 1, size(gaslist)
          ! Find the index of the current gas in the constituents array
          call ccpp_constituent_index(trim(gaslist(gas_idx)), const_idx, errcode, errmsg)
          if (errcode /= 0) then
             return
          end if
          if (const_idx /= int_unassigned) then
             rad_const_array(:,:,gas_idx) = const_array(:,:,const_idx)
          end if
       end do

   end subroutine rrtmgp_constituents_run

end module rrtmgp_constituents
