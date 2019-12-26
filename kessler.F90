module kessler_mod

  use ccpp_kinds, only: kind_phys

  implicit none


contains

  !> \section arg_table_kessler_run  Argument Table
  !! \htmlinclude arg_table_kessler_run.html

  subroutine kessler_run()

    use ppgrid,            only: pcols, pver, pverp
    use physics_types_cam7,only: physics_state

    use physics_types_cam7, only: state, tend, physics_type_alloc
    use constituents,       only: pcnst, ix_qv, ix_qc, ix_qr

    use cam7_kessler_ccpp_cap,     only: cam7_kessler_ccpp_physics_initialize
    use cam7_kessler_ccpp_cap,     only: cam7_kessler_ccpp_physics_timestep_initial
    use cam7_kessler_ccpp_cap,     only: cam7_kessler_ccpp_physics_run
    use cam7_kessler_ccpp_cap,     only: cam7_kessler_ccpp_physics_timestep_final
    use cam7_kessler_ccpp_cap,     only: cam7_kessler_ccpp_physics_finalize
    use cam7_kessler_ccpp_cap,     only: ccpp_physics_suite_list
    use cam7_kessler_ccpp_cap,     only: ccpp_physics_suite_part_list

    integer,            parameter   :: begchunk = 33 ! Not needed for CAM7
    integer,            parameter   :: endchunk = 33 ! Not needed for CAM7
    integer,            parameter   :: ncols = 1
    integer,            parameter   :: ntimes = 3

    integer                         :: i, j, k, rk
    integer                         :: ierr
    integer                         :: col_start, col_end
    integer                         :: ncol, nwrite, pver_in, nwrite_in, nstep
    real(kind_phys)                 :: ztodt
    real(kind_phys)                 :: precl(pcols)
    real(kind_phys)                 :: scratch(pcols,pver)
    character(len=20)               :: string
    character(len=512)              :: errmsg
    character(len=128), allocatable :: part_names(:)
    integer                         :: errflg
    real(kind_phys)                 :: pmiddry_top2bot(pcols,pver)
    real(kind_phys)                 :: pint_top2bot(pcols,pverp)
    real(kind_phys)                 :: pmid_top2bot(pcols,pver)
    real(kind_phys)                 :: pdel_top2bot(pcols,pver)
    real(kind_phys)                 :: rpdel_top2bot(pcols,pver)
    real(kind_phys)                 :: pdeldry_top2bot(pcols,pver)
    real(kind_phys)                 :: zi_top2bot(pcols,pverp)
    real(kind_phys)                 :: zm_top2bot(pcols,pver)
    real(kind_phys)                 :: exner_top2bot(pcols,pver)
    real(kind_phys)                 :: t_top2bot(pcols,pver)
    real(kind_phys)                 :: s_top2bot(pcols,pver)
    real(kind_phys)                 :: q_top2bot(pcols,pver,pcnst)
    real(kind_phys)                 :: lnpint_top2bot(pcols,pverp)
    real(kind_phys)                 :: lnpmid_top2bot(pcols,pver)
    real(kind_phys)                 :: ttend_top2bot(pcols,pver)

    ! Allocate the host variables
    call physics_type_alloc(state, tend, pcols)

    ! Use the suite information to setup the run
    call cam7_kessler_ccpp_physics_initialize('cam_kessler_test',  &
           state, tend, precl, ztodt, errmsg, errflg)
    if (errflg /= 0) then
       write(6, *) trim(errmsg)
       stop
    end if

    ! loop over all time steps
    nstep=0
    open (60, file='fort.60')
    do j = 1, ntimes

       ncol = pcols
       read(60,fmt='(a10,i4)') string(1:8),nwrite_in
       read(60,fmt='(a20,2i4,f20.13)') string(1:19),ncol, pver_in, ztodt
       read(60,fmt='(a20,(e25.18))') string(1:20),pmiddry_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e25.18))') string(1:20),pint_top2bot(:ncol,:pverp)
       read(60,fmt='(a20,(e25.18))') string(1:20),pmid_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e25.18))') string,pdel_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e25.18))') string,rpdel_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e25.18))') string,pdeldry_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e25.18))') string,exner_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e25.18))') string,t_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e25.18))') string,s_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e25.18))') string,q_top2bot(:ncol,:pver_in,1)
       read(60,fmt='(a20,(e25.18))') string,q_top2bot(:ncol,:pver_in,2)
       read(60,fmt='(a20,(e25.18))') string,q_top2bot(:ncol,:pver_in,3)
       read(60,fmt='(a20,(e25.18))') string,zi_top2bot(:ncol,:pverp)
       read(60,fmt='(a20,(e25.18))') string,zm_top2bot(:ncol,:pver)
       read(60,fmt='(a20,(e25.18))') string,lnpint_top2bot(:ncol,:pverp)
       read(60,fmt='(a20,(e25.18))') string,lnpmid_top2bot(:ncol,:pver)
       read(60,fmt='(a20,(e25.18))') string,precl(:ncol)
       read(60,fmt='(a20,(e25.18))') string,state%phis(:ncol)
       read(60,fmt='(a20,(e25.18))') string,ttend_top2bot(:ncol,:pver_in)

       ! Need to swap the bottom and top
       do k=1,pver
         rk= pver - k +1 
         state%pmiddry(:ncol,rk) = pmiddry_top2bot(:ncol,k)
         state%pmid(:ncol,rk)    = pmid_top2bot(:ncol,k)
         state%pdel(:ncol,rk)    = pdel_top2bot(:ncol,k)
         state%rpdel(:ncol,rk)   = rpdel_top2bot(:ncol,k)
         state%pdeldry(:ncol,rk) = pdeldry_top2bot(:ncol,k)
         state%exner(:ncol,rk)   = exner_top2bot(:ncol,k)
         state%t(:ncol,rk)       = t_top2bot(:ncol,k)
         state%s(:ncol,rk)       = s_top2bot(:ncol,k)
         state%zm(:ncol,rk)      = zm_top2bot(:ncol,k)
         state%lnpmid(:ncol,rk)  = lnpmid_top2bot(:ncol,k)
         state%q(:ncol,rk,1)     = q_top2bot(:ncol,k,1)
         state%q(:ncol,rk,2)     = q_top2bot(:ncol,k,2)
         state%q(:ncol,rk,3)     = q_top2bot(:ncol,k,3)
       end do

       do k=1,pverp
         rk= pverp - k +1 
         state%pint(:ncol,rk)    = pint_top2bot(:ncol,k)
         state%lnpint(:ncol,rk)  = lnpint_top2bot(:ncol,k)
         state%zi(:ncol,rk)      = zi_top2bot(:ncol,k)
       end do

       ! Initialize the timestep
!       call cam7_kessler_ccpp_physics_timestep_initial('cam_kessler_test', col_start, col_end, &
!           ncol, state, tend, precl, ztodt, errmsg, errflg)
       call cam7_kessler_ccpp_physics_timestep_initial('cam_kessler_test', &
            state, tend, precl, ztodt, errmsg, errflg)
       col_start = 1
       col_end = ncol

       ! Initialize the total tendency after the timestep initialization
       do k=1,pver
         rk= pver - k +1 
         tend%dtdt(:ncol,rk)     = ttend_top2bot(:ncol,k)
       end do

       call cam7_kessler_ccpp_physics_run('cam_kessler_test', 'physics', col_start, col_end, &
           ncol, state, tend, precl, ztodt, errmsg, errflg)
       if (errflg /= 0) then
          write(6, *) trim(errmsg)
          call ccpp_physics_suite_part_list('cam_kessler_test', part_names, errmsg, errflg)
          write(6, *) 'Available suite parts are:'
          do nwrite = 1, size(part_names)
             write(6, *) trim(part_names(nwrite))
          end do
          stop
       end if

       write(6,*) 'At time step', j, 'in host model Temperature =', state%T(8, :pver)


       call cam7_kessler_ccpp_physics_timestep_final('cam_kessler_test', &
           state, tend, precl, ztodt, errmsg, errflg)

         write(61,'(a10,i4)') 'nstep=',nstep
         write(61,'(a20,2i4,f20.13)') 'ncol, pver, ztodt=',ncol, pver, ztodt
         write(61,fmt='(a10,(e25.18))') 'pmiddry=',state%pmiddry(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e25.18))') 'pint=',state%pint(:ncol,pverp:1:-1)
         write(61,fmt='(a10,(e25.18))') 'pmid=',state%pmid(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e25.18))') 'pdel=',state%pdel(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e25.18))') 'rpdel=',state%rpdel(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e25.18))') 'pdeldry=',state%pdeldry(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e25.18))') 'exner=',state%exner(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e25.18))') 'state%t=',state%t(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e25.18))') 'state%s=',state%s(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e22.15))') 'qv=',state%q(:ncol,pver:1:-1,1)
         write(61,fmt='(a10,(e25.18))') 'qc=',state%q(:ncol,pver:1:-1,2)
         write(61,fmt='(a10,(e25.18))') 'qr=',state%q(:ncol,pver:1:-1,3)
         write(61,fmt='(a20,(e25.18))') 'zi=',state%zi(:ncol,pverp:1:-1)
         write(61,fmt='(a20,(e25.18))') 'zm=',state%zm(:ncol,pver:1:-1)
         write(61,fmt='(a20,(e25.18))') 'lnpint=',state%lnpint(:ncol,pverp:1:-1)
         write(61,fmt='(a20,(e25.18))') 'lnpmid=',state%lnpmid(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e25.18))') 'precl=',precl(:ncol)
         write(61,fmt='(a10,(e25.18))') 'phis=',state%phis(:ncol)
         write(61,fmt='(a10,(e25.18))') 'tend%dtdt=',tend%dtdt(:ncol,pver:1:-1)

       nstep = nstep + 1

    end do


    call cam7_kessler_ccpp_physics_finalize('cam_kessler_test', &
           state, tend, precl, ztodt, errmsg, errflg)
    if (errflg /= 0) then
       write(6, *) trim(errmsg)
       write(6,'(a)') 'An error occurred in ccpp_timestep_final, Exiting...'
       stop
    end if

  end subroutine kessler_run

end module kessler_mod
