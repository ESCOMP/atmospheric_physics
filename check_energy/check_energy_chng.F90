module check_energy_chng

	use ccpp_kinds, only: kind_phys

	implicit none
	private

	public :: check_energy_chng_run

contains

  subroutine check_energy_chng_run( &
  	   state, tend, name, nstep, ztodt,        &
       flx_vap, flx_cnd, flx_ice, flx_sen)
    use cam_thermo,      only: get_hydrostatic_energy
    use dyn_tests_utils, only: vc_physics, vc_dycore, vc_height, vc_dry_pressure
    use cam_abortutils,  only: endrun
    use physics_types,   only: phys_te_idx, dyn_te_idx
!-----------------------------------------------------------------------
! Check that the energy and water change matches the boundary fluxes
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state)    , intent(inout) :: state
    type(physics_tend )    , intent(inout) :: tend
    character*(*),intent(in) :: name               ! parameterization name for fluxes
    integer , intent(in   ) :: nstep               ! current timestep number
    real(r8), intent(in   ) :: ztodt               ! 2 delta t (model time increment)
    real(r8), intent(in   ) :: flx_vap(:)          ! (pcols) - boundary flux of vapor         (kg/m2/s)
    real(r8), intent(in   ) :: flx_cnd(:)          ! (pcols) -boundary flux of liquid+ice    (m/s) (precip?)
    real(r8), intent(in   ) :: flx_ice(:)          ! (pcols) -boundary flux of ice           (m/s) (snow?)
    real(r8), intent(in   ) :: flx_sen(:)          ! (pcols) -boundary flux of sensible heat (w/m2)

!******************** BAB ******************************************************
!******* Note that the precip and ice fluxes are in precip units (m/s). ********
!******* I would prefer to have kg/m2/s.                                ********
!******* I would also prefer liquid (not total) and ice fluxes          ********
!*******************************************************************************

!---------------------------Local storage-------------------------------

    real(r8) :: te_xpd(state%ncol)                 ! expected value (f0 + dt*boundary_flux)
    real(r8) :: te_dif(state%ncol)                 ! energy of input state - original energy
    real(r8) :: te_tnd(state%ncol)                 ! tendency from last process
    real(r8) :: te_rer(state%ncol)                 ! relative error in energy column

    real(r8) :: tw_xpd(state%ncol)                 ! expected value (w0 + dt*boundary_flux)
    real(r8) :: tw_dif(state%ncol)                 ! tw_inp - original water
    real(r8) :: tw_tnd(state%ncol)                 ! tendency from last process
    real(r8) :: tw_rer(state%ncol)                 ! relative error in water column

    real(r8) :: te(state%ncol)                     ! vertical integral of total energy
    real(r8) :: tw(state%ncol)                     ! vertical integral of total water
    real(r8) :: cp_or_cv(state%psetcols,pver)      ! cp or cv depending on vcoord
    real(r8) :: scaling(state%psetcols,pver)       ! scaling for conversion of temperature increment
    real(r8) :: temp(state%ncol,pver)              ! temperature

    real(r8) :: se(state%ncol)                     ! enthalpy or internal energy (J/m2)
    real(r8) :: po(state%ncol)                     ! surface potential or potential energy (J/m2)
    real(r8) :: ke(state%ncol)                     ! kinetic energy    (J/m2)
    real(r8) :: wv(state%ncol)                     ! column integrated vapor       (kg/m2)
    real(r8) :: liq(state%ncol)                    ! column integrated liquid      (kg/m2)
    real(r8) :: ice(state%ncol)                    ! column integrated ice         (kg/m2)

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i                                     ! column index
!-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol

    ! If psetcols == pcols, cpairv is the correct size and just copy into cp_or_cv
    ! If psetcols > pcols and all cpairv match cpair, then assign the constant cpair

    if (state%psetcols == pcols) then
       cp_or_cv(:,:) = cpairv(:,:,lchnk)
    else if (state%psetcols > pcols .and. all(cpairv(:,:,:) == cpair)) then
       cp_or_cv(:,:) = cpair
    else
       call endrun('check_energy_chng: cpairv is not allowed to vary when subcolumns are turned on')
    end if

    call get_hydrostatic_energy(state%q(1:ncol,1:pver,1:pcnst),.true.,               &
         state%pdel(1:ncol,1:pver), cp_or_cv(1:ncol,1:pver),                         &
         state%u(1:ncol,1:pver), state%v(1:ncol,1:pver), state%T(1:ncol,1:pver),     &
         vc_physics, ptop=state%pintdry(1:ncol,1), phis = state%phis(1:ncol),        &
         te = te(1:ncol), H2O = tw(1:ncol), se=se(1:ncol),po=po(1:ncol),             &
         ke=ke(1:ncol),wv=wv(1:ncol),liq=liq(1:ncol),ice=ice(1:ncol))
    ! compute expected values and tendencies
    do i = 1, ncol
       ! change in static energy and total water
       te_dif(i) = te(i) - state%te_cur(i,phys_te_idx)
       tw_dif(i) = tw(i) - state%tw_cur(i,phys_te_idx)

       ! expected tendencies from boundary fluxes for last process
       te_tnd(i) = flx_vap(i)*(latvap+latice) - (flx_cnd(i) - flx_ice(i))*1000._r8*latice + flx_sen(i)
       tw_tnd(i) = flx_vap(i) - flx_cnd(i) *1000._r8

       ! cummulative tendencies from boundary fluxes
       tend%te_tnd(i) = tend%te_tnd(i) + te_tnd(i)
       tend%tw_tnd(i) = tend%tw_tnd(i) + tw_tnd(i)

       ! expected new values from previous state plus boundary fluxes
       te_xpd(i) = state%te_cur(i,phys_te_idx) + te_tnd(i)*ztodt
       tw_xpd(i) = state%tw_cur(i,phys_te_idx) + tw_tnd(i)*ztodt

       ! relative error, expected value - input state / previous state
       te_rer(i) = (te_xpd(i) - te(i)) / state%te_cur(i,phys_te_idx)
    end do

    ! relative error for total water (allow for dry atmosphere)
    tw_rer = 0._r8
    where (state%tw_cur(:ncol,phys_te_idx) > 0._r8)
       tw_rer(:ncol) = (tw_xpd(:ncol) - tw(:ncol)) / state%tw_cur(:ncol,1)
    end where

    ! error checking
    if (print_energy_errors) then
       if (any(abs(te_rer(1:ncol)) > 1.E-14_r8 .or. abs(tw_rer(1:ncol)) > 1.E-10_r8)) then
          do i = 1, ncol
             ! the relative error threshold for the water budget has been reduced to 1.e-10
             ! to avoid messages generated by QNEG3 calls
             ! PJR- change to identify if error in energy or water
             if (abs(te_rer(i)) > 1.E-14_r8 ) then
                state%count = state%count + 1
                write(iulog,*) "significant energy conservation error after ", name,        &
                      " count", state%count, " nstep", nstep, "chunk", lchnk, "col", i
                write(iulog,*) te(i),te_xpd(i),te_dif(i),tend%te_tnd(i)*ztodt,  &
                      te_tnd(i)*ztodt,te_rer(i)
             endif
             if ( abs(tw_rer(i)) > 1.E-10_r8) then
                state%count = state%count + 1
                write(iulog,*) "significant water conservation error after ", name,        &
                      " count", state%count, " nstep", nstep, "chunk", lchnk, "col", i
                write(iulog,*) tw(i),tw_xpd(i),tw_dif(i),tend%tw_tnd(i)*ztodt,  &
                      tw_tnd(i)*ztodt,tw_rer(i)
             end if
          end do
       end if
    end if

    ! copy new value to state

    do i = 1, ncol
      state%te_cur(i,phys_te_idx) = te(i)
      state%tw_cur(i,phys_te_idx) = tw(i)
    end do

    !
    ! Dynamical core total energy
    !
    if (vc_dycore == vc_height) then
      !
      ! compute cv if vertical coordinate is height: cv = cp - R
      !
      ! Note: cp_or_cv set above for pressure coordinate
      if (state%psetcols == pcols) then
        cp_or_cv(:ncol,:) = cp_or_cv_dycore(:ncol,:,lchnk)
      else
        cp_or_cv(:ncol,:) = cpair-rair
      endif
      scaling(:,:)   = cpairv(:,:,lchnk)/cp_or_cv(:,:) !cp/cv scaling
      temp(1:ncol,:) = state%temp_ini(1:ncol,:)+scaling(1:ncol,:)*(state%T(1:ncol,:)-state%temp_ini(1:ncol,:))
      call get_hydrostatic_energy(state%q(1:ncol,1:pver,1:pcnst),.true.,               &
           state%pdel(1:ncol,1:pver), cp_or_cv(1:ncol,1:pver),                         &
           state%u(1:ncol,1:pver), state%v(1:ncol,1:pver), temp(1:ncol,1:pver),        &
           vc_dycore, ptop=state%pintdry(1:ncol,1), phis = state%phis(1:ncol),         &
           z_mid = state%z_ini(1:ncol,:),                                              &
           te = state%te_cur(1:ncol,dyn_te_idx), H2O = state%tw_cur(1:ncol,dyn_te_idx))
    else if (vc_dycore == vc_dry_pressure) then
      !
      ! SE specific hydrostatic energy
      !
      if (state%psetcols == pcols) then
        cp_or_cv(:ncol,:) = cp_or_cv_dycore(:ncol,:,lchnk)
        scaling(:ncol,:)  = cpairv(:ncol,:,lchnk)/cp_or_cv_dycore(:ncol,:,lchnk)
      else
        cp_or_cv(:ncol,:) = cpair
        scaling(:ncol,:)  = 1.0_r8
      endif
      !
      ! enthalpy scaling for energy consistency
      !
      temp(1:ncol,:)   = state%temp_ini(1:ncol,:)+scaling(1:ncol,:)*(state%T(1:ncol,:)-state%temp_ini(1:ncol,:))
      call get_hydrostatic_energy(state%q(1:ncol,1:pver,1:pcnst),.true.,               &
           state%pdel(1:ncol,1:pver), cp_or_cv(1:ncol,1:pver),                         &
           state%u(1:ncol,1:pver), state%v(1:ncol,1:pver), temp(1:ncol,1:pver),        &
           vc_dry_pressure, ptop=state%pintdry(1:ncol,1), phis = state%phis(1:ncol),   &
           te = state%te_cur(1:ncol,dyn_te_idx), H2O = state%tw_cur(1:ncol,dyn_te_idx))
    else
      state%te_cur(1:ncol,dyn_te_idx) = te(1:ncol)
      state%tw_cur(1:ncol,dyn_te_idx) = tw(1:ncol)
    end if
  end subroutine check_energy_chng

end module check_energy_chng