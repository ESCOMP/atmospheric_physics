!< \section arg_table_held_suarez_update
!! \htmlinclude held_suarez_update.html
module held_suarez_update

  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  public held_suarez_update_run

!=======================================================================
contains
!=======================================================================

!> \section arg_table_held_suarez_update_run Argument Table
!! \htmlinclude held_suarez_update_run.html
  subroutine held_suarez_update_run(dt,cpair,du,dv,s,u,v,t,subnam,errmsg,     &
    errflg)

    real(kind_phys),   intent(in)   :: dt        ! time step
    real(kind_phys),   intent(in)   :: cpair
    real(kind_phys),   intent(in)   :: du(:,:)   ! Zonal wind tend
    real(kind_phys),   intent(in)   :: dv(:,:)   ! Meridional wind tend
    real(kind_phys),   intent(in)   :: s(:,:)    ! Heating rate
    real(kind_phys),   intent(inout):: u(:,:)    ! Zonal wind (m/s)
    real(kind_phys),   intent(inout):: v(:,:)    ! Meridional wind (m/s)
    real(kind_phys),   intent(inout):: t(:,:)    ! Temperature (K)
    character(len=64), intent(out)  :: subnam
    character(len=512),intent(out)  :: errmsg
    integer,           intent(out)  :: errflg

    errmsg = ' '
    errflg = 0
    subnam = "HELD_SUAREZ_1994"

    ! Add the tendencies to the state variables
    u(:,:) = u(:,:) + du(:,:) * dt
    v(:,:) = v(:,:) + dv(:,:) * dt
    t(:,:) = t(:,:) + s(:,:) * dt/cpair

  end subroutine held_suarez_update_run

end module held_suarez_update
