module dadadj_apply_qv_tendency

   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: dadadj_apply_qv_tendency_run

CONTAINS

   !> \section arg_table_dadadj_apply_qv_tendency_run  Argument Table
   !! \htmlinclude dadadj_apply_qv_tendency_run.html
   subroutine dadadj_apply_qv_tendency_run(q_tend, state_q, dt, errmsg, errcode)

     ! update the constituent state.
     ! Replace this with standard constitutent update function.

     ! Dummy arguments
      real(kind_phys),    intent(in)    :: q_tend(:,:)     ! water vapor tendency
      real(kind_phys),    intent(inout) :: state_q(:,:)    ! water vapor
      real(kind_phys),    intent(in)    :: dt              ! physics time step
      character(len=512), intent(out)   :: errmsg
      integer,            intent(out)   :: errcode

      errcode = 0
      errmsg = ''

      state_q = state_q + (q_tend * dt)

   end subroutine dadadj_apply_qv_tendency_run

end module dadadj_apply_qv_tendency
