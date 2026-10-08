! Ungather the Zhang-McFarlane deep convective cloud top from the gathered
! convective-column arrays (jt, ideep, lengath) onto the full column list.
! Columns without deep convection get the index of the lowest layer (pver),
! matching the CLDTOP field CAM's convect_deep_tend fills.
module zm_conv_cloud_top
  implicit none
  private

  public :: zm_conv_cloud_top_run

contains

!> \section arg_table_zm_conv_cloud_top_run Argument Table
!! \htmlinclude zm_conv_cloud_top_run.html
  subroutine zm_conv_cloud_top_run(ncol, pver, lengath, jt, ideep, jctop, errmsg, errflg)

    integer,            intent(in)  :: ncol
    integer,            intent(in)  :: pver
    integer,            intent(in)  :: lengath     ! number of convective columns [index]
    integer,            intent(in)  :: jt(:)       ! top layer of deep convection, gathered [index]
    integer,            intent(in)  :: ideep(:)    ! column index of each convective column [index]
    integer,            intent(out) :: jctop(:)    ! top layer of deep convection, all columns [index]
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer :: i

    errmsg = ''
    errflg = 0

    jctop(:ncol) = pver
    do i = 1, lengath
       jctop(ideep(i)) = jt(i)
    end do

  end subroutine zm_conv_cloud_top_run

end module zm_conv_cloud_top
