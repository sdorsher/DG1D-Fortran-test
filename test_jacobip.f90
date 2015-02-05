program test_jacobip
! n from zero to sizex, alpha and beta =0 or alpha and beta =1
!tested for sizex =4, unrealistically evenly spaced points

  use kinds
  use tools

  implicit none
  real(wp), dimension(4) :: myx,myy
  real(wp) :: myalpha, mybeta
  integer(ip) :: myn
  integer :: i
  integer :: sizex =4
  real :: step, neglim=-1.0, poslim=1.0
! sizex=4 is the dimension of x


  myalpha=1
  mybeta=1
  myn=3

  step = (poslim-neglim)/(sizex-1.0);

  do i=1,4
     myx(i)=neglim+step*(i-1)
  enddo
  myy= JacobiP(myx,myalpha,mybeta,myn)
  do i=1,sizex
     print*, myy(i)
  end do
end program test_jacobip

