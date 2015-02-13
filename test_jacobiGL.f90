program test_jacobigl

  use kinds
  use tools
  use structures
  
  implicit none
  real(wp) :: myalpha= 0.0, mybeta=0.0
  integer(ip) :: myn=5
  integer :: i
  real(wp), dimension(6):: myx

  myx= JacobiGL(myalpha, mybeta, myn)
  do i=1,myn+1
     print*,myx(i)
  enddo
end program
