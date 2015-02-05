program test_jacobigq

  use kinds
  use tools
  use structures
  
  implicit none
  real(wp), dimension(4):: myx, myw
  real(wp) :: myalpha= 0.0, mybeta=0.0
  integer(ip) :: myn=3
  integer :: sizex = 4
  integer :: i
  
  
  call JacobiGQ(myalpha, mybeta, myn, myx, myw)
  print*, "x: "
  do i=1,sizex
     print*, myx(i)
  end do
  print*, "-----"
  print*, "w: "
  do i=1,sizex
     print*, myw(i)
  end do

end program test_jacobigq
