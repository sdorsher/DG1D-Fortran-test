module kinds

  implicit none

  integer, parameter :: sp = selected_real_kind(5,30)
  integer, parameter :: dp = selected_real_kind(9,99)
  integer, parameter :: qp = selected_real_kind(20,199)

  integer, parameter :: wp = dp
  integer, parameter :: ip = selected_int_kind(8)

! These empty arrays are used to initialize variables to either the min or
! max possible number of kind wp or integer.
  real(wp), dimension(2:1) :: empty
  integer(ip), dimension(2:1) :: iempty

end module kinds
