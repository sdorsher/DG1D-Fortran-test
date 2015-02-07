module schw_physics

  use kinds
  use structures

  implicit none
  real(wp) :: mass, r0, sigma
  integer(ip) :: ll
  type(element), pointer :: extract, extract_scri
  real(wp) :: extract_radius 
  
  real(wp), parameter :: sqrt2 = sqrt(2.0_wp)
  real(wp), parameter :: invsqrt2 = 1.0_wp / sqrt2
  real(wp), parameter :: one = 1.0_wp
  real(wp), parameter :: two = 2.0_wp
end module schw_physics
