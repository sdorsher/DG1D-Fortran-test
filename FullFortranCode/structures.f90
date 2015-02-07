module structures

  use kinds

  implicit none

  type element
    integer(ip) :: n                              ! Order of this element
    integer(ip) :: ind                            ! Global address of this
                                                  ! element
    type(element), pointer :: left                ! Pointer to left neighbour
    type(element), pointer :: right               ! Pointer to right neighbour 
    real(wp), dimension(2) :: vx                  ! Coordinates of the
                                                  ! vertices of this element
    real(wp), dimension(:), allocatable :: xk     ! Physical coordinates of the
                                                  ! nodes in this element
    real(wp), dimension(:), allocatable :: r      ! Node location within this
                                                  ! element
    real(wp), dimension(:,:), allocatable :: v    ! Vandermonde matrix for this
                                                  ! element
    real(wp), dimension(:,:), allocatable :: dr   ! Derivative matrix for this
                                                  ! element
    real(wp), dimension(:), allocatable :: rx     ! Jacobian (dr/dx) for this
                                                  ! element
    real(wp), dimension(:,:), allocatable :: lift ! Lift matrix for this element
    real(wp), dimension(:,:), allocatable :: u    ! Data for this element
    real(wp), dimension(:,:), allocatable :: urhs ! RHS for this element
    real(wp), dimension(:,:), allocatable :: resu ! Temporary storage for RK4
                                                  ! for this element.
    real(wp), dimension(:,:), allocatable :: c    ! Position dependent
                                                  ! coefficients
    real(wp), dimension(2,2) :: lambda            ! Position dependent
                                                  ! eigenvalues
    real(wp), dimension(2,2,2) :: s, sinv         ! Position dependent
                                                  ! diagonalization matrices
  end type element

  integer(ip), dimension(:), allocatable :: order
  integer(ip) :: nvar
  real(wp) :: time, dtime, localtime
  real(wp) :: update_resu_factor, update_u_factor
  real(wp) :: xmin
  real(wp) :: a, alpha
end module structures
