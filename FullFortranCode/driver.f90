program driver

  use kinds
  use structures
  use tools

  implicit none

  real(wp), parameter :: pi = 3.1415926535897932385_wp
  type(element), pointer :: first
  integer(ip) :: k, n, i, j
  real(wp), dimension(:), allocatable :: pvec, nup, emat, ematt, nupt, rhs
  real(wp) :: finaltime, cfl
  integer(ip) :: nsteps, output_every

  nvar = 1
  a = 2.0_wp*pi
  alpha = 1.0_wp
  k = 8
  allocate(order(k))
  order = 3
!  order(4) = 6
!  order(5) = 6
!  order(6) = 6
!  order(7) = 6
  call MeshGen1D ( 0.0_wp, 2.0_wp*pi, k, first )
  xmin = minval (empty)
  call loop ( first, setup_element )
!  print*,'xmin = ', xmin

  call loop ( first, set_initial )
!  call loop ( first, printall )
   
  time = 0.0_wp
  finaltime = 1000.0_wp
  cfl = 0.75_wp
  dtime = cfl/(2*pi)*xmin
  dtime = 0.25_wp*dtime

  nsteps = ceiling(finaltime/dtime)
  dtime = finaltime/nsteps
  dtime = 0.01_wp
  nsteps = 400
  output_every = nsteps/200
  output_every = 1

  open(1,file='output.dat',status='replace',action='write')
  write(1,*), '"x-label x'
  write(1,*), '"y-label u'
  write(1,*)
  write(1,*)
  open(2,file='error.dat',status='replace',action='write')
  write(2,*), '"x-label x'
  write(2,*), '"y-label u-error'
  write(2,*)
  write(2,*)

  do i = 0, nsteps-1
    if (mod(i,output_every) == 0) then
      write(1,*), '"Time = ', time
      write(2,*), '"Time = ', time
      call loop ( first, output_element )
      call loop ( first, output_error )
      write(1,*) 
      write(1,*) 
      write(2,*) 
      write(2,*) 
    end if
    call rk4 ( first, advectrhs1d )
!    print*,'u = ', first%u
!    print*,'urhs = ', first%urhs 
!    print*,'err = ', first%u(:,1)-sin(first%xk-a*time)
  end do

  write(1,*), '"Time = ', time
  write(2,*), '"Time = ', time
  call loop ( first, output_element )
  call loop ( first, output_error )

  close(1)
  close(2)

!  allocate(pvec(order(1)+1),nup(order(1)+1),emat(order(1)+1),ematt(order(1)+1),nupt(order(1)+1))
!
!  do i=1,order(1)+1
!    pvec(i:i) = JacobiP(first%r(4:4),0.0_wp,0.0_wp,i-1)
!  end do
!  print*,'pvec = ', pvec
!  
!  nup = 0.0_wp
!  do i = 1, order(1)+1
!    do j = 1, order(1)+1
!      nup(i) = nup(i) + first%v(i,j)*pvec(j)
!    end do
!  end do
!
!  print*,'nup = ', nup
!
!  emat = 0.0_wp
!  emat(4) = 1.0_wp
!  ematt = 0.0
!  do i = 1, order(1)+1
!    do j = 1, order(1)+1
!      ematt(i) = ematt(i) + first%v(j,i)*emat(j)
!    end do
!  end do
!
!  print*,'ematt = ', ematt
  contains
    subroutine set_initial ( work )
      use structures

      type(element), pointer, intent(inout) :: work

      allocate(work%u(work%n+1,nvar),work%urhs(work%n+1,nvar), &
                                     work%resu(work%n+1,nvar))

      work%u(:,1) = sin(work%xk+0.5_wp*pi)
      work%urhs = 0.0_wp
    end subroutine set_initial

    subroutine advectrhs1d ( work )
      use structures

      type(element), pointer, intent(inout) :: work
      real(wp), dimension(2) :: uint, uext, nx, uavg, ujump, du
      

      nx(1) = -1.0_wp
      nx(2) = 1.0_wp
      uint(1) = work%u(1,1)
      uint(2) = work%u(work%n+1,1)
      if ( associated(work%left) ) then
        uext(1) = work%left%u(work%left%n+1,1)
      else
        uext(1) = -sin(a*localtime)
      end if
      if ( associated(work%right) ) then
        uext(2) = work%right%u(1,1)
      else
        uext(2) = 0.0_wp
      end if
      du = 0.5_wp*(uint-uext)*(nx*a - (1-alpha)*abs(nx*a))
      if ( .not. associated(work%right) ) then
        du(2) = 0.0_wp
      end if
      work%urhs(:,1) = -a * work%rx * matmul ( work%dr, work%u(:,1) ) &
                    + matmul ( work%lift, work%rx(1)*du ) 
      if ( work%ind == 0 ) then
        print*, 'u = ', work%u(:,1)
        print*, 'uint = ', uint
        print*, 'uext = ', uext
        print*, 'du = ', du
        print*, 'dudx = ', work%rx * matmul ( work%dr, work%u(:,1) )
        print*,'boundary = ', matmul ( work%lift, work%rx(1)*du )
        print*, 'urhs = ', work%urhs(:,1)
        print*
        print*
      end if
!      if (.not. associated(work%left) ) then
!        print*,'localtime = ', localtime
!        print*,'du = ', du
!        print*,'u = ', work%u
!        print*,'tmp = ',  -a * work%rx * matmul ( work%dr, work%u(:,1) )
!        print*,'tmp2 = ', matmul ( work%lift, work%rx(1)*du )
!        print*,'rhs = ', work%urhs
!        print*
!      end if
    end subroutine advectrhs1d

    subroutine output_error ( work )
      use structures

      type(element), pointer, intent(inout) :: work
      integer(ip) :: i

      do i = 1, work%n+1
        write(2,*) work%xk(i), work%u(i,1)-sin(work%xk(i)+0.5_wp*pi-a*time)
      end do
    end subroutine output_error
end program driver
