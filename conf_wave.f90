program conf_wave

  use kinds
  use structures
  use tools
  use schw_physics

  implicit none

  type(element), pointer :: first
  integer(ip) :: k, n, i, j
  real(wp), dimension(:), allocatable :: pvec, nup, emat, ematt, nupt, rhs
  real(wp) :: finaltime, cfl, rbh, rscri, deltar, rmin
  integer(ip) :: nsteps, output_every

  nvar = 3
  mass = 1.0_wp
  r0 = 0.7_wp
  sigma = 0.1_wp
  rbh = 0.5_wp*(sqrt(5.0_wp)-1.0_wp)
  rscri = 1.0_wp
  ll = 0
!  k = 982
  k = 21
  deltar = (rscri-rbh)/(k-1)
  rmin = rbh-deltar
  extract_radius = rscri-deltar
!  k = 250
  allocate(order(k))
  order = 5
  print*,'rmin = ', rmin
  print*,'rscri = ', rscri
  print*,'deltar = ', deltar
  call MeshGen1D ( rmin, rscri, k, first )
!  call MeshGen1D ( 1.8_wp, 1025.8_wp, k, first )
!  call MeshGen1D ( 1.8_wp, 33.8_wp, k, first )
  xmin = minval (empty)
  call loop ( first, setup_element )
!  print*,'xmin = ', xmin

  call loop ( first, set_initial )
!  call loop ( first, check_deriv )
!  call loop ( first, printall )
   
  call loop ( first, set_coeffs )
!  stop
  time = 0.0_wp
  finaltime = 2000.0_wp
  cfl = 50.0_wp
  dtime = cfl*xmin
!  dtime = 0.5_wp*dtime

  nsteps = ceiling(finaltime/dtime)
  print*,'nsteps = ', nsteps
!  nsteps = 100000
  dtime = finaltime/nsteps
!  nsteps = 200
  output_every = max(1,nsteps/10000)
!  output_every = max(4,nsteps/1000)
!  output_every = nsteps
!  output_every = 1

!  open(1,file='psi.dat',status='replace',action='write')
!  write(1,*), '"x-label x'
!  write(1,*), '"y-label psi'
!  write(1,*)
!  write(1,*)
!  open(2,file='rho.dat',status='replace',action='write')
!  write(2,*), '"x-label x'
!  write(2,*), '"y-label rho'
!  write(2,*)
!  write(2,*)
!  open(3,file='phi.dat',status='replace',action='write')
!  write(3,*), '"x-label x'
!  write(3,*), '"y-label phi'
!  write(3,*)
!  write(3,*)
!  open(7,file='rhs.dat',status='replace',action='write')
!  write(7,*), '"x-label x'
!  write(7,*), '"y-label rhs(rho)'
!  write(7,*)
!  write(7,*)
!  open(8,file='rhs2.dat',status='replace',action='write')
!  write(8,*), '"x-label x'
!  write(8,*), '"y-label rhs(phi)'
!  write(8,*)
!  write(8,*)
  open(4,file='extract_psi.dat',status='replace',action='write')
  open(5,file='extract_rho.dat',status='replace',action='write')
  open(6,file='extract_phi.dat',status='replace',action='write')
  call loop ( first, conf_wave_rhs )
  do i = 0, nsteps-1
    if (mod(i,output_every) == 0) then
      print*,'Time = ', time
!      write(1,*), '"Time = ', time
!      write(2,*), '"Time = ', time
!      write(3,*), '"Time = ', time
!      call loop ( first, output_element )
!      write(1,*) 
!      write(1,*) 
!      write(2,*) 
!      write(2,*) 
!      write(3,*) 
!      write(3,*) 
!      write(7,*) 
!      write(7,*) 
!      write(8,*) 
!      write(8,*) 
      write(4,10) time, extract%u(1,1), extract_scri%u(extract_scri%n+1,1)
      write(5,10) time, extract%u(1,2), extract_scri%u(extract_scri%n+1,2)
      write(6,10) time, extract%u(1,3), extract_scri%u(extract_scri%n+1,3)
    end if
    10 format(3g32.24)
    call rk4 ( first, conf_wave_rhs )
!    print*,'u = ', first%u
!    print*,'urhs = ', first%urhs 
!    print*,'err = ', first%u(:,1)-sin(first%xk-a*time)
  end do

!  write(1,*), '"Time = ', time
!  write(2,*), '"Time = ', time
!  write(3,*), '"Time = ', time
  write(4,10) time, extract%u(1,1), extract_scri%u(extract_scri%n+1,1)
  write(5,10) time, extract%u(2,2), extract_scri%u(extract_scri%n+1,2)
  write(6,10) time, extract%u(3,3), extract_scri%u(extract_scri%n+1,3)
!  call loop ( first, output_element )

  close(1)
  close(2)
  close(3)
  close(4)
  close(5)
  close(6)
  close(7)

  contains
    subroutine set_initial ( work )
      use structures
      use schw_physics

      type(element), pointer, intent(inout) :: work

      allocate(work%u(work%n+1,nvar),work%urhs(work%n+1,nvar), &
                                     work%resu(work%n+1,nvar))

      work%u(:,1) = 0.0_wp
      where ( work%xk >= r0-0.5_wp*sigma .and. work%xk <= r0+0.5_wp*sigma )
        work%u(:,2) = 2.0_wp*work%xk/(1.0_wp-work%xk**2) &
                      *exp(-1.0_wp/abs(work%xk-r0+0.5_wp*sigma) &
                           -1.0_wp/abs(work%xk-r0-0.5_wp*sigma) &
                           + 4.0_wp/sigma)
      elsewhere
        work%u(:,2) = 0.0_wp
      end where  
      work%u(:,3) = 0.0_wp
      work%urhs = 0.0_wp
      if ( abs(work%xk(1)-extract_radius) <= 1.0e-12_wp ) then
        extract => work
      end if
      if ( abs(work%xk(work%n+1)-1.0_wp) <= 1.0e-12_wp ) then
        extract_scri => work
      end if
    end subroutine set_initial

    subroutine set_coeffs ( work )
      use structures
      use schw_physics

      type(element), pointer, intent(inout) :: work
      integer(ip) :: i, np
      real(wp) :: rad
      integer(ip), dimension(2) :: ind

      np = work%n+1
      ind(1) = 1
      ind(2) = np

      allocate(work%c(np,6))

      do i = 1, np
        rad = work%xk(i)
        work%c(i,1) = (1.0_wp + 7.0_wp*rad - 16.0_wp*rad**3 - 7.0_wp*rad**4 &
                   + 7.0_wp*rad**5 + 4.0_wp*rad**6) &
                / (1.0_wp + 14.0_wp*rad + 52.0_wp*rad**2 + 23.0_wp*rad**3 &
                   - 24.0_wp*rad**4 - 16.0_wp*rad**5)
        work%c(i,2) = -((-1.0_wp + rad)**2*(1.0_wp + rad)**3* &
                       (-1.0_wp + rad + rad**2)) &
                / (4.0_wp*(-1.0_wp - 5.0_wp*rad + 4.0_wp*rad**2) &
                      *(1.0_wp + 9.0_wp*rad + 11.0_wp*rad**2 + 4.0_wp*rad**3))
        work%c(i,3) = -((1.0_wp + rad)*(1.0_wp + 12.0_wp*rad - 19.0_wp*rad**2  &
                                  - 40.0_wp*rad**3 + 9.0_wp*rad**4 &
                                  + 4.0_wp*rad**5 + 9.0_wp*rad**6  &
                                  + 16.0_wp*rad**7)) &
                /(2.0_wp*rad*(1.0_wp + rad**2) &
                             *(-1.0_wp - 5.0_wp*rad + 4.0_wp*rad**2) &
                             *(1.0_wp + 9.0_wp*rad + 11.0_wp*rad**2 &
                                                   + 4.0_wp*rad**3))
        work%c(i,4) = -((-1.0_wp + rad)*(1.0_wp + rad)**2 &
                    *(1.0_wp - 2.0_wp*rad - 9.0_wp*rad**2 + 6.0_wp*rad**3 &
                          + 3.0_wp*rad**4 + 4.0_wp*rad**5 + 5.0_wp*rad**6)) &
                / (4.0_wp*rad*(1.0_wp + rad**2) &
                 *(-1.0_wp - 5.0_wp*rad + 4.0_wp*rad**2) &
                 *(1.0_wp + 9.0_wp*rad + 11.0_wp*rad**2 + 4.0_wp*rad**3))
        work%c(i,5) = -((-1.0_wp + rad)*(1.0_wp + rad)**3 &
                   *(-2.0_wp + 5.0_wp*rad - rad**2 + 2.0_wp*rad**3)) &
                / (2.0_wp*(1.0_wp + rad**2) &
                 *(-1.0_wp - 5.0_wp*rad + 4.0_wp*rad**2) &
                 *(1.0_wp + 9.0_wp*rad + 11.0_wp*rad**2 + 4.0_wp*rad**3))
        work%c(i,6) = -((1.0_wp + rad)*(1.0_wp + rad**2)**2) &
                 /(4.0_wp*rad*(-1.0_wp - 5.0_wp*rad + 4.0_wp*rad**2) &
                             *(1.0_wp + 9.0_wp*rad + 11.0_wp*rad**2 + &
                                                      4.0_wp*rad**3))
      end do
      do i = 1, 2
        rad = work%xk(ind(i))
        work%lambda(i,1) = ((-1.0_wp + rad)**2*(1.0_wp + rad)) &
                    / (-2.0_wp - 10.0_wp*rad + 8.0_wp*rad**2)
        work%lambda(i,2) = ((1.0_wp + rad)**2*(-1.0_wp + rad + rad**2)) &
                    / (2.0_wp + 18*rad + 22.0_wp*rad**2 + 8.0_wp*rad**3)
        work%s(i,1,1) = -((-1.0_wp + rad)**2*(1.0_wp + rad)) &
                 / (2.0_wp*(-1.0_wp - 5.0_wp*rad + 4.0_wp*rad**2) &
                        *sqrt(1.0_wp + ((-1.0_wp + rad)**4*(1.0_wp + rad)**2) &
                                       /(2 + 10.0_wp*rad - 8.0_wp*rad**2)**2))
        work%s(i,1,2) = -(((1.0_wp + rad)**2*(-1.0_wp + rad + rad**2)) &
                 / sqrt(5.0_wp + 74.0_wp*rad + 409.0_wp*rad**2 &
                        + 814.0_wp*rad**3 + 768.0_wp*rad**4 &
                         +362.0_wp*rad**5 + 77.0_wp*rad**6 &
                         +6.0_wp*rad**7 + rad**8))
        work%s(i,2,1) = 1.0_wp/sqrt(1.0_wp + ((-1.0_wp + rad)**4*(1.0_wp + rad)**2) &
                                        /(2 + 10.0_wp*rad - 8.0_wp*rad**2)**2)
        work%s(i,2,2) = 1.0_wp/sqrt(1.0_wp + ((1.0_wp + rad)**4 &
                                        *(-1.0_wp + rad + rad**2)**2) &
                                      / (2.0_wp + 18.0_wp*rad &
                                        + 22.0_wp*rad**2 + 8.0_wp*rad**3)**2)
        work%sinv(i,1,1) = -(((-1.0_wp - 5.0_wp*rad + 4.0_wp*rad**2) &
                        *(1.0_wp + 9.0_wp*rad + 11.0_wp*rad**2  &
                                               + 4.0_wp*rad**3) &
                        *sqrt(1.0_wp + ((-1.0_wp + rad)**4 &
                                       *(1.0_wp + rad)**2) &
                                       /( 2.0_wp + 10.0_wp*rad  &
                                        - 8.0_wp*rad**2)**2)) &
                       /(rad*(1.0_wp + rad + rad**2 + rad**3)))
        work%sinv(i,1,2) = -((1.0_wp + rad)*(-1.0_wp + rad + rad**2) &
                                      *(-1.0_wp - 5.0_wp*rad + 4.0_wp*rad**2) &
                                      *sqrt(1.0_wp + ((-1.0_wp + rad)**4 &
                                                     *(1.0_wp + rad)**2) &
                                                    /( 2.0_wp + 10.0_wp*rad &
                                                     - 8.0_wp*rad**2)**2)) &
                       /(2.0_wp*(rad + rad**3))
        work%sinv(i,2,1) = ((-1.0_wp - 5.0_wp*rad + 4.0_wp*rad**2) &
                      *(1.0_wp + 9.0_wp*rad + 11.0_wp*rad**2 + 4.0_wp*rad**3) &
                      *sqrt(1.0_wp + ((1.0_wp + rad)**4 &
                                     *(-1.0_wp + rad + rad**2)**2) &
                                     /(2.0_wp + 18.0_wp*rad + 22.0_wp*rad**2 &
                                      + 8.0_wp*rad**3)**2)) &
                      /(rad*(1.0_wp + rad + rad**2 + rad**3))
        work%sinv(i,2,2) = ((-1.0_wp + rad)**2*(1.0_wp + 9.0_wp*rad &
                                         + 11.0_wp*rad**2 + 4.0_wp*rad**3) &
                                         *sqrt(1.0_wp + ((1.0_wp + rad)**4 &
                                                       *(-1.0_wp + rad &
                                                        + rad**2)**2) &
                                                        /(2.0_wp + 18.0_wp*rad &
                                                         + 22.0_wp*rad**2 &
                                                         + 8.0_wp*rad**3)**2)) &
                      /(2.0_wp*(rad + rad**3))
      end do
    end subroutine set_coeffs

    subroutine conf_wave_rhs ( work )
      use structures
      use schw_physics

      type(element), pointer, intent(inout) :: work
      real(wp), dimension(2) :: nx
      real(wp), dimension(2,2) :: lambdaminus, lambdaplus
      real(wp), dimension(2,2) :: uint, uext, du, nflux
      
      integer(ip) :: np, i, j
      real(wp) :: tmp, rad, rad2, mass2
      integer(ip), dimension(2) :: ind

      np = work%n+1
      ind(1) = 1
      ind(2) = np
      do i =1, 2
        uint(i,1) = work%u(ind(i),2)
        uint(i,2) = work%u(ind(i),3)
      end do

      nx(1) = -one
      nx(2) = one

      if ( associated(work%left) ) then
        uext(1,1) = work%left%u(work%left%n+1,2)
        uext(1,2) = work%left%u(work%left%n+1,3)
      else
        uext(1,1) = 0.0_wp
        uext(1,2) = 0.0_wp
      end if
      if ( associated(work%right) ) then
        uext(2,1) = work%right%u(1,2)
        uext(2,2) = work%right%u(1,3)
      else
        uext(2,1) = 0.0_wp
        uext(2,2) = 0.0_wp
      end if

      do i = 1, 2
        lambdaminus = 0.0_wp
        lambdaplus = 0.0_wp
        do j=1, 2
          if ( nx(i) * work%lambda(i,j) <= 0.0_wp ) then
            lambdaminus(j,j) = nx(i)*work%lambda(i,j)
          else
            lambdaplus(j,j) = nx(i)*work%lambda(i,j)
          end if
        end do
        nflux(i,:) = matmul(lambdaplus,matmul(work%sinv(i,:,:),uint(i,:)))
!        print*,'i = ', i
!        print*,'lambdaplus  = ', lambdaplus
!        print*,'lambdaminus  = ', lambdaminus
!        print*,'sinv = ', sinv(i,:,:)
!        print*,'s = ', s(i,:,:)
!        print*,'uint = ', uint(i,:)
!        print*,'sinv*uint = ', matmul(sinv(i,:,:),uint(i,:))
!        print*,'lambdaplus*sinv*uint = ', matmul(lambdaplus,matmul(sinv(i,:,:),uint(i,:)))
        nflux(i,:) = nflux(i,:) + matmul(lambdaminus,matmul(work%sinv(i,:,:),uext(i,:)))
!        print*,'uext  = ', uext(i,:)
!        print*,'sinv*uxt = ', matmul(sinv(i,:,:),uext(i,:))
!        print*,'lambdaminus*sinv*uint = ', matmul(lambdaminus,matmul(sinv(i,:,:),uext(i,:)))
        nflux(i,:) = matmul(work%s(i,:,:),nflux(i,:))
      end do

!     The factors fac1, fac2 and 1 are minus the factors needed for the
!     Flux, since the expression for the RHS has moved to the other side.
      do i = 1, 2
        du(i,1) = -nx(i)*( work%c(ind(i),1)*uint(i,1) &
                         + work%c(ind(i),2)*uint(i,2) ) &
                         - nflux(i,1)
        du(i,2) = -nx(i)*( uint(i,1)) - nflux(i,2)
      end do
!      print*,'xk = ', work%xk(1), work%xk(np)
!      print*,'fac1 = ', fac1(ind(:))
!      print*,'fac2 = ', fac2(ind(:))
!      print*,'uint(:,1) = ', uint(:,1)
!      print*,'uint(:,2) = ', uint(:,2)
!      print*,'uext(:,1) = ', uext(:,1)
!      print*,'uext(:,2) = ', uext(:,2)
!      print*,'flux(:,1) = ', nx(:)*( -fac1(ind(:))*uint(:,1) - fac2(ind(:))*uint(:,2) )
!      print*,'flux(:,2) = ', nx(:)*(-uint(:,1))
!      print*,'nflux(:,1) = ', nflux(:,1)
!      print*,'nflux(:,2) = ', nflux(:,2)
!      print*,'du(:,1) = ', du(:,1)
!      print*,'du(:,2) = ', du(:,2)
     
      work%urhs(:,1) = work%u(:,2)
!      work%urhs(:,2) = work%rx * ( fac2 * matmul ( work%dr, work%u(:,3) ) )
      work%urhs(:,2) = work%rx * ( work%c(:,1)*matmul(work%dr, work%u(:,2)) &
                                 + work%c(:,2)*matmul(work%dr, work%u(:,3)) ) &
                       + work%c(:,3) * work%u(:,2) + work%c(:,4) * work%u(:,3) &
                       + (work%c(:,5) - (ll*(ll+1))*work%c(:,6)) *work%u(:,1) &
                       + matmul ( work%lift, work%rx(1)*du(:,1) )
      work%urhs(:,3) = work%rx * matmul ( work%dr, work%u(:,2) ) &
                       + matmul ( work%lift, work%rx(1)*du(:,2) )
                      
!      work%urhs(:,1) = -a * work%rx * matmul ( work%dr, work%u(:,1) ) &
!                    + matmul ( work%lift, work%rx(1)*du ) 
!      if (.not. associated(work%left) ) then
!        print*,'localtime = ', localtime
!        print*,'du = ', du
!        print*,'u = ', work%u
!        print*,'tmp = ',  -a * work%rx * matmul ( work%dr, work%u(:,1) )
!        print*,'tmp2 = ', matmul ( work%lift, work%rx(1)*du )
!        print*,'rhs = ', work%urhs
!        print*
!      end if
    end subroutine conf_wave_rhs

    subroutine check_deriv ( work )
      use structures
      use schw_physics

      implicit none

      type(element), pointer, intent(inout) :: work
      integer(ip) :: i

      real(wp), dimension(:), allocatable :: dudr

      allocate(dudr(work%n+1))
      dudr = work%rx * matmul ( work%dr, work%u(:,2) )+ 2.0_wp*(work%xk-r0)/sigma**2*exp(-(work%xk-r0)**2/sigma**2)
!      write(*,*) 'element(',work%ind,')'
      write(*,*) '#error(du/dr) = '
      do i = 1, work%n+1
        write(*,100) work%xk(i), dudr(i)
      end do
      100 format(2g32.24)
    end subroutine check_deriv
end program conf_wave
