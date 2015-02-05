module tools

  use kinds
  use structures

  implicit none

  contains

    subroutine MeshGen1D ( xmin, xmax, k, first )
      real(wp), intent(in) :: xmin, xmax
      integer(ip), intent(in) :: k
      type(element), pointer, intent(out) :: first

      integer(ip) :: i
      real(wp) :: nk
      type(element), pointer :: work

      nk = real(k,wp)
      
      allocate(work)
      nullify(work%left)
      do i=1,k
        work%vx(1) = (xmax-xmin)*real(i-1,wp)/nk + xmin
        work%vx(2) = (xmax-xmin)*real(i,wp)/nk + xmin
        work%ind = i
        if (i==1) then
          first => work
        end if
        if (i<k) then
          allocate(work%right)
          work%right%left => work
          work => work%right
        else
!          nullify(work%right)
          first%left => work
          work%right => first
        end if
      end do
    end subroutine MeshGen1D

    subroutine printall ( work )
      use structures

      type(element), pointer, intent(inout) :: work

      print*, 'element(',work%ind,') '
      print*, 'n (order) = ', work%n
      print*, 'vx = ', work%vx
      print*, 'r = ', work%r
      print*, 'v = ', work%v
      print*, 'dr = ', work%dr
      print*, 'xk = ', work%xk
      print*, 'rx = ', work%rx
      print*, 'lift = ', work%lift
      print*, 'u = ', work%u
      print*
    end subroutine printall
 
    subroutine setup_element ( work )
      use structures

      type(element), pointer, intent(inout) :: work

      work%n = order(work%ind)
      allocate ( work%r(work%n+1), work%v(work%n+1,work%n+1) )
      allocate ( work%dr(work%n+1,work%n+1) )
      work%r = JacobiGL ( 0.0_wp, 0.0_wp, work%n )
      work%v = VanderMonde1D ( work%n, work%r )
      work%dr = Dmatrix1D ( work%n, work%r, work%v )
      allocate ( work%xk(work%n+1), work%rx(work%n+1) )
      call Jacobian ( work%n, work%vx, work%r, work%dr, work%xk, work%rx )
      xmin = min ( xmin, abs(work%xk(1)-work%xk(2)) )
      allocate (work%lift(work%n+1,2) )
      work%lift = Lift1D ( work%n, work%v ) 
    end subroutine setup_element

    subroutine loop ( first, workfunction )
      type(element), pointer, intent(inout) :: first
      interface
        subroutine workfunction(work)
          use structures

          type(element), pointer, intent(inout) :: work
        end subroutine workfunction
      end interface
      type(element), pointer :: work

      if ( associated(first) ) then
        work => first
        call workfunction(work)
        do
!          if (.not. associated(work%right) ) exit
          if (associated(work%right, first) ) exit
          call workfunction(work%right)
          work => work%right
        end do
      else
        print*,'loop error: first element not associated'
        stop
      end if
    end subroutine loop

    subroutine JacobiGQ(alpha, beta, n, x, w)
      real(wp), intent(in) :: alpha, beta
      integer(ip), intent(in) :: n
      real(wp), dimension(n+1), intent(out) :: x, w

      real(wp), dimension(1:n+1) :: h1, jdiag1, ireal
      real(wp), dimension(1:n) :: jdiag2, d1, d2, d3
      real(wp), dimension(max(1,2*n)) :: work
      real(wp), dimension(1:n+1,1:n+1) :: vect
      real(wp) :: lngammaab, lngammaa, lngammab
      integer(ip) :: info, i

      if (n==0) then
        x(1) = ( alpha - beta ) / ( alpha + beta + 2.0_wp )
        w(1) = 2.0_wp
        return
      end if

      forall(i=0:n) ireal(i+1) = real(i,wp)
      h1 = 2.0_wp*ireal+alpha+beta
      where (h1>0.0_wp)
        jdiag1 = -(alpha**2-beta**2)/(h1*(h1+2.0_wp))
      elsewhere
        jdiag1 = 0.0_wp
      end where
      d1 = 2.0_wp/(h1(1:n)+2.0_wp)
      d2 = ireal(2:n+1)*(ireal(2:n+1)+alpha+beta)* &
                        (ireal(2:n+1)+alpha)*(ireal(2:n+1)+beta)
      d3 = 1.0_wp/((h1(1:n)+1)*(h1(1:n)+3))
      jdiag2 = d1*sqrt(d2*d3)
      if (wp == dp ) then
        call dsteqr('I', n+1, jdiag1, jdiag2, vect, n+1, work, info )
      else if (wp == qp ) then
        call qsteqr('I', n+1, jdiag1, jdiag2, vect, n+1, work, info )
      end if
      if (info <0) then
        print*,'Parameter ', i, ' in call to dsteqr has illegal value'
        stop
      end if
      if (info >0) then
        print*, i, ' off-diagonal elements have not converged to zero in call to dsteqr'
        stop
      end if
      x = jdiag1
      lngammaab = log_gamma(alpha+beta+1.0_wp)
      lngammaa  = log_gamma(alpha+1.0_wp)
      lngammab  = log_gamma(beta+1.0_wp)
      w = vect(1,:)**2*2.0_wp**(alpha+beta+1)/(alpha+beta+1)* &
                       exp(lngammaa+lngammab-lngammaab)
      return
    end subroutine JacobiGQ

    function JacobiGL(alpha, beta, n)
      real(wp), intent(in) :: alpha, beta
      integer(ip), intent(in) :: n
      real(wp), dimension(n+1) :: JacobiGL

      real(wp), dimension(n-1) :: w

      if ( n<=0 ) then
        print*,'JacobiGL called with n<=0. Aborting'
        stop
      end if
      JacobiGL = 0.0_wp
      if ( n==1 ) then
        JacobiGL(1) = -1.0_wp
        JacobiGL(2) = 1.0_wp
        return
      end if

      call JacobiGQ ( alpha+1.0_wp, beta+1.0_wp, n-2, JacobiGL(2:n), w )
      JacobiGL(1) = -1.0_wp
      JacobiGL(n+1) = 1.0_wp
      return
    end function JacobiGL

!    pure function JacobiP(x, alpha, beta, n)
    function JacobiP(x, alpha, beta, n)
      real(wp), dimension(:), intent(in) :: x
      real(wp), intent(in) :: alpha, beta
      integer(ip), intent(in) :: n
      real(wp), dimension(size(x)) :: JacobiP

      real(wp), dimension(n+1,size(x)) :: pl
      real(wp) :: lngammaab, lngammaa, lngammab, invsqgamma0, gamma0, gamma1
      real(wp) :: fac1, fac2
      real(wp) :: aold, anew, bnew, h1, ireal, irealp1
      integer(ip) :: i

      pl = 0.0_wp

      lngammaab = log_gamma(alpha+beta+1.0_wp)
      lngammaa  = log_gamma(alpha+1.0_wp)
      lngammab  = log_gamma(beta+1.0_wp)

      invsqgamma0 = 2.0_wp**(alpha+beta+1.0_wp)/(alpha+beta+1.0_wp)* &
                             exp(lngammaa+lngammab-lngammaab)
      gamma0 = 1.0_wp/sqrt(invsqgamma0)

      pl(1,:) = gamma0

      if (n==0) then
        JacobiP(:) = pl(1,:)
        return
      end if

      gamma1 = 0.5_wp*sqrt((alpha+beta+3.0_wp) &
                           /((alpha+1.0_wp)*(beta+1.0_wp)))*gamma0
      fac1 = (alpha+beta+2.0_wp)
      fac2 = (alpha-beta)
!      print*, " gamma 1 ", gamma1, " fac1 ", fac1, " fac2 ", fac2 
! checks out up to here
      pl(2,:) = gamma1 * ( fac1*x(:) + fac2 )

      if (n==1) then
        JacobiP(:) = pl(2,:)
        return
      end if

      aold = 2.0_wp / (2.0_wp+alpha+beta) * &
             sqrt ( (1.0_wp+alpha)*(1.0_wp+beta) / (3.0_wp+alpha+beta) )

      do i = 1, n-1
        ireal = real(i,wp)
        irealp1 = ireal+1.0_wp
        h1 = 2.0_wp*ireal+alpha+beta
        anew = 2.0_wp/(h1+2.0_wp)* &
               sqrt( irealp1*(irealp1+alpha+beta) * &
                     (irealp1+alpha) * (irealp1+beta) / &
                     (h1+1.0_wp)/(h1+3.0_wp) )
        bnew = - (alpha**2-beta**2) / (h1*(h1+2.0_wp) )
        pl(i+2,:) = 1.0_wp / anew * ( -aold*pl(i,:) + (x(:)-bnew)*pl(i+1,:) )
        aold = anew
      end do
      JacobiP(:) = pl(n+1,:)
      return
    end function JacobiP

!    pure function Vandermonde1D ( n, x )
    function Vandermonde1D ( n, x )
      integer(ip), intent(in) :: n
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(size(x),n+1) :: Vandermonde1D

      integer :: j

!      forall(j=1:n+1) Vandermonde1D(:,j) = JacobiP(x, 0.0_wp, 0.0_wp, j-1)
      do j=1,n+1
        Vandermonde1D(:,j) = JacobiP(x, 0.0_wp, 0.0_wp, j-1)
      end do
!      do j = 1, n+1
!        print*,'n = ', j
!        print*,Vandermonde1D(j,:)
!        print*
!      end do

      return
    end function VanderMonde1D

!    pure function GradJacobiP ( x, alpha, beta, n )
    function GradJacobiP ( x, alpha, beta, n )
      real(wp), dimension(:), intent(in) :: x
      real(wp), intent(in) :: alpha, beta
      integer(ip), intent(in) :: n
      real(wp), dimension(size(x)) :: GradJacobiP

      if (n==0) then
        GradJacobiP = 0.0_wp
        return
      end if
      GradJacobiP = sqrt(n*(n+alpha+beta+1))* &
                      JacobiP(x,alpha+1.0_wp,beta+1.0_wp,n-1)
      return
    end function GradJacobiP

!    pure function GradVandermonde1D ( n, x )
    function GradVandermonde1D ( n, x )
      integer(ip), intent(in) :: n
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(size(x),n+1) :: GradVandermonde1D

      integer :: j

!      forall(j=0:n) GradVandermonde1D(:,j+1) = GradJacobiP(x, 0.0_wp, 0.0_wp, j)
      do j=0,n
        GradVandermonde1D(:,j+1) = GradJacobiP(x, 0.0_wp, 0.0_wp, j)
      end do

      return
    end function GradVandermonde1D

    function Dmatrix1D ( n, x, v )
      integer(ip), intent(in) :: n
      real(wp), dimension(n+1), intent(in) :: x
      real(wp), dimension(n+1,n+1), intent(in) :: v
      real(wp), dimension(n+1,n+1) :: Dmatrix1D

      real(wp), dimension(n+1,n+1) :: vr, vrt, vt, drt
      integer(ip), dimension(n+1) :: ipiv
      integer(ip) :: i, j, info

      vt = transpose(v)
      vr = GradVandermonde1D ( n, x )
      vrt = transpose(vr)

      if ( wp == dp ) then
        call dgesv ( n+1, n+1, vt, n+1, ipiv, vrt, n+1, info )
      else if ( wp == qp ) then
        call qgesv ( n+1, n+1, vt, n+1, ipiv, vrt, n+1, info )
      endif

      if (info <0) then
        print*,'Parameter ', i, ' in call to dgesv has illegal value'
        stop
      end if
      if (info >0) then
        print*, 'Matrix in call to dgesv is singular'
        stop
      end if

      Dmatrix1D = transpose(vrt)
    end function Dmatrix1D

!    pure function Lift1D ( n, v )
    function Lift1D ( n, v )
      integer(ip), intent(in) :: n
      real(wp), dimension(n+1,n+1), intent(in) :: v
      real(wp), dimension(n+1,2) :: Lift1D

      real(wp), dimension(n+1,2) :: emat
      integer(ip) :: i, j

      emat = 0.0_wp
      emat(1,1) = 1.0_wp
      emat(n+1,2) = 1.0_wp

      Lift1D = matmul(v,matmul(transpose(v),emat))
!      print*,'Lift1D = '
!      print*,Lift1D(:,1)
!      print*
!      print*,Lift1D(:,2)
!      print*
!      print*
!      stop
      return
    end function Lift1D

    subroutine Jacobian ( n, vx, r, dr, x, rx )
      integer(ip), intent(in) :: n
      real(wp), dimension(2), intent(in) :: vx
      real(wp), dimension(n+1), intent(in) :: r
      real(wp), dimension(:,:), intent(in) :: dr
      real(wp), dimension(n+1), intent(out) :: x
      real(wp), dimension(n+1), intent(out) :: rx
      real(wp), dimension(n+1) :: xr
      integer(ip) :: i,j
      
      do i = 1, n+1
        x(i) = vx(1)+0.5_wp*(1.0_wp+r(i))*(vx(2)-vx(1))
      end do
      xr = 0.0_wp
      do i = 1, n+1
        do j = 1, n+1
          xr(i) = xr(i) + dr(i,j) * x(j)
        end do
      end do
      rx = 1.0_wp / xr
    end subroutine Jacobian

    subroutine rk4 (first, rhsroutine )
      use structures

      type(element), pointer, intent(inout) :: first
      interface
        subroutine rhsroutine(work)
          use structures
          type(element), pointer, intent(inout) :: work
        end subroutine rhsroutine
      end interface

      integer(ip), parameter :: nsteps = 5
      real(wp), dimension(nsteps), parameter :: &
        rk4a = (/ 0.0_wp, &
                  -567301805773.0_wp/1357537059087.0_wp, &
                  -2404267990393.0_wp/2016746695238.0_wp, &
                  -3550918686646.0_wp/2091501179385.0_wp, &
                  -1275806237668.0_wp/842570457699.0_wp /)
      real(wp), dimension(nsteps), parameter :: &
        rk4b = (/ 1432997174477.0_wp/9575080441755.0_wp, &
                  5161836677717.0_wp/13612068292357.0_wp, &
                  1720146321549.0_wp/2090206949498.0_wp, &
                  3134564353537.0_wp/4481467310338.0_wp, &
                  2277821191437.0_wp/14882151754819.0_wp /)
      real(wp), dimension(nsteps), parameter :: &
        rk4c = (/ 0.0_wp, &
                  1432997174477.0_wp/9575080441755.0_wp, &
                  2526269341429.0_wp/6820363962896.0_wp, &
                  2006345519317.0_wp/3224310063776.0_wp, &
                  2802321613138.0_wp/2924317926251.0_wp /)
      integer(ip) :: i

!      print*,'rk4a = ', rk4a
!      print*,'rk4b = ', rk4b
!      print*,'rk4c = ', rk4c
!      stop
      call loop (first, init_resu )

      do i = 1, nsteps
        localtime = time + rk4c(i)*dtime
!        print*,'localtime = ', localtime
        call loop ( first, rhsroutine )
        update_resu_factor = rk4a(i)
        call loop ( first, update_resu )
!        resu = rk4a(i) * resu + dtime * rhsroutine(u,localtime)
!        print*,'resu = '
!        do j=1,n+1
!          print*,resu(j,:)
!        end do
        update_u_factor = rk4b(i)
        call loop ( first, update_u )
!        u = u + rk4b(i) * resu
!        print*,'u = '
!        do j=1,n+1
!          print*,u(j,:)
!        end do
      end do
      time = time + dtime
    end subroutine rk4

    subroutine init_resu ( work )
      use structures

      type(element), pointer, intent(inout) :: work

      work%resu = 0.0_wp
    end subroutine init_resu

    subroutine update_resu ( work )
      use structures

      type(element), pointer, intent(inout) :: work
 
      work%resu = update_resu_factor * work%resu + dtime*work%urhs
    end subroutine update_resu

    subroutine update_u ( work )
      use structures

      type(element), pointer, intent(inout) :: work

      work%u = work%u + update_u_factor * work%resu
    end subroutine update_u
      

    subroutine output_element ( work )
      use structures

      type(element), pointer, intent(inout) :: work
      integer(ip) :: i, j

      do i = 1, work%n+1 
        do j = 1, nvar
          write(j,10) work%xk(i), work%u(i,j)
        end do
        write(7,10) work%xk(i), work%urhs(i,2)
        write(8,10) work%xk(i), work%urhs(i,3)
!        10 format(2f20.15)
        10 format(2g32.24)
      end do
    end subroutine output_element
end module tools
