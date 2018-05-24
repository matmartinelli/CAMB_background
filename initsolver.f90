
module initsolver
use precision
use constants
use ModelParams

implicit none

!global variable for solver
real(dl)                           :: initial_z                   !starting scale factor for ODE
real(dl)                           :: final_z                     !final scale factor for ODE
integer                            :: nsteps    = 10000           !number of ODE integration steps
real(dl), dimension(:),allocatable :: z_ode,sol1, sol2, solH      !rhoDE(z) and H(z)
real(dl), dimension(:),allocatable :: tempz, temp1, temp2, temp3  !temporary arrays
real(dl), dimension(:),allocatable :: b1, c1, d1, b2, c2, d2      !coefficient of polynomial for interpolation
real(dl), dimension(:),allocatable :: bh, ch, dh

logical                            :: debugging = .true.         !if T prints some files to check solver



contains

subroutine getH(a,out_hubble)
!this subroutine just returns rho_DE at any(ish) redshift
!CHANGE ME: need to check which solution is rho_DE
real(dl), intent(in)  :: a
real(dl)              :: z
real(dl), intent(out) :: out_hubble
integer               :: i
real(dl) :: xi_dhost

z=-1+1/a

if ((z.ge.initial_z).and.(z.le.final_z)) then
   out_hubble = ispline(z, z_ode, solH, bh, ch, dh, nsteps)
else
   out_hubble = ispline(final_z, z_ode, solH, bh, ch, dh, nsteps)
end if

end subroutine getH

subroutine deinterface(CP)
      Type(CAMBparams) CP
      integer, parameter      :: n = 1                          !number of dependent variables (-1)
      real, dimension(0:n)    :: x                              !dependent variables: rho_m, rho_v
      real                    :: h                              !step size
      real(dl)                :: rhoc_init, rhov_init
      integer                 :: i,j,k
      !debugging stuff
      real(dl)                :: debug_a, debug_c, debug_v, first_a_debug
      real(dl)                :: debug_q
      real(dl)                :: zp, zm, dz, fplus, fminus, integral
      real(dl)                :: xi_dhost



      !initializing global ODE solver parameters from CAMB
      initial_z = 0._dl
      final_z   = 10._dl

      !allocating arrays
      if (allocated(z_ode) .eqv. .false.) allocate(z_ode(nsteps), sol1(nsteps), sol2(nsteps),solH(nsteps))
      if (allocated(b1) .eqv. .false.) allocate(b1(nsteps), c1(nsteps), d1(nsteps))
      if (allocated(b2) .eqv. .false.) allocate(b2(nsteps), c2(nsteps), d2(nsteps))
      if (allocated(bH) .eqv. .false.) allocate(bH(nsteps), cH(nsteps), dH(nsteps))
      if (allocated(tempz) .eqv. .false.) allocate(tempz(nsteps), temp1(nsteps), temp2(nsteps),temp3(nsteps))


      !setting initial conditions for rho_c and rho_v at z=0
      xi_dhost = (3*CP%c3_dhost-sqrt(9*CP%c3_dhost**2.-6*CP%c2_dhost*(3*CP%beta_dhost+8*CP%c4_dhost)))/(9*CP%beta_dhost+24*CP%c4_dhost)
      if (debugging) write(*,*) 'xi_dhost=',xi_dhost
      x = (/ 3*(1-CP%omegav), xi_dhost/CP%H0 /)                                     !initial conditions: CHANGE ME!
      h = (final_z - initial_z)/(nsteps-1)                         !step size for runge-kutta

      call rk4sys(CP,n,h,x)

      write(*,*) 'solution done: computed f(z)'

      !inverting order: interpolation routine works with (x,y) table
      !WARNING: x needs to be strictly increasing
!      do i=1,nsteps !if needed
!         tempz(nsteps-(i-1)) = z_ode(i)
!         temp1(nsteps-(i-1)) = sol1(i)
!         temp2(nsteps-(i-1)) = sol2(i)
!         temp3(nsteps-(i-1)) = solH(i)
!      end do
      
!      z_ode(:)   = tempz(:)
!      sol1(:)    = temp1(:)
!      sol2(:)    = temp2(:)
!      solH(:)    = temp3(:)
      deallocate(tempz,temp1,temp2,temp3)
      
      if (debugging) write(*,*) 'computed H(z)'
      if (debugging) write(*,*) 'H(z=0)=',solH(1)
      if (debugging) write(*,*) 'H0    =',CP%H0

      !getting everything ready to interpolate
      call newspline(z_ode, sol1, b1, c1, d1, nsteps)
      call newspline(z_ode, sol2, b2, c2, d2, nsteps)
      call newspline(z_ode, solH, bh, ch, dh, nsteps)


      if (debugging) then
         open(656, file='test_solution.dat')
         do i=1,nsteps
            write(656,*) z_ode(i), sol1(i)*CP%H0/(3*solH(i)), solH(i)
         end do
         close(656)
         stop
      end if



end subroutine deinterface




!differential equations utilities

!
! Numerical Mathematics and Computing, Fifth Edition
! Ward Cheney & David Kincaid
! Brooks/Cole Publ. Co.
! (c) 2003
!
! Section 11.1
!
! File: rk4sys.f90
!
! Runge-Kutta method of order 4 for a system of ode's (rk4sys,xpsys)


subroutine xpsys(CP,n,k,h,x,f)
      Type(CAMBparams) CP
      real, dimension (0:n) ::  x, f
      integer n,k
      real :: h      !stepsize
      real(dl) :: redshift, derivative1, derivative2


      !Gets redshift for this step 
      redshift = initial_z + (k-1)*h

      !Equations are insane, thus I compute the r.h.s. in a separate routine
      !to make things less messy
      call get_derivative(CP,redshift,n,x,solH(k),derivative1, derivative2)
      solH(k) = solH(k)*CP%H0

      !These are the actual derivatives CHANGE ME
      !x'(0) = f(0)---> x(0) = phi
      !x'(1) = f(1)---> x(1) = psi = phi'

      f(0) = derivative1
      f(1) = derivative2

end subroutine xpsys

subroutine rk4sys(CP,n,h,x)
      Type(CAMBparams) CP
      real ::  x(0:n)
      real, allocatable :: y(:), f(:,:)
      integer :: i, k, n
      real :: h


      z_ode(1) = initial_z
      sol1(1) = x(0)
      sol2(1) = x(1)

      allocate (y(0:n), f(0:n,4))
out:  do k = 1,nsteps
        call xpsys(CP,n,k,h,x,f(0,1))
in1:    do i = 0,n
          y(i) = x(i) + 0.5*h*f(i,1)
        end do in1
        call xpsys(CP,n,k,h,y,f(0,2))
in2:    do i = 0,n
          y(i) = x(i) + 0.5*h*f(i,2)
        end do in2
        call xpsys(CP,n,k,h,y,f(0,3))
in3:    do i = 0,n
          y(i) = x(i) + h*f(i,3)
        end do in3
        call xpsys(CP,n,k,h,y,f(0,4))
in4:    do i = 0,n
          x(i) = x(i) + (h/6.0)* (f(i,1) + 2.0*(f(i,2) + f(i,3)) + f(i,4))
        end do in4
!        print *, k, x
        !storing functions at each step
        z_ode(k+1) = initial_z + k*h
        sol1(k+1)  = x(0)
        sol2(k+1)  = x(1)
      end do out
end subroutine rk4sys

!INTERPOLATION ROUTINES-----------------------------------------------
   subroutine newspline (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
implicit none
integer n
real(dl) x(n), y(n), b(n), c(n), d(n)
integer i, j, gap
real(dl) h

gap = n-1
! check input
if ( n < 2 ) return
if ( n < 3 ) then
  b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
  c(1) = 0.
  d(1) = 0.
  b(2) = b(1)
  c(2) = 0.
  d(2) = 0.
  return
end if
!
! step 1: preparation
!
d(1) = x(2) - x(1)
c(2) = (y(2) - y(1))/d(1)
do i = 2, gap
  d(i) = x(i+1) - x(i)
  b(i) = 2.0*(d(i-1) + d(i))
  c(i+1) = (y(i+1) - y(i))/d(i)
  c(i) = c(i+1) - c(i)
end do
!
! step 2: end conditions
!
b(1) = -d(1)
b(n) = -d(n-1)
c(1) = 0.0
c(n) = 0.0
if(n /= 3) then
  c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
  c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
  c(1) = c(1)*d(1)**2/(x(4)-x(1))
  c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
end if
!
! step 3: forward elimination
!
do i = 2, n
  h = d(i-1)/b(i-1)
  b(i) = b(i) - h*d(i-1)
  c(i) = c(i) - h*c(i-1)
end do
!
! step 4: back substitution
!
c(n) = c(n)/b(n)
do j = 1, gap
  i = n-j
  c(i) = (c(i) - d(i)*c(i+1))/b(i)
end do
!
! step 5: compute spline coefficients
!
b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
do i = 1, gap
  b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
  d(i) = (c(i+1) - c(i))/d(i)
  c(i) = 3.*c(i)
end do
c(n) = 3.0*c(n)
d(n) = d(n-1)
end subroutine newspline

  function ispline(u, x, y, b, c, d, n)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================
implicit none
double precision ispline
integer n
double precision  u, x(n), y(n), b(n), c(n), d(n)
integer i, j, k
double precision dx

! if u is ouside the x() interval take a boundary value (left or right)
if(u <= x(1)) then
  ispline = y(1)
  return
end if
if(u >= x(n)) then
  ispline = y(n)
  return
end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
i = 1
j = n+1
do while (j > i+1)
  k = (i+j)/2
  if(u < x(k)) then
    j=k
    else
    i=k
   end if
end do
!*
!  evaluate spline interpolation
!*
dx = u - x(i)
ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function ispline


subroutine get_derivative(CP,z,n,x,hubble,derivative1,derivative2)
!This subroutine obtains the r.h.s. of the phi'' equation
!Taken from Jorgos' notebook
Type(CAMBparams) CP
real, dimension (0:n), intent(in) ::  x
real(dl), intent(in) :: z
integer n,k
real(dl)              :: rhom !dimensionless rho_matter
real(dl), intent(out) :: hubble
real(dl), intent(out) :: derivative1, derivative2

rhom = x(0)!3*(1-CP%omegav)*(1+z)**3. !This is the adimensional rho_m(z)



!Equations coming from notebook
hubble = (2*(CP%c2_dhost*(2 + (3*CP%beta_dhost + 20*CP%c4_dhost)*x(1)**4)*(6*CP%c3_dhost*x(1)**3*(1 + 2*CP%c4_dhost*x(1)**4)* &
     &   (12 + 4*(3*CP%beta_dhost - 10*CP%c4_dhost)*x(1)**4 - 3*(3*CP%beta_dhost**2 + 16*CP%c4_dhost**2)*x(1)**8 + &
     &   2*CP%c4_dhost*(9*CP%beta_dhost**2 + 72*CP%beta_dhost*CP%c4_dhost + 80*CP%c4_dhost**2)*x(1)**12) +& 
     &   Sqrt(3.)*(4 + 12*CP%beta_dhost*x(1)**4 + 3*(3*CP%beta_dhost**2 + 32*CP%beta_dhost*CP%c4_dhost + 80*CP%c4_dhost**2)*x(1)**8)*&
     &   Sqrt((1 + 2*CP%c4_dhost*x(1)**4)**3*(-(CP%c2_dhost*x(1)**2*(4 + 12*(CP%beta_dhost + 4*CP%c4_dhost)*x(1)**4 + &
     &   (9*CP%beta_dhost**2 + 72*CP%beta_dhost*CP%c4_dhost + 80*CP%c4_dhost**2)*x(1)**8)) + &
     &   4*(1 + 2*CP%c4_dhost*x(1)**4)*(3*CP%c3_dhost**2*x(1)**6 + (2 + (3*CP%beta_dhost + 20*CP%c4_dhost)*x(1)**4)*rhom)))) + & 
     &   x(1)*(-72*CP%c3_dhost**3*x(1)**6*(1 + 2*CP%c4_dhost*x(1)**4)**2*(10 - 3*(3*CP%beta_dhost + 8*CP%c4_dhost)*x(1)**4 + &
     &   2*CP%c4_dhost*(3*CP%beta_dhost + 20*CP%c4_dhost)*x(1)**8) + &
     &   96*CP%c3_dhost*(1 + 2*CP%c4_dhost*x(1)**4)**2*(-2 + (3*CP%beta_dhost + 8*CP%c4_dhost)*x(1)**4 + (9*CP%beta_dhost**2 + &
     &   102*CP%beta_dhost*CP%c4_dhost + 280*CP%c4_dhost**2)*x(1)**8)*rhom - &
     &   12*Sqrt(3.)*CP%c3_dhost**2*x(1)**3*(10 - 3*(3*CP%beta_dhost + 8*CP%c4_dhost)*x(1)**4 + 2*CP%c4_dhost*(3*CP%beta_dhost + &
     &   20*CP%c4_dhost)*x(1)**8)*Sqrt((1 + 2*CP%c4_dhost*x(1)**4)**3*(-(CP%c2_dhost*x(1)**2*(4 + 12*(CP%beta_dhost + 4*CP%c4_dhost)*x(1)**4 + &
     &   (9*CP%beta_dhost**2 + 72*CP%beta_dhost*CP%c4_dhost + 80*CP%c4_dhost**2)*x(1)**8)) + &
     &   4*(1 + 2*CP%c4_dhost*x(1)**4)*(3*CP%c3_dhost**2*x(1)**6 + (2 + (3*CP%beta_dhost + 20*CP%c4_dhost)*x(1)**4)*rhom))) + &
     &   Sqrt(3.)*x(1)*(-6400*CP%c4_dhost**3*x(1)**8 - 80*CP%c4_dhost**2*x(1)**4*(-16 + 39*CP%beta_dhost*x(1)**4) + &
     &   9*CP%beta_dhost*(4 + 4*CP%beta_dhost*x(1)**4 - 3*CP%beta_dhost**2*x(1)**8) - 24*CP%c4_dhost*(-8 - 18*CP%beta_dhost*x(1)**4 +& 
     &   21*CP%beta_dhost**2*x(1)**8))*rhom* &
     &   Sqrt((1 + 2*CP%c4_dhost*x(1)**4)**3*(-(CP%c2_dhost*x(1)**2*(4 + 12*(CP%beta_dhost + 4*CP%c4_dhost)*x(1)**4 + &
     &   (9*CP%beta_dhost**2 + 72*CP%beta_dhost*CP%c4_dhost + 80*CP%c4_dhost**2)*x(1)**8)) + &
     &   4*(1 + 2*CP%c4_dhost*x(1)**4)*(3*CP%c3_dhost**2*x(1)**6 + (2 + (3*CP%beta_dhost + 20*CP%c4_dhost)*x(1)**4)*rhom))))))/&
     &   (3.*(1 + z)*(1 + 2*CP%c4_dhost*x(1)**4)*(2 + (3*CP%beta_dhost + 20*CP%c4_dhost)*x(1)**4)*(CP%c2_dhost*&
     &   (16 + 96*(CP%beta_dhost + 2*CP%c4_dhost)*x(1)**4 + 8*(27*CP%beta_dhost**2 + 204*CP%beta_dhost*CP%c4_dhost + &
     &   160*CP%c4_dhost**2)*x(1)**8 + 24*(9*CP%beta_dhost**3 + 126*CP%beta_dhost**2*CP%c4_dhost + 512*CP%beta_dhost*CP%c4_dhost**2 + 480*CP%c4_dhost**3)*x(1)**12 + &
     &   3*(3*CP%beta_dhost + 20*CP%c4_dhost)**2*(3*CP%beta_dhost**2 + 16*CP%beta_dhost*CP%c4_dhost + 16*CP%c4_dhost**2)*x(1)**16) -& 
     &   2*x(1)*(6*CP%c3_dhost**2*x(1)**3*(1 + 2*CP%c4_dhost*x(1)**4)*(12 + 9*CP%beta_dhost**2*x(1)**8 + 240*CP%c4_dhost**2*x(1)**8 + &
     &   16*CP%c4_dhost*(x(1)**4 + 6*CP%beta_dhost*x(1)**8)) + &
     &   x(1)*(1 + 2*CP%c4_dhost*x(1)**4)*(6400*CP%c4_dhost**3*x(1)**8 + 80*CP%c4_dhost**2*x(1)**4*(-16 + 39*CP%beta_dhost*x(1)**4) +&
     &   9*CP%beta_dhost*(-4 - 4*CP%beta_dhost*x(1)**4 + 3*CP%beta_dhost**2*x(1)**8) + 24*CP%c4_dhost*(-8 - 18*CP%beta_dhost*x(1)**4 + 21*CP%beta_dhost**2*x(1)**8))*rhom - &
     &   4*Sqrt(3.)*CP%c3_dhost*(-2 + (3*CP%beta_dhost + 20*CP%c4_dhost)*x(1)**4)*Sqrt((1 + 2*CP%c4_dhost*x(1)**4)**3*&
     &   (-(CP%c2_dhost*x(1)**2*(4 + 12*(CP%beta_dhost + 4*CP%c4_dhost)*x(1)**4 + (9*CP%beta_dhost**2 + 72*CP%beta_dhost*CP%c4_dhost + 80*CP%c4_dhost**2)*x(1)**8)) + &
     &   4*(1 + 2*CP%c4_dhost*x(1)**4)*(3*CP%c3_dhost**2*x(1)**6 + (2 + (3*CP%beta_dhost + 20*CP%c4_dhost)*x(1)**4)*rhom))))))


derivative1 = -3._dl*hubble*x(0)

derivative2 = (3*(1 + z)**3*x(1)*(2 + (3*CP%beta_dhost + 20*CP%c4_dhost)*x(1)**4)*(-12*CP%c3_dhost*x(1)*(1 +2*CP%c4_dhost*x(1)**4)* &
     &       (6*CP%c3_dhost**2*x(1)**6*(1 + 2*CP%c4_dhost*x(1)**4)**2 + (1 + 2*CP%c4_dhost*x(1)**4)**2*(2 + (3*CP%beta_dhost + 20*CP%c4_dhost)*x(1)**4)*rhom + &
     &         Sqrt(3.)*CP%c3_dhost*x(1)**3*Sqrt((1 + 2*CP%c4_dhost*x(1)**4)**3*(-(CP%c2_dhost*x(1)**2*(4 + 12*(CP%beta_dhost + &
     &       4*CP%c4_dhost)*x(1)**4 + (9*CP%beta_dhost**2 + 72*CP%beta_dhost*CP%c4_dhost + 80*CP%c4_dhost**2)*x(1)**8)) + &
     &       4*(1 + 2*CP%c4_dhost*x(1)**4)*(3*CP%c3_dhost**2*x(1)**6 + (2 + (3*CP%beta_dhost + 20*CP%c4_dhost)*x(1)**4)*rhom)))) + &
     &       CP%c2_dhost*(2 + (3*CP%beta_dhost + 20*CP%c4_dhost)*x(1)**4)*(6*CP%c3_dhost*x(1)**3*(1 + 2*CP%c4_dhost*x(1)**4)**2*(2 + (3*CP%beta_dhost + 4*CP%c4_dhost)*x(1)**4) + &
     &       Sqrt(3.)*(2 + 3*(CP%beta_dhost + 4*CP%c4_dhost)*x(1)**4)*Sqrt((1 + 2*CP%c4_dhost*x(1)**4)**3*(-(CP%c2_dhost*x(1)**2*(4 +&
     &       12*(CP%beta_dhost + 4*CP%c4_dhost)*x(1)**4 + (9*CP%beta_dhost**2 + 72*CP%beta_dhost*CP%c4_dhost + 80*CP%c4_dhost**2)*x(1)**8)) + &
     &       4*(1 + 2*CP%c4_dhost*x(1)**4)*(3*CP%c3_dhost**2*x(1)**6 + (2 + (3*CP%beta_dhost + 20*CP%c4_dhost)*x(1)**4)*rhom))))))/&
     &       (CP%c2_dhost*(2 + (3*CP%beta_dhost + 20*CP%c4_dhost)*x(1)**4)*(6*CP%c3_dhost*x(1)**3*(1 + 2*CP%c4_dhost*x(1)**4)*&
     &       (12 + 4*(3*CP%beta_dhost - 10*CP%c4_dhost)*x(1)**4 - 3*(3*CP%beta_dhost**2 + 16*CP%c4_dhost**2)*x(1)**8 + &
     &       2*CP%c4_dhost*(9*CP%beta_dhost**2 + 72*CP%beta_dhost*CP%c4_dhost + 80*CP%c4_dhost**2)*x(1)**12) + &
     &       Sqrt(3.)*(4 + 12*CP%beta_dhost*x(1)**4 + 3*(3*CP%beta_dhost**2 + 32*CP%beta_dhost*CP%c4_dhost + 80*CP%c4_dhost**2)*x(1)**8)*&
     &       Sqrt((1 + 2*CP%c4_dhost*x(1)**4)**3*(-(CP%c2_dhost*x(1)**2*(4 + 12*(CP%beta_dhost + 4*CP%c4_dhost)*x(1)**4 + &
     &       (9*CP%beta_dhost**2 + 72*CP%beta_dhost*CP%c4_dhost + 80*CP%c4_dhost**2)*x(1)**8)) + &
     &       4*(1 + 2*CP%c4_dhost*x(1)**4)*(3*CP%c3_dhost**2*x(1)**6 + (2 + (3*CP%beta_dhost + 20*CP%c4_dhost)*x(1)**4)*rhom)))) + &
     &       x(1)*(-72*CP%c3_dhost**3*x(1)**6*(1 + 2*CP%c4_dhost*x(1)**4)**2*(10 - 3*(3*CP%beta_dhost + 8*CP%c4_dhost)*x(1)**4 + 2*CP%c4_dhost*(3*CP%beta_dhost + 20*CP%c4_dhost)*x(1)**8) + &
     &       96*CP%c3_dhost*(1 + 2*CP%c4_dhost*x(1)**4)**2*(-2 + (3*CP%beta_dhost + 8*CP%c4_dhost)*x(1)**4 + (9*CP%beta_dhost**2 + 102*CP%beta_dhost*CP%c4_dhost + 280*CP%c4_dhost**2)*x(1)**8)*rhom - &
     &       12*Sqrt(3.)*CP%c3_dhost**2*x(1)**3*(10 - 3*(3*CP%beta_dhost + 8*CP%c4_dhost)*x(1)**4 + 2*CP%c4_dhost*(3*CP%beta_dhost + 20*CP%c4_dhost)*x(1)**8)*&
     &       Sqrt((1 + 2*CP%c4_dhost*x(1)**4)**3*(-(CP%c2_dhost*x(1)**2*(4 + 12*(CP%beta_dhost + 4*CP%c4_dhost)*x(1)**4 + (9*CP%beta_dhost**2 + 72*CP%beta_dhost*CP%c4_dhost + 80*CP%c4_dhost**2)*x(1)**8)) + &
     &       4*(1 + 2*CP%c4_dhost*x(1)**4)*(3*CP%c3_dhost**2*x(1)**6 + (2 + (3*CP%beta_dhost + 20*CP%c4_dhost)*x(1)**4)*rhom))) + &
     &       Sqrt(3.)*x(1)*(-6400*CP%c4_dhost**3*x(1)**8 - 80*CP%c4_dhost**2*x(1)**4*(-16 + 39*CP%beta_dhost*x(1)**4) + &
     &       9*CP%beta_dhost*(4 + 4*CP%beta_dhost*x(1)**4 - 3*CP%beta_dhost**2*x(1)**8) - 24*CP%c4_dhost*(-8 - 18*CP%beta_dhost*x(1)**4 + 21*CP%beta_dhost**2*x(1)**8))*rhom*&
     &       Sqrt((1 + 2*CP%c4_dhost*x(1)**4)**3*(-(CP%c2_dhost*x(1)**2*(4 + 12*(CP%beta_dhost + 4*CP%c4_dhost)*x(1)**4 + (9*CP%beta_dhost**2 + 72*CP%beta_dhost*CP%c4_dhost + 80*CP%c4_dhost**2)*x(1)**8)) + &
     &       4*(1 + 2*CP%c4_dhost*x(1)**4)*(3*CP%c3_dhost**2*x(1)**6 + (2 + (3*CP%beta_dhost + 20*CP%c4_dhost)*x(1)**4)*rhom)))))



end subroutine get_derivative



end module initsolver
