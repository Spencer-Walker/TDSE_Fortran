! NOTE if itterating over i,j and performing operations on M(i,j)
! having j as the outside loop will give a significant performance boost
! this is because Fortran is column major which means that for some matrix
! [v1,v2,v3] with vj as column vectors that elements within one vector vj
! are closest to one another in memory. C is the opposite!
module prec_def
  implicit none
  integer, parameter :: dp = kind(1.d0)
end module prec_def

program main
  use prec_def
  implicit none
! Declarations
!--------------------------------------------------------------------
  real(dp), parameter :: pi = 4*atan(1.d0)
  real(dp), parameter :: rmin = -100.d0
  real(dp), parameter :: rmax = 100.d0
  integer, parameter :: nr0 = -2000
  integer, parameter :: nr = 2000
  real(dp) :: r(nr0:nr),hr,ihr2,V0(nr0:nr)
  integer, parameter :: nl = 0
  real(dp), parameter :: dt = 0.1d0
  real(dp) :: T0
  real(dp), parameter :: cycles = 20.d0
  real(dp), parameter :: w = 0.05
  real(dp), parameter :: E0 = 0.05
  complex(dp) :: psi(nr0:nr,0:nl),rhs(nr0:nr,0:nl),diag(nr0:nr,0:nl),off_diag
  real(dp) :: ground_state(nr0:nr)
  real(dp) :: E,t
  integer :: ir,l,it,nt
  real(dp), allocatable :: d(:)
! Precomputation
!--------------------------------------------------------------------
  T0 = 2*cycles*pi/w  
  nt = ceiling(T0/dt)
  print*, T0
  allocate(d(nt))
  hr = (rmax-rmin)/real(nr-nr0+1,dp)
  do ir =nr0,nr
    r(ir) = ir*hr
    V0(ir) = -1.d0/sqrt(2.d0+r(ir)**2)
  end do

  call readdble1d(ground_state,nr-nr0+1,'ground_state.vec')

  psi = 0.d0
  psi(:,0) = ground_state

! Time evolution of the wavefunction
!--------------------------------------------------------------------
  ! Split operator
  do it = 0,nt-1
    t = it*dt
    write(*,*) t,E(t,w,E0)
  ! CN is an implicit scheme therefore we need to solve LHS = A*psi = RHS
  ! Here RHS = (1 - 0.5idtH(t))*psi_old
    rhs = psi
    ihr2 = 1.d0/hr**2
    off_diag = (0.d0,0.25d0)*dt*ihr2
    do l = 0,nl
      do ir = nr0,nr
        diag(ir,l) = (1.d0,0.d0) + (0.d0,-0.5d0)*dt*(ihr2 + V0(ir)+r(ir)*E(t,w,E0))
      end do
    end do
    call zmultri_special(off_diag,diag(:,0),rhs(:,0),nr0,nr)
    ! We then solve for psi using our tridiagonal solver
    ! LHS = (1 + 0.5idtH(t))*psi
    off_diag = (0.d0,-0.25d0)*dt*ihr2
    do l = 0,nl
      do ir = nr0,nr
        diag(ir,l) = (1.d0,0.d0) + (0.d0,0.5d0)*dt*(ihr2 + V0(ir)+r(ir)*E(t,w,E0))
      end do
    end do
    call zsoltri_special(off_diag,diag,rhs,psi,nr0,nr)

    d(it+1) = dot_product(psi(:,0),r*psi(:,0)) 
  end do

  call printdble1d(real(psi),nr-nr0+1,'psi_real.txt')
  call printdble1d(aimag(psi),nr-nr0+1,'psi_imag.txt')
  call printdble1d(V0,nr-nr0+1,'V0.txt')
  call printdble1d(r,nr-nr0+1,'r.txt')
  call printdble1d(d,nt,'d.txt')
end program main


! Evolves one l block
!--------------------------------------------------------------------


! Evolves one r block
!--------------------------------------------------------------------

! Multiplies a vector by a tridiagonal matrix
subroutine zmultri(a,b,c,x,n0,n1)
  use prec_def
  implicit none
  integer, intent(in) :: n0,n1
  complex(dp), intent(in) :: a(n0:n1-1),b(n0:n1),c(n0:n1-1)
  complex(dp), intent(inout) :: x(n0:n1)
  complex(dp) :: t0,t1
  integer :: j

  t0 = x(n0)
  x(n0) = b(n0)*t0 + c(n0)*x(n0+1)
  do j = n0+1,n1-1
    t1 = x(j)
    x(j) = a(j-1)*t0+b(j)*t1+c(j)*x(j+1)
    t0 = t1
  end do
  t1 = x(n1)
  x(n1) = a(n1-1)*t0 + b(n1)*t1

end subroutine zmultri

! Solves tridiagonal system
subroutine zsoltri(a,b,c,d,x,n0,n1)
  use prec_def
  implicit none
  integer, intent(in) :: n0,n1
  complex(dp) :: a(n0:n1-1),b(n0:n1),c(n0:n1-1),d(n0:n1),x(n0:n1)
  integer :: j
  complex(dp) :: w

  do j = n0+1,n1
    w = a(j-1)/b(j-1)
    b(j) = b(j) - w*c(j-1)
    d(j) = d(j) - w*d(j-1)
  end do
  x(n1) = d(n1)/b(n1)
  do j = n1-1,n0,-1
    x(j) = (d(j)-c(j)*x(j+1))/b(j)
  end do

end subroutine zsoltri

subroutine zsoltri_special(a,b,d,x,n0,n1)
  use prec_def
  implicit none
  integer, intent(in) :: n0,n1
  complex(dp) :: a,b(n0:n1),d(n0:n1),x(n0:n1)
  integer :: j
  complex(dp) :: w

  do j = n0+1,n1
    w = a/b(j-1)
    b(j) = b(j) - w*a
    d(j) = d(j) - w*d(j-1)
  end do
  x(n1) = d(n1)/b(n1)
  do j = n1-1,n0,-1
    x(j) = (d(j)-a*x(j+1))/b(j)
  end do
end subroutine zsoltri_special

subroutine zmultri_special(a,b,x,n0,n1)
  use prec_def
  implicit none
  integer, intent(in) :: n0,n1
  complex(dp), intent(in) :: a,b(n0:n1)
  complex(dp), intent(inout) :: x(n0:n1)
  complex(dp) :: t0,t1
  integer :: j

  t0 = x(n0)
  x(n0) = b(n0)*t0 + a*x(n0+1)
  do j = n0+1,n1-1
    t1 = x(j)
    x(j) = a*t0+b(j)*t1+a*x(j+1)
    t0 = t1
  end do
  t1 = x(n1)
  x(n1) = a*t0 + b(n1)*t1
end subroutine zmultri_special

subroutine printdble1d(u,nxl,str)
  use prec_def
  implicit none
  integer, intent(in) :: nxl
  real(dp), intent(in) :: u(nxl)
  character(len=*), intent(in) :: str
  integer :: i
  open(2,file=trim(str),status='unknown')
  do i=1,nxl,1
    if (abs(u(i)) .gt. 1d-99) then
      write(2,fmt='(E24.16)',advance='no') u(i)
    else
      write(2,fmt='(E24.16)',advance='no') 0.d0
    end if
  end do
  close(2)
end subroutine printdble1d

subroutine printdble2d(u,nx1,nx2,ny1,ny2,str)
  use prec_def
  implicit none
  integer, intent(in) :: nx1,nx2,ny1,ny2
  real(DP), intent(in) :: u(nx1:nx2,ny1:ny2)
  character(len=*), intent(in) :: str
  integer :: i,j
  open(2,file=trim(str),status='unknown')
  do j=ny1,ny2,1
    do i=nx1,nx2,1
      if(abs(u(i,j)) .lt. 1e-40) then
        write(2,fmt='(E24.16)',advance='no') 0.d0
      else
        write(2,fmt='(E24.16)',advance='no') u(i,j)
      end if
    end do
    write(2,'()')
  end do
  close(2)
end subroutine printdble2d

subroutine readdble2d(u,nx1,nx2,ny1,ny2,str)
  use prec_def
  implicit none
  integer, intent(in) :: nx1,nx2,ny1,ny2
  real(DP), intent(inout) :: u(nx1:nx2,ny1:ny2)
  character(len=*), intent(in) :: str
  integer :: j
  open(2,file=trim(str),status='old',action='read')
  do j=ny1,ny2,1
    read(2,*) u(:,j)
  end do
  close(2)
end subroutine readdble2d

subroutine readdble1d(u,nx,str)
  use prec_def
  implicit none
  integer, intent(in) :: nx
  real(dp), intent(inout) :: u(nx)
  character(len=*), intent(in) :: str
  open(2,file=trim(str),status='old',action='read')
  read(2,*) u
  close(2)
end subroutine readdble1d

function E(t,w,E0) 
  use prec_def
  implicit none
  real(dp) :: E,E0,w,t
  
  E = E0*sin(w*t)
end function E