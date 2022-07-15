!pgf90 -acc -cuda -cudalib -Minfo=accel mainno2.f90 -o mainno2
real function leastsquares (x, y, mi, mj)
!Acuracia - quadrados minimos
	real, dimension(mi,mj) :: x
	real(2), dimension(mi,mj) :: y
	integer mi, mj
	
	real, allocatable, dimension(:,:) :: z
	real :: soma
	integer :: i, j
	allocate(z(mi,mj))
	z = (x - y)**2
	soma = 0.0d0
	do j = 1, mj
		do i = 1, mi
		    soma = soma + abs(z(i, j))
		end do
	end do
	leastsquares = soma / (mi * mj)
	return 
end function
	 
program main
    use cutensorex
	use cudadevice
	IMPLICIT NONE
	real, external :: leastsquares
    integer, parameter :: ni=1024, nj=1024
	integer :: nt, i, j, ntimes=100
	real :: alpha = 0.25
    real(2), allocatable, dimension(:,:) :: a, d
	real, allocatable, dimension(:,:) :: a1, e
	real :: t1, t2, flops
    allocate(a(ni,nj),d(ni,nj))
	allocate(a1(ni,nj),e(ni,nj))
	a = 0.
	d = 0.
	a1 = 0.0d0
	e = 0.0d0
	!Inicializacao - lateral esquerda com 1.0
	do j = 1, nj
		a(j,1) = 1.0
		d(j,1) = 1.0
		a1(j,1) = 1.0
		e(j,1) = 1.0
	end do
	! print *, "a1(1,1) = ", a1(1,1)
	! print *, "a(1,1) = ", a(1,1)
	!print *, "a1 - a = ", leastsquares(a1, a, ni, nj)

      print *,"mainno2.f90 ", ni, " x ", nj, " ntimes = ", ntimes
      call cpu_time(t1)
	  !$acc data copy(a(1:ni,1:nj)) copyin(d(1:ni,1:nj))
	do nt = 1, ntimes
	!$acc parallel loop collapse(2)
      do j = 2, nj-1
         do i = 2, ni-1
			d(i,j) = alpha * (a(i,j-1) + a(i,j+1) + a(i-1,j) + a(i+1,j))
         end do
      end do
	!$acc end parallel
	!$acc kernels
	 a = d
	 !$acc end kernels
	end do
	!$acc end data
      call cpu_time(t2)
	  print  "(8(f12.5))", ((d(i,j) , j= 1, 4 ), i= 1, 4)
	  print *
	  
	! CPU para validacao
	print *,"CPU"
	do nt = 1, ntimes
	do j = 2, nj-1
		do i = 2, ni-1
			e(i,j) = alpha * (a1(i,j-1) + a1(i,j+1) + a1(i-1,j) + a1(i+1,j))
		end do
	end do
	a1 = e
	end do
	print  "(4(f12.5))", ((e(i,j) , j= 1, 4 ), i= 1, 4)
	print *
	print *, "e - a = ", leastsquares(e, a, ni, nj)
	
      flops = 2.0*ni*nj*4
      flops = flops*ntimes
      print *,"times",t2,t1,t2-t1
      print *,"GFlops",flops/(t2-t1)/1.e9
	  end program