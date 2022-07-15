!pgf90 -acc -cuda -cudalib -Minfo=accel mainn2k.f90 -o mainn2k
real function leastsquares (x, y, mi, mj, mx, my)
!Acuracia - quadrados minimos
	real, dimension(mi,mj) :: x
	real(2), dimension(mi,mj) :: y
	integer mi, mj
	real, intent(out) :: mx
	real(2), intent(out) :: my
	real, allocatable, dimension(:,:) :: w, z
	real :: soma, diff
	integer :: i, j
	allocate(w(mi,mj))
	allocate(z(mi,mj))
	w = (x - y)
	z = w**2
	soma = 0.0d0
	diff = 0.0d0
	do j = 1, mj
		do i = 1, mi
		    soma = soma + abs(z(i, j))
			if (abs(w(i,j)) > diff) then
				diff = abs(w(i,j))
				mx = x(i,j)
				my = y(i,j)
			end if
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
    integer, parameter :: ni=5120, nj=5120
	integer :: nt, i, j, ntimes=500
	real :: alpha = 0.25
    real(2), allocatable, dimension(:,:) :: a, d
	real, allocatable, dimension(:,:) :: a1, e
	real(2) :: soma, c, y, t !Kahan
	real(2) :: my
	real :: t1, t2, flops, mx
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
	!print *, "a1(1,1) = ", a1(1,1)
	!print *, "a(1,1) = ", a(1,1)
	!print *, "a1 - a = ", leastsquares(a1, a, ni, nj, mx, my)
	! print *, "a1 (max) = ", mx
	! print *, "a (max) = ", my

      print *,"mainn2k.f90 base 4", ni, " x ", nj, " ntimes = ", ntimes
      call cpu_time(t1)
	do nt = 1, ntimes
	!$acc kernels
      do j = 2, nj-1
         do i = 2, ni-1
			!d(i,j) = alpha * (a(i,j-1) + a(i,j+1) + a(i-1,j) + a(i+1,j))
			
			!Kahan Summation Algorithm
			soma = 0.
			c = 0.
            !do k = 1, nk
			   y = (a(i,j-1) + a(i,j+1)) - c
			   t = soma + y
			   c = (t - soma) - y
			   soma = t
			   
			   y = (c + a(i-1,j)) - c
			   t = soma + y
			   c = (t - soma) - y
			   soma = t
			   
			   y = (c + a(i+1,j)) - c
			   t = soma + y
			   c = (t - soma) - y
			   soma = t
            !end do
			d(i,j) = alpha * soma
			
         end do
      end do
	 a = d
	!$acc end kernels
	end do
      call cpu_time(t2)
	  print  "(4(f12.5))", ((d(i,j) , j= 1, 4 ), i= 1, 4)
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
	print *, "e - d (real(2)) = ", leastsquares(e, d, ni, nj, mx, my)
	print *, "e (max) = ", mx
	print *, "d (max) = ", my
	
      flops = 2.0*ni*nj*4
      flops = flops*ntimes
      print *,"times",t2,t1,t2-t1
      print *,"GFlops",flops/(t2-t1)/1.e9
end program