!pgf90 -acc -cuda -cudalib -Minfo=accel mainn4d.f90 -o mainn4d
real(8) function leastsquares (x, y, mi, mj, mx, my)
!Acuracia - quadrados minimos
	real(8), dimension(mi,mj) :: x
	real, dimension(mi,mj) :: y
	integer mi, mj
	real(8), intent(out) :: mx
	real, intent(out) :: my
	real(8), allocatable, dimension(:,:) :: w, z
	real(8) :: soma, diff
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
			if (w(i,j) > diff) then
				diff = w(i,j)
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
	IMPLICIT NONE
	real, external :: leastsquares
    integer, parameter :: ni=5120, nj=5120
	integer :: nt, i, j, ntimes=500
	real :: alpha = 0.25
    real, allocatable, dimension(:,:) :: a, d
	real(8), allocatable, dimension(:,:) :: a1, e
	real :: t1, t2, flops, my
	real(8) :: mx
    allocate(a(ni,nj),d(ni,nj))
	allocate(a1(ni,nj),e(ni,nj))
	a = 0.0d0
	d = 0.0d0
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
	! print *, "a1 - a = ", leastsquares(a1, a, ni, nj, mx, my)
	! print *, "a1 (max) = ", mx
	! print *, "a (max) = ", my

      print *,"mainnd4.f90 (base real(8)) ", ni, " x ", nj, " ntimes = ", ntimes
      call cpu_time(t1)
	do nt = 1, ntimes
	!$acc kernels
      do j = 2, nj-1
         do i = 2, ni-1
			d(i,j) = alpha * (a(i,j-1) + a(i,j+1) + a(i-1,j) + a(i+1,j))
         end do
      end do
	 a = d
	!$acc end kernels
	end do
      call cpu_time(t2)
	  print  "(4(f12.5))", ((d(i,j) , j= 1, 4 ), i= 1, 4)
	  print *
	  
	! CPU para validacao (double)
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
	print *, "e - d (real) = ", leastsquares(e, d, ni, nj, mx, my)
	 print *, "e (max) = ", mx
	 print *, "d (max) = ", my
	
      flops = 2.0*ni*nj*4
      flops = flops*ntimes
      print *,"times", t2, t1, t2-t1
      print *,"GFlops",flops/(t2-t1)/1.e9
      end program
