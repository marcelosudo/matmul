!pgf90 -acc -cuda -cudalib -Minfo=accel mainmm8.f90 -o mainmm8
real(8) function leastsquares (x, y, mi, mj, mx, my)
!Acuracia - quadrados minimos
	real(8), dimension(mi,mj) :: x
	real(8), dimension(mi,mj) :: y
	integer mi, mj
	real(8), intent(out) :: mx
	real(8), intent(out) :: my
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
	IMPLICIT NONE
	real(8), external :: leastsquares
    integer, parameter :: ni=1024, nj=1024, nk=1024
	integer :: nt, i, j, k, n = 1, ntimes=1
	integer, allocatable, dimension(:) :: state
    real(8), allocatable, dimension(:,:) :: a, b, a1, b1, d, e
	real :: t1, t2, flops
	real(8) :: mx, my
    allocate(a(ni,nk),b(nk,nj),a1(ni,nk),b1(nk,nj),d(ni,nj),e(ni,nj))
	allocate(state(ni))
	state = 20220315
	call random_seed(put = state)
	call random_number(a1)
	a = a1
	print *, "a1 = ", a1(1,1)
	print *, "a = ", a(1,1)
	print *, "a1 - a = ", leastsquares(a1, a, ni, nj, mx, my)
	print *, "a1 (max) = ", mx
	print *, "a (max) = ", my
	call random_seed(put = state)
    call random_number(b1)
	b = b1
	print *, "b1 = ", b1(1,1)
	print *, "b = ", b(1,1)
	print *, "b1 - b = ", leastsquares(b1, b, ni, nj, mx, my)
	print *, "b1 (max) = ", mx
	print *, "b (max) = ", my
    d = 0.
	e = 0.0d0
 
      print *,"mainmm8.f90 ", ni, " x ", nj, " ntimes = ", ntimes
      call cpu_time(t1)
	  !$acc data copyin(a,b) copyout(d)
	  !$acc host_data use_device(a, b, d)
      do nt = 1, ntimes
        d = d + matmul(a,b)
      end do
	  !$acc end host_data
	  !$acc end data
      call cpu_time(t2)
	  print  "(8(f12.5))", ((d(i,j) , j= 1, 4 ), i= 1, 4)
	  print *
	  
	! CPU para validacao (double)
	 print *,"CPU"
	 do nt = 1, ntimes
	 do j = 1, nj
		 do i = 1, ni
			 do k = 1, nk
				 e(i,j)=e(i,j)+a1(i,k)*b1(k,j)
			 end do
		 end do
	 end do
	 end do
	 print  "(8(f12.5))", ((e(i,j) , j= 1, 4 ), i= 1, 4)
	 print *
	 print *, "e - d (real(8)) = ", leastsquares(e, d, ni, nj, mx, my)
	 print *, "e (max) = ", mx
	 print *, "d (max) = ", my
	
      flops = 2.0*ni*nj*nk
      flops = flops*ntimes
      print *,"times",t2,t1,t2-t1
      print *,"GFlops",flops/(t2-t1)/1.e9
      end program
