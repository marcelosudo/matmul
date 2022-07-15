!pgf90 -acc -cuda -cudalib main.f90
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
    integer, parameter :: ni=5120, nj=5120, nk=5120
	integer :: nt, i, j, k, n = 1, ntimes=1
	integer, allocatable, dimension(:) :: state
    real(2), allocatable, dimension(:,:) :: a, b, d
	real, allocatable, dimension(:,:) :: a1, b1, e, f
	real(2) :: my
	real :: t1, t2, flops, soma, mx, myp
	real(2) :: tmp
    allocate(a(ni,nk),b(nk,nj),d(ni,nj))
	allocate(a1(ni,nk),b1(nk,nj),e(ni,nj),f(ni,nj))
	allocate(state(ni))
	state = 20220315
	call random_seed(put = state)
	call random_number(a1)
	a = a1
	print *, "a1 = ", a1(1,1)
	print *, "a = ", a(1,1)
	print *, "a1 - a = ", leastsquares(a1, a, ni, nj, mx, my)
	print *, "a1 (max) = ", mx
	myp = my
	print *, "a (max) = ", myp
	call random_seed(put = state)
    call random_number(b1)
	b = b1
	print *, "b1 = ", b1(1,1)
	print *, "b = ", b(1,1)
	print *, "b1 - b = ", leastsquares(b1, b, ni, nj, mx, my)
	print *, "b1 (max) = ", mx
	myp = my
	print *, "b (max) = ", myp
    d = 0.
	e = 0.0d0
	f = 0.0d0
 
      print *,"mainno2.f90"
      call cpu_time(t1)

	!$acc enter data copyin(a,b) create(d)  	
	!$acc parallel loop collapse(2) private(tmp)
	do nt = 1, ntimes
	do j = 1, nj
		do i = 1, ni
			!tmp = 0.
			!do k = 1, nk
			!	tmp = tmp+a(i,k)*b(k,j)
			!end do
			!d(i,j) = d(i,j) + tmp
			do k = 1, nk
				d(i,j) = d(i,j) + a(i,k)*b(k,j)!tentar fazer a conta em real, tentar verificar se tem cast OBJETIVO: tentar diminuir o erro quadrático médio para real(2)
			end do
		end do
	end do
	end do
	!$acc end parallel
	!$acc exit data copyout(d)
      call cpu_time(t2)

	  !converte alguns valores para real(8) para imprimir
	  do j = 1, 4
         do i = 1, 4
            f(i,j) = d(i,j)
         end do
      end do
	  print  "(8(f12.5))", ((f(i,j) , j= 1, 4 ), i= 1, 4)
	  print *
	  
	! CPU para validacao (double)
	print *,"CPU"
	!do nt = 1, ntimes
	do j = 1, nj
		do i = 1, ni
			do k = 1, nk
				e(i,j)=e(i,j)+a1(i,k)*b1(k,j)
			end do
			e(i,j) = e(i,j) * NTIMES
		end do
	end do
	!end do
	print  "(8(f12.5))", ((e(i,j) , j= 1, 4 ), i= 1, 4)
	print *
	print *, "e - d (real) = ", leastsquares(e, d, ni, nj, mx, my)
	print *, "e (max) = ", mx
	myp = my
	print *, "d (max) = ", myp
	
      flops = 2.0*ni*nj*nk
      flops = flops*ntimes
      print *,"times",t2,t1,t2-t1
      print *,"GFlops",flops/(t2-t1)/1.e9
      end program