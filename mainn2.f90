!pgf90 -acc -cuda -cudalib -Minfo=accel mainn2.f90 -o mainn2
!real(8) function difference (x, y, mi, mj)
real(8) function leastsquares (x, y, mi, mj)
!Acuracia - quadrados minimos
	real(8), dimension(mi,mj) :: x
	real(2), dimension(mi,mj) :: y
	integer mi, mj
	
	real(8), allocatable, dimension(:,:) :: z
	real(8) :: soma
	integer :: i, j
	allocate(z(mi,mj))
	!z = (x - y)
	z = (x - y)**2
	soma = 0.0d0
	do j = 1, mj
		do i = 1, mi
		    soma = soma + abs(z(i, j))
		end do
	end do
	!difference = soma
	leastsquares = soma / (mi * mj)
	return 
end function
	 
program main
    use cutensorex
	IMPLICIT NONE
	!real(8), external :: difference
	real(8), external :: leastsquares
    integer, parameter :: ni=5120, nj=5120, nk=5120
	integer :: nt, i, j, k, n = 1, ntimes=10
	integer, allocatable, dimension(:) :: state
    real(2), allocatable, dimension(:,:) :: a, b, d
	real(8), allocatable, dimension(:,:) :: a1, b1, e, f
	real(8) :: t1, t2, flops
    allocate(a(ni,nk),b(nk,nj),d(ni,nj))
	allocate(a1(ni,nk),b1(nk,nj),e(ni,nj),f(ni,nj))
	allocate(state(ni))
	state = 20220315
	call random_seed(put = state)
	call random_number(a1)
	a = a1
	print *, "a1 = ", a1(1,1)
	print *, "a = ", a(1,1)
	! verificar magnitude de diferença nos dados de entrada, entra a e a1, e entre b e b1
	! ignorar quadrados mínimos, e fazer diferença simples, modular
	!print *, "a1 - a = ", difference(a1, a, ni, nj)
	print *, "a1 - a = ", leastsquares(a1, a, ni, nj)
	call random_seed(put = state)
    call random_number(b1)
	b = b1
	print *, "b1 = ", b1(1,1)
	print *, "b = ", b(1,1)
	!print *, "b1 - b = ", difference(b1, b, ni, nj)
	print *, "b1 - b = ", leastsquares(b1, b, ni, nj)
    d = 0.
	e = 0.0d0
	f = 0.0d0

      print *,"mainn2.f90"
      call cpu_time(t1)
	do nt = 1, ntimes
	!$acc kernels
      do j = 1, nj
         do i = 1, ni
            do k = 1, nk
               d(i,j) = d(i,j) + a(i,k) * b(k,j)
            end do
         end do
      end do
	!$acc end kernels
	end do
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
	
	!print *, "e - d = ", difference(e, d, ni, nj)
	print *, "e - d = ", leastsquares(e, d, ni, nj)
	
      flops = 2.0*ni*nj*nk
      flops = flops*ntimes
      print *,"times",t2,t1,t2-t1
      print *,"GFlops",flops/(t2-t1)/1.e9
end program