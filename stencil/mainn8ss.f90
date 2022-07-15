!pgf90 -acc -cuda -cudalib -Minfo=accel mainn8ss.f90 -o mainn8ss
program main
    use cutensorex
	IMPLICIT NONE
    integer, parameter :: ni=5120, nj=5120
	integer :: nt, it, i, j, ntimes=1000000
	real(8) :: alpha = 0.25
	real(8) :: tmp
    real(8), allocatable, dimension(:,:) :: a, c, d
	real :: t1, t2, flops
    allocate(a(ni,nj),c(ni,nj),d(ni,nj))
	a = 0.0d0
	d = 0.0d0
	!Inicializacao - lateral esquerda com 1.0
	do j = 1, nj
		a(j,1) = 1.0
		d(j,1) = 1.0
	end do
 
      print *,"mainn8ss.f90 ", ni, " x ", nj, " ntimes = ", ntimes
      call cpu_time(t1)
	do nt = 1, ntimes
	!$acc kernels
      do j = 2, nj-1
         do i = 2, ni-1
			d(i,j) = alpha * (a(i,j-1) + a(i,j+1) + a(i-1,j) + a(i+1,j))
         end do
      end do
	 c = d - a !maior elemento de c for 10e-7 para 
	 do j = 2, nj-1
       do i = 2, ni-1
		if (tmp < c(i,j)) tmp = c(i,j)
       end do
    end do
	 a = d
	 it = nt
	!$acc end kernels
	if (tmp < 0.000001) EXIT
	if (nt < ntimes) tmp = 0.0d0
	end do
      call cpu_time(t2)
	  print  "(4(f12.5))", ((d(i,j) , j= 1, 4 ), i= 1, 4)
	  print *, "steady point ", it, " / max value ", tmp
	  print *
	
      flops = 2.0*ni*nj*4
      flops = flops*ntimes
      print *,"times",t2,t1,t2-t1
      print *,"GFlops",flops/(t2-t1)/1.e9
      end program
