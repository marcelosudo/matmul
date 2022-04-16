# matmul
Experimentos em Multiplicação de Matrizes com FORTRAN, com variações de versões naive (ingênua), naive otimizada (com diretivas OpenACC) e com função intrínseca matmul. Além disso há variações entre os tipos de dados sendo real(2), real(4) e real(8).

Tests with Matrix Multiplication in FORTRAN, with variations of naive, optimized naive (OpenACC directives) and with intrinsic function matmul. Also variations with data types real(2), real(4) and real(8).

Forma de compilação:
pgf90 -acc -cuda -cudalib -Minfo=accel [program].f90 -o [program]
