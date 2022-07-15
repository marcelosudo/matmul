# matmul
Experimentos em Multiplicação de Matrizes com FORTRAN, com variações de versões naive (ingênua), naive otimizada (com diretivas OpenACC) e com função intrínseca matmul. Além disso há variações entre os tipos de dados sendo real(2), real(4) e real(8).

Tests with Matrix Multiplication in FORTRAN, with variations of naive, optimized naive (OpenACC directives) and with intrinsic function matmul. Also variations with data types real(2), real(4) and real(8).

Forma de compilação:
pgf90 -acc -cuda -cudalib -Minfo=accel [program].f90 -o [program]

Nome dos arquivos (racional):
- mainmm[x].f90 (chamada a matmul)
- mainn[x].f90 (versão naive)
- mainno[x].f90 (versão naive otimizada)
Sendo [x]:
- 2 - real(2) com matriz de referência real(4)
- 2d - real(2) com matriz de referência real(8)
- 4 - real(4) com matriz de referência real(4)
- 4d - real(4) com matriz de referência real(8)
- 8 - real(8)
- 2a - real(2) com acumulador em real(4)
- 2ad - real(2) com acumulador em real(8)
- 2i - real(2) com funções intrínsecas half2float e float2half, referência real(4)
- 2id - real(2) com funções intrínsecas half2float e float2half, referência real(8)
- 2k - real(2) com algoritmo de Kahan, referência real(4)
- 2kd - real(2) com algoritmo de Kahan, referência real(8)
