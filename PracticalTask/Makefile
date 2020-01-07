main: main.cpp
	mpixlcxx_r -qsmp=omp main.cpp -o main

mainT: main.cpp
	mpixlcxx_r -qsmp=omp:noopt main.cpp -o main

runT:
	mpisubmit.bg -n 128 -w 00:15:00 -m SMP main --  8192 8192 32 4 20 0.00000001


run:
	mpisubmit.bg -n 512 --stdout 16384_512_3_noOpt.out -w 00:05:00 -m SMP -e "OMP_NUM_THREADS=3" main -- 16384 16384 32 16 20 0.0000000007

