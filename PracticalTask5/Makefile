main: main.cpp Matrix.h
	mpicxx -O5 main.cpp -o main

genMatrix: genMatrix.cpp
	g++ -O3 genMatrix.cpp -o genMatrix

genMap: genMap.cpp
	g++ -O3 genMap.cpp -o genMap

printMtr: printMtr.cpp
	g++ -O3 printMtr.cpp -o printMtr

test8: genMatrix main
	mpisubmit.bg -n 1 -w 00:01:00 -m SMP main -- A B C1 1
	mpisubmit.bg -n 8 -w 00:01:00 -m SMP main -- A B C8 2

test9: genMatrix main
	mpisubmit.bg -n 1 -w 00:01:00 -m SMP main -- A9 B9 C1 1
	mpisubmit.bg -n 27 -w 00:01:00 -m SMP main -- A9 B9 C9 3

test16: genMatrix main
	mpisubmit.bg -n 1 -w 00:02:00 -m SMP main -- A16 B16 C1 1
	mpisubmit.bg -n 8 -w 00:02:00 -m SMP main -- A16 B16 C16 2


test512: genMatrix main
	mpisubmit.bg -n 512 -w 00:05:00 -m SMP main -- A512 B512 C512 8

test1024: genMatrix main
	mpisubmit.bg -n 512 -w 00:04:00 -m SMP main -- A1024 B1024 C1024 8

test2048: genMatrix main
	mpisubmit.bg -n 512 -w 00:05:00 -m SMP main -- A2048 B2048 C2048 8

test4096: genMatrix main
	mpisubmit.bg -n 64 -w 00:05:00 -m SMP main -- A4096 B4096 C4096 4

