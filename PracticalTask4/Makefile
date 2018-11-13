main: main.cpp Matrix.h
	mpicxx main.cpp -o main

genMatrix: genMatrix.cpp
	g++ genMatrix.cpp -o genMatrix

printMtr: printMtr.cpp
	g++ printMtr.cpp -o printMtr


test: main 
	mpisubmit.bg -n 1 -w 00:05:00 -m SMP main -- A1 B1 C
	sleep 60
	mpisubmit.bg -n 2 -w 00:05:00 -m SMP main -- A1 B1 C

	
testRow: main 
	mpisubmit.bg -n 1 -w 00:05:00 -m SMP main -- A B C
	sleep 60
	mpisubmit.bg -n 2 -w 00:05:00 -m SMP main -- A B C

	
report512: main
	./genMatrix d 512 512 A
	./genMatrix d 512 1 B
	mpisubmit.bg -n 1 -w 00:05:00 -m SMP main -- A B C
	sleep 60
	mpisubmit.bg -n 32 -w 00:05:00 -m SMP main -- A B C
	sleep 60
	mpisubmit.bg -n 64 -w 00:05:00 -m SMP main -- A B C
	sleep 60
	mpisubmit.bg -n 128 -w 00:05:00 -m SMP main -- A B C
	sleep 60
	mpisubmit.bg -n 256 -w 00:05:00 -m SMP main -- A B C
	sleep 60
	mv result result512

report1024: main
	./genMatrix d 1024 1024 A
	./genMatrix d 1024 1 B
	mpisubmit.bg -n 1 -w 00:05:00 -m SMP main -- A B C
	sleep 60
	mpisubmit.bg -n 32 -w 00:05:00 -m SMP main -- A B C
	sleep 60
	mpisubmit.bg -n 64 -w 00:05:00 -m SMP main -- A B C
	sleep 60
	mpisubmit.bg -n 128 -w 00:05:00 -m SMP main -- A B C
	sleep 60
	mpisubmit.bg -n 256 -w 00:05:00 -m SMP main -- A B C
	sleep 60
	mv result result1024

report2048: main
	./genMatrix d 2048 2048 A
	./genMatrix d 2048 1 B
	mpisubmit.bg -n 1 -w 00:05:00 -m SMP main -- A B C
	sleep 60
	mpisubmit.bg -n 32 -w 00:05:00 -m SMP main -- A B C
	sleep 60
	mpisubmit.bg -n 64 -w 00:05:00 -m SMP main -- A B C
	sleep 60
	mpisubmit.bg -n 128 -w 00:05:00 -m SMP main -- A B C
	sleep 60
	mpisubmit.bg -n 256 -w 00:05:00 -m SMP main -- A B C
	sleep 60
	mv result result2048

report4096: main
	./genMatrix d 4096 4096 A
	./genMatrix d 4096 1 B
	mpisubmit.bg -n 1 -w 00:05:00 -m SMP main -- A B C
	sleep 100
	mpisubmit.bg -n 32 -w 00:05:00 -m SMP main -- A B C
	sleep 100
	mpisubmit.bg -n 64 -w 00:05:00 -m SMP main -- A B C
	sleep 100
	mpisubmit.bg -n 128 -w 00:05:00 -m SMP main -- A B C
	sleep 100
	mpisubmit.bg -n 256 -w 00:05:00 -m SMP main -- A B C
	sleep 100
	mv result result4096

report4096Row: main
	./genMatrix d 4096 1024 A
	./genMatrix d 1024 1 B
	mpisubmit.bg -n 1 -w 00:05:00 -m SMP main -- A B C
	sleep 100
	mpisubmit.bg -n 32 -w 00:05:00 -m SMP main -- A B C
	sleep 100
	mpisubmit.bg -n 64 -w 00:05:00 -m SMP main -- A B C
	sleep 100
	mpisubmit.bg -n 128 -w 00:05:00 -m SMP main -- A B C
	sleep 100
	mpisubmit.bg -n 256 -w 00:05:00 -m SMP main -- A B C
	sleep 100
	mv result result4096Row

report4096Col: main
	./genMatrix d 1024 4096 A
	./genMatrix d 4096 1 B
	mpisubmit.bg -n 1 -w 00:05:00 -m SMP main -- A B C
	sleep 100
	mpisubmit.bg -n 32 -w 00:05:00 -m SMP main -- A B C
	sleep 100
	mpisubmit.bg -n 64 -w 00:05:00 -m SMP main -- A B C
	sleep 100
	mpisubmit.bg -n 128 -w 00:05:00 -m SMP main -- A B C
	sleep 100
	mpisubmit.bg -n 256 -w 00:05:00 -m SMP main -- A B C
	sleep 100
	mv result result4096Col




	