main: main.cpp Matrix.h
	mpicxx -O3 main.cpp -o main

main_cppLint: main_cppLint.cpp Matrix.h
	mpicxx -O3 main_cppLint.cpp -o main_cppLint

genMatrix: genMatrix.cpp
	g++ -O3 genMatrix.cpp -o genMatrix

genMap: genMap.cpp
	g++ -O3 genMap.cpp -o genMap

printMtr: printMtr.cpp
	g++ -O3 printMtr.cpp -o printMtr

testCpplint: genMatrix main_cppLint
	./genMatrix d 4 4 Alint
	sleep 5
	./genMatrix d 4 1 Blint
	mpisubmit.bg -n 1 -w 00:05:00 -m SMP main_cppLint -- Alint Blint C
	sleep 30
	mpisubmit.bg -n 2 -w 00:05:00 -m SMP main_cppLint -- Alint Blint C

testRow: genMatrix main
	./genMatrix d 4 3 Ar
	sleep 5
	./genMatrix d 3 1 Br
	mpisubmit.bg -n 1 -w 00:05:00 -m SMP main -- Ar Br C
	sleep 30
	mpisubmit.bg -n 2 -w 00:05:00 -m SMP main -- Ar Br C
	
testCol: genMatrix main
	./genMatrix d 3 4 Ac
	sleep 5
	./genMatrix d 4 1 Bc
	mpisubmit.bg -n 1 -w 00:05:00 -m SMP main -- Ac Bc C
	sleep 30
	mpisubmit.bg -n 2 -w 00:05:00 -m SMP main -- Ac Bc C

	
report512: main
	./genMatrix d 512 512 A
	sleep 5
	./genMatrix d 512 1 B
	mpisubmit.bg -n 1 -w 00:05:00 -m SMP main -- A B C
	sleep 70
	mpisubmit.bg -n 32 -w 00:05:00 -m SMP main -- A B C
	sleep 70
	mpisubmit.bg -n 64 -w 00:05:00 -m SMP main -- A B C
	sleep 70
	mpisubmit.bg -n 128 -w 00:05:00 -m SMP main -- A B C
	sleep 70
	mpisubmit.bg -n 256 -w 00:05:00 -m SMP main -- A B C
	sleep 70
	mpisubmit.bg -n 512 -w 00:05:00 -m SMP main -- A B C
	sleep 70
	llsubmit j512.jcf
	sleep 100
	mv result result512

report1024: main
	./genMatrix d 1024 1024 A
	./genMatrix d 1024 1 B
	mpisubmit.bg -n 1 -w 00:05:00 -m SMP main -- A B C
	sleep 70
	mpisubmit.bg -n 32 -w 00:05:00 -m SMP main -- A B C
	sleep 70
	mpisubmit.bg -n 64 -w 00:05:00 -m SMP main -- A B C
	sleep 70
	mpisubmit.bg -n 128 -w 00:05:00 -m SMP main -- A B C
	sleep 70
	mpisubmit.bg -n 256 -w 00:05:00 -m SMP main -- A B C
	sleep 70
	mpisubmit.bg -n 512 -w 00:05:00 -m SMP main -- A B C
	sleep 70
	llsubmit j512.jcf
	sleep 100
	mv result result1024

report2048: main
	./genMatrix d 2048 2048 A
	./genMatrix d 2048 1 B
	mpisubmit.bg -n 1 -w 00:05:00 -m SMP main -- A B C
	sleep 70
	mpisubmit.bg -n 32 -w 00:05:00 -m SMP main -- A B C
	sleep 70
	mpisubmit.bg -n 64 -w 00:05:00 -m SMP main -- A B C
	sleep 70
	mpisubmit.bg -n 128 -w 00:05:00 -m SMP main -- A B C
	sleep 70
	mpisubmit.bg -n 256 -w 00:05:00 -m SMP main -- A B C
	sleep 70
	mpisubmit.bg -n 512 -w 00:05:00 -m SMP main -- A B C
	sleep 70
	llsubmit j512.jcf
	sleep 100
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
	mpisubmit.bg -n 512 -w 00:05:00 -m SMP main -- A B C
	sleep 100
	llsubmit j512.jcf
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
	mpisubmit.bg -n 512 -w 00:05:00 -m SMP main -- A B C
	sleep 100
	llsubmit j512.jcf
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
	mpisubmit.bg -n 512 -w 00:05:00 -m SMP main -- A B C
	sleep 100
	llsubmit j512.jcf
	sleep 100
	mv result result4096Col




	