main: main.cpp
	mpicxx main.cpp -o main

test: main 
	mpisubmit.bg -n 17 -w 00:05:00 -m dual main -- 2910 2980 test_out.txt
	sleep 10
	mpisubmit.bg -n 5 -w 00:05:00 -m dual main -- 100 220 test_out.txt
	sleep 10
	mpisubmit.bg -n 2 -w 00:05:00 -m dual main -- 1 99 test_out.txt
	sleep 10
	mpisubmit.bg -n 3 -w 00:05:00 -m dual main -- 10010 10325 test_out.txt
	sleep 10
	mpisubmit.bg -n 8 -w 00:05:00 -m dual main -- 115 217 test_out.txt

report: main max_time.sh total_time.sh
	mpisubmit.bg -n 1 -w 00:05:00 -m dual main -- 1 10000000 out.txt
	sleep 60
	mpisubmit.bg -n 2 -w 00:05:00 -m dual main -- 1 10000000 out.txt
	sleep 60
	mpisubmit.bg -n 4 -w 00:05:00 -m dual main -- 1 10000000 out.txt
	sleep 60
	mpisubmit.bg -n 8 -w 00:05:00 -m dual main -- 1 10000000 out.txt
	sleep 60
	mpisubmit.bg -n 16 -w 00:05:00 -m dual main -- 1 10000000 out.txt
	sleep 100
	./max_time.sh
	./total_time.sh
