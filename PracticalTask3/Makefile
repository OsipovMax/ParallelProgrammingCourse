main: main.cpp
	mpicxx main.cpp -o main

test: main 
	mpisubmit.bg -n 3 -w 00:05:00 -m dual main -- 1 100 test_out.txt

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
