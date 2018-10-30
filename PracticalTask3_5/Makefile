main: main.cpp
	g++ -lpthread -fopenmp main.cpp -o main

test: main
	./main 2910 2980 out_test.txt 17
	./main 100  220 out_test.txt 5
	./main 1 99 out_test.txt 2
	./main 10010 10325 out_test.txt 3
	./main 115 217 out_tst.txt 8

report: main
	./main 1 10000000 out.txt 1
	./main 1 10000000 out.txt 2
	./main 1 10000000 out.txt 4
	./main 1 10000000 out.txt 6
	./main 1 10000000 out.txt 8
	./main 1 10000000 out.txt 12
	./main 1 10000000 out.txt 16
	./script.sh
