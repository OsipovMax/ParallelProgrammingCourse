main: main.cpp
	g++ -lpthread -fopenmp main.cpp

test: main
	./main 1 100 out_test.txt 3

report: main
	./main 1 10000000 out.txt 1
	./main 1 10000000 out.txt 2
	./main 1 10000000 out.txt 4
	./main 1 10000000 out.txt 6
	./main 1 10000000 out.txt 8
