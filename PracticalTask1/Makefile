main: main.cpp
	g++ main.cpp -o main

compare: compare.cpp
	g++ compare.cpp -o _comp

generate:
	g++ genMatrix.cpp -o _getMatrix

print_matrix:
	g++ prinMtr.cpp -o _printMtr

test: main compare
	./main A B C 0
	./_comp C C0
	./main A B C 1
	./_comp C C0
	./main A B C 2
	./_comp C C0
	./main A B C 3
	./_comp C C0
	./main A B C 4
	./_comp C C0
	./main A B C 5
	./_comp C C0

report: generate main plotScript.sh
	./_getMatrix d 700 700 _rA
	./_getMatrix d 700 700 _rB
	./plotScript.sh
