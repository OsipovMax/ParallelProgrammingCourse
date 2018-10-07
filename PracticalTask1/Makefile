 main: main.cpp Matrix.h
	g++ main.cpp -o main

_comp: compare.cpp
	g++ compare.cpp -o _comp

_getMatrix: genMatrix.cpp
	g++ genMatrix.cpp -o _getMatrix

_printMtr: printMtr.cpp
	g++ printMtr.cpp -o _printMtr

test: main _comp
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

report: _getMatrix main plotScript.sh
	./_getMatrix d 700 700 _rA
	./_getMatrix d 700 700 _rB
	./plotScript.sh
