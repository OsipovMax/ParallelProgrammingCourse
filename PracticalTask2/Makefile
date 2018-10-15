PPT_prac: PPT_prac2.cpp Matrix.h
	g++ -std=c++11 -I/usr/include - O3 PPT_prac2.cpp /usr/lib64/libpapi.so -o PPT_prac

compare: compare.cpp
	g++ compare.cpp -o compare

genMatrix: genMatrix.cpp
	g++ genMatrix.cpp -o genMatrix

printMtr: printMtr.cpp
	g++ printMtr.cpp -o printMtr

test: PPT_prac compare
	./PPT_prac A B C 0 0
	./compare C C0
	./PPT_prac A B C 1 0
	./compare C C0
	./PPT_prac A B C 2 0
	./compare C C0

report: plotScript.sh gnuScript.sh
	./plotScript.sh
	./gnuScript.sh
