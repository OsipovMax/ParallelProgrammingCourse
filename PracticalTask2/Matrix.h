#pragma once
#include <vector>
#include <iostream>

template <typename T> 
class Matrix
{
	public:
	Matrix();
	Matrix(int _rows, int _columns, std::vector < std::vector<T> > _mtr);
	~Matrix();

	int getColumns() const; 
	int getRows() const;
	//char getTypeVal() const;

	void setColumns(int);
	void setRows(int);
	//void setTypeVal(char);

	void printMatrix();

	std::vector < std::vector<T> > mtr; 
private:
	int columns;
	int rows;
	char typeVal;
};

template<typename T>
Matrix<T>::Matrix()
{

	/*columns = 3;
	rows = 3;
	mtr.reserve(3);
	mtr.resize(3);

	for (std::vector<T> &vec : mtr) {
		vec.resize(3);
		for (T &i : vec)
			i = i+1;
	}*/
	/*for (int i = 0; i < mtr.size(); ++i) {
		mtr[i].resize(3);
		for (int j = 0; j < mtr[i].size(); ++j)
			mtr[i][j] = 1;
	}*/
}
template<typename T>
Matrix<T>::Matrix(int _rows, int _columns,std::vector<std::vector<T > > _mtr) {
	rows = _rows;
	columns = _columns;
	mtr = _mtr;
}
template<typename T>
Matrix<T>::~Matrix()
{
}

template<typename T>
void Matrix<T>::printMatrix() {
	
	/*for (const std::vector<T> &vec : mtr) {
		for (const T i : vec)
			std::cout << i;
		std::cout<<std::endl;
	}*/
	for (int i = 0; i < mtr.size(); ++i) {
		for (int j = 0; j < mtr[i].size(); ++j)
			std::cout << mtr[i][j];
		std::cout << std::endl;
	}
}
template<typename T>
int Matrix<T>::getColumns() const {
	return columns;
}

template<typename T>
int Matrix<T>::getRows() const {
	return rows;
}


template<typename T>
void Matrix<T>::setColumns(int columns) {
	this->columns = columns;
}

template<typename T>
void Matrix<T>::setRows(int rows) {
	this->rows = rows;
}
//template<typename T>
//std::vector < std::vector<T> > Matrix<T>::mul(const Matrix& rightMatrix, const int mode) {
//
//}
