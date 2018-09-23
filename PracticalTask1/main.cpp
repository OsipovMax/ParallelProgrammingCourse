#include "Matrix.h"
#include <fstream>
#include <string>
#include <time.h>
#include <typeinfo>
#include <chrono>
template <typename T> 
Matrix<T> mulMatrix(const Matrix<T> &m1,  const Matrix<T> &m2, const int mode, unsigned int &res_time ) {

	Matrix <T> resMtr;
	resMtr.mtr.resize(m1.getRows());
	for (int i = 0; i < resMtr.mtr.size(); ++i) {
		resMtr.mtr[i].resize(m2.getColumns());
	}
	resMtr.setRows(m1.getRows());
	resMtr.setColumns(m2.getColumns());
	if (mode == 0) {
		time_t start = time(NULL);
		for (int i = 0; i < m1.getRows(); ++i)
			for (int j = 0; j < m2.getColumns(); ++j)
				for (int k = 0; k < m1.getColumns(); ++k)
					resMtr.mtr[i][j] += m1.mtr[i][k] * m2.mtr[k][j];
	       	time_t end = time(NULL);
		res_time = difftime(end,start);
		//res_time = (double)(end - start) / CLOCKS_PER_SEC;
	}
	else if (mode == 1) {
		time_t start = time(NULL);
		for (int i = 0; i < m1.getRows(); ++i)
			for (int k = 0; k < m1.getColumns(); ++k)
				for (int j = 0; j < m2.getColumns(); ++j)
					resMtr.mtr[i][j] += m1.mtr[i][k] * m2.mtr[k][j];
		time_t end = time(NULL);
		res_time = difftime(end,start);
		//res_time = (double)(end - start) / CLOCKS_PER_SEC;
	}
	else if (mode == 2) {
		time_t start = time(NULL);
		for (int k = 0; k < m1.getColumns(); ++k)
			for (int i = 0; i < m1.getRows(); ++i)
				for (int j = 0; j < m2.getColumns(); ++j)
					resMtr.mtr[i][j] += m1.mtr[i][k] * m2.mtr[k][j];
		time_t end = time(NULL);
		res_time = difftime(end,start);
		//res_time = (double)(end - start) / CLOCKS_PER_SEC;
	}
	else if (mode == 3) {
		int p;
		time_t start = time(NULL);
		for (int j = 0; j < m2.getColumns(); ++j)
			for (int i = 0; i < m1.getRows(); ++i)
				for (int k = 0; k < m1.getColumns(); ++k)
					resMtr.mtr[i][j] += m1.mtr[i][k] * m2.mtr[k][j];
		time_t end = time(NULL);
		res_time = difftime(end,start);
	}
	else if (mode == 4) {
		time_t start = time(NULL);
		for (int j = 0; j < m2.getColumns(); ++j)
			for (int k = 0; k < m1.getColumns(); ++k)
				for (int i = 0; i < m1.getRows(); ++i)
					resMtr.mtr[i][j] += m1.mtr[i][k] * m2.mtr[k][j];
		time_t end = time(NULL);
		res_time = difftime(end,start);
		//res_time = (double)(end - start) / CLOCKS_PER_SEC;
	}
	else if (mode == 5) {
		time_t start = time(NULL);
		//clock_t start = clock();
		for (int k = 0; k < m1.getColumns(); ++k)
			for (int j = 0; j < m2.getColumns(); ++j)
				for (int i = 0; i < m1.getRows(); ++i)
					resMtr.mtr[i][j] += m1.mtr[i][k] * m2.mtr[k][j];
		time_t end = time(NULL);
		res_time = difftime(end,start);
		//clock_t end = clock();
		//res_time = (double)(end - start) / CLOCKS_PER_SEC;
	}

	return resMtr;
}


template<typename T>
Matrix<T> binReader(const std::string &strr) {

	int modeVal, n, m;
	T val;
	char type_val;
	Matrix<T> tempObj;


	std::fstream f1(strr, std::ios::binary | std::ios::in);
	f1.read((char*)&type_val, sizeof(char));
	f1.read((char*)&n, sizeof(int));
	f1.read((char*)&m, sizeof(int));
	tempObj.setRows(n);
	tempObj.setColumns(m);

	tempObj.mtr.resize(n);
	for (int i = 0; i < n; ++i) {
		tempObj.mtr[i].resize(m);
		for (int j = 0; j < m; j++) {
			f1.read((char*)&val, sizeof(T));
			tempObj.mtr[i][j] = val;
		}
	}
	f1.close();
	return tempObj;
}


int main(int argc, char *argv[])
//int main()
{
	
	int modeVal, n, m,val;
	char type_val;
	std::string fName1, fName2,fResult;
	unsigned int res_time;

	/*fName1 = "mtr1";
	fName2 = "mtr2";*/
	// Разбор командной строки 

	if (argc != 5) {
		std::cout << "invalid number of arguments" << std::endl;
		exit(1);
	}
	else {
		fName1 = argv[1];
		fName2 = argv[2];
		fResult = argv[3];
		modeVal = atoi(argv[4]);
		if (modeVal < 0 || modeVal>5) {
			std::cout << "invalid mode" << std::endl;
			exit(1);
		}
	}


	std::fstream f1(fName1, std::ios::binary | std::ios::in);
	f1.read((char*)&type_val, sizeof(char));
	f1.close();

	modeVal = atoi(argv[4]);

	if (type_val == 'f') {
		//std::cout << "lets f"<< std::endl;
		Matrix<float> Lmatr;
		Lmatr = binReader<float>(fName1);
		Matrix<float> Rmatr;
		Rmatr = binReader<float>(fName2);
		Matrix<float> Ematr;
		Ematr = mulMatrix<float>(Lmatr, Rmatr, modeVal,res_time);
		std::ofstream fout(argv[3], std::ios::binary | std::ios::out);
		if (fout.is_open()) {
			float temp_val;
			int newRow = Ematr.getRows();
			int newColmn = Ematr.getColumns();
			fout.write((char*)&type_val, sizeof(char));
			fout.write((char*)&newRow, sizeof(int));
			fout.write((char*)&newColmn, sizeof(int));
			for (int i = 0; i< newRow; ++i) {
				for (int j = 0; j<newColmn; ++j) {
					temp_val = Ematr.mtr[i][j];
					fout.write((char*)&temp_val, sizeof(float)); //*n*m
				}
			}
		}
		fout.close();
	}

	else if (type_val == 'd') {
		//std::cout<<"lets d"<<std::endl;
		Matrix<double> Lmatr;
		Lmatr = binReader<double>(fName1);
		Matrix<double> Rmatr;
		Rmatr = binReader<double>(fName2);
		Matrix<double> Ematr;
		Ematr = mulMatrix<double>(Lmatr, Rmatr, modeVal,res_time);
		std::ofstream fout(argv[3], std::ios::binary | std::ios::out);
		if (fout.is_open()) {
			int newRow = Ematr.getRows();
			int newColmn = Ematr.getColumns();
			fout.write((char*)&type_val, sizeof(char));
			fout.write((char*)&newRow, sizeof(int));
			fout.write((char*)&newColmn, sizeof(int));
			for (int i = 0; i< newRow; ++i) {
				for (int j = 0; j< newColmn; ++j) {
					fout.write((char*)&Ematr.mtr[i][j], sizeof(double)); //*n*m
				}
			}
		}
	}
	
	std::cout << res_time << std::endl;
    return 0;
}
