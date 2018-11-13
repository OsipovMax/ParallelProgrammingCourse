#include "Matrix.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <string>
#include <mpi.h>



using namespace std;
template <typename T>
double* mult(Matrix<T> &m, double* v ) {
	double *result;
	result = new double[m.getRows()];
	for (int i = 0; i<m.getRows(); ++i){
		for (int j = 0; j<m.getColumns(); ++j){
			result[i] += m.mtr[i][j] * v[j];
		}
	}

	
	return result;
	
}


template <typename T>
Matrix<T> binReader(const std::string &strr) {

	int modeVal, n, m;
	double val;
	char type_val;
	Matrix<T> tempObj;

	const char* cstr = strr.c_str();
	std::fstream f1(cstr, std::ios::binary | std::ios::in);
	f1.read((char*)&type_val, sizeof(char));
	f1.read((char*)&n, sizeof(int));
	f1.read((char*)&m, sizeof(int));
	tempObj.setRows(n);
	tempObj.setColumns(m);

	tempObj.mtr.resize(n);
	for (int i = 0; i < n; ++i) {
		tempObj.mtr[i].resize(m);
		for (int j = 0; j < m; j++) {
			f1.read((char*)&val, sizeof(double));
			tempObj.mtr[i][j] = val;
		}
	}
	f1.close();
	return tempObj;
}


template <typename T>
double* toBuf ( Matrix<T> &m , int trig ) {
	
	int inner_counter = 0;
	int numproces;
	double* buf = new double[m.getColumns() * m.getRows()];
	if (trig == 0) {
		for (int i = 0; i < m.getRows(); ++i){
			for (int j =0; j< m.getColumns(); ++j){
				buf[inner_counter] = m.mtr[i][j];
				inner_counter++;
			}
		}
	}
	else {
		MPI_Comm_size(MPI_COMM_WORLD, &numproces);
		int temp_size = (int)(m.getColumns()/numproces);
		cout << "temp_size =" << temp_size <<endl;
		for (int k = 0; k < numproces; ++k){
			for (int i = 0; i < m.getRows();++i){
				for (int j =0; j < temp_size; ++j){
					buf[inner_counter] = m.mtr[i][j+k*temp_size];
					inner_counter++;
				}
			}
		}
	}
	/*for (int i =0; i <m.getColumns() * m.getRows(); i++)
		cout << buf[i] << "  ";
	cout << endl;*/
	return buf;
}

template <typename T>
Matrix<T> toMatr(int r_size, int c_size, double *buf ) {
	int inner_counter = 0;

	Matrix<T> tempObj;
	tempObj.setRows(r_size);
	tempObj.setColumns(c_size);
	tempObj.mtr.resize(r_size);
	for ( int i = 0; i < r_size; i++ ){
		tempObj.mtr[i].resize(c_size);
		for ( int j = 0; j < c_size; j++ ){
			tempObj.mtr[i][j] = buf[inner_counter];
			inner_counter++;
		}
	}
	return tempObj;
}

int main(int argc, char *argv[])
{
	int triger;
	int retval;
	int left_limit, right_limit;
	int count_rows, count_columns;
	int n_size, m_size;
	int numproces,rank;
	double time_total = 0.0, time_max = 0.0;
	double time_begin = 0.0, time_end = 0.0;
	string file_A, file_B,file_C;


	file_A = argv[1];
	file_B = argv[2];
	file_C = argv[3];

	Matrix<double> matr;
	Matrix<double> vect;
	
	if ((retval = MPI_Init(&argc, &argv)) != MPI_SUCCESS) {
		cout << " MPI_Init error " << endl;
		MPI_Abort(MPI_COMM_WORLD, retval);
	}

	MPI_Comm_size(MPI_COMM_WORLD, &numproces);
	//cout << "numproc = " << numproces << endl;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	double *matr_buf, *vector_buf;
	double *proc_matr_buf, *proc_vect_buf;
	double *rb;
	double *mull_res_vect;
	
	if (rank == 0 ) {
		
		matr = binReader<double>(file_A);
		vect = binReader<double>(file_B);
		
		n_size = matr.getRows();
		m_size = matr.getColumns();

		/*	cout << " matr" << endl;
			cout << "rows_size" << matr.mtr.size() << endl;
			cout << "colm_size" << matr.mtr[0].size()<< endl;

			for (int p = 0 ; p < matr.mtr.size(); p++){
				for (int s = 0; s< matr.mtr[p].size(); s++){
					cout << matr.mtr[p][s] << "  ";
				}
				cout << endl;
			}
			cout << " vect" << endl;
			cout << "rows_size" << vect.mtr.size() << endl;
			cout << "colm_size" << vect.mtr[0].size()<< endl;

			for (int p = 0 ; p< 2; p++){
				for (int s = 0; s< 1; s++){
					cout << vect.mtr[p][s] << "  ";
				}
				cout << endl;
			} */

		
		if (matr.getRows() > matr.getColumns()){
			triger = 0;
			left_limit = rank*n_size/numproces;
			right_limit = (rank+1)*n_size/numproces - 1;
			count_rows = right_limit - left_limit +1 ;
		}
		else {
			triger = 1;
			left_limit = rank*m_size/numproces;
			right_limit = (rank+1)*m_size/numproces - 1;
			count_columns = right_limit - left_limit +1 ;
		}
		
	}
		MPI_Bcast(&triger, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		
		if ( triger == 0){
			cout << " zashel " << endl;
			MPI_Bcast(&n_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&m_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&count_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);

			vector_buf = new double[m_size];
			
			if ( rank == 0 ){
				//matr_buf = new double[m_size*n_size];

				matr_buf = toBuf<double>(matr,0);
				vector_buf = toBuf<double>(vect,0);
			}

			MPI_Bcast(vector_buf, m_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

			MPI_Barrier(MPI_COMM_WORLD);
			
			proc_matr_buf = new double [m_size*count_rows];

			MPI_Scatter(matr_buf, m_size*count_rows, MPI_DOUBLE, proc_matr_buf, m_size*count_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

			Matrix<double> mull_matr = toMatr<double>(count_rows,m_size,proc_matr_buf); 

			mull_res_vect = new double[count_rows]; 

			time_begin = MPI_Wtime();

			mull_res_vect = mult(mull_matr,vector_buf);  
		
			time_end = MPI_Wtime() - time_begin;

		
			if( rank!=0){
			MPI_Send(&time_end, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
			}
			
			delete proc_matr_buf;
			delete vector_buf;

			if (rank == 0) {
				delete matr_buf;
			}


			MPI_Status status;

			if ( rank == 0 ){
				time_total+=time_end;
				rb = new double [count_rows];
				if( time_end > time_max)
					time_max = time_end;
				for (int i = 0; i< numproces -1; ++i){
				MPI_Recv(&time_end, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				if( time_end > time_max)
					time_max = time_end;
				time_total+=time_end;
				}

			}
			MPI_Gather(mull_res_vect, count_rows , MPI_DOUBLE, rb, count_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			
			delete mull_res_vect;
			cout << "vishel" << endl;
		}
		else {
			cout<<"numproces = " << numproces << endl;

			MPI_Bcast(&n_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&m_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&count_columns, 1, MPI_INT, 0, MPI_COMM_WORLD);
			
			vector_buf= new double[m_size];
			if ( rank == 0 ){
				matr_buf = toBuf<double>(matr,1);
				vector_buf = toBuf<double>(vect,0);
			}
			
			proc_matr_buf = new double [n_size*count_columns];
			proc_vect_buf = new double [count_columns];

			MPI_Scatter(matr_buf, n_size*count_columns, MPI_DOUBLE, proc_matr_buf, n_size*count_columns, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Scatter(vector_buf, count_columns, MPI_DOUBLE, proc_vect_buf, count_columns, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			
			Matrix<double> mull_matr = toMatr<double>(n_size,count_columns, proc_matr_buf);

			time_begin = MPI_Wtime();
			double *mull_res_vect = mult(mull_matr,proc_vect_buf);
			time_end = MPI_Wtime() - time_begin;
			
			if( rank!=0){
			MPI_Send(&time_end, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
			}
			
			delete proc_matr_buf;
			delete proc_vect_buf;

			MPI_Status status;
			if ( rank == 0 ){
				time_total+=time_end;
				if (time_end > time_max)
					time_max = time_end;
				rb = new double [count_rows];
				for (int i = 0; i< numproces -1; ++i){
				MPI_Recv(&time_end, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				if( time_end > time_max)
					time_max = time_end;
				time_total += time_end;
				}
	
			}
			MPI_Reduce(mull_res_vect, rb, n_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	

	/*MPI_File fh;
	MPI_Info i;
	MPI_Info_create(&i);
	MPI_Status status;
		if ( rank == 0 ){
		cout<<"pes"<<endl;
		MPI_File_open(MPI_COMM_WORLD, "C.txt", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
		double pp = 3.3;
		MPI_File_write(fh, &pp, 1, MPI_DOUBLE, &status);
		MPI_File_close(&fh);
		}*/
	if (rank ==0){
		if (numproces == 1) {
			FILE *f = fopen( "rproc1", "w");
			fprintf(f, "%f %f %d\n", time_max, time_total, numproces);
			fclose(f);
		}
		else {
			float tm;
			double speedup, efficiency;
			cout << "mmm" << endl;
			FILE *f1 = fopen( "rproc1", "r");
			cout << "m333" << endl;
			fscanf(f1,"%f",&tm);
			fclose(f1);
			cout << "tm = " << tm << endl;
			speedup = (double)tm/time_max;
			efficiency = speedup/numproces;
			FILE *ss = fopen("result", "a");
			fprintf(ss, "%d %f %f %f %f \n", numproces, time_total, time_max, speedup, efficiency);
			fclose(ss);
		}
		cout << "result_vect" <<endl;
		for (int i = 0; i < 3; i++)
			cout << rb[i] << "  " ;
		
	}
	MPI_Finalize();
	return 0;
}
