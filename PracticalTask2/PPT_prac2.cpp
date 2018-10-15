#include <iostream>
#include "Matrix.h"
#include <fstream>
#include <papi.h>
#include <cstdlib>
#include <cmath>
#include <string>

#define NUM_EVENTS 2
const int cachesize = 4096;

using namespace std;
template <typename T>
Matrix<T> mulMatrix(const Matrix<T> &m1, const Matrix<T> &m2, const int mode, int f) {
	Matrix <T> resMtr;
	
	long long start_cycles, end_cycles;

	int counters = 0;

	int retval,rt, EventSet = PAPI_NULL;
	int Events[NUM_EVENTS] = { PAPI_L1_DCM, PAPI_L2_DCM};
	long long values[NUM_EVENTS];
	float real_time, proc_time,mflops;
	float ireal_time, iproc_time, imflops;
	long long iflpops;
	long long flpops;


	resMtr.mtr.resize(m1.getRows());
	for (int i = 0; i < resMtr.mtr.size(); ++i) {
		resMtr.mtr[i].resize(m2.getColumns());
	}
	resMtr.setRows(m1.getRows());
	resMtr.setColumns(m2.getColumns());
	
	if (f == 0) {
		start_cycles=PAPI_get_real_cyc();
		if( PAPI_create_eventset(&EventSet) != PAPI_OK){
			cout << "error -1 "<< endl;
			exit(1);
		}
		if (PAPI_add_events(EventSet,Events,NUM_EVENTS)!= PAPI_OK){
			cout << "error-2 " <<endl;
			exit(1);
		}
		if (PAPI_start(EventSet) != PAPI_OK) {
			cout << "error-3" <<endl;
			exit(1);
		}
	}
	

	
	if ( f == 1) {
		if (( rt = PAPI_flops( &ireal_time, &iproc_time, &iflpops, &imflops )) < PAPI_OK){
			printf("error__");
			exit(1);
		}
	}

	if (mode == 0) { // 32x32
		cout << "mode = " << mode << endl;
		int blocksize = 32;
		for (int i = 0; i < m1.getRows(); i += blocksize) {
			for (int j = 0; j < m1.getColumns(); j += blocksize) {
				for (int k = 0; k < m2.getColumns(); k += blocksize) {

					for (int i1 = i; i1 < i + blocksize && i1 < m1.getRows(); i1++) {
						for (int j1 = j; j1 < j + blocksize && j1 < m1.getColumns(); j1++) {
							for (int k1 = k; k1 < k + blocksize && k1 < m2.getColumns(); k1++) {
								resMtr.mtr[i1][j1] += m1.mtr[i1][k1] * m2.mtr[k1][j1];
							}
						}
					}
		
				}
			}
		}
	}
	else if (mode == 1) {
		int blocksize = 32;
		cout << "mode = " << mode <<endl;
		for (int i = 0; i < m1.getRows(); i += blocksize) {
			for (int k = 0; k < m2.getColumns(); k += blocksize) {
				for (int j = 0; j < m1.getColumns(); j += blocksize) {

					for (int i1 = i; i1 < i + blocksize && i1 < m1.getRows(); i1++) {
						for (int k1 = k; k1 < k + blocksize && k1 < m2.getColumns(); k1++) {
							for (int j1 = j; j1 < j + blocksize && j1 < m1.getColumns(); j1++) {
								resMtr.mtr[i1][j1] += m1.mtr[i1][k1] * m2.mtr[k1][j1];
							}
						}
					}
	
				}
			}
		}
	}
	
	if (mode == 2) {
		int lib = PAPI_library_init(PAPI_VER_CURRENT);
		if (lib == 0) 
			cout << "lib(";
		const PAPI_hw_info_t *info = PAPI_get_hardware_info();

		int blocksize = sqrt(cachesize/(3*sizeof(float)));
		
		cout <<"blocksize = "<< blocksize << endl;
		for(int i = 0; i < m1.getRows(); i += blocksize) {
			for (int k = 0; k < m2.getColumns(); k += blocksize) {
				for (int j = 0; j < m1.getColumns(); j += blocksize) {

					for (int i1 = i; i1 < i + blocksize && i1 < m1.getRows(); i1++) {
						for (int k1 = k; k1 < k + blocksize && k1 < m2.getColumns(); k1++) {
							for (int j1 = j; j1 < j + blocksize && j1 < m1.getColumns(); j1++) {
								resMtr.mtr[i1][j1] += m1.mtr[i1][k1] * m2.mtr[k1][j1];
							}
						}
					}

				}
			}
		}

	}
	
	if ( f == 1) {
		if (( rt = PAPI_flops(&real_time,&proc_time,&flpops,&mflops))<PAPI_OK){
			printf("retval: %d\n", retval);
			exit(1);
		}
		string num = to_string(mode);
		string RealTime = "RealTime" + num;
		string ProcTime = "ProcTime" + num;
		string MFLOPS = "Mflops" + num;
		ofstream fout (RealTime, ios::app);
		fout<< m1.getRows() << " " << real_time<< endl;
		fout.close(); 
		ofstream _fout (ProcTime,ios::app);
		_fout<< m1.getRows() << " " << proc_time<< endl;
		_fout.close();
		ofstream fout_ (MFLOPS,ios::app);
		fout_<< m1.getRows() << " " << mflops<< endl;
		fout_.close();

		//cout << "Real_time :" << real_time << " Proc_time :" << proc_time  <<" Total flpops : " << flpops << " MFLOPS : " << mflops <<endl;
	}

	if (f == 0) {
		end_cycles = PAPI_get_real_cyc();
		string num = to_string(mode);
		retval = PAPI_stop(EventSet,values);
		string L1_DCM = "L1_DCM" + num;
		string L2_DCM = "L2_DCM" + num;
		string Cycles = "Cycles" + num;
		ofstream fout (L1_DCM, ios::app);
		fout<< m1.getRows() << " " << values[0]<< endl;
		fout.close(); 
		ofstream _fout (L2_DCM,ios::app);
		_fout<< m1.getRows() << " " << values[1]<< endl;
		_fout.close();
		ofstream fcycles (Cycles,ios::app);
		fcycles << m1.getRows() << " " << end_cycles - start_cycles<< endl;
		fcycles.close();
		//cout << "L1_DCM = " << values[0] << endl;
		//cout << "L2_DCM = " << values[1] << endl;
	}

	//end_cycles = PAPI_get_real_usec();
	
	//printf("L1_DCM = %lld L2_DCM = %lld \n", values[0],values[1]);
	//printf("Real_time: %f Proc_time: %f Total flpops: %lld MFLOPS: %f\n,
		//real_time,proc_time,flpops,mflops);

	//cout << "Real_time :" << real_time << " Proc_time :" << proc_time 
	//<<" Total flpops : " << flpops << " MFLOPS : " << mflops <<endl; 
	
	//printf("Wallclock cycles: %lld\n", end_cycles - start_cycles);
	//cout << "Time_result = " << end_cycles - start_cycles<<endl;
	
	return resMtr;
}

template<typename T>
Matrix<T> binReader(const std::string &strr) {

	int modeVal, n, m;
	float val;
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
			f1.read((char*)&val, sizeof(float));
			tempObj.mtr[i][j] = val;
		}
	}
	f1.close();
	return tempObj;
}

int main(int argc, char *argv[])
{
	int modeVal, n, m, val, retval , flag;
	char type_val;
	std::string fName1, fName2, fResult;
	unsigned int res_time;

	/*fName1 = "mtr1";
	fName2 = "mtr2";*/
	flag = atoi(argv[5]);
	
	if (( retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT){
		cout << "ret: " <<retval << "\n curr ver: " << PAPI_VER_CURRENT<<endl;
		//fprint (stderr, "PAPI_library_init\n");
	}
	if (argc != 6) {
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

	const char* cstr_m = fName1.c_str();
	std::fstream f1(cstr_m, std::ios::binary | std::ios::in);
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
		Ematr = mulMatrix<float>(Lmatr, Rmatr, modeVal,flag/*, res_time*/);
		std::ofstream fout(argv[3], std::ios::binary | std::ios::out);
		if (fout.is_open()) {
			float temp_val;
			int newRow = Ematr.getRows();
			int newColmn = Ematr.getColumns();
			cout << "newRow =" << newRow <<endl;
			cout << "newCol =" << newColmn <<endl;
			fout.write((char*)&type_val, sizeof(char));
			fout.write((char*)&newRow, sizeof(int));
			fout.write((char*)&newColmn, sizeof(int));
			for (int i = 0; i < newRow; ++i) {
				for (int j = 0; j < newColmn; ++j) {
					temp_val = Ematr.mtr[i][j];
					fout.write((char*)&temp_val, sizeof(float)); //*n*m
				}
			}
		}
		fout.close();
	}
}
