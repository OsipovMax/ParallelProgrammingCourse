//#include <mpi.h>
#include <omp.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <pthread.h>
#include <sys/types.h>
#include <cmath>
#include <time.h>

pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t condition = PTHREAD_COND_INITIALIZER;

using namespace std;

struct thread_info {
	pthread_t thread_id;
	int thread_num;
	char **argv_string;
};

struct result_struct {
	int res;
	double res_time;
};

void *sieve_function (void *arg){
	int lval, rval, rank;
	int temp_lval, temp_rval;
	int retval;
	int sqrt_rval;
	int counter = 0;
	struct thread_info *tinfo = (thread_info*)arg;
	result_struct *rStruct = new result_struct;
	double time_begin, time_end, time_max = 0.0, time_total = 0.0;
	double buf[2];
	int NumThread = atoi(tinfo->argv_string[4]);
	lval = atoi(tinfo->argv_string[1]);
	rval = atoi(tinfo->argv_string[2]);
	
	//if (lval > rval) {
		//cout << "left value less than right value" << endl;
		//exit(1);
	//}

	sqrt_rval = sqrt(rval);

	vector<int> sieve(sqrt_rval + 1, 1);
	sieve[0] = 0;
	sieve[1] = 0;

	temp_lval = tinfo->thread_num * rval / NumThread + 1;
	temp_rval = (tinfo->thread_num + 1)*rval / NumThread;
	//pthread_mutex_lock(&lock);
	//cout << "thread_lval:"<< tinfo->thread_num<< "---" << temp_lval<<endl;
	//cout << "thread_rval:"<< tinfo->thread_num<< "---" << temp_rval << endl;
	//pthread_mutex_unlock(&lock);
	vector<int> mass(temp_rval - temp_lval + 1, 1);
	//mass[0] = 0;

	//time_begin = MPI_Wtime();
	time_begin = omp_get_wtime();
	for (int i = 2; i <= sqrt_rval; ++i) {
		if (sieve[i] == 1)
			for (int j = 2 * i; j <= sqrt_rval; j += i)
				sieve[j] = 0;
	}

	int temp_val, cycle;

	for (int i = 2; i <= sqrt_rval; i++) {
		if (sieve[i] == 1) {
			temp_val = temp_lval % i;
			if (temp_val == 0)
				temp_val = i;
			for (int j = i - temp_val; j <= temp_rval - temp_lval; j += i) {
					if( temp_lval + j >= sqrt_rval +1){
						mass[j] = 0;
					}
					else {
						mass[j] = sieve[temp_lval + j];
					}
			}
		}
	}
	if (tinfo->thread_num == 0)
		mass[0]=0;
	vector<int>().swap(sieve);
	for (int i = 0; i <= temp_rval - temp_lval; ++i) {
		if (mass[i] == 1) {
			sieve.push_back(temp_lval + i);
		}
	}

	//time_end = MPI_Wtime() - time_begin;
	time_end = omp_get_wtime() - time_begin;
	
	pthread_mutex_lock(&lock);
	ofstream fout(tinfo->argv_string[3], ios::app);
	for (int i = 0; i < sieve.size(); ++i) {
		fout << sieve[i] << " ";
	}
	fout.close();
	pthread_mutex_unlock(&lock);

	rStruct->res = sieve.size();
	rStruct->res_time = time_end;
	vector<int>().swap(mass);
	vector<int>().swap(sieve);

	return rStruct;
}


int main(int argc, char **argv)
{
	int counter = 0;
	double counter_time;
	int NumThread = atoi(argv[4]);
	thread_info *tdInfo = new thread_info[NumThread];
	void *rtStruct;
	struct result_struct* res_struct;
	for (int num = 0; num < NumThread; num++){
		tdInfo[num].thread_num = num;
		tdInfo[num].argv_string = argv;
		pthread_create(&tdInfo[num].thread_id, NULL, (void*(*)(void*))sieve_function, &tdInfo[num]);
	}
	for (int num = 0; num < NumThread; num++){
		pthread_join(tdInfo[num].thread_id, &rtStruct);
		res_struct =(result_struct*)rtStruct;
		counter+=res_struct->res;
		counter_time+=res_struct->res_time;
	}
	cout << "number of primes = " << counter << endl;
	//cout << "result time = " << counter_time << endl;
	ofstream fout_("result", ios::app);
		fout_ << NumThread << " " << counter_time << endl;
	delete tdInfo;
	delete res_struct;
	return 0;
}
