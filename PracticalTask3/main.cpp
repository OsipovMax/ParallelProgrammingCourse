#include <mpi.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>


using namespace std;

int main (int argc, char* argv[]) {
	
	int lval, rval, rank;
	int temp_lval, temp_rval;
	int retval,numproces;
	int sqrt_rval;
	int counter = 0;
	
	double time_begin, time_end,time_max = 0.0, time_total =0.0;
	double buf[2];

	lval = atoi(argv[1]);
	rval = atoi(argv[2]);
	
	if (lval > rval) {
		cout <<"left value less than right value" << endl;
		exit(1);
	}
				
	if ( (  retval = MPI_Init(&argc,&argv)) != MPI_SUCCESS) {
		cout << " MPI_Init error " <<endl;
		MPI_Abort(MPI_COMM_WORLD,retval);
	}
	
	MPI_Comm_size(MPI_COMM_WORLD, &numproces);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	
	sqrt_rval = sqrt(rval);
	
	vector<int> sieve (rval + 1 , 1);
	sieve[0] = 0;
	sieve[1] = 0;

	temp_lval = rank*rval/numproces +1;
	temp_rval = (rank+1)*rval/numproces;
	
	vector<int> mass;
	
	time_begin = MPI_Wtime();
	
	for (int i = 2 ; i <= sqrt_rval; ++i){
		if(sieve[i] == 1)	
			for ( int j = 2*i ; j <= rval; j+=i) 
				sieve[j]=0;
	}


	for (int p = temp_lval ; p<=temp_rval; ++p){
		if (sieve[p] == 1){
			mass.push_back(p);
		}
	}

	if ( rank == 0 ) {
		
		time_end = MPI_Wtime() - time_begin;
		
		if ( time_end > time_max) 
			time_max = time_end;
		
		time_total += time_end;
	}	

	if ( rank!=0) {
		time_end = MPI_Wtime() - time_begin;
		buf[0] = time_end;
		buf[1] = (double)mass.size();
		MPI_Send(buf, 2, MPI_DOUBLE, 0 , rank,MPI_COMM_WORLD);
	}

	ofstream fout (argv[3], ios::app);
	for (int i = 0; i < mass.size(); ++i){
		fout<< mass[i] << " ";
	}
	fout.close(); 		
			
		
	if ( rank == 0 ) {
		
		counter = mass.size();
		
		MPI_Status status;
		for (int i = 0 ; i < numproces-1; ++i){
			MPI_Recv(buf,2,MPI_DOUBLE,MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			time_end = buf[0];
			if (time_end > time_max)
				time_max = time_end;
			time_total +=time_end;
			
			counter += buf[1];
		}
		cout << counter <<endl;
		ofstream fout_ ("result",ios::app);
		fout_ << numproces << " " << time_total << " " << time_max <<endl;		
	}
	vector<int>().swap(mass);
	vector<int>().swap(sieve);
	MPI_Finalize();
	return 0;
	
}