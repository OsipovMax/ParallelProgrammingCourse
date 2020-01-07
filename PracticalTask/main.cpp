#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <mpi.h>

using namespace std;
MPI_Comm comm;
int* part;
int numprocess;
int rank;

const double PI = 3.14159265358979323;
const double Lx = 3.14159265358979323;
const double Ly = 3.14159265358979323;

double phi1Function(int i, int j, double h_x, double h_y) {
	double x = i * h_x;
	double y = j * h_y;
	return sin(x) * sin(y);
}

double u1Function(int i, int j, double h_x, double h_y, int n, double tau) {
	double x = i * h_x;
	double y = j * h_y;
	return (n * tau * n * tau + 1) *  sin(x) * sin(y);
}

double phi2Function(int i, int j, double h_x, double h_y) {
	double x = i * h_x;
	double y = j * h_y;
	return cos(2 * x) * cos(y);
}

double u2Function(int i, int j, double h_x, double h_y, int n, double tau) {
	double x = i * h_x;
	double y = j * h_y;
	return (n * tau * n * tau + 1) *  cos(2 * x) * cos(y);
}


double getError1(vector<vector<double> >& v, vector<vector<double> >& u, double h_x, double h_y, int n, double tau) {
	double temp = 0.0;
	double max = -1.0;
	int tempI = 0;
	int tempJ = 0;
	for (int i = 1; i < v.size() - 1; ++i) {
		for (int j = 1; j < v[i].size() - 1; ++j) {
			temp = fabs(v[i][j] - u[i][j]);
			if (max < temp) {
				max = temp;
			}
		}
	}
	return max;
}

double laplasCalc(vector<vector<double> > &v, int x, int y, double h_x, double h_y){ 
	return ((v[x - 1][y] - 2 * v[x][y] + v[x + 1][y]) / (h_x * h_x)) + 
			((v[x][y - 1] - 2 * v[x][y] + v[x][y + 1]) / (h_y * h_y));
}

double grad1(vector<vector<double> > &v1,vector<vector<double> > &v2,
		int x, int y, double h_x, double h_y){
	return ((v1[x - 1][y] - 2 * v1[x][y] + v1[x + 1][y]) / (h_x * h_x)) + 
			((v2[x + 1][y + 1] - v2[x + 1][y - 1] - v2[x - 1][y + 1] + v2[x - 1][y - 1]) / (4 * h_x * h_y));
}

double grad2(vector<vector<double> > &v1,vector<vector<double> > &v2,
		int x, int y, double h_x, double h_y){
	return ((v2[x][y - 1] - 2 * v2[x][y] + v2[x][y + 1]) / (h_y * h_y)) +
			((v1[x + 1][y + 1] - v1[x + 1][y - 1] - v1[x - 1][y + 1] + v1[x - 1][y - 1]) / (4 * h_x * h_y));		
}

//Determining the direction
vector<int> selectionDirrect(int coords[2], int procPerXAxis, int procPerYAxis){
	//dirrection[0] - up	dirrection[4] - top-left
	//dirrection[1] - down	dirrection[5] - top-right
	//dirrection[2] - left	dirrection[6] - bottom left
	//dirrection[3] - right	dirrection[7] - bottom-right
	// if dirrection[ind] == 0, the transfer is not performed
	vector<int> dirrection(8,1);
	
	if (coords[1] - 1 < 0){
		dirrection[0] = 0;
	}
	
	if (coords[1] + 1 == procPerYAxis){
		dirrection[1] = 0;
	}
	
	if (coords[0] - 1 < 0){
		dirrection[2] = 0;
	}
	
	if (coords[0] + 1 == procPerXAxis){
		dirrection[3] = 0;
	}
	
	if (dirrection[0] == 0){
		dirrection[4] = 0;
		dirrection[5] = 0;
	}

	if (dirrection[1] == 0){
		dirrection[6] = 0;
		dirrection[7] = 0;
	}
	
	if (dirrection[2] == 0){
		dirrection[4] = 0;
		dirrection[6] = 0;
	}
	
	if (dirrection[3] == 0){
		dirrection[5] = 0;
		dirrection[7] = 0;
	}
	return dirrection;
}
// The distribution of grid nodes in the process
vector<vector<double> > gridDistribution(int xSizeGlobalGrid, int ySizeGlobalGrid, 
	int processorPerXAxis, int processorPerYAxis, double h_x, double h_y, int n, double tau, int flag) {

	int numCells = xSizeGlobalGrid * ySizeGlobalGrid;
	int coordRank;
	int cnt = 0;
	int processorCoordinate[2];
	int procRank;
	int coords[2];

	MPI_Comm_size(comm, &numprocess);
	MPI_Comm_rank(comm, &procRank);

	int ib, ie, jb, je;
	MPI_Cart_coords(comm, procRank, 2, coords);
	
	ib = coords[0] * (xSizeGlobalGrid / processorPerXAxis);
	ie = (coords[0] + 1) * (xSizeGlobalGrid / processorPerXAxis);

	jb = coords[1] * (ySizeGlobalGrid / processorPerYAxis);
	je = (coords[1] + 1) * (ySizeGlobalGrid / processorPerYAxis);
	
	//adding elements on XAxis 

	if ((coords[0] + 1 == processorPerXAxis) && (xSizeGlobalGrid % processorPerXAxis != 0)){
		ib = coords[0] * (xSizeGlobalGrid / processorPerXAxis);
		ie = xSizeGlobalGrid;
		if ((coords[1] + 1 == processorPerYAxis) && (ySizeGlobalGrid % processorPerYAxis != 0)){
			jb = coords[1] * (ySizeGlobalGrid / processorPerYAxis);
			je = ySizeGlobalGrid;
		}
	}
	
	//adding elements on YAxis
	
	if ((coords[1] + 1 == processorPerYAxis) && (ySizeGlobalGrid % processorPerYAxis != 0)){
		jb = coords[1] * (ySizeGlobalGrid / processorPerYAxis);
	    	je = ySizeGlobalGrid;
		if ((coords[0] + 1 == processorPerXAxis) && (xSizeGlobalGrid % processorPerXAxis != 0)){
			ib = coords[0] * (xSizeGlobalGrid / processorPerXAxis);
			ie = xSizeGlobalGrid;
		}
	}

	vector<vector<double> > res;
	res.resize(je - jb + 2);
	for (int i = 0; i < res.size(); ++i){
		res[i].resize(ie - ib + 2);
		for (int j = 0; j < res[i].size(); ++j){
			res[i][j] = -1.0;
		}
	}
	
	// v1 distribution flag = 0
	// v2 distribution flag = 1
	// u1 distribution flag = 2
	// u2 distribution flag = 3
	int iInd = jb - 1;
	int jInd = ib - 1;
	if (flag == 0){
		for (int i = 1; i < je - jb + 1; ++i){
			iInd++;
			for (int j = 1; j < ie - ib + 1; ++j){
				jInd++;	
				res[i][j] = phi1Function(iInd, jInd, h_x, h_y);
			}
			jInd = ib - 1;
		}
	}
	if (flag == 1){
		for (int i = 1; i < je - jb + 1; ++i){
			iInd++;
			for (int j = 1; j < ie - ib + 1; ++j){
				jInd++;	
				res[i][j] = phi2Function(iInd, jInd, h_x, h_y);
			}
			jInd = ib - 1;
		}
	}
	if (flag == 2){
		for (int i = 1; i < je - jb + 1; ++i){
			iInd++;
			for (int j = 1; j < ie - ib + 1; ++j){
				jInd++;	
				res[i][j] = u1Function(iInd, jInd, h_x, h_y, n, tau);
			}
			jInd = ib - 1;
		}
		return res;
	}
	if (flag == 3){
		for (int i = 1; i < je - jb + 1; ++i){
			iInd++;
			for (int j = 1; j < ie - ib + 1; ++j){
				jInd++;	
				res[i][j] = u2Function(iInd, jInd, h_x, h_y, n, tau);
			}
			jInd = ib - 1;
		}
	}
	return res;
}

//-------------------Identify the different transfer of boundary values-----------------------------------
void sendBottomSide(vector<vector<double> > &v, int coords[2], int procPerXAxis, int procPerYAxis){
	int nbCoords[2];
	int nbRank;
	MPI_Request sendReq[1];
	int bufferSize = v[0].size() - 2;
	double sendBuf[bufferSize];
	nbCoords[0] = coords[0];
	nbCoords[1] = 0;
	MPI_Cart_rank(comm, nbCoords, &nbRank);
	for (int i = 0; i < bufferSize; ++i){
			sendBuf[i] = v[v.size() - 3][i + 1]; 
		}
	MPI_Isend(sendBuf, bufferSize, MPI_DOUBLE, nbRank, 0, comm, &sendReq[0]);
	MPI_Wait(&sendReq[0], MPI_STATUSES_IGNORE);
	if (coords[0] == 0){
		nbCoords[0] = coords[0] + 1;
		nbCoords[1] = 0;
		MPI_Cart_rank(comm, nbCoords, &nbRank);
		MPI_Isend(&v[v.size() - 3][v[v.size() - 3].size() - 2], 1, MPI_DOUBLE, nbRank, 0, comm, &sendReq[0]);
		MPI_Wait(&sendReq[0], MPI_STATUSES_IGNORE);
		return;
	}
	if (coords[0] == procPerXAxis - 1){
		nbCoords[0] = coords[0]  - 1;
		nbCoords[1] = 0;
		MPI_Cart_rank(comm, nbCoords, &nbRank);
		MPI_Isend(&v[v.size() - 3][1], 1, MPI_DOUBLE, nbRank, 0, comm, &sendReq[0]);
		MPI_Wait(&sendReq[0], MPI_STATUSES_IGNORE);
		return;	
	}
	nbCoords[0] = coords[0] - 1;
	nbCoords[1] = 0;
	MPI_Cart_rank(comm, nbCoords, &nbRank);
    	MPI_Isend(&v[v.size() - 3][1], 1, MPI_DOUBLE, nbRank, 0, comm, &sendReq[0]);
	MPI_Wait(&sendReq[0], MPI_STATUSES_IGNORE);
	
	nbCoords[0] = coords[0] + 1;
	nbCoords[1] = 0;
	MPI_Cart_rank(comm, nbCoords, &nbRank);
    	MPI_Isend(&v[v.size() - 3][v[v.size() - 3].size() - 2], 1, MPI_DOUBLE, nbRank, 0, comm, &sendReq[0]);
	MPI_Wait(&sendReq[0], MPI_STATUSES_IGNORE);
}

void recvBottomSide(vector<vector<double> > &v, int coords[2], int procPerXAxis, int procPerYAxis){
	int nbCoords[2];
	int nbRank;
	MPI_Request recvReq[1];
	int bufferSize = v[0].size() - 2;
	double recvBuf[bufferSize];
	nbCoords[0] = coords[0];
	nbCoords[1] = 0;
	MPI_Cart_rank(comm, nbCoords, &nbRank);
	MPI_Irecv(recvBuf, bufferSize, MPI_DOUBLE, nbRank, 0, comm, &recvReq[0]);
	MPI_Wait(&recvReq[0], MPI_STATUSES_IGNORE);
	for (int i = 0; i < bufferSize; ++i){
		v[v.size() - 2][i + 1] = recvBuf[i];
	}
}

void sendTopSide(vector<vector<double> > &v, int coords[2], int procPerXAxis, int procPerYAxis){
	int nbCoords[2];
	int nbRank;
	MPI_Request sendReq[1];
	int bufferSize = v[0].size() - 2;
	double sendBuf[bufferSize];
	nbCoords[0] = coords[0];
	nbCoords[1] = procPerYAxis - 1;
	MPI_Cart_rank(comm, nbCoords, &nbRank);
	for (int i = 0; i < bufferSize; ++i){
		sendBuf[i] = v[1][i + 1];
	}
	MPI_Isend(sendBuf, bufferSize, MPI_DOUBLE, nbRank, 0, comm, &sendReq[0]);
	MPI_Wait(&sendReq[0], MPI_STATUSES_IGNORE);
}

void recvTopSide(vector<vector<double> > &v, int coords[2], int procPerXAxis, int procPerYAxis){
	int nbCoords[2];
	int nbRank;
	MPI_Request recvReq[1];
	int bufferSize = v[0].size() - 2;
	double recvBuf[bufferSize];
	nbCoords[0] = coords[0];
	nbCoords[1] = procPerYAxis - 1;
	MPI_Cart_rank(comm, nbCoords, &nbRank);
	MPI_Irecv(recvBuf, bufferSize, MPI_DOUBLE, nbRank, 0, comm, &recvReq[0]);
	MPI_Wait(&recvReq[0], MPI_STATUSES_IGNORE);
	for (int i = 0; i < bufferSize; ++i){
		v[0][i + 1] = recvBuf[i];
	}
	if (coords[0] == 0){
		nbCoords[0] = coords[0] + 1;
		nbCoords[1] = procPerYAxis - 1;
		MPI_Cart_rank(comm, nbCoords, &nbRank);
		MPI_Irecv(&v[0][v[0].size() - 1], 1, MPI_DOUBLE, nbRank, 0, comm, &recvReq[0]);
		MPI_Wait(&recvReq[0], MPI_STATUSES_IGNORE);
		return;
	}
	if (coords[0] == procPerXAxis - 1){
		nbCoords[0] = coords[0]  - 1;
		nbCoords[1] = procPerYAxis - 1;
		MPI_Cart_rank(comm, nbCoords, &nbRank);
		MPI_Irecv(&v[0][0], 1, MPI_DOUBLE, nbRank, 0, comm, &recvReq[0]);
		MPI_Wait(&recvReq[0], MPI_STATUSES_IGNORE);
		return;	
	}
	nbCoords[0] = coords[0] - 1;
	nbCoords[1] = procPerYAxis - 1;
	MPI_Cart_rank(comm, nbCoords, &nbRank);
        MPI_Irecv(&v[0][0], 1, MPI_DOUBLE, nbRank, 0, comm, &recvReq[0]);
	MPI_Wait(&recvReq[0], MPI_STATUSES_IGNORE);
	
	nbCoords[0] = coords[0] + 1;
	nbCoords[1] = procPerYAxis - 1;
	MPI_Cart_rank(comm, nbCoords, &nbRank);
        MPI_Irecv(&v[0][v[0].size() - 1], 1, MPI_DOUBLE, nbRank, 0, comm, &recvReq[0]);
	MPI_Wait(&recvReq[0], MPI_STATUSES_IGNORE);
}
//-----------------Calculation of values on the borders------------------------------------------------
void boundaryConditionV2 (vector<vector<double> > &v1,vector<vector<double> > &v2, 
	int coords[2], int procPerXAxis, int procPerYAxis, double h_x, double h_y, double tau) {
	vector<int> dirrections;
	dirrections = selectionDirrect(coords, procPerXAxis, procPerYAxis);
	vector<vector<double> > next2;
	next2 = v2;
	
	if (dirrections[0] == 0 && dirrections[1] == 1
		&& dirrections[2] == 1 && dirrections[3] == 1){
		recvTopSide(v2, coords, procPerXAxis, procPerYAxis);
		recvTopSide(v1, coords, procPerXAxis, procPerYAxis);
		#pragma omp parallel for
		for (int i = 1; i < next2[1].size() - 1; ++i){ 
			next2[1][i] = tau * (laplasCalc(v2, 1, i, h_x, h_y) + grad2(v1, v2, 1, i, h_x, h_y)) + v2[1][i];
		}
		sendTopSide(next2, coords, procPerXAxis, procPerYAxis);
		v2 = next2;
		return;
	}

	if (dirrections[0] == 1 && dirrections[1] == 0
		&& dirrections[2] == 1 && dirrections[3] == 1){
		sendBottomSide(v2, coords, procPerXAxis, procPerYAxis);
		sendBottomSide(v1, coords, procPerXAxis, procPerYAxis);
		recvBottomSide(next2, coords, procPerXAxis, procPerYAxis);
		v2 = next2;
		return;
	}
	
	if (dirrections[0] == 1 && dirrections[1] == 1
		&& dirrections[2] == 0 && dirrections[3] == 1){
		#pragma omp parallel for
		for (int i = 1; i < next2.size() - 1; ++i){
			next2[i][1] = (4.0 * next2[i][2] - next2[i][3]) / 3.0;
		}
		v2 = next2;
		return;
	}

	if (dirrections[0] == 1 && dirrections[1] == 1
		&& dirrections[2] == 1 && dirrections[3] == 0){
		#pragma omp parallel for
		for (int i = 1; i < next2.size() - 1; ++i){
			next2[i][next2[i].size() - 2] = ( 4.0 * next2[i][next2[i].size() - 3] - next2[i][next2[i].size() - 4]) / 3.0;
		}
		v2 = next2;
		return;
	}

	if (dirrections[0] == 0 && dirrections[1] == 1
		&& dirrections[2] == 0 && dirrections[3] == 1){
		recvTopSide(v2, coords, procPerXAxis, procPerYAxis);
		recvTopSide(v1, coords, procPerXAxis, procPerYAxis);
		#pragma omp parallel for
		for (int i = 2; i < next2[1].size() - 1; ++i){
			next2[1][i] = tau * (laplasCalc(v2, 1, i, h_x, h_y) + grad2(v1, v2, 1, i, h_x, h_y)) + v2[1][i];
		}
		sendTopSide(next2, coords, procPerXAxis, procPerYAxis);
		#pragma omp parallel for
		for (int i = 1; i < next2.size() - 1; ++i){
			next2[i][1] = (4.0 * next2[i][2] - next2[i][3]) / 3.0;
		}
		v2 = next2;
		return;
	}

	if (dirrections[0] == 0 && dirrections[1] == 1
		&& dirrections[2] == 1 && dirrections[3] == 0){
		recvTopSide(v2, coords, procPerXAxis, procPerYAxis);
		recvTopSide(v1, coords, procPerXAxis, procPerYAxis);
		#pragma omp parallel for
		for (int i = 1; i < next2[1].size() - 2; ++i){
			next2[1][i] = tau * (laplasCalc(v2, 1, i, h_x, h_y) + grad2(v1, v2, 1, i, h_x, h_y)) + v2[1][i];
		}
		sendTopSide(next2, coords, procPerXAxis, procPerYAxis);
		#pragma omp parallel for
		for (int i = 1; i < next2.size() - 1; ++i){
			next2[i][next2[i].size() - 2] = ( 4.0 * next2[i][next2[i].size() - 3] - next2[i][next2[i].size() - 4]) / 3.0;
		}
		v2 = next2;
		return;
	}

	if (dirrections[0] == 1 && dirrections[1] == 0
		&& dirrections[2] == 0 && dirrections[3] == 1){
		sendBottomSide(v2, coords, procPerXAxis, procPerYAxis);
		sendBottomSide(v1, coords, procPerXAxis, procPerYAxis);
		recvBottomSide(next2, coords, procPerXAxis, procPerYAxis);
		#pragma omp parallel for
		for (int i = 1; i < next2.size() - 1; ++i){
			next2[i][1] = (4.0 * next2[i][2] - next2[i][3]) / 3.0;
		}
		v2 = next2;
		return;
	}

	if (dirrections[0] == 1 && dirrections[1] == 0
		&& dirrections[2] == 1 && dirrections[3] == 0){
		sendBottomSide(v2, coords, procPerXAxis, procPerYAxis);
		sendBottomSide(v1, coords, procPerXAxis, procPerYAxis);
		recvBottomSide(next2, coords, procPerXAxis, procPerYAxis);
		#pragma omp parallel for
		for (int i = 1; i < next2.size() - 1; ++i){
			next2[i][next2[i].size() - 2] = ( 4.0 * next2[i][next2[i].size() - 3] - next2[i][next2[i].size() - 4]) / 3.0;
		}
		v2 = next2;
		return;
	}
	

	if (dirrections[0] == 0 && dirrections[1] == 1
		&& dirrections[2] == 0 && dirrections[3] == 0){
		recvTopSide(v2, coords, procPerXAxis, procPerYAxis);
		recvTopSide(v1, coords, procPerXAxis, procPerYAxis);
		#pragma omp parallel for
		for (int i = 2; i < next2[1].size() - 2; ++i){
			next2[1][i] = tau * (laplasCalc(v2, 1, i, h_x, h_y) + grad2(v1, v2, 1, i, h_x, h_y)) + v2[1][i];
		}
		sendTopSide(next2, coords, procPerXAxis, procPerYAxis);
		#pragma omp parallel for
		for (int i = 1; i < next2.size() - 1; ++i){
			next2[i][1] = (4.0 * next2[i][2] - next2[i][3]) / 3.0;
			next2[i][next2[i].size() - 2] = ( 4.0 * next2[i][next2[i].size() - 3] - next2[i][next2[i].size() - 4]) / 3.0;
		}
		v2 = next2;
		return;
	}

	if (dirrections[0] == 1 && dirrections[1] == 0
		&& dirrections[2] == 0 && dirrections[3] == 0){
		sendBottomSide(v2, coords, procPerXAxis, procPerYAxis);
		sendBottomSide(v1, coords, procPerXAxis, procPerYAxis);
		recvBottomSide(next2, coords, procPerXAxis, procPerYAxis);
		#pragma omp parallel for
		for (int i = 1; i < next2.size() - 1; ++i){
			next2[i][1] = (4.0 * next2[i][2] - next2[i][3]) / 3.0;
			next2[i][next2[i].size() - 2] = ( 4.0 * next2[i][next2[i].size() - 3] - next2[i][next2[i].size() - 4]) / 3.0;
		}
		v2 = next2;
		return;
	}

	if (dirrections[0] == 0 && dirrections[1] == 0
		&& dirrections[2] == 0 && dirrections[3] == 1){
		sendBottomSide(v2, coords, procPerXAxis, procPerYAxis);
		sendBottomSide(v1, coords, procPerXAxis, procPerYAxis);
		recvTopSide(v2, coords, procPerXAxis, procPerYAxis);
		recvTopSide(v1, coords, procPerXAxis, procPerYAxis);
		#pragma omp parallel for
		for (int i = 2; i < next2[1].size() - 1; ++i){
			next2[1][i] = tau * (laplasCalc(v2, 1, i, h_x, h_y) + grad2(v1, v2, 1, i, h_x, h_y)) + v2[1][i];
		}
		sendTopSide(next2, coords, procPerXAxis, procPerYAxis);
		recvBottomSide(next2, coords, procPerXAxis, procPerYAxis);
		#pragma omp parallel for
		for (int i = 1; i < next2.size() - 1; ++i){
			next2[i][1] = (4.0 * next2[i][2] - next2[i][3]) / 3.0;
		}
		v2 = next2;
		return;
	}

	if (dirrections[0] == 0 && dirrections[1] == 0
		&& dirrections[2] == 1 && dirrections[3] == 0){
		sendBottomSide(v2, coords, procPerXAxis, procPerYAxis);
		sendBottomSide(v1, coords, procPerXAxis, procPerYAxis);
		recvTopSide(v2, coords, procPerXAxis, procPerYAxis);
		recvTopSide(v1, coords, procPerXAxis, procPerYAxis);
		#pragma omp parallel for
		for (int i = 1; i < next2[1].size() - 2; ++i){
			next2[1][i] = tau * (laplasCalc(v2, 1, i, h_x, h_y) + grad2(v1, v2, 1, i, h_x, h_y)) + v2[1][i];
		}
		sendTopSide(next2, coords, procPerXAxis, procPerYAxis);
		recvBottomSide(next2, coords, procPerXAxis, procPerYAxis);
		#pragma omp parallel for
		for (int i = 1; i < next2.size() - 1; ++i){
			next2[i][next2[i].size() - 2] = ( 4.0 * next2[i][next2[i].size() - 3] - next2[i][next2[i].size() - 4]) / 3.0;
		}
		v2 = next2;
		return;
	}

	if (dirrections[0] == 0 && dirrections[1] == 0
		&& dirrections[2] == 1 && dirrections[3] == 1){
		sendBottomSide(v2, coords, procPerXAxis, procPerYAxis);
		sendBottomSide(v1, coords, procPerXAxis, procPerYAxis);
		recvTopSide(v2, coords, procPerXAxis, procPerYAxis);
		recvTopSide(v1, coords, procPerXAxis, procPerYAxis);
		#pragma omp parallel for
		for (int i = 1; i < next2[1].size() - 1; ++i){
			next2[1][i] = tau * (laplasCalc(v2, 1, i, h_x, h_y) + grad2(v1, v2, 1, i, h_x, h_y)) + v2[1][i];
		}
		sendTopSide(next2, coords, procPerXAxis, procPerYAxis);
		recvBottomSide(next2, coords, procPerXAxis, procPerYAxis);
		v2 = next2;
		return;
	}

	if (dirrections[0] == 1 && dirrections[1] == 1
		&& dirrections[2] == 0 && dirrections[3] == 0){
		#pragma omp parallel for
		for (int i = 1; i < next2.size() - 1; ++i){
			next2[i][1] = (4.0 * next2[i][2] - next2[i][3]) / 3.0;
			next2[i][next2[i].size() - 2] = ( 4.0 * next2[i][next2[i].size() - 3] - next2[i][next2[i].size() - 4]) / 3.0;
		}
		v2 = next2;
		return;
	}

	if (dirrections[0] == 1 && dirrections[1] == 1
		&& dirrections[2] == 1 && dirrections[3] == 1){
		return;
	}
		
	if (dirrections[0] == 0 && dirrections[1] == 0 
		&& dirrections[2] == 0 && dirrections[3] == 0){
		for (int i = 1; i < v2.size() - 1; ++i){
			v2[0][i] = v2[v2.size() - 3][i];
			v1[0][i] = v1[v1.size() - 3][i];
		}
		for (int i = 2; i < next2[1].size() - 2; ++i){
			next2[1][i] = tau * (laplasCalc(v2, 1, i, h_x, h_y) + grad2(v1, v2, 1, i, h_x, h_y)) + v2[1][i];
		}
		for (int i = 2; i < next2[1].size() - 2; ++i){
			next2[next2.size() - 2][i] = next2[1][i];
		}
		for (int i = 1; i < next2.size() - 1; ++i){
			next2[i][1] = (4.0 * next2[i][2] - next2[i][3]) / 3.0;
			next2[i][next2[i].size() - 2] = ( 4.0 * next2[i][next2[i].size() - 3] - next2[i][next2[i].size() - 4]) / 3.0;			
		}

		v2 = next2;
		return;
	}
}
//-------------------------Updating halo cells-------------------------
void ExchangeHalo(int procRank, vector<vector<double> > &localFragment, int procPerXAxis, int procPerYAxis){
	
	int coords[2];
	vector<int> dirrections;
	MPI_Cart_coords(comm, procRank, 2, coords);
	dirrections = selectionDirrect(coords, procPerXAxis, procPerYAxis);

	//top-left
	if (dirrections[4] == 1){
		int nbCoords[2];
		int nbRank;
		MPI_Request sendReq[1], recvReq[1];
		nbCoords[0] = coords[0] - 1;
		nbCoords[1] = coords[1] - 1;
		MPI_Cart_rank(comm, nbCoords, &nbRank);
		MPI_Isend(&localFragment[1][1], 1, MPI_DOUBLE, nbRank, 0, comm, &sendReq[0]);
		MPI_Irecv(&localFragment[0][0], 1, MPI_DOUBLE, nbRank, 0, comm, &recvReq[0]);
		MPI_Wait(&sendReq[0], MPI_STATUSES_IGNORE);
		MPI_Wait(&recvReq[0], MPI_STATUSES_IGNORE);
	}
	
	//bottom - right
	if (dirrections[7] == 1){
		int nbCoords[2];
		int nbRank;
		MPI_Request sendReq[1], recvReq[1];
		nbCoords[0] = coords[0] + 1;
		nbCoords[1] = coords[1] + 1;
		MPI_Cart_rank(comm, nbCoords, &nbRank);
		MPI_Isend(&localFragment[localFragment.size() - 2][localFragment[localFragment.size() - 2].size() - 2], 
			1, MPI_DOUBLE, nbRank, 0, comm, &sendReq[0]);
		MPI_Irecv(&localFragment[localFragment.size() - 1][localFragment[localFragment.size() - 1].size() - 1], 
			1, MPI_DOUBLE, nbRank, 0, comm, &recvReq[0]);
		MPI_Wait(&sendReq[0], MPI_STATUSES_IGNORE);
		MPI_Wait(&recvReq[0], MPI_STATUSES_IGNORE);
	}

	//top - right
	if (dirrections[5] == 1){
		int nbCoords[2];
		int nbRank;
		MPI_Request sendReq[1], recvReq[1];
		nbCoords[0] = coords[0] + 1;
		nbCoords[1] = coords[1] - 1;
		MPI_Cart_rank(comm, nbCoords, &nbRank);
		MPI_Isend(&localFragment[1][localFragment[1].size() - 2], 1, MPI_DOUBLE, nbRank, 0, comm, &sendReq[0]);
		MPI_Irecv(&localFragment[0][localFragment[1].size() - 1], 1, MPI_DOUBLE, nbRank, 0, comm, &recvReq[0]);
		MPI_Wait(&sendReq[0], MPI_STATUSES_IGNORE);
		MPI_Wait(&recvReq[0], MPI_STATUSES_IGNORE);
	}

	//bot - left
	if (dirrections[6] == 1){
		int nbCoords[2];
		int nbRank;
		MPI_Request sendReq[1], recvReq[1];
		nbCoords[0] = coords[0] - 1;
		nbCoords[1] = coords[1] + 1;
		MPI_Cart_rank(comm, nbCoords, &nbRank);
		MPI_Isend(&localFragment[localFragment.size() - 2][1], 1, MPI_DOUBLE, nbRank, 0, comm, &sendReq[0]);
		MPI_Irecv(&localFragment[localFragment.size() - 1][0], 1, MPI_DOUBLE, nbRank, 0, comm, &recvReq[0]);
		MPI_Wait(&sendReq[0], MPI_STATUSES_IGNORE);
		MPI_Wait(&recvReq[0], MPI_STATUSES_IGNORE);
	}
	
	//down
	if (dirrections[1] == 1){
		int nbCoords[2];
		int nbRank;
		MPI_Request sendReq[1], recvReq[1];
		nbCoords[0] = coords[0];
		nbCoords[1] = coords[1] + 1;
		MPI_Cart_rank(comm, nbCoords, &nbRank);
		int buffersSize = localFragment[0].size() - 2;
		double sendBuf_down[buffersSize];
		double recvBuf_down[buffersSize];

		for (int i = 0; i < buffersSize; ++i){
			sendBuf_down[i] = localFragment[localFragment.size() - 2][i + 1];
		}
		MPI_Isend(sendBuf_down, buffersSize, MPI_DOUBLE, nbRank, 0, comm, &sendReq[0]);
		MPI_Irecv(recvBuf_down, buffersSize, MPI_DOUBLE, nbRank, 0, comm, &recvReq[0]);
		
		MPI_Wait(&sendReq[0], MPI_STATUSES_IGNORE);
		MPI_Wait(&recvReq[0], MPI_STATUSES_IGNORE);

		for (int i = 0; i < buffersSize; ++i){
			localFragment[localFragment.size() - 1][i + 1] = recvBuf_down[i];
		}
	}
	
	//up
	if (dirrections[0] == 1){
		int nbCoords[2];
		int nbRank;
		MPI_Request sendReq[1], recvReq[1];
		nbCoords[0] = coords[0];
		nbCoords[1] = coords[1] - 1;
		MPI_Cart_rank(comm, nbCoords, &nbRank);
		int buffersSize = localFragment[0].size() - 2;
		double sendBuf_up[buffersSize];
		double recvBuf_up[buffersSize];

		for (int i = 0; i < buffersSize; ++i){
			sendBuf_up[i] = localFragment[1][i + 1];
		}
		MPI_Isend(sendBuf_up, buffersSize, MPI_DOUBLE, nbRank, 0, comm, &sendReq[0]);
		MPI_Irecv(recvBuf_up, buffersSize, MPI_DOUBLE, nbRank, 0, comm, &recvReq[0]);

		MPI_Wait(&sendReq[0], MPI_STATUSES_IGNORE);
		MPI_Wait(&recvReq[0], MPI_STATUSES_IGNORE);

		for (int i = 0; i < buffersSize; ++i){
			localFragment[0][i + 1] = recvBuf_up[i];
		}
	}

	//left
	if (dirrections[2] == 1){
		int nbCoords[2];
		int nbRank;
		MPI_Request sendReq[1], recvReq[1];
		nbCoords[0] = coords[0] - 1;
		nbCoords[1] = coords[1];
		MPI_Cart_rank(comm, nbCoords, &nbRank);
		int buffersSize = localFragment.size() - 2;
		double sendBuf_left[buffersSize];
		double recvBuf_left[buffersSize];

		for (int i = 0; i < buffersSize; ++i){
			sendBuf_left[i] = localFragment[i + 1][1];
		}
		MPI_Isend(sendBuf_left, buffersSize, MPI_DOUBLE, nbRank, 0, comm, &sendReq[0]);
		MPI_Irecv(recvBuf_left, buffersSize, MPI_DOUBLE, nbRank, 0, comm, &recvReq[0]);


		MPI_Wait(&sendReq[0], MPI_STATUSES_IGNORE);
		MPI_Wait(&recvReq[0], MPI_STATUSES_IGNORE);

		for (int i = 0; i < buffersSize; ++i){
			localFragment[i + 1][0] = recvBuf_left[i];
		}
	}

	//right
	if (dirrections[3] == 1){
		int nbCoords[2];
		int nbRank;
		MPI_Request sendReq[1], recvReq[1];
		nbCoords[0] = coords[0] + 1;
		nbCoords[1] = coords[1];
		MPI_Cart_rank(comm, nbCoords, &nbRank);
		int buffersSize = localFragment.size() - 2;
		double sendBuf_right[buffersSize];
		double recvBuf_right[buffersSize];

		for (int i = 0; i < buffersSize; ++i){
			sendBuf_right[i] = localFragment[i + 1][localFragment[0].size() - 2];
		}
		MPI_Isend(sendBuf_right, buffersSize, MPI_DOUBLE, nbRank, 0, comm, &sendReq[0]);
		MPI_Irecv(recvBuf_right, buffersSize, MPI_DOUBLE, nbRank, 0, comm, &recvReq[0]);

		MPI_Wait(&sendReq[0], MPI_STATUSES_IGNORE);
		MPI_Wait(&recvReq[0], MPI_STATUSES_IGNORE);

		for (int i = 0; i < buffersSize; ++i){
			localFragment[i + 1][localFragment[0].size() - 1] = recvBuf_right[i];
		}
	}

}

//-----------------------------------Calculation of values in the Central nodes of the grid-------------------
void excuteCalculationNew (vector<vector<double> > &v1, vector<vector<double> > &v2,
			vector<vector<double> > &next1, vector<vector<double> > &next2,
			int coords[2], int procPerXAxis, int procPerYAxis, double h_x, double h_y, double tau) {
	vector<int> dirrections;
	dirrections = selectionDirrect(coords, procPerXAxis, procPerYAxis);
	next1 = v1;
	next2 = v2;

	if (dirrections[0] == 0 && dirrections[1] == 1
		&& dirrections[2] == 1 && dirrections[3] == 1){
		#pragma omp parallel for
		for (int i = 2; i < v1.size() - 1; ++i){
			for (int j = 1; j < v1[i].size() - 1; ++j){
				next1[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad1(v1, v2, i, j, h_x, h_y)) + v1[i][j];
				next2[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad2(v1, v2, i, j, h_x, h_y)) + v2[i][j];
			}
		}
		return;
	}

	if (dirrections[0] == 1 && dirrections[1] == 0
		&& dirrections[2] == 1 && dirrections[3] == 1){
		#pragma omp parallel for
		for (int i = 1; i < v1.size() - 2; ++i){
			for (int j = 1; j < v1[i].size() - 1; ++j){
				next1[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad1(v1, v2, i, j, h_x, h_y)) + v1[i][j];
				next2[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad2(v1, v2, i, j, h_x, h_y)) + v2[i][j];
			}
		}
		return;
	}

	if (dirrections[0] == 1 && dirrections[1] == 1
		&& dirrections[2] == 0 && dirrections[3] == 1){
		#pragma omp parallel for
		for (int i = 1; i < v1.size() - 1; ++i){
			for (int j = 2; j < v1[i].size() - 1; ++j){
				next1[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad1(v1, v2, i, j, h_x, h_y)) + v1[i][j];
				next2[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad2(v1, v2, i, j, h_x, h_y)) + v2[i][j];
			}
		}
		return;
	}

	if (dirrections[0] == 1 && dirrections[1] == 1
		&& dirrections[2] == 1 && dirrections[3] == 0){
		#pragma omp parallel for
		for (int i = 1; i < v1.size() - 1; ++i){
			for (int j = 1; j < v1[i].size() - 2; ++j){
				next1[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad1(v1, v2, i, j, h_x, h_y)) + v1[i][j];
				next2[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad2(v1, v2, i, j, h_x, h_y)) + v2[i][j];
			}
		}
		return;	
	}

	if (dirrections[0] == 0 && dirrections[1] == 1
		&& dirrections[2] == 0 && dirrections[3] == 1){
		#pragma omp parallel for
		for (int i = 2; i < v1.size() - 1; ++i){
			for (int j = 2; j < v1[i].size() - 1; ++j){
				next1[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad1(v1, v2, i, j, h_x, h_y)) + v1[i][j];
				next2[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad2(v1, v2, i, j, h_x, h_y)) + v2[i][j];
			}
		}
		return;	
	}

	if (dirrections[0] == 0 && dirrections[1] == 1
		&& dirrections[2] == 1 && dirrections[3] == 0){
		#pragma omp parallel for
		for (int i = 2; i < v1.size() - 1; ++i){
			for (int j = 1; j < v1[i].size() - 2; ++j){
				next1[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad1(v1, v2, i, j, h_x, h_y)) + v1[i][j];
				next2[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad2(v1, v2, i, j, h_x, h_y)) + v2[i][j];
			}
		}
		return;
	}

	if (dirrections[0] == 1 && dirrections[1] == 0
		&& dirrections[2] == 0 && dirrections[3] == 1){
		#pragma omp parallel for
		for (int i = 1; i < v1.size() - 2; ++i){
			for (int j = 2; j < v1[i].size() - 1; ++j){
				next1[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad1(v1, v2, i, j, h_x, h_y)) + v1[i][j];
				next2[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad2(v1, v2, i, j, h_x, h_y)) + v2[i][j];
			}
		}
		return;
	}

	if (dirrections[0] == 1 && dirrections[1] == 0
		&& dirrections[2] == 1 && dirrections[3] == 0){
		#pragma omp parallel for
		for (int i = 1; i < v1.size() - 2; ++i){
			for (int j = 1; j < v1[i].size() - 2; ++j){
				next1[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad1(v1, v2, i, j, h_x, h_y)) + v1[i][j];
				next2[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad2(v1, v2, i, j, h_x, h_y)) + v2[i][j];
			}
		}
		return;
	}
	
	if (dirrections[0] == 0 && dirrections[1] == 1
		&& dirrections[2] == 0 && dirrections[3] == 0){
		#pragma omp parallel for
		for (int i = 2; i < v1.size() - 1; ++i){
			for (int j = 2; j < v1[i].size() - 2; ++j){
				next1[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad1(v1, v2, i, j, h_x, h_y)) + v1[i][j];
				next2[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad2(v1, v2, i, j, h_x, h_y)) + v2[i][j];
			}
		}
		return;
	}

	if (dirrections[0] == 1 && dirrections[1] == 0
		&& dirrections[2] == 0 && dirrections[3] == 0){
		#pragma omp parallel for
		for (int i = 1; i < v1.size() - 2; ++i){
			for (int j = 2; j < v1[i].size() - 2; ++j){
				next1[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad1(v1, v2, i, j, h_x, h_y)) + v1[i][j];
				next2[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad2(v1, v2, i, j, h_x, h_y)) + v2[i][j];
			}
		}
		return;
	}

	if (dirrections[0] == 0 && dirrections[1] == 0
		&& dirrections[2] == 0 && dirrections[3] == 1){
		#pragma omp parallel for
		for (int i = 2; i < v1.size() - 2; ++i){
			for (int j = 2; j < v1[i].size() - 1; ++j){
				next1[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad1(v1, v2, i, j, h_x, h_y)) + v1[i][j];
				next2[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad2(v1, v2, i, j, h_x, h_y)) + v2[i][j];
			}
		}
		return;
	}

	if (dirrections[0] == 0 && dirrections[1] == 0
		&& dirrections[2] == 1 && dirrections[3] == 0){
		#pragma omp parallel for
		for (int i = 2; i < v1.size() - 2; ++i){
			for (int j = 1; j < v1[i].size() - 2; ++j){
				next1[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad1(v1, v2, i, j, h_x, h_y)) + v1[i][j];
				next2[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad2(v1, v2, i, j, h_x, h_y)) + v2[i][j];
			}
		}
		return;
	}

	if (dirrections[0] == 0 && dirrections[1] == 0
		&& dirrections[2] == 1 && dirrections[3] == 1){
		#pragma omp parallel for
		for (int i = 2; i < v1.size() - 2; ++i){
			for (int j = 1; j < v1[i].size() - 1; ++j){
				next1[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad1(v1, v2, i, j, h_x, h_y)) + v1[i][j];
				next2[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad2(v1, v2, i, j, h_x, h_y)) + v2[i][j];
			}
		}	
	}

	if (dirrections[0] == 1 && dirrections[1] == 1
		&& dirrections[2] == 0 && dirrections[3] == 0){
		#pragma omp parallel for
		for (int i = 1; i < v1.size() - 1; ++i){
			for (int j = 2; j < v1[i].size() - 2; ++j){
				next1[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad1(v1, v2, i, j, h_x, h_y)) + v1[i][j];
				next2[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad2(v1, v2, i, j, h_x, h_y)) + v2[i][j];
			}
		}
	}

	if (dirrections[0] == 1 && dirrections[1] == 1
		&& dirrections[2] == 1 && dirrections[3] == 1){
		#pragma omp parallel for
		for (int i = 1; i < v1.size() - 1; ++i){
			for (int j = 1; j < v1[i].size() - 1; ++j){
				next1[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad1(v1, v2, i, j, h_x, h_y)) + v1[i][j];
				next2[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad2(v1, v2, i, j, h_x, h_y)) + v2[i][j];
			}
		}
	}	

	if (dirrections[0] == 0 && dirrections[1] == 0 
		&& dirrections[2] == 0 && dirrections[3] == 0){
		#pragma omp parallel for
		for (int i = 2; i < v1.size() - 2; ++i){
			for (int j = 2; j < v1[i].size() - 2; ++j){
				next1[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad1(v1, v2, i, j, h_x, h_y)) + v1[i][j];
				next2[i][j] = tau * (laplasCalc(v1, i, j, h_x, h_y) + grad2(v1, v2, i, j, h_x, h_y)) + v2[i][j];
			}
		}
		return;
	}
}


int main(int argc, char **argv) {
   MPI_Init(&argc, &argv);
   int numprocess;
   vector<vector<double> > localRes1;
   vector<vector<double> > localRes2;
   vector<vector<double> > localU1;
   vector<vector<double> > localU2;

   int xGlobalValue = atoi(argv[1]);
   int yGlobalValue = atoi(argv[2]);
   int procPerXAxis = atoi(argv[3]);
   int procPerYAxis = atoi(argv[4]);
   int numIteration = atoi(argv[5]);
   double h_x = Lx / (xGlobalValue - 1); // Lx/Nx
   double h_y = Ly / (yGlobalValue - 1); // Ly/Ny
   double tau = atof(argv[6]);
   int coords[2];
   

   MPI_Comm_size(MPI_COMM_WORLD, &numprocess);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   int periods[] = { 0, 0 };
   int dims[] = { procPerXAxis, procPerYAxis };
   MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm);

   localRes1 = gridDistribution(xGlobalValue, yGlobalValue, procPerXAxis, procPerYAxis, h_x, h_y, numIteration, tau, 0);
   localRes2 = gridDistribution(xGlobalValue, yGlobalValue, procPerXAxis, procPerYAxis, h_x, h_y, numIteration, tau, 1);
   localU1 = gridDistribution(xGlobalValue, yGlobalValue, procPerXAxis, procPerYAxis, h_x, h_y, numIteration, tau, 2);
   localU2 = gridDistribution(xGlobalValue, yGlobalValue, procPerXAxis, procPerYAxis, h_x, h_y, numIteration, tau, 3);
	
   int tempRank;
   MPI_Comm_rank(comm, &tempRank);
   vector<vector<double> > nextLvl1;
   vector<vector<double> > nextLvl2;
   MPI_Cart_coords(comm, tempRank, 2, coords);
   long double localTime = 0.0;
   localTime = MPI_Wtime();

   for (int it = 0; it < numIteration; ++it){
	ExchangeHalo(tempRank, localRes1, procPerXAxis, procPerYAxis);
	ExchangeHalo(tempRank, localRes2, procPerXAxis, procPerYAxis);
	excuteCalculationNew(localRes1, localRes2, nextLvl1, nextLvl2, coords, procPerXAxis, 
		procPerYAxis, h_x, h_y, tau);
	boundaryConditionV2 (nextLvl1, nextLvl2, coords, procPerXAxis, procPerYAxis, h_x, h_y, tau);
	localRes1 = nextLvl1;
	localRes2 = nextLvl2;
   }

   localTime = MPI_Wtime() - localTime;


   if (rank != 0) {
	long double discrepancy[3];
	discrepancy[0] = getError1(localRes1, localU1, h_x, h_y, numIteration, tau);
	discrepancy[1] = getError1(localRes2, localU2, h_x, h_y, numIteration, tau);
	discrepancy[2] = localTime;
	MPI_Send(discrepancy, 3, MPI_LONG_DOUBLE, 0, rank, comm);
   }

   if (rank == 0) {
	MPI_Status status;
	long double discrepancy[3];
	long double max1 = getError1(localRes1, localU1, h_x, h_y, numIteration, tau);
	long double max2 = getError1(localRes2, localU2, h_x, h_y, numIteration, tau);
	long double maxTime = localTime;
	for (int i = 0; i < numprocess - 1; ++i) {
		MPI_Recv(discrepancy, 3, MPI_LONG_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);
		if (max1 < discrepancy[0]){
			max1 = discrepancy[0];
		}
		if (max2 < discrepancy[1]){
			max2 = discrepancy[1];
		}
		if (maxTime < discrepancy[2]){
			maxTime = discrepancy[2];
		}
	}
	cout << "maxTime = " << maxTime << endl;
	cout << "discrepancy V1 = " << max1 << endl;
	cout << "discrepancy V2 = " << max2 << endl;
   }

   MPI_Barrier(comm);
   MPI_Finalize();
}