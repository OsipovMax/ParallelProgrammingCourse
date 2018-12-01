#include <mpi.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "Matrix.h"

using namespace std;

int coord[3];

MPI_Comm ij_comm;
int numproces_axis;
int N_block_size, part_size;
int rz;

template <typename T>
Matrix<T> mult(Matrix<T> &m1, Matrix<T> &m2) {
  vector<vector<T> > temp_vv(m1.getRows(), vector<T>(m2.getColumns(), 0));
  Matrix<T> temp(m1.getRows(), m2.getColumns(), temp_vv);
  int blocksize = 36;
  for (int i = 0; i < m1.getRows(); i += blocksize) {
    for (int k = 0; k < m2.getColumns(); k += blocksize) {
      for (int j = 0; j < m1.getColumns(); j += blocksize) {
        for (int i1 = i; i1 < i + blocksize && i1 < m1.getRows(); i1++) {
          for (int k1 = k; k1 < k + blocksize && k1 < m2.getColumns(); k1++) {
            for (int j1 = j; j1 < j + blocksize && j1 < m1.getColumns(); j1++) {
              temp.mtr[i1][j1] += m1.mtr[i1][k1] * m2.mtr[k1][j1];
            }
          }
        }
      }
    }
  }
  return temp;
}

template <typename T>
Matrix<T> toMatr(int r_size, int c_size, double *buf) {
  int inner_counter = 0;
  Matrix<T> tempObj;
  tempObj.setRows(r_size);
  tempObj.setColumns(c_size);
  tempObj.mtr.resize(r_size);
  for (int i = 0; i < r_size; i++) {
    tempObj.mtr[i].resize(c_size);
    for (int j = 0; j < c_size; j++) {
      tempObj.mtr[i][j] = buf[inner_counter];
      inner_counter++;
    }
  }

  return tempObj;
}

template <typename T>
double *toBuf(Matrix<T> &m) {
  int inner_counter = 0;
  int numproces;
  double *buf = new double[m.getColumns() * m.getRows()];
  for (int i = 0; i < m.getRows(); ++i) {
    for (int j = 0; j < m.getColumns(); ++j) {
      buf[inner_counter] = m.mtr[i][j];
      inner_counter++;
    }
  }
  return buf;
}

template <typename T>
Matrix<T> binReader(const std::string &strr) {
  int modeVal, n, m;
  double val;
  char type_val;
  Matrix<T> tempObj;

  const char *cstr = strr.c_str();
  std::fstream f1(cstr, std::ios::binary | std::ios::in);
  f1.read((char *)&type_val, sizeof(char));
  f1.read((char *)&n, sizeof(int));
  f1.read((char *)&m, sizeof(int));
  tempObj.setRows(n);
  tempObj.setColumns(m);

  tempObj.mtr.resize(n);
  for (int i = 0; i < n; ++i) {
    tempObj.mtr[i].resize(m);
    for (int j = 0; j < m; j++) {
      f1.read((char *)&val, sizeof(double));
      tempObj.mtr[i][j] = val;
    }
  }
  f1.close();
  return tempObj;
}

/*template <typename T>
double *getBlock(char *file_name) {

  int matrixSize;
  char tp;
  double *buf, mp;

  MPI_Status status;
  MPI_File fl;

  MPI_File_open(ij_comm, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fl);
  MPI_File_read_at(fl, 0, &tp, 1, MPI_CHAR, &status);
  MPI_File_read_at(fl, 1, &matrixSize, 1, MPI_INT, &status);
  MPI_File_read_at(fl, 5, &matrixSize, 1, MPI_INT, &status);

  rz = matrixSize;

  N_block_size = matrixSize / numproces_axis;
  part_size = N_block_size * N_block_size;
  buf = new double[part_size];

  int offset = sizeof(double) *
                  (coord[1] * N_block_size + N_block_size * matrixSize * coord[0]) +
              2 * sizeof(int) + sizeof(char);

  MPI_File_seek(fl, offset, MPI_SEEK_SET);
  for (int i = 0; i < N_block_size; i++) {
    for (int j = 0; j < N_block_size; j++) {

      MPI_File_read(fl, &mp, 1, MPI_DOUBLE, &status);
      buf[i * N_block_size + j] = mp;
    }
    offset = sizeof(double) * (matrixSize - N_block_size);
    MPI_File_seek(fl, offset, MPI_SEEK_CUR);
  }
   return buf;
}*/


template <typename T>
double *getBlock (char *file_name){
	
   ifstream fin(file_name, ios::binary | ios::in);
   int matrixSize;
   char tp;
   double *buf, mp;
   fin.read((char*)&tp, sizeof(char));
   fin.read((char*)&matrixSize, sizeof(int));
   fin.read((char*)&matrixSize, sizeof(int));

   rz = matrixSize;

   N_block_size = matrixSize/numproces_axis;
   part_size = N_block_size * N_block_size;
   buf = new double [part_size];

   int offset = sizeof(double) * (coord[1] * N_block_size + N_block_size * matrixSize * coord[0]) + 2*sizeof(int) + sizeof(char);

   fin.seekg(offset, fin.beg);
   for ( int i = 0; i < N_block_size; i++ ){
      for ( int j = 0; j < N_block_size; j++ ){
         fin.read((char*)&mp, sizeof(double));
         buf[i*N_block_size + j] = mp;
      }
      offset = sizeof(double)*(matrixSize - N_block_size);
      fin.seekg(offset, fin.cur);
   }
   return buf;
}

void writeBlock(double *buf, char *file_name) {
  if (coord[2] != 0) {
    return;
  }
  char cha = 'd';
  MPI_Status status;
  int matrixSize = N_block_size * numproces_axis;

  MPI_File fh;
  MPI_File_open(ij_comm, file_name, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL,
                &fh);
  if (coord[0] == coord[1] == coord[2] == 0) {
    MPI_File_write_at(fh, 0, &cha, 1, MPI_CHAR, &status);
    MPI_File_write_at(fh, 1, &matrixSize, 1, MPI_INT, &status);
    MPI_File_write_at(fh, 5, &matrixSize, 1, MPI_INT, &status);
  }


  int offset = sizeof(double) *
                   (coord[0] * matrixSize * N_block_size + coord[1] * N_block_size) +
               2 * sizeof(int) + sizeof(char);
  for (int i = 0; i < part_size; i++) {
    int blockOffset =
        sizeof(double) * ((i / N_block_size) * matrixSize + i % N_block_size);
    MPI_File_write_at(fh, offset + blockOffset, &buf[i], 1, MPI_DOUBLE,
                      &status);
  }

  return;
}

 
int main(int argc, char **argv) {
  int numproces, grid_rank, kk;
  double time_total = 0.0, time_max = 0.0;
  double time_begin = 0.0, time_end = 0.0;
  double *proc_bufA, *proc_bufB, *proc_bufC;

  MPI_Comm grid_comm, i_comm, j_comm, k_comm;
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numproces);

  numproces_axis = atoi(argv[4]);
  if (numproces == 1) {
    Matrix<double> a = binReader<double>(argv[1]);
    Matrix<double> b = binReader<double>(argv[2]);
    time_begin = MPI_Wtime();
    Matrix<double> c = mult(a, b);
    time_end = MPI_Wtime() - time_begin;

    time_max = time_end;
    time_total = time_end;
    ofstream fout1p("rproc1", ios::out);
    fout1p << time_max << " " << time_total << " " << numproces << endl;
    fout1p.close();

    std::ofstream ft(argv[3], std::ios::binary | std::ios::out);
    if (ft.is_open()) {
      double temp_val;
      char type_val = 'd';
      int newRow = c.getRows();
      int newColmn = c.getColumns();
      cout << "newRow =" << newRow << endl;
      cout << "newCol =" << newColmn << endl;
      ft.write((char *)&type_val, sizeof(char));
      ft.write((char *)&newRow, sizeof(int));
      ft.write((char *)&newColmn, sizeof(int));
      for (int i = 0; i < newRow; ++i) {
        for (int j = 0; j < newColmn; ++j) {
          temp_val = c.mtr[i][j];
          ft.write((char *)&temp_val, sizeof(double));
        }
      }
    }
    ft.close();
    MPI_Finalize();
    return 0;
  }

  int dims[3] = {numproces_axis, numproces_axis, numproces_axis};
  int period[3] = {0, 0, 0};
  int reorder = 0;

  MPI_Cart_create(MPI_COMM_WORLD, 3, dims, period, reorder, &grid_comm);

  int remain_dims_ij[3] = {1, 1, 0};

  int remain_dims_i[3] = {1, 0, 0};
  int remain_dims_j[3] = {0, 1, 0};
  int remain_dims_k[3] = {0, 0, 1};

  MPI_Cart_create(MPI_COMM_WORLD, 3, dims, period, reorder, &grid_comm);
  MPI_Cart_sub(grid_comm, remain_dims_ij, &ij_comm);
  MPI_Cart_sub(grid_comm, remain_dims_i, &i_comm);
  MPI_Cart_sub(grid_comm, remain_dims_j, &j_comm);
  MPI_Cart_sub(grid_comm, remain_dims_k, &k_comm);

  MPI_Comm_rank(grid_comm, &grid_rank);


  MPI_Cart_coords(grid_comm, grid_rank, 3, coord);

  if (coord[2] == 0) {
    proc_bufA = getBlock<double>(argv[1]);
  }
  if (coord[2] == 0) {
       proc_bufB = getBlock<double>(argv[2]);
  }

  MPI_Barrier(grid_comm);
  MPI_Bcast(&part_size, 1, MPI_DOUBLE, 0, k_comm);

  MPI_Barrier(grid_comm);
  MPI_Bcast(&rz, 1, MPI_INT, 0, k_comm);

  if (coord[2] != 0) {
    proc_bufA = new double[part_size];
    proc_bufB = new double[part_size];
  }

  int sendCoord[3], recvCoord[3], bcoord[1];

  int i = coord[0], j = coord[1], k = coord[2];
  if (k == 0 && j != 0) {
    int rrank;
    recvCoord[0] = i;
    recvCoord[1] = j;
    recvCoord[2] = j;
    MPI_Cart_rank(grid_comm, recvCoord, &rrank);
    MPI_Send(proc_bufA, part_size, MPI_DOUBLE, rrank, grid_rank, grid_comm);
  }

  if (j == k && j != 0) {
    MPI_Status status;
    int srank;
    sendCoord[0] = i;
    sendCoord[1] = j;
    sendCoord[2] = 0;
    MPI_Cart_rank(grid_comm, sendCoord, &srank);
    MPI_Recv(proc_bufA, part_size, MPI_DOUBLE, srank, MPI_ANY_TAG, grid_comm,
             &status);
  }

  int rank;

  bcoord[0] = coord[2];
  MPI_Cart_rank(j_comm, bcoord, &rank);
  MPI_Bcast(proc_bufA, part_size, MPI_DOUBLE, rank, j_comm);

  if (k == 0 && i != 0) {
    int rrank;
    recvCoord[0] = i;
    recvCoord[1] = j;
    recvCoord[2] = i;
    MPI_Cart_rank(grid_comm, recvCoord, &rrank);
    MPI_Send(proc_bufB, part_size, MPI_DOUBLE, rrank, grid_rank, grid_comm);
  }
  if (i == k && i != 0) {
    MPI_Status status;
    int srank;
    sendCoord[0] = i;
    sendCoord[1] = j;
    sendCoord[2] = 0;
    MPI_Cart_rank(grid_comm, sendCoord, &srank);
    MPI_Recv(proc_bufB, part_size, MPI_DOUBLE, srank, MPI_ANY_TAG, grid_comm,
             &status);
  }

  bcoord[0] = coord[2];
  MPI_Cart_rank(i_comm, bcoord, &rank);
  MPI_Bcast(proc_bufB, part_size, MPI_DOUBLE, rank, i_comm);

  Matrix<double> a =
      toMatr<double>(rz / numproces_axis, rz / numproces_axis, proc_bufA);
  Matrix<double> b =
      toMatr<double>(rz / numproces_axis, rz / numproces_axis, proc_bufB);
 
  time_begin = MPI_Wtime();
  Matrix<double> c = mult(a, b);
  time_end = MPI_Wtime() - time_begin;

  if (grid_rank != 0) {
    MPI_Send(&time_end, 1, MPI_DOUBLE, 0, rank, grid_comm);
  }


  if (grid_rank == 0) {
    time_total += time_end;
    if (time_end > time_max) time_max = time_end;
    for (int i = 0; i < numproces - 1; ++i) {
      MPI_Recv(&time_end, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, grid_comm,
               &status);
      if (time_end > time_max) time_max = time_end;
      time_total += time_end;
    }
  }

  proc_bufC = toBuf<double>(c);

  double *sumBuffer;

  sumBuffer = new double[part_size];
  for (int i = 0; i < part_size; i++) sumBuffer[i] = 0;
  MPI_Reduce(proc_bufC, sumBuffer, part_size, MPI_DOUBLE, MPI_SUM, 0, k_comm);

  double t_b = 0.0, t_e = 0.0;
  
  t_b = MPI_Wtime();
  writeBlock(sumBuffer, argv[3]);
  t_e = MPI_Wtime() - t_b;

  double t_t;
  if (grid_rank != 0) {
    MPI_Send(&t_e, 1, MPI_DOUBLE, 0, rank, grid_comm);
  }

  if (grid_rank == 0) {
    t_t += t_e;
    for (int i = 0; i < numproces - 1; ++i) {
      MPI_Recv(&t_e, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, grid_comm,
               &status);
      t_t += t_e;
    }
  }


  if (grid_rank == 0) {
    int rcol = 1;
    char tp = 'd';
    if (numproces == 1) {
      ofstream fout1p("rproc1", ios::out);
      fout1p << time_max << " " << time_total << " " << numproces << endl;
      fout1p.close();
    } else {
      float tm;
      double speedup, efficiency;
      ifstream fin("rproc1", ios::in);
      fin >> tm;
      fin.close();
      speedup = (double)tm / time_max;
      efficiency = speedup / numproces;
      ofstream fout("result", ios::app);
      fout << numproces << " " << time_total << " " << time_max << " "
           << speedup << " " << efficiency << endl;
      fout.close();
      ofstream file_out("write_block", ios::app);
      file_out << numproces << " " << t_t << endl;
      file_out.close();
    }
  }

  delete[] proc_bufA;
  delete[] proc_bufB;
  delete[] sumBuffer;

  MPI_Finalize();
  return 0;
}