#include <vector>
#include <fstream>
#include <iostream>
#include <time.h>

int main(int argc, char* argv[])
{


	int n = atoi(argv[2]);
	int m = atoi(argv[3]);

	if (*argv[1] == 'f') {
		std::vector<std::vector< float >> matr;
		matr.resize(n);
		srand(time(0));
		for (int i = 0; i< n; ++i) {
			matr[i].resize(m);
			for (int j = 0; j<m; ++j) {
				matr[i][j] = (float)(rand()) / RAND_MAX;
			}
		}
		std::ofstream fout(argv[4], std::ios::binary | std::ios::out);
		if (fout.is_open()) {
			fout.write(argv[1], sizeof(char));
			fout.write((char*)&n, sizeof(int));
			fout.write((char*)&m, sizeof(int));
			for (int i = 0; i< n; ++i) {
				for (int j = 0; j<m; ++j) {
					fout.write((char*)&matr[i][j], sizeof(float)); //*n*m
				}
			}
		}
		fout.close();
	}

	if (*argv[1] == 'd') {
		std::vector<std::vector< double >> matr;
		matr.resize(n);
		srand(time(0));
		for (int i = 0; i<n; ++i) {
			matr[i].resize(m);
			for (int j = 0; j<m; ++j) {
				matr[i][j] = (double)(rand()) / RAND_MAX;
			}
		}
		std::ofstream fout(argv[4], std::ios::binary | std::ios::out);
		if (fout.is_open()) {
			fout.write(argv[1], sizeof(char));
			fout.write((char*)&n, sizeof(int));
			fout.write((char*)&m, sizeof(int));
			for (int i = 0; i<n; ++i) {
				for (int j = 0; j<m; ++j) {
					fout.write((char*)&matr[i][j], sizeof(double));
				}
			}
		}
		fout.close();
	}
	return 0;
}
