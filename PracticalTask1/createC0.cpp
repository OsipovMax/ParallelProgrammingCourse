#include <vector>
#include <fstream>
#include <iostream>
#include <boost\numeric\ublas\matrix.hpp>

int main(int argc,char* argv[]) {
	boost::numeric::ublas::matrix<double> Ad;
	boost::numeric::ublas::matrix<float> Af;
	boost::numeric::ublas::matrix<double> Bd;
	boost::numeric::ublas::matrix<float> Bf;

	int n, m;
	char type_val;
	std::fstream fin(argv[1], std::ios::binary | std::ios::in);
	if (fin.is_open()) {

		fin.read((char*)&type_val, sizeof(char));
		fin.read((char*)&n, sizeof(int));
		fin.read((char*)&m, sizeof(int));
		std::cout << type_val << std::endl;
		std::cout << n << std::endl;
		std::cout << m << std::endl;

		if (type_val == 'f') {
			Af.resize(n, m);
			std::vector<std::vector< float >> temp;
			float val;;
			for (int i = 0; i<Af.size1(); ++i) {
				temp[i].resize(m);
				for (int j = 0; j<Af.size2(); ++j) {
					fin.read((char*)&val, sizeof(float));
					temp[i][j] = val;
					Af(i, j) = val;
					std::cout << Af(i, j) << " ";
				}
				std::cout << std::endl;
			}
		}

		else if (type_val == 'd') {
			std::vector < std::vector < double >> temp;
			temp.resize(n);
			double val;
			Ad.resize(n, m);
			for (int i = 0; i<Ad.size1(); ++i) {
				temp[i].resize(m);
				for (int j = 0; j<Ad.size2(); ++j) {
					fin.read((char*)&val, sizeof(double));
					temp[i][j] = val;
					Ad(i, j) = val;
					std::cout << Ad(i, j) << " ";
					//std::cout << temp[i][j] << "  ";
				}
				std::cout << std::endl;
			}
		}
		fin.close();
	}
	std::fstream fin1(argv[2], std::ios::binary | std::ios::in);
	if (fin1.is_open()) {

		fin1.read((char*)&type_val, sizeof(char));
		fin1.read((char*)&n, sizeof(int));
		fin1.read((char*)&m, sizeof(int));
		std::cout << type_val << std::endl;
		std::cout << n << std::endl;
		std::cout << m << std::endl;
		if (type_val == 'f') {
			std::vector<std::vector< float >> temp;
			float val;
			Bf.resize(n, m);
			temp.resize(n);
			for (int i = 0; i<Bf.size1(); ++i) {
				temp[i].resize(m);
				for (int j = 0; j<Bf.size2(); ++j) {
					fin.read((char*)&val, sizeof(float));
					temp[i][j] = val;
					Bf(i, j) = val;
					std::cout << Bf(i, j) << " ";
				}
				std::cout << std::endl;
			}
		}
		else if (type_val == 'd') {
			std::vector<std::vector< float >> temp;
			temp.resize(n);
			double val;
			Bd.resize(n, m);
			for (int i = 0; i<Bd.size1(); ++i) {
				temp[i].resize(m);
				for (int j = 0; j<Bd.size2(); ++j) {
					fin1.read((char*)&val, sizeof(double));
					temp[i][j] = val;
					Bd(i, j) = val;
					std::cout << Bd(i, j) << " ";
				}
				std::cout << std::endl;
			}
		}
		fin.close();
	}

	if (type_val = 'f') {
	boost::numeric::ublas::matrix<float> Cf;
	Cf = prod(Af,Bf);
	std::ofstream foutC0 ("C0", std::ios::binary | std::ios::out);
	if (foutC0.is_open()) {
		foutC0.write((char*)&type_val, sizeof(char));
		foutC0.write((char*)&n, sizeof(int));
		foutC0.write((char*)&m, sizeof(int));
		for (int i = 0; i< Cf.size1(); ++i) {
			for (int j = 0; Cf.size2(); ++j) {
				foutC0.write((char*)&Cf(i,j), sizeof(float)); //*n*m
			}
		}
	}
	foutC0.close();
	}

	if (type_val = 'd') {
		boost::numeric::ublas::matrix<double> Cd;
		Cd = prod(Ad, Bd);
		std::ofstream foutC0("C0", std::ios::binary | std::ios::out);
		if (foutC0.is_open()) {
			foutC0.write((char*)&type_val, sizeof(char));
			foutC0.write((char*)&n, sizeof(int));
			foutC0.write((char*)&m, sizeof(int));
			for (int i = 0; i< Cd.size1(); ++i) {
				for (int j = 0; j<Cd.size2(); ++j) {
					foutC0.write((char*)&Cd(i, j), sizeof(double)); //*n*m
				}
			}
		}
		foutC0.close();
	}
	std::cout << std::endl;
	boost::numeric::ublas::matrix<double> Cd;
	Cd = prod(Ad, Bd);
	for (int i = 0; i<Cd.size1(); ++i) {
		for (int j = 0; j<Bd.size2(); ++j) {
			std::cout << Cd(i, j) << " ";
		}
		std::cout << std::endl;
	}

	return 0;
}
