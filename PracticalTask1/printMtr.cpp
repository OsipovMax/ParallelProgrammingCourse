
#include <fstream>
#include <iostream>
#include <vector>


int main(int argc, char* argv[])
{
	char type_val;
	int n, m;

	std::fstream fin(argv[1], std::ios::binary | std::ios::in);
	if (fin.is_open()) {

	fin.read((char*)&type_val, sizeof(char));
	fin.read((char*)&n, sizeof(int));
	fin.read((char*)&m, sizeof(int));
	std::cout << type_val << std::endl;
	std::cout << n << std::endl;
	std::cout << m << std::endl;
	if (type_val == 'f') {
		std::vector<std::vector< float >> temp;
		float val;
		temp.resize(n);
		for (int i = 0; i<n; ++i) {
			temp[i].resize(m);
			for (int j = 0; j< m; ++j) {
				fin.read((char*)&val, sizeof(float));
				std::cout << val <<"  ";
			}
			std::cout << std::endl;
		}
	}
	else if (type_val == 'd') {
		double val;
		for (int i = 0; i<n; ++i) {
			for (int j = 0; j<m; ++j) {
				fin.read((char*)&val, sizeof(double));
				std::cout << val << "  ";
			}
			std::cout << std::endl;
		}
	}
	fin.close();
	}
    else  std::cout<<"error"<<std::endl;
	return 0;
}




