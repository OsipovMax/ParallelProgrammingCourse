#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>
#include <limits>



int main(int argc, char* argv[])
{

    char type_valC0,type_valC;
    int nC,mC,nC0,mC0;



    std::fstream fC ( argv[1], std::ios::binary | std::ios::in);
    if (!fC.is_open()){
        std::cout << "file not open" <<std::endl;
        exit(1);
    }
    fC.read((char*)&type_valC, sizeof(char));
    fC.read((char*)&nC, sizeof(int));
    fC.read((char*)&mC, sizeof(int));


    std::fstream fC0 ( argv[2], std::ios::binary | std::ios::in);
    if (!fC0.is_open()){
        std::cout << "file not open" <<std::endl;
        exit(1);
    }
    fC0.read((char*)&type_valC0, sizeof(char));
    fC0.read((char*)&nC0, sizeof(int));
    fC0.read((char*)&mC0, sizeof(int));

    if (type_valC != type_valC0){
        std::cout << " matrix type does not match " << std::endl;
        exit(1);
    }

    if (nC != nC0 || mC != mC0){
        std::cout << "matrix size does not match" << std::endl;
        exit(1);
    }

    if (type_valC == 'f'){
        float tempC,tempC0;
        for (int i = 0; i<nC ; ++i){
            for (int j= 0; j<mC; ++j){
                fC.read((char*)&tempC, sizeof(float));
                fC0.read((char*)&tempC0, sizeof(float));
                if (abs(tempC-tempC0) > std::numeric_limits<float>::epsilon()) {
                    std::cout << "matrix are not equal" << std::endl;
                    exit(1);
                }
            }
        }
    }
    if (type_valC == 'd'){
        double tempC,tempC0;
        for (int i = 0; i<nC ; ++i){
            for (int j= 0; j<mC; ++j){
                fC.read((char*)&tempC, sizeof(double));
                fC0.read((char*)&tempC0, sizeof(double));
                if (abs(tempC-tempC0) > std::numeric_limits<double>::epsilon()) {
                    std::cout << "matrix are not equal" << std::endl;
                    exit(1);
                }
            }
        }
    }
    std::cout << " matrix are equal" << std:: endl;

    return 0;
}
