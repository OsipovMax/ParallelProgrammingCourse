#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <algorithm>

using namespace std;

struct coord {
	int x;
	int y;
	int z;
	int t;
};
int main()
{
	vector<coord>v;
	coord c;

	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			for (int k = 0; k < 8; k++) {
				c.x = i;
				c.y = j;
				c.z = k;
				c.t = 0;
				v.push_back(c);
			}
		}
	}
	srand(unsigned(time(0)));
	random_shuffle(v.begin(), v.end());
	
	ofstream fout("myMap.map", ios::out);
	for (int i = 0; i < v.size(); i++) {
		fout << v[i].x << " " << v[i].y << " " << v[i].z << " " << v[i].t << endl;
	}
	fout.close();

	return 0;
}