#include <iostream>
using namespace std;

class Water {
	public:
		Water();
 		~Water();
		Water(double** coords);
		void getOW(double xyz[3]);
		void getHW1(double xyz[3]);
		void getHW2(double xyz[3]);
		void print();
	private:
		double** watCoords;

};
