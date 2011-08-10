#include "Water.h"

using namespace std;

Water::Water()
:watCoords(0)
{

}

Water::Water(double** coords)
:watCoords(coords)
{

}

Water::~Water(){
	for(int natoms = 0; natoms < 3; natoms++){
		delete [] watCoords[natoms];
	}
	delete [] watCoords;
}

void Water::getOW(double xyz[3]){
	//int oxy_ind = 3*watind;
	xyz[0] = watCoords[0][0];
	xyz[1] = watCoords[0][1];
	xyz[2] = watCoords[0][2];
}
void Water::getHW1(double xyz[3]){
	//int oxy_ind = 3*watind + 1;
	xyz[0] = watCoords[1][0];
	xyz[1] = watCoords[1][1];
	xyz[2] = watCoords[1][2];
}

void Water::getHW2(double xyz[3]){
	//int oxy_ind = 3*watind + 2;
	xyz[0] = watCoords[2][0];
	xyz[1] = watCoords[2][1];
	xyz[2] = watCoords[2][2];
}

void Water::print(){
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			cout<<watCoords[i][j];
		}
		cout<<endl;
	}
}

