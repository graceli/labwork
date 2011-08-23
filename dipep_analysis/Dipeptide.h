#ifndef _DIPEPTIDE_H
#define _DIPEPTIDE_H

#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#define PI 3.14159265
//#define DEBUG
using namespace std;

class Dipeptide {

	public:
		Dipeptide(double** dipepCoords);
		//~Dipeptide();
		void debug();
		double getNHdist(int groupNum);
		double getDonorAngle(int groupNum, double accepter[3]);	//0<=groupNum<=1
		bool isHbonded(int groupNum, double accepter[3]);
		double getDonorAcceptorDist(int groupNum, double acceptor[3]);
		void getCO(int groupNum, double Ocoord[3]);
		double getAcceptorHDist(int groupNum,double acceptor[3]);
	private:
		double norm(double v[3]);
		double dot(double v1[3], double v2[3]);
		double** dipep;
		vector<string> atomNameIndex;		//map of atom name to atom index
		string title;		//grofile title
		string numAtoms;	//number of atoms in file
};
#endif

Dipeptide::Dipeptide(double** dipepCoords)
:dipep(dipepCoords) {

#ifdef DEBUG
	cerr<<"Initialized Dipeptide object with coordinates:"<<endl;

	for(int i=0; i<10; i++){
		for(int j=0; j<3; j++){
			cerr<<dipep[i][j]<<"\t";
		}
		cerr<<endl;
	}
#endif
}

//returns the coordinates of the O atom of CO groups 
void Dipeptide::getCO(int groupNum, double Ocoord[3]){
	if(groupNum == 0){
		for(int i=0; i<3; i++){
			Ocoord[i]=dipep[1][i];
		}
	}else{
		for(int i=0; i<3; i++){
			Ocoord[i]=dipep[7][i];
		}
	}
}

//checks if the hydrogen bond criteria (used by stephanie) is satisfied
//i.e. angle DHA >=120 degrees, dist(D,A) <= 0.35 nm, dist(H,A)<=0.25nm
bool Dipeptide::isHbonded(int groupNum, double acceptor[3]){
	double DAdist = getDonorAcceptorDist(groupNum,acceptor);
	double DAangle = getDonorAngle(groupNum, acceptor);
	double HAdist = getAcceptorHDist(groupNum,acceptor);

#ifdef DEBUG
	cerr<<"isHbonded(...)"<<endl;
	cerr<<"group "<<groupNum<<": DAdist = "<<DAdist<<endl;
	cerr<<"group "<<groupNum<<": DAangle = "<<DAangle<<endl;
#endif

	if(DAdist <= 0.351 && DAangle >= 120 && HAdist <=0.25){
		return true;
	}
	return false;
}

//returns the distance between Acceptor atom and Hydrogen atom of the donor; dist(H,A)
double Dipeptide::getAcceptorHDist(int groupNum, double acceptor[3]){
	double dist;
	if(groupNum == 0){
		dist=sqrt( pow(acceptor[0]-dipep[3][0] , 2) + pow(acceptor[1]-dipep[3][1] , 2) + pow(acceptor[2]-dipep[3][2] , 2) );
#ifdef DEBUG
		cerr<<"group "<<groupNum<<": Acceptor - DonorH dist="<<dist<<endl;
#endif
		return dist;
	}else{
		dist=sqrt( pow(acceptor[0] - dipep[9][0] , 2) + pow(acceptor[1] - dipep[9][1] , 2) + pow(acceptor[2] - dipep[9][2] , 2) );
#ifdef DEBUG
		cerr<<"group "<<groupNum<<": Acceptor - DonorH dist="<<dist<<endl;
#endif
		return dist;
	}
	
}
double Dipeptide::getNHdist(int groupNum){
	double dist=0;
	if(groupNum == 0 ){
		for(int i=0; i<3; i++){
			dist+=pow(dipep[3][i] - dipep[2][i],2);
		}
	}else{
		for(int i=0; i<3; i++){
			dist+=pow(dipep[9][i] - dipep[8][i],2);
		}
	}
	return sqrt(dist);
}

//returns the distance between the Acceptor atom and the Donor atom
double Dipeptide::getDonorAcceptorDist(int groupNum, double acceptor[3]) {
	double dist;
	if(groupNum == 0){
		dist=sqrt( pow(acceptor[0]-dipep[2][0] , 2) + pow(acceptor[1]-dipep[2][1] , 2) + pow(acceptor[2]-dipep[2][2] , 2) );
#ifdef DEBUG
		cerr<<"group "<<groupNum<<": Acceptor-Donor dist="<<dist<<endl;
#endif
		return dist;
	}else{
		dist=sqrt( pow(acceptor[0] - dipep[8][0] , 2) + pow(acceptor[1] - dipep[8][1] , 2) + pow(acceptor[2] - dipep[8][2] , 2) );
#ifdef DEBUG
		cerr<<"group "<<groupNum<<": Acceptor-Donor dist="<<dist<<endl;
#endif
		return dist;
	}
}

//returns the angle H-D-A; previously used with the Gromacs hydrogen bonding criteria
/*double Dipeptide::getDonorAngle(int donorGroup, double acceptor[3]){
	//returns the H-D-A angle in degrees
	double vHN[3];
	double vAN[3];
	double alpha; //the angle between vHN and vAN
	for(int i=0; i<3; i++){
		if(donorGroup == 0) {
			//vector HN = point H - point N
			vHN[i]=dipep[3][i] - dipep[2][i];
			//cout<<vHN[i]<<endl;
			vAN[i]=acceptor[i] - dipep[2][i];
		}else{
			vHN[i]=dipep[9][i] - dipep[8][i];
			vAN[i]=acceptor[i] - dipep[8][i];
		}
	}

#ifdef DEBUG 
	cerr<<"group "<<donorGroup<<":vHN="<<vHN[0]<<" "<<vHN[1]<<" "<<vHN[2]<<endl;
	cerr<<"group "<<donorGroup<<":vAN="<<vAN[0]<<" "<<vAN[1]<<" "<<vAN[2]<<endl;
#endif

	alpha=acos(dot(vHN,vAN)/(norm(vHN)*norm(vAN)));

#ifdef DEBUG
	cerr<<"group"<<donorGroup<<": angle(vHN,vAN) = "<<alpha*180.0/PI<<endl;
#endif
	return alpha*180.0/PI;
}*/

//returns the angle D-H-A in degrees
double Dipeptide::getDonorAngle(int donorGroup, double acceptor[3]){
	double vNH[3];	//vector N - H
	double vAH[3];	//vector A - H
	double alpha; //the angle between vHN and vAN
	for(int i=0; i<3; i++){
		if(donorGroup == 0) {
			//vector NH = point N - point H
			vNH[i]=dipep[2][i] - dipep[3][i];
			//vector AH = point A - point H
			vAH[i]=acceptor[i] - dipep[3][i];
		}else{
			vNH[i]=dipep[8][i] - dipep[9][i];
			vAH[i]=acceptor[i] - dipep[9][i];
		}
	}
	alpha=acos(dot(vNH,vAH)/(norm(vNH)*norm(vAH)));

	return alpha*180.0/PI;
}
//returns norm of a vector
double Dipeptide::norm(double v[3]){
	return sqrt(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));
}

//returns the dot product of two vectors
double Dipeptide::dot(double v1[3], double v2[3]){
	return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}

void Dipeptide::debug(){
	cerr<<"atomNameIndex map"<<endl;

	for(int k=0; k<atomNameIndex.size(); k++){
		cerr<<k<<"\t"<<atomNameIndex[k]<<endl;
	}
	cerr<<endl;
	for(int i=0; i<10; i++){
		for(int j=0; j<4; j++){
			cerr<<dipep[i][j]<<"\t";
		}
		cerr<<endl;
	}
}

