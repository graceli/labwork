#ifndef _INOSITOL_H_
#define _INOSITOL_H_

#include <cmath>
#define PI 3.14159265

// #define DEBUG
enum bidstates {NONE,
		CO0CO2_12,CO0NH3_12,NH1CO2_12,NH1NH3_12,
		CO0CO2_13,CO0NH3_13,NH1CO2_13,NH1NH3_13,
		CO0CO2_14,CO0NH3_14,NH1CO2_14,NH1NH3_14,
		MONO_CO0,MONO_NH1,MONO_CO2,MONO_NH3,MONO_TWO}; 	//combinations of peptide groups

class Inos{
	public:
		Inos(double** inos);
		double getDonorAngle(int groupNum, double accepter[3]);
		bool isHbonded(int groupNum, double accepter[3]);
		double getDonorAcceptorDist(int groupNum, double acceptor[3]);
		double norm(double v[3]);
		double dot(double v1[3], double v2[3]);
		void getOH(int groupNum, double Ocoords[3]);
		void setHbondedInos(int inosID);
		double getAcceptorHDist(int groupNum,double acceptor[3]);
		bool isInosHbonded(int inosID);
		bool isHbondedPep();
		void setHbondedPep();
		void setBindMode(bidstates mode);
		bidstates getBindMode();
	private:
		double** inosCoords;
		int hbondedInos[4];	//vector of inositol ids that *this is hbonded to
		bool hbonded2Pep;
		bidstates binding_mode;
};
#endif

Inos::Inos(double** inos)
:inosCoords(inos),hbonded2Pep(false), binding_mode(NONE)
{
#ifdef DEBUG
	for(int i=0; i<12; i++){
		for(int j=0; j<3; j++){
			cerr<<inosCoords[i][j]<<"\t";
		}
		cerr<<endl;
	}
#endif
}

void Inos::setBindMode(bidstates mode){
	binding_mode = mode;
}
bidstates Inos::getBindMode(){
	return binding_mode;
}
void Inos::setHbondedPep(){
	hbonded2Pep=true;
}

bool Inos::isHbondedPep(){
	return hbonded2Pep;
}

//records whether this is hbonded
void Inos::setHbondedInos(int inosID){
	hbondedInos[inosID]=1;
}

bool Inos::isInosHbonded(int inosID){
	return (hbondedInos[inosID]==1);
}

void Inos::getOH(int groupNum, double Ocoords[3]){
	int groupIndex=2*groupNum;
	for(int i=0; i<3; i++){
		Ocoords[i]=inosCoords[groupIndex][i];
	}
}

bool Inos::isHbonded(int groupNum, double acceptor[3]){
	if( getDonorAcceptorDist(groupNum,acceptor) <= 0.351 && getDonorAngle(groupNum, acceptor) >=120 && getAcceptorHDist(groupNum, acceptor) <= 0.25 ) {
		return true;
	}
	return false;
}

//pre: 0<=groupNum <=5
//returns the distance of Acceptor atom to Hydrogen atom of OH groupNum
double Inos::getAcceptorHDist(int groupNum,double acceptor[3]){
	int groupIndex=2*groupNum+1;	//index to the H atom of OH group number groupNum
	double dist = sqrt( pow(acceptor[0] - inosCoords[groupIndex][0] , 2) + pow(acceptor[1] - inosCoords[groupIndex][1] , 2) + pow(acceptor[2] - inosCoords[groupIndex][2] , 2) );
	return dist;
}

//pre: 0<=groupNum <=5
//returns the distance from Donor to acceptor
double Inos::getDonorAcceptorDist(int groupNum, double acceptor[3]) {
	int groupIndex=2*groupNum;
	double dist = sqrt( pow(acceptor[0] - inosCoords[groupIndex][0] , 2) + pow(acceptor[1] - inosCoords[groupIndex][1] , 2) + pow(acceptor[2] - inosCoords[groupIndex][2] , 2) );

#ifdef DEBUG
	cerr<<"Acceptor-Donor dist="<<dist<<endl;
#endif

	return dist;
}
/*
//returns the H-Donor-Acceptor angle in degrees, where the donor is an OH of inositol
double Inos::getDonorAngle(int donorGroup, double acceptor[3]){
	//returns the H-D-A angle in degrees
	double vHO[3];
	double vAO[3];
	double alpha; //the angle between vHN and vAN
	int groupIndex=2*donorGroup;
	for(int i=0; i<3; i++){
		//vector HN = point H - point N
		vHO[i]=inosCoords[groupIndex+1][i] - inosCoords[groupIndex][i];
		vAO[i]=acceptor[i] - inosCoords[groupIndex][i];
	} 
#ifdef DEBUG
	cerr<<"vHO="<<vHO[0]<<" "<<vHO[1]<<" "<<vHO[2]<<endl;
	cerr<<"vAO="<<vAO[0]<<" "<<vAO[1]<<" "<<vAO[2]<<endl;
#endif
	alpha=acos(dot(vHO,vAO)/(norm(vHO)*norm(vAO)));

#ifdef DEBUG
	cerr<<"angle(vHO,vAO) = "<<alpha*180.0/PI<<endl;
#endif

	return alpha*180.0/PI;
}
*/

//returns the H-Donor-Acceptor angle in degrees, where the donor is an OH of inositol
double Inos::getDonorAngle(int donorGroup, double acceptor[3]){
	//returns the H-D-A angle in degrees
	double vOH[3]; //vector O - H
	double vAH[3]; //vector A - H
	double alpha; //the angle between vHN and vAN
	int groupIndex=2*donorGroup;
	for(int i=0; i<3; i++){
		//vector HO = point O - point H
		vOH[i]=inosCoords[groupIndex][i] - inosCoords[groupIndex+1][i];
		//vector HA = point A - point H
		vAH[i]=acceptor[i] - inosCoords[groupIndex+1][i];
	} 
	alpha=acos(dot(vOH,vAH)/(norm(vOH)*norm(vAH)));

	return alpha*180.0/PI;
}

double Inos::norm(double v[3]){
	return sqrt(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));
}

double Inos::dot(double v1[3], double v2[3]){
	return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}

