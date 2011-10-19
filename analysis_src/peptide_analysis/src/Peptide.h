#ifndef _PEPTIDE_H_
#define _PEPTIDE_H_

#include <string>
#include <iostream>
#include <cmath>
#include "HB.h"
#include "PepGroup.h"

using namespace std;

class Peptide {
	public:
		Peptide();	//default constructor
		~Peptide();
		Peptide(int pId, const string& seq, double** bbCoords);
		void getNH(int resId, double NHcoords[3]);  //get coord of backbone Nitrogen atom
		void getO(int resId, double Ocoords[3]);  //get coord of backbone Oxygen atom
		void getC(int resId, double Ccoords[3]);  //get coord of backbone Oxygen atom
		void getCbeta(int resId, double Cbcoords[3]); // get coord of cbeta atom of methyl sidechains
		void getHofNH(int resId, double Hcoords[3]); //get backbone H
		bool isHbondedNH(int resId, double acceptor[3], double acceptorH[3], double box[3]);
		bool isHbondedCO(int resId, double donorHeavy[3], double donorH[3], double box[3]);
		int numIntraHB(double box[3]);
		void setBoundGroup(int bbgroup, PepGroup* bgroup);
		PepGroup* getBoundGroup(int bbgroup);
		void print();

	private:
		int pepId; 	//id of the peptide because a system might have multiple strands
		string sequence; //the amino acid sequence of the peptide
		double** backboneCoords;
		HB hb;
		vector<PepGroup*> boundGroups;
};
#endif


