#include "Peptide.h"
//#define DEBUG_HB

Peptide::Peptide()
:pepId(0), sequence("GAGAGAGA"), backboneCoords(0), hb(0.35,0.25, 120), boundGroups(16)
{
}

//destructor; free up the memory taken up by the peptide backbone coordinates
//when the object goes out of memory
Peptide::~Peptide(){
	//NOTE: the number of atoms changed for the peptide -- find out the correct number
	for(int i=0; i<44; i++){
		delete [] backboneCoords[i];
	}
	delete [] backboneCoords;
}

Peptide::Peptide(int pId, const string& seq, double** bbCoords)
:pepId(pId), sequence(seq), backboneCoords(bbCoords), hb(0.35,0.25, 120),boundGroups(16)
{
}

//returns the coordinates of N of the NH group of residue resId
void Peptide::getNH(int resId, double NHcoords[3]){
	//the index of atoms in peptide had to be corrected because the amount of data input 
	//changed -- this reflects poor software engineering 
	int index = resId*5 + int(floor(resId/2.0));	

	NHcoords[0]=backboneCoords[index][0];
	NHcoords[1]=backboneCoords[index][1];
	NHcoords[2]=backboneCoords[index][2];
}

void Peptide::getHofNH(int resId, double Hcoords[3]){
	int index = resId*5+1+int(floor(resId/2.0));

	Hcoords[0]=backboneCoords[index][0];
	Hcoords[1]=backboneCoords[index][1];
	Hcoords[2]=backboneCoords[index][2];
}


//returns the coordinates of O of the CO group of residue resId
void Peptide::getO(int resId, double Ocoords[3]){
	int index = 3+resId*5+1 + int(ceil(resId/2.0));
	Ocoords[0]=backboneCoords[index][0];
	Ocoords[1]=backboneCoords[index][1];
	Ocoords[2]=backboneCoords[index][2];
}

//returns the coordinates of C of the CO group of residue resId
void Peptide::getC(int resId, double Ccoords[3]){
	int index = 3+resId*5 + int(ceil(resId/2.0));

	Ccoords[0]=backboneCoords[index][0];
	Ccoords[1]=backboneCoords[index][1];
	Ccoords[2]=backboneCoords[index][2];
}

void Peptide::getCbeta(int resId, double Cbcoords[3]){
	//pre: resId is an alanine 
	int index = 3+resId*5 + int(floor(resId/2.0));
	Cbcoords[0]=backboneCoords[index][0];
	Cbcoords[1]=backboneCoords[index][1];
	Cbcoords[2]=backboneCoords[index][2];
}

// char Peptide::getResidueName(int resId){
// 	return sequence[resId];
// }

void Peptide::print(){
	int numAtoms = 5*sequence.size();
	for(int natoms=0; natoms < numAtoms; natoms++){
		for(int i=0; i<3; i++){
			cout<<backboneCoords[natoms][i]<<"  ";
		}
		cout<<endl;
	}
}

//determines if the N-H group of resId is hbonded to "heavy"
//returns true if N-H is hbonded, else false
bool Peptide::isHbondedNH(int resId, double acceptor[3], double acceptorH[3], double box[3]){
	double N[3]; //Nitrogen of the NH group
	double H[3]; //Hydrogen of the NH group
	getNH(resId, N);
	getHofNH(resId,H);

#ifdef DEBUG_HB
	cout<<"isHbondedNH("<<resId<<")"<<endl;
	cout<<"dist(N,heavy)="<<dist(N,heavy,box)<<endl;
	cout<<"getNHAtomAngle(..)="<<getNHAtomAngle(resId, heavy)<<endl;
	cout<<"dist(H,heavy)="<<dist(H,heavy,box)<<endl;
#endif

	if(hb.isHbonded(N, H, acceptor, acceptorH, box))
		return true;
	return false;
}


//determines if CO group of resId is hbonded to "donorHeavy"
//return true if hbonded, else false
bool Peptide::isHbondedCO(int resId, double donorHeavy[3], double donorH[3], double box[3]) {

	double O[3],C[3];
	getO(resId,O);
	getC(resId,C);


#ifdef DEBUG_HB
	cout<<"isHbondedCO("<<resId<<")"<<endl;
	cout<<"dist(D,O)="<<dist(O,donorHeavy,box)<<endl;
	cout<<"getAtomHCOAngle(..)="<<getAtomHCOAngle(resId, donorHeavy, donorH)<<endl;
	cout<<"dist(O,H)="<<dist(O,donorH,box)<<endl;
#endif

	if(hb.isHbonded(donorHeavy,donorH,O,C,box))
		return true;
	return false;
}

//calculates the number of intra hydrogen bonds for this peptide
int Peptide::numIntraHB(double box[3]){
	//ofstream out("intra.hb");
	int numIntra=0;
	int numRes = sequence.size();
	for(int n=0; n < numRes; n++){
		double NH_pep[3], H_pep[3], O_pep[3], C_pep[3];
		getNH(n, NH_pep);
		getHofNH(n, H_pep);
		getO(n, O_pep);
		getC(n, C_pep);
		for(int m=n+1; m < numRes; m++){
			if(isHbondedCO(m, NH_pep, H_pep, box)){
				numIntra++;
			}
			if(isHbondedNH(m, O_pep, C_pep, box)){
				numIntra++;
			}
		}
	}
	
	return numIntra;
}

//sets the backbone hydrogen bonding group to it's bound groups
//bbgroup -- index of the backbone hydrogen bonding group
void Peptide::setBoundGroup(int bbgroup, PepGroup* bgroup) {
	boundGroups[bbgroup] = bgroup;
}

PepGroup* Peptide::getBoundGroup(int bbgroup){
	return boundGroups.at(bbgroup);
}
