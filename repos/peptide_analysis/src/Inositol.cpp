#include "Inositol.h"

Inositol::Inositol(int id, double** coords)
:inosId(id), OHcoords(coords),boundGroups(6)
{
	//note boundGroups(6) is eg of proper way to initialize vectors in constructors of classes
	//in theory each OH can bind 3 peptide groups ( accepts 2, donates 1 -- correct?)
}

//Destructor; free up memory taken up by the coordinates of the inositol
Inositol::~Inositol(){
	//free up coordinates mem
	for(int i=0; i<12; i++){
		delete [] OHcoords[i];
	}
	delete [] OHcoords;

	//free up bound groups mem
	for(int b=0; b<boundGroups.size(); b++){
		delete boundGroups.at(b);
	}
	
}

//gets the xyz coordinates of the Oxygen atom of an OH group at groupId
//groupId is in [0,5]
void Inositol::getOCoords(int groupId, double O[3]) {
	//groupId -> OHcoords O index as 2*groupId
	O[0] = OHcoords[2*groupId][0];
	O[1] = OHcoords[2*groupId][1];
	O[2] = OHcoords[2*groupId][2];
}

//gets the xyz coordinates of the H atom of an OH group at groupId
//groupId is in [0,5]
void Inositol::getHCoords(int groupId, double H[3]) {
	//groupId -> OHcoords H index as 2*groupId+1
	H[0] = OHcoords[2*groupId+1][0];
	H[1] = OHcoords[2*groupId+1][1];
	H[2] = OHcoords[2*groupId+1][2];
}

//sets the group(s) bound to OH group at OHgroupId of an inositol
//note that we didn't take into account whether it is the H or the O that the group is bound to 
void Inositol::setBoundGroup(int OHgroupId, PepGroup* aGroup) {
	boundGroups[OHgroupId] = aGroup;
}

//gets the group(s) bound to OH group at OHgroupId
PepGroup* Inositol::getBoundGroup(int OHgroupId) {
	return boundGroups[OHgroupId];
}

//prints the coordinates of each atoms of OH groups stored
//on separate lines to standard error
void Inositol::print() {
	for(int natoms=0; natoms<12; natoms++) {
		for(int i=0; i<3; i++) {
			cout<<OHcoords[natoms][i]<<"  ";
		}
		cout<<endl;
	}

}


//returns the total number of hydrogen bonds made by this inositol
//I do not expect this number to be greater than 18 (for a single OH group: O accepts 2, H donates 1)
int Inositol::totalNumHBs(){
	int total=0;
	for(int i=0; i<boundGroups.size(); i++){
		//the total number of HBs made by one inositol is 
		//the sum of number of bound groups at each OH group of the inositol over all OH groups of the inositol
		total+=boundGroups.at(i)->numGroups();
	}
	return total;
}

//returns the total number of bound groups of inositol
//note that this corresponds to the size of the boundGroups vector
int Inositol::totalBoundGroups(){
	int total = 0;
	for(int i=0; i<boundGroups.size();i++){
		if(boundGroups.at(i)->numGroups()){
			total++;
		}
	}
	return total;
}

//returns the number of chains the inositols is bound to
int Inositol::numChainsBound(int numChains){
	//eg num_boundtochain[0] is the number of inositol groups bound to chain 0
	//int num_boundtochain[4] = {0,0,0,0}; - this variable needs to be a vector; needs to change with # of chains!
	vector<int> num_boundtochain(numChains,0);
	for(int i=0; i<boundGroups.size(); i++){
		PepGroup* bg = boundGroups.at(i);
		if(bg->numGroups()!=0){
			num_boundtochain[bg->getPepId(0)]++;
		}
	}
	//count of number of chains the inositol is bound to
	int numchains=0;
	for(int chainid=0; chainid < numChains; chainid++){
		if(num_boundtochain[chainid]){
			numchains++;
		}
	}
	return numchains;
}



