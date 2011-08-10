//#define DEBUG_PRINT
//#define DEBUG_HB
//#define DEBUG_CHAIN_BOUND
//#define DEBUG_SEG_FAULT
//#define DEBUG

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "Inositol.h"
#include "Peptide.h"
#include "Water.h"

using namespace std;

const int numRes = 8; 

//reads the comment section and the number of atoms (what I call the header of a gro snapshot
bool read_header(ifstream &gro, string& comments, int& numAtoms){
	string atomnum;
	getline(gro,comments);
	if(gro.eof())
		return true;

	getline(gro,atomnum);
	numAtoms = atoi(atomnum.c_str());
	
	return false;
}

//reads the box size (assuming using rectangular boxes, no angle terms)
void read_box_dims(ifstream &gro, double box[3]){
	//string holding the box dimensions
	string boxdims;
	getline(gro, boxdims);
	//use a stringstream to parse boxdims into doubles
	istringstream ist(boxdims);
	ist>>box[0]>>box[1]>>box[2];
}

//reads the peptide portion of the trajectory
void readPeptides(ifstream &gro, vector<Peptide*> &peptides, int numPeptides){
	string line;
	string resName, atomName;
	int atomNum;
	//read in body of coordinates
	//read in peptide coordinates
	for(int npep = 0; npep < numPeptides; npep++){		
		int numPepAtoms = 5*numRes+4;		//number of residue backbone atoms in a peptide
		int numCterm = 2;			//number of atoms of the Cterminal methyl 
		int numNterm = 3;			//number of atoms of the N terminal NH2

		//initialize array to hold coordinates of a residue coords of a single peptide
		double** pepCoord = new double*[numPepAtoms];
		for(int i=0; i<numPepAtoms; i++){
			pepCoord[i]=new double[3];
		}
		
		//read and discard C-terminus
		for(int cterm=0; cterm < numCterm; cterm++){
			getline(gro, line);
		}

		//read peptide backbone coordinates + CB (for non polar contact analysis)
		for(int i=0; i<numPepAtoms; i++){
			getline(gro,line);
			istringstream ist(line);
			ist >> resName >> atomName >> atomNum;
			ist >> pepCoord[i][0] >> pepCoord[i][1] >> pepCoord[i][2];
			//cerr << pepCoord[i][0] << pepCoord[i][1] << pepCoord[i][2];
		}

		//read and discard N-terminus
		for(int nterm=0; nterm<numNterm; nterm++){
			getline(gro,line);
		}
		Peptide* aPep = new Peptide(npep, "GAGAGAGA", pepCoord);
		peptides.push_back(aPep);
	}
}

//reads the inositol portion of the trajectory
void readInositols(ifstream&gro, vector<Inositol*> &inositols, int numIns){
	string line;
	string resName, atomName;
	int atomNum;
	//read in inositol coordinates
	for(int nins = 0; nins < numIns; nins++){
		//initialize the coordinate array to hold coordinates of inositols
		double** inosCoords = new double*[12];
		for(int nOH=0; nOH<12; nOH++){
			inosCoords[nOH] = new double[3];
			getline(gro,line);
			istringstream ist(line);
			ist>>resName>>atomName>>atomNum;
			ist>>inosCoords[nOH][0] >> inosCoords[nOH][1] >> inosCoords[nOH][2];
		}
		Inositol* aInos=new Inositol(nins, inosCoords);
		inositols.push_back(aInos);
	}
}

//void readWaters(ifstream&gro, vector<Water*> &waters, int numWaters){
void readWaters(ifstream&gro, vector<Water*> &waters, double box[3]){
	string line;
	string resName, atomName;
	int atomNum;
	int resNum;
	//read in inositol coordinates
	//for(int nwat = 0; nwat < numWaters; nwat++){
	bool isWater=true;
	while(true) {
		//initialize the coordinate array to hold coordinates of waters
		double** waterCoords = new double*[3];
		for(int n=0; n<3; n++){
			waterCoords[n] = new double[3];
			getline(gro,line);
			istringstream ist(line);
			ist>>resNum>>resName>>atomName>>atomNum;

			if(resName != "SOL"){
				isWater = false;
				break;
			}

			ist>>waterCoords[n][0] >> waterCoords[n][1] >> waterCoords[n][2];
		}

		if(!isWater){
			break;
		}

		Water* aWater=new Water(waterCoords);
		waters.push_back(aWater);
	}
	//if we broke out of the while loop, this means the last line read in is the box dimensions
	istringstream ist(line);
	ist>>box[0]>>box[1]>>box[2];

	//cerr<<"Number of waters read = "<<waters.size()<<endl;
	//cerr<<box[0]<<" "<<box[1]<<" "<<box[2]<<endl;
}

//reads in a single snapshot that is part of a multi-snapshot gro file
//pre: peptides, and inositol are empty vectors, gro is the file stream containing the coordinate gro file
//post: one snapshot from the gro file stream is read, and the coordinates are parsed into vectors peptides, and inositols
//	returns true if end of file encountered, false if not end of file
bool readGroFile(ifstream& gro, vector<Peptide*>& peptides, vector<Inositol*>& inositols, vector<Water*> &waters, int numPeptides, int numIns, double boxDims[3]){
	string comments;	//dummy holder for comment of the gro snapshot
	int numAtoms=0;		//number of atoms in the snapshot
	string line;		//dummy holder for lines of the file

	//begin parsing gro file
	//read header of gro snapshot
	bool isEnd = read_header(gro, comments, numAtoms);
	
	if(isEnd){
		return true;
	}

	readPeptides(gro, peptides, numPeptides);
	readInositols(gro,inositols, numIns);
	//readWaters(gro,waters, numWaters);
	readWaters(gro,waters,boxDims);
	//read box dimensions
	//read_box_dims(gro,boxDims);

	return false;
}

//cleans up the space allocated for storing peptides, and inositol coords etc
void delete_vectors (vector<Peptide*> &peptides, vector<Inositol*> &inositols, vector<Water*> &waters){
	// *********delete memory allocated for snapshot*********
	for(int npep=1; npep <=peptides.size(); npep++){
		delete peptides.at(npep-1);
	}

	for(int nins=1; nins <=inositols.size(); nins++){
		delete inositols.at(nins-1);
	}
	for(int nwat=1; nwat <=waters.size(); nwat++){
		delete waters.at(nwat-1);
	}
	peptides.resize(0);
	inositols.resize(0);
	waters.resize(0);
}
