#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "Inositol.h"
#include "Peptide.h"
//#include "Water.h"
#include "FileRead.h"


int main(int argc, char* argv[]) {
	
	if(argc < 5){
		cerr<<"usage: peptide_analysis_testing_simple <gro-file> <num-pep> <num-ins> <base>"<<endl;
		return 0;
	}
	
	ifstream gro(argv[1]);
	int numPeptides = atoi(argv[2]);	//number of peptides in the snapshot
	int numIns = atoi(argv[3]);		//number of inositols in the snapshot
	string base(argv[4]);

	vector<Peptide*> peptides;
	vector<Inositol*> inositols;
	vector<Water*> waters;

	int time = 0;
	double boxDims[3];


	ofstream out("count");
	while(!readGroFile(gro, peptides, inositols, waters, numPeptides, numIns, boxDims)) {
		//we only have 1 peptide in the system
		Peptide* pep = peptides.at(0);
		for(int nres = 0; nres < numRes; nres++) {
			PepGroup* bgroupNH = new PepGroup;
			PepGroup* bgroupCO = new PepGroup;
			for(int nins = 0; nins < numIns; nins++) {
				Inositol* aInos = inositols.at(nins);
				for(int noh = 0; noh < 6; noh++) {
					double inosO[3], inosH[3];
					aInos->getOCoords(noh,inosO);
					aInos->getHCoords(noh,inosH);
					if (pep->isHbondedNH(nres, inosO, inosH, boxDims)) {
						//cout<<nres<<" has "<<nins*6+noh<<" NH bound"<<endl;
						bgroupNH->addPepGroup(nins,nins*6+noh,"OH");
					}
			
					if (pep->isHbondedCO(nres, inosO, inosH, boxDims)) {
						//cout<<nres<<" has "<<nins*6+noh<<" CO bound"<<endl;
						bgroupCO->addPepGroup(nins,nins*6+noh,"OH");
					}
				}
			}
			pep->setBoundGroup(2*nres, bgroupNH);
			pep->setBoundGroup(2*nres+1, bgroupCO);
			//delete bgroupNH;
			//delete bgroupCO;
		}

		//output the data computed above
		for(int nbb = 0; nbb < 16; nbb++) {
			PepGroup* bgroup=pep->getBoundGroup(nbb);
			int numBoundGroups = bgroup->numGroups();
			if(numBoundGroups){
				for(int nbgroup=0; nbgroup<numBoundGroups; nbgroup++){
					out<<bgroup->getResId(nbgroup);
					if(nbgroup<numBoundGroups-1){
						out<<" ";
					}
				}
			}else{
				out<<"-";
			}
			if(nbb%2 == 0 && nbb < 16){
				out<<";";
			}
			if(nbb%2 && nbb < 15){
				out<<"|";
			}

				//cout<<nbb<<" has " << bgroup->numGroups()<< " bound groups"<<endl;

		}
		out<<endl;

		time++;
		delete_vectors(peptides, inositols, waters);
	}
	out.close();
}
