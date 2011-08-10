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



//count and return the total number inter- peptide strand hydrogen bonds in a snapshot
int countInterpeptideHB(vector<Peptide*> &peptides, double box[3]){
	int totalNumInter=0;
	int totalNumPep = peptides.size();

	for(int i=0; i < totalNumPep; i++) {
		Peptide* pepA = peptides.at(i);
		for(int n=0; n < numRes; n++) {
			double NH_pepA[3], H_pepA[3], O_pepA[3], C_pepA[3];
			pepA->getNH(n,NH_pepA);
			pepA->getHofNH(n,H_pepA);
			pepA->getO(n,O_pepA);
			pepA->getC(n,C_pepA);
			for(int j = i+1; j < totalNumPep; j++) {
				Peptide* pepB = peptides.at(j);
				for(int m=0; m < numRes; m++) {
					if(pepB->isHbondedNH(m,O_pepA, C_pepA, box)){
						totalNumInter++;
					}
					if(pepB->isHbondedCO(m,NH_pepA, H_pepA, box)){
						totalNumInter++;
					}
				}
			}
		}
	}
	return totalNumInter;
}



//prints each entry of vector v separated by a single space, on a single line
void print_vector(ostream &out, const vector<int> &v){
	for(int i=1; i<=v.size(); i++){
		out<<v[i-1]<<" ";
	}
}

void print_data(const vector<Peptide*> &peptides, const vector<Inositol*> &inositols, const vector<Water*> &waters){
	cout<<"peptides"<<endl;
	for(int i=0; i<peptides.size(); i++){
		peptides.at(i)->print();
	}
	cout<<"inositols"<<endl;
	for(int i=0; i<inositols.size(); i++){
		inositols.at(i)->print();
	}
	cout<<"water"<<endl;
	for(int i=0; i<waters.size(); i++){
		waters.at(i)->print();
	}
}

//computes the (vector) distribution of number of inositols over #HB formed
void insHBdistr(vector<int> &distr, vector<int> &count){
	for(int i=1; i<=count.size(); i++){
		distr[ count.at(i-1) ]++;
	}
}

void waterContacts(const vector<Peptide*> &peptides, const vector<Inositol*> &inositols, 
		const vector<Water*> &waters, int* water_bound_res, double boxDims[3]) {
		double OWcoord[3],HW1coord[3],HW2coord[3];
		double N[3],H[3];
		int numPeptides = peptides.size();
		int numWaters = waters.size();

		for(int npep = 0; npep < numPeptides; npep++) {
			for(int nres = 0; nres < numRes; nres++) {
				for(int nw=0; nw < numWaters; nw++) {
					waters.at(nw)->getOW(OWcoord);
					waters.at(nw)->getHW1(HW1coord);
					waters.at(nw)->getHW2(HW2coord);

					if( peptides.at(npep)->isHbondedNH(nres, OWcoord , HW1coord, boxDims) ||
					peptides.at(npep)->isHbondedNH(nres, OWcoord , HW2coord, boxDims) ) {
#ifdef DEBUG_BBWATER
						cout<<"NH "<<npep+1<<" "<<nres+1<<" "<<resid+1<<" "<<nw+43<<endl;
						for(int i=0; i<3; i++){
							cout<<OWcoord[i]<<" "<<HW1coord[i]<<" "<<HW2coord[i]<<endl;
						}
#endif

						peptides.at(npep)->getNH(nres,N);
						peptides.at(npep)->getHofNH(nres,H);

#ifdef DEBUG_BBWATER
						cout<<"N "<<endl;
						for(int i=0; i<3; i++){
							cout<<N[i]<<" ";
						}
						cout<<endl;
						cout<<"H "<<endl;
						for(int i=0; i<3; i++){
							cout<<H[i]<<" ";
						}
						cout<<endl;
						cout<<"water bound = "<<nw+43<<endl;
#endif
						water_bound_res[2*nres]++;
					}

					if( peptides.at(npep)->isHbondedCO(nres, OWcoord, HW1coord, boxDims) ||
					peptides.at(npep)->isHbondedCO(nres, OWcoord, HW2coord, boxDims) ){
#ifdef DEBUG_BBWATER						
						cout<<"CO "<<npep+1<<" "<<nres+1<<" "<<resid+1<<" "<<nw+43<<endl;
						
						for(int i=0; i<3; i++){
							cout<<OWcoord[i]<<" "<<HW1coord[i]<<" "<<HW2coord[i]<<endl;
						}
#endif
						double C[3],O[3];
						peptides.at(npep)->getO(nres,O);
						peptides.at(npep)->getC(nres,C);
#ifdef DEBUG_BBWATER
						cout<<"C "<<endl;			
						for(int i=0; i<3; i++){
							cout<<C[i]<<" ";
						}
						cout<<endl;
						cout<<"O "<<endl;
						for(int i=0; i<3; i++){
							cout<<O[i]<<" ";
						}
						cout<<endl;
						cout<<"water bound = "<<nw + 43<<endl;
#endif
						water_bound_res[2*nres+1]++;
					}
				}
			}
		}

#ifdef DEBUG_BBWATER	
		for(int i=0; i < total_bound_groups; i++){
			if(water_bound_res[i]){
				cout<<i+2<<" "<<water_bound_res[i]<<endl;
			}
		}
		cout<<endl;
#endif
}

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
	
	int total_bound_groups = numPeptides*16;

	//assumes that max number of HB that can be made by an inositol to be 6
	vector<int> distr_inos_numHB_tot(7, 0); //distribution of number of inositols over #HB
	ofstream contactwat((base + "_water_contact").c_str());
	//ofstream nonpolar((base + "_np_contact").c_str());
	ofstream cmap((base + "_hb_contact_map").c_str());
	ofstream inos((base + "_inos_contact").c_str());
	//read each snapshot of the partial gro file
	//build vector of peptides, inositols, detect hydrogen bonds,collect statistics
	//delete memory allocated for the snapshots
	while(!readGroFile(gro, peptides, inositols, waters, numPeptides, numIns, boxDims)){
#ifdef DEBUG
		cout<<"################frame # "<<totalNumSnap<<" #####################"<<endl;
#endif
		int* water_bound_res = new int[total_bound_groups];
		//initialize all to 0
		for(int i=0; i<total_bound_groups; i++){
			water_bound_res[i]=0;
		}
		//calculate water contacts
		waterContacts(peptides, inositols, waters, water_bound_res, boxDims);
		//output to file
		for(int i=0; i < total_bound_groups; i++){
			contactwat<<water_bound_res[i]<<" ";
		}
		contactwat<<endl;
		//de-allocate memory
		delete [] water_bound_res;

		Peptide* pep = peptides.at(0);
		// Contact maps calculations
		for(int n=0; n < numRes; n++){
			double NH_pep[3], H_pep[3], O_pep[3], C_pep[3];
			pep->getNH(n, NH_pep);
			pep->getHofNH(n, H_pep);
			pep->getO(n, O_pep);
			pep->getC(n, C_pep);
			for(int m=n+1; m < numRes; m++){
				if(pep->isHbondedCO(m, NH_pep, H_pep, boxDims)){
					cmap << n << " " << m <<endl;
				}
				if(pep->isHbondedNH(m, O_pep, C_pep, boxDims)){
					cmap << m << " "<< n <<endl;
				}
			}
		}
		
		//inositol peptide contact calculations
		//we only have 1 peptide in the system
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
						//cinos<<nres<<" has "<<nins*6+noh<<" NH bound"<<endl;
						bgroupNH->addPepGroup(nins,nins*6+noh,"OH");
					}
			
					if (pep->isHbondedCO(nres, inosO, inosH, boxDims)) {
						//cinos<<nres<<" has "<<nins*6+noh<<" CO bound"<<endl;
						bgroupCO->addPepGroup(nins,nins*6+noh,"OH");
					}
				}
			}
			pep->setBoundGroup(2*nres, bgroupNH);
			pep->setBoundGroup(2*nres+1, bgroupCO);
			//delete bgroupNH;
			//delete bgroupCO;
		}
	
		//inosput the data computed above
		for(int nbb = 0; nbb < 16; nbb++) {
			PepGroup* bgroup=pep->getBoundGroup(nbb);
			int numBoundGroups = bgroup->numGroups();
			if(numBoundGroups){
				for(int nbgroup=0; nbgroup<numBoundGroups; nbgroup++){
					inos<<bgroup->getResId(nbgroup);
					if(nbgroup<numBoundGroups-1){
						inos<<" ";
					}
				}
			}else{
				inos<<"-";
			}
			if(nbb%2 == 0 && nbb < 16){
				inos<<";";
			}
			if(nbb%2 && nbb < 15){
				inos<<"|";
			}
	
				//cinos<<nbb<<" has " << bgroup->numGroups()<< " bound groups"<<endl;
	
		}
		inos<<endl;
	
		time++;
		delete_vectors(peptides, inositols, waters);
	}

	contactwat.close();
	cmap.close();
	//nonpolar.close();
	inos.close();
	return 0;
}

