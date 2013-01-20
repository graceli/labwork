#define CPLUSPLUS

using namespace std;

//c++ libraries
#include <map>
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>

//gromacs c libraries
extern "C" {
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>
#include <stdlib.h>
#include "sysstuff.h"
#include "string.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "xvgr.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "index.h"
#include "pdbio.h"
#include "confio.h"
#include "rmpbc.h"
#include "names.h"
#include "atomprop.h"
#include "physics.h"
#include "gstat.h"
#include "tpxio.h"
}

//#define DEBUG_INOS

class PheMolecule {
	// This class represents the phenyl group of phenylalanine (with the Cbeta atom attached)
	public:
		PheMolecule(int residueNumber, const string &residueName);
//		~PheMolecule();
		void addAtom(real* atom);
		void center_of_geometry(real* com);
//		void addResidueNumber(int residueNumber);
//		void addResidueName(const string &residueName);

	private:
		vector<real*> atoms;
		int residueNumber;
		string residueName;
};

PheMolecule::PheMolecule(int resNumber, const string &resName) {
	residueNumber = resNumber;
	residueName.assign(resName);
}

void PheMolecule::addAtom(real* atom) {
	atoms.push_back(atom);
}

void PheMolecule::center_of_geometry(real* com) {
	com[0] = 0;
	com[1] = 0;
	com[2] = 0;

	for(int i = 0; i < atoms.size(); i++) {
		for(int d = 0; d < DIM; d++) {
			com[d] += atoms.at(i)[d];
		}
	}

	for(int d = 0; d < DIM; d++) {
		com[d] = com[d] / atoms.size();
	}
}

//void PheMolecule::addResidueNumber(int residueNumber) {
//	residueNumber = residueNumber;
//}
//
//void PheMolecule::addResidueName(const string &residueName) {
//	residueName(residueName);
//}


//takes a topology structure
//and outputs atom index, atom name, and the residue the atom is in
//do not explicitly exclude residues (such as lysine or glutamate from calculations)
void dump_top ( t_topology *top ) {
	int totalNumAtoms = top->atoms.nr;
	int totalNumResidues = top->atoms.nres;

#ifdef DEBUG
	cerr<<"Number of atoms = " << totalNumAtoms <<endl;
	cerr<<"Number of residues = " << totalNumResidues <<endl;
#endif

	for(int i = 0; i < totalNumAtoms; i++){
		int atomResidueNum = top->atoms.atom[i].resnr;  //residue number corresponding
		int atomIndex = i+1;
		cout << atomIndex << " " << *(top->atoms.atomname[ atomIndex-1 ]) << " " << *(top->atoms.resname[ atomResidueNum ]) <<endl;
	}
}

void dump_index ( t_topology *top, atom_id ** index, int* isize, int ngrps ) {
	for(int g=0; g < ngrps; g++){
		for(int gIndex=0; gIndex < isize[g]; gIndex ++){
			cout<<g<< " " << index[g][gIndex]<< " " <<*(top->atoms.atomname[ gIndex ])<<" "<<*(top->atoms.resname[ top->atoms.atom[gIndex].resnr ])<<endl;
		}
		cout<<endl;
	}
}

//int calculate_com ( t_topology *top, atom_id** index, rvec* x, int groupNum, int groupSize, real* com) {
//	com[0] = 0;
//	com[1] = 0;
//	com[2] = 0;
//
//	int num_atoms_in_residue = 0;
//	int atomIndex=0;
//	for(int i = 0; i < groupSize; i++) {
//		atomIndex = index[groupNum][i];
//
//		inositol_com.push_back(com);
//
//		for(int d = 0; d < DIM); d++) {
//			com[d] += x[atomIndex][d];
//		}
//	}
//
//	for(int d = 0; d < DIM; d++) {
//		com[d] = com[d] / num_atoms_in_residue;
//	}
//
//
//#ifdef DEBUG_COM
//	cerr << "current = " << currentResnum << " " << num_atoms_in_residue << endl;
//	cerr << com[0]*10 << " " << com[1]*10 << " " << com[2]*10 << endl;
//#endif
//}

int calculate_com_for_inositols ( t_topology *top, atom_id** index, rvec* x, int start_group_num, int num_inositols,
								  int group_size, vector<real*> &inositol_com, vector<int> &inositol_residue_indices ) {

	real* com = new real[3];
	com[0] = 0;
	com[1] = 0;
	com[2] = 0;

	int resnum = 0;
	int atomIndex = 0;

	for(int ins_idx = 0; ins_idx < num_inositols; ins_idx++) {
		for(int i = 0; i < group_size; i++ ) {
			atomIndex = index[ins_idx][i];
			for(int d = 0; d < DIM; d++) {
				com[d] += x[atomIndex][d] / group_size;
			}
			resnum = top->atoms.atom[atomIndex].resnr;
		}
		inositol_com.push_back(com);
		inositol_residue_indices.push_back(resnum);
	}
}

void delete_vector(vector<real*> &v) {
	for(int i = 0; i < v.size(); i++) {
		delete v[i];
	}
}

void print_map(map<string,int>& m, ofstream& f_residue_table) {
	map<string,int>::iterator iter;
	for(iter=m.begin(); iter!=m.end(); iter++) {
		f_residue_table<<(*iter).first<<" "<<(*iter).second<<endl;
	}
}

//initialize a map with the same keys as a reference map
void initialize_maps(map<string, int>& uninit, map<string, int>&ref) {
	map<string,int>::iterator ref_iter;
	for(ref_iter=ref.begin(); ref_iter!=ref.end(); ref_iter++) {
		uninit[ref_iter->first]=0;
	}
}

void initialize_vectors(vector<int> &uninit, map<string,int> &ref){
	uninit.resize(ref.size(),0);
}

bool is_in_contact(t_pbc* pbc, rvec p1, rvec p2, real cutoff, real &dist) {
	rvec dx;
	pbc_dx(pbc, p1, p2, dx);
	dist = norm(dx);
	if(dist < cutoff) {
		return true;
	}

	return false;
}

int main(int argc,char *argv[]) {
	static char *desc[] = {
			" This program calculates the number of contacts for each atom in each of the",
			" residue of the first group to the center of mass of each of the residues in",
			" the second group",
			" Specifically, it is designed for computing the contacts per residue to the C.O.M of",
			" inositol molecules in the system",
			" v2 outputs the number of contacts for each residue in the protein input group (atom)",
			" to C.O.M of inositol"
	};

	t_topology *top=NULL;
	rvec *xtop;
	real t,cut2,dist2;
	rvec *x=NULL,*v=NULL,dx;
	matrix box;
	int status;
	int g=0,d=0,i=0,j=0,res=0,teller=0;
	atom_id aid;

	atom_id **max;   /* the index for the atom numbers */
	real    *mass;
	FILE    *fp=NULL;
	t_pbc   pbc;

	const char  *leg[4] = { "|d|","d\\sx\\N","d\\sy\\N","d\\sz\\N" };

	// App (non-GROMACS) variables
	rvec    *com;
	rvec    *com_phe;
	//static char *f_contact="residue_nonpolar_contact.dat";
	static char* per_phe_stacking_fname = "per_phe_stacking.dat";
	static char* per_inositol_stacking_fname = "per_inositol_stacking.dat";
	static char* stacking_info = "stacking_info.dat";
	static int numInositols = 0;
	static int NPHE = 0;
	real CUTOFF = 0.45;     /*the np cbeta-inositol cutoff in nm*/

	static t_pargs pa[] = {
		{ "-dist",      FALSE, etREAL, {&CUTOFF},
		"the cutoff"},

		{"-per_phe_stacking_fname", FALSE, etSTR, {&per_phe_stacking_fname},
		"Dataset describing whether the Phe residues are stacked, bound, or not for each time frame.  "
		"Stacked = 2; Bound = 1; Unbound = 0"},

		{"-per_inositol_phe_contacts", FALSE, etSTR, {&per_inositol_stacking_fname},
		"Dataset describing whether each inositol molecule is stacked, bound, or not for each time frame. "
		"Stacked = 2; Bound = 1; Unbound = 0"},

		{"-stacking_info", FALSE, etSTR, {&stacking_info},
		"Dataset with the residue numbers of (phe, inosito) pairs which are stacked. This data file is useful for "
		"debugging purposes."},

		{"-num_inositols", FALSE, etINT, {&numInositols},
			"number of inositols to read in from index"
		}
	};

	t_filenm fnm[] = {
			{ efTPS, NULL, NULL, ffREAD },
			{ efTRX, "-f", NULL, ffREAD },
			{ efNDX, NULL, NULL, ffOPTRD }
			/* { efXVG, NULL, "dist", ffOPTWR },*/
	};

#define NPA asize(pa)
#define NFILE asize(fnm)

	int ngrps = 1 + numInositols;     /* the number of index groups */
	atom_id **index = new atom_id*[ngrps];
	int *isize = new int[ngrps];    /* the size of each group */
	char **grpname = new char*[ngrps]; /* the name of each group */

    // Index of the protein group
    // By setting this index to 0, this forces the user to
    // specify the protein group first when running this tool
    // Here we assume that there is only one protein group
    const int PROTEIN_GROUP_START_IDX = 0;

    // Start index of the inositol group
    const int INOSITOL_GROUP_START_IDX = 1;

	// Number of atoms in the protein group
    const int NUM_ATOMS_PROTEIN = isize[PROTEIN_GROUP_START_IDX];

    // Number of atoms in the inositol group
    const int NUM_ATOMS_INOSITOL = isize[INOSITOL_GROUP_START_IDX];

	char title[STRLEN];

	CopyRight(stderr, argv[0]);
	parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_BE_NICE, NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL);

	int ePBC = -1;
	top = read_top(ftp2fn(efTPS, NFILE, fnm), &ePBC);
	get_index(&top->atoms, ftp2fn(efNDX, NFILE, fnm), ngrps, isize, index, grpname);
	int natoms = read_first_x(&status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);

	string residue_name;
	int residue_id;
	bool first_time = true;

//	ofstream f_inos_phe_contacts(perInositolPheContacts);
//	vector<int> inositol_residue_indices;
//	do {
//		set_pbc(&pbc, -1, box);
//		// GMX v 4.0.5 I think this is needed to make system whole
//		rm_pbc(&top->idef, ePBC, natoms, box, x, x);
//
//		vector<real*> inositol_com;
//		vector<real*> phe_atoms;
//		real* phe_com = new real[3];
//		real dist;
//
//		// Calculate the COM of all the inositols in the system and populate a list of these vectors
//		int num_inositols_found = calculate_com_for_inositols(top, index, x, INOSITOL_GROUP_START_IDX, numInositols,
//																	 NUM_ATOMS_INOSITOL, inositol_com, inositol_residue_indices);
//
//		// TODO: Parse out the coordinates of all PHE residues in a list of PHE molecules
//
//		// TODO: Calculate the COM of the phe residues -- populate a list of phe center of masses
//
//		for(int i = 0; i < inositol_com.size(); i++) {
//			for(int protein_atom_num = 0; protein_atom_num < NUM_ATOMS_PROTEIN; protein_atom_num++) {
//				// TODO: get the residue name
//	        	int protein_atom_idx = index[PROTEIN_GROUP_START_IDX][protein_atom_num];
//	            residue_id = top->atoms.atom[protein_atom_idx].resnr;
//	            residue_name = *(top->atoms.resname[residue_id]);
//
//				if(residue_name == "PHE") {
//					phe_atoms.push_back(phe_atom);
//				}
//
//				if(phe_atoms.size() == 6) {
//					calculate_com(top, index, phe_atoms, phe_com);
//
//					if(is_in_contact(&pbc, inositol_com[i], phe_com, CUTOFF, dist)) {
//						// TODO: Calculate the planar angle between this inositol molecule and the PHE residue that
//						// it is in contact with
//
//#ifdef DEBUG_COM
//						cout <<t<<" "<< inositol_residue_indices[i] << " ";
//						for(int d=0; d<3;d++){
//							cout<<inositol_com[i][d]<<" ";
//						}
//
//						cout << phe_table[j] << " ";
//						for(int d=0; d<3; d++){
//							cout<<phe_com[j][d]<<" ";
//						}
//						cout<<dist<<endl;
//#endif
//
//						// TODO: Perform various result outputs
//					}
//				}
//			}
//		}
//		delete_vector(inositol_com);
//	} while (read_next_x(status, &t, natoms, x, box));
}
