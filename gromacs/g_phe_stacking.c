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
		~PheMolecule();
		void add_atom(int atom_index, real* atom);
		void center_of_geometry(real* com);
		void get_angle_index_group(atom_id* index);
		int get_resid();
		const string& get_resname();
        int get_num_atoms();
        void print_info();
 
	private:
		vector<real*> atoms;
		vector<int> atom_indices;
		int residueNumber;
		string residueName;
};

PheMolecule::PheMolecule(int resNumber, const string &resName) {
	residueNumber = resNumber;
	residueName.assign(resName);
}

PheMolecule::~PheMolecule() {
	for(int i = 0; i<atoms.size(); i++) {
		delete atoms.at(i);
	}
}

void PheMolecule::print_info() {
    cout << "resName:" << residueName << endl;
    cout << "resNumber:" << residueNumber << endl;
    cout << "atom indies: ";
    for(int i = 0; i < atom_indices.size(); i++) {
        cout << atom_indices[i] << " ";
    }
    cout << endl;

    atom_id* index = new atom_id[3];
    get_angle_index_group(index); 
    for(int i = 0; i < 3; i++) {
        cout << index[i] << " ";
    }
    cout << endl;
}

int PheMolecule::get_num_atoms() {
    return atoms.size();
}

const string& PheMolecule::get_resname() {
	return residueName;
}

int PheMolecule::get_resid() {
	return residueNumber;
}

void PheMolecule::add_atom(int atom_index, real* atom) {
	atom_indices.push_back(atom_index);
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

void PheMolecule::get_angle_index_group(atom_id* index) {
	// Pick every 2nd atom to use for the planar group for plane-plane angle calculation
	index[0] = atom_indices.at(1);
	index[1] = atom_indices.at(3);
	index[2] = atom_indices.at(5);
}

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
		cout << atomIndex << " " << *(top->atoms.atomname[ atomIndex-1 ]) << " " << atomResidueNum << " " << *(top->atoms.resname[ atomResidueNum ]) <<endl;
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

/* 	GROMACS SGANGLE START */
static void calculate_normal(atom_id index[],rvec x[],rvec result,rvec center) {
	rvec c1,c2;
	int i;

	/* calculate centroid of triangle spanned by the three points */
	for(i=0;i<3;i++)
		center[i] = (x[index[0]][i] + x[index[1]][i] + x[index[2]][i])/3;

	/* use P1P2 x P1P3 to calculate normal, given three points P1-P3 */
	rvec_sub(x[index[1]],x[index[0]],c1);    /* find two vectors */
	rvec_sub(x[index[2]],x[index[0]],c2);

	cprod(c1,c2,result);                    /* take crossproduct between these */
}

/* calculate the angle and distance between the two groups */
static void calc_angle(int ePBC,matrix box,rvec x[], atom_id index1[],
		atom_id index2[], int gnx1, int gnx2,
		real *angle,      real *distance,
		real *distance1,  real *distance2)

/* distance is distance between centers, distance 1 between center of plane
   and atom one of vector, distance 2 same for atom two
 */
{
	rvec
	normal1,normal2,  	/* normals on planes of interest */
	center1,center2,  	/* center of triangle of points given to define plane,*/
	/* or center of vector if a vector is given */
	h1,h2,h3,h4,h5;  	/* temp. vectors */
	t_pbc pbc;

	set_pbc(&pbc,ePBC,box);

	switch(gnx1)
	{
	case 3:           /* group 1 defines plane */
		calculate_normal(index1,x,normal1,center1);
		break;
	default:          /* group 1 does none of the above */
		gmx_fatal(FARGS,"Something wrong with contents of index file.\n");
	}

	switch(gnx2)
	{
	case 3:          /* group 2 defines plane */
		calculate_normal(index2,x,normal2,center2);
		break;
	default:         /* group 2 does none of the above */
		gmx_fatal(FARGS,"Something wrong with contents of index file.\n");
	}

	*angle = cos_angle(normal1,normal2);

	if (box)
		pbc_dx(&pbc,center1,center2,h3);
	else
		rvec_sub(center1,center2,h3);
	*distance = norm(h3);

	// Case 3
	*distance1 = 0; *distance2 = 0;
}
/* GROMACS SGANGLE END */

void populate_inositol_coordinates(t_topology *top, atom_id** index, rvec* x, int start_group_num, int num_inositols,
								int group_size, vector<PheMolecule*> &inositol_molecules) {
	int atomIndex = 0;
	PheMolecule* aInositol;
	for(int ins_idx = start_group_num; ins_idx < start_group_num + num_inositols; ins_idx++) {
		int resid = top->atoms.atom[index[ins_idx][1]].resnr;
		aInositol = new PheMolecule(resid, "INS");
		for(int i = 0; i < group_size; i++) {
			atomIndex = index[ins_idx][i];
			aInositol->add_atom(atomIndex, x[atomIndex]);
		}
		inositol_molecules.push_back(aInositol);
	}
}

void populate_phe_coordinates(t_topology *top, atom_id** index, rvec* x, int group_start_num, int group_size,
							  vector<PheMolecule*> &phe_molecules) {
	// Parse out the coordinates of all PHE residues in a list of PHE molecules
	bool firstAtom = true;
	PheMolecule* aPheMolecule = 0; 
	int residue_id;
	string residue_name;
	for(int protein_atom_num = 0; protein_atom_num < group_size; protein_atom_num++) {
		// Get the residue name corresponding to this protein atom number
		int protein_atom_idx = index[group_start_num][protein_atom_num];
		residue_id = top->atoms.atom[protein_atom_idx].resnr;
		residue_name = *(top->atoms.resname[residue_id]);
		if(residue_name.compare("PHE") == 0) {
			if(firstAtom) {
				aPheMolecule = new PheMolecule(residue_id, residue_name);
				aPheMolecule->add_atom(protein_atom_idx, x[protein_atom_idx]);
				firstAtom = false;
			} else {
				// Store both atom index and its coordinates
				aPheMolecule->add_atom(protein_atom_idx, x[protein_atom_idx]);
			}
		}

        if(aPheMolecule != 0 && aPheMolecule->get_num_atoms() == 7) {
            phe_molecules.push_back(aPheMolecule);
            firstAtom = true;
            aPheMolecule = 0;
        }
	}
}

void delete_vector(vector<PheMolecule*> &v) {
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
	static char* per_phe_bound_fname = "per_phe_bound.dat";
	static char* per_inositol_bound_fname = "per_inositol_bound.dat";
	static char* stacking_info = "stacking_info.dat";
	static int numInositols = 0;
	static int NPHE = 0;
	real CUTOFF = 0.45;     /*the np cbeta-inositol cutoff in nm*/

	static t_pargs pa[] = {
		{ "-dist",      FALSE, etREAL, {&CUTOFF},
		  "cutoff used for com-com calculations"},

		{ "-per_phe_stacking_fname", FALSE, etSTR, {&per_phe_stacking_fname},
		  "Number of times each Phe residues are stacked"},

		{ "-per_phe_bound_fname", FALSE, etSTR, {&per_phe_stacking_fname},
		  "Number of times each phe is bound"},

		{ "-per_inositol_stacking_fname", FALSE, etSTR, {&per_inositol_stacking_fname},
		  "Dataset describing the number of times each inositol molecule is stacked"},

		{ "-per_inositol_bound_fname", FALSE, etSTR, {&per_inositol_bound_fname},
		  "Dataset describing the number of times inositol molecule is stacked."},

		{ "-stacking_info", FALSE, etSTR, {&stacking_info},
		  "Dataset with the residue numbers of (phe, inosito) pairs which are stacked. This data file is useful for "
		  "debugging purposes."},

		{ "-num_inositols", FALSE, etINT, {&numInositols},
		  "number of inositols to read in from index"}
	};

	t_filenm fnm[] = {
			{ efTPS, NULL, NULL, ffREAD },
			{ efTRX, "-f", NULL, ffREAD },
			{ efNDX, NULL, NULL, ffOPTRD }
			/* { efXVG, NULL, "dist", ffOPTWR },*/
	};

#define NPA asize(pa)
#define NFILE asize(fnm)

	char title[STRLEN];

	CopyRight(stderr, argv[0]);
	parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_BE_NICE, NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL);

	int ngrps = 1 + numInositols;     /* the number of index groups */
	atom_id **index = new atom_id*[ngrps];
	int *isize = new int[ngrps];    /* the size of each group */
	char **grpname = new char*[ngrps]; /* the name of each group */

	int ePBC = -1;
	top = read_top(ftp2fn(efTPS, NFILE, fnm), &ePBC);
	get_index(&top->atoms, ftp2fn(efNDX, NFILE, fnm), ngrps, isize, index, grpname);

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

	int natoms = read_first_x(&status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);

	string residue_name;
	int residue_id;
	bool first_time = true;

	ofstream f_inos_stacking(per_inositol_stacking_fname);
	ofstream f_inos_bound(per_inositol_bound_fname);
	ofstream f_phe_stacking(per_phe_stacking_fname);
	ofstream f_phe_bound(per_phe_bound_fname);

	cerr << endl;

	do {
		set_pbc(&pbc, -1, box);
		// GMX v 4.0.5 I think this is needed to make system whole
		rm_pbc(&top->idef, ePBC, natoms, box, x, x);

		vector<PheMolecule*> inositol_molecules;
		vector<PheMolecule*> phe_molecules;

		real* phe_com = new real[3];
		real* inositol_com = new real[3];
		real dist;

		// Parse inositol coordinates
		populate_inositol_coordinates(top, index, x, INOSITOL_GROUP_START_IDX, numInositols,
									  NUM_ATOMS_INOSITOL, inositol_molecules);

		// Parse phe coordinates
		populate_phe_coordinates(top, index, x, PROTEIN_GROUP_START_IDX, NUM_ATOMS_PROTEIN, phe_molecules);

		vector<int> phe_stacking(phe_molecules.size());
		vector<int> phe_bound(phe_molecules.size());

		vector<int> inos_stacking(numInositols);
		vector<int> inos_bound(numInositols);

		atom_id* phe_angle_index = new atom_id[3];
		atom_id* inositol_angle_index = new atom_id[3];
		real angle, distance, distance1, distance2;

		// for(int ins_num = 0; ins_num < inositol_molecules.size(); ins_num++) {
        //   inositol_molecules.at(ins_num)->print_info(); 
        //}
 
		//for(int phe_num = 0; phe_num < phe_molecules.size(); phe_num++) {
        //   phe_molecules.at(phe_num)->print_info(); 
        //}

		for(int ins_num = 0; ins_num < inositol_molecules.size(); ins_num++) {
			for(int phe_num = 0; phe_num < phe_molecules.size(); phe_num++) {
				phe_molecules.at(phe_num)->center_of_geometry(phe_com);
				inositol_molecules.at(ins_num)->center_of_geometry(inositol_com);
                // Calculate the planar angle between this inositol molecule and the PHE residue that
                // it is in contact with.

                // Gets the index of the atoms for the phe molecules which are used in the angle calculations
                phe_molecules.at(phe_num)->get_angle_index_group(phe_angle_index);

                // Gets the index of the atoms for the phe molecules which are used
                inositol_molecules.at(ins_num)->get_angle_index_group(inositol_angle_index);

                calc_angle(ePBC, box, x, phe_angle_index, inositol_angle_index, 3, 3,
                           &angle, &distance, &distance1, &distance2);

                // fprintf(sg_angle, "%12g  %12g  %12g\n", t, angle, acos(angle)*180.0 / M_PI);

                double angle_degrees = acos(angle)*180.0 / M_PI;
                bool in_contact = is_in_contact(&pbc, inositol_com, phe_com, CUTOFF, dist);

                if(dist < 0.8) {
                    cerr << "t=" << t << " (" << phe_molecules.at(phe_num)->get_resname()
                    					      << phe_molecules.at(phe_num)->get_resid() << ","
                    						  << inositol_molecules.at(ins_num)->get_resname()
                    						  << inositol_molecules.at(ins_num)->get_resid() << ") "
                    						  << "angle=" << angle_degrees << " dist=" << dist << " ";
                }

				if(in_contact) {
                    if(angle_degrees < 15.0 || (180 - angle_degrees) < 15) {
                        cerr << "STACKED " << endl;
                        phe_stacking[phe_num]++;
                        inos_stacking[ins_num]++;
                    	phe_bound[phe_num]++;
                    	inos_bound[ins_num]++;
                    } else {
                    	cerr << "BOUND " << endl;
                    	phe_bound[phe_num]++;
                    	inos_bound[ins_num]++;
                    }
				}
			}
		}

		f_phe_stacking << t << " ";
		f_inos_stacking << t << " ";
		f_phe_bound << t << " ";
		f_inos_bound << t << " ";

		// Output stacking count for phes
		for(int i = 0; i < phe_stacking.size(); i++) {
			f_phe_stacking << phe_stacking[i] << " ";
		}

		// Output bound count for phes
		for(int i = 0; i < phe_bound.size(); i++) {
			f_phe_bound << phe_bound[i] << " ";
		}

		// Output stacking count for phes
		for(int i = 0; i < inos_stacking.size(); i++) {
			f_inos_stacking << inos_stacking[i] << " ";
		}

		// Output bound count for inositol
		for(int i = 0; i < inos_bound.size(); i++) {
			f_inos_bound << inos_bound[i] << " ";
		}

		f_phe_stacking << endl;
		f_inos_stacking << endl;
		f_phe_bound << endl;
		f_inos_bound << endl;

//		delete_vector(inositol_molecules);
//		delete_vector(phe_molecules);
	} while (read_next_x(status, &t, natoms, x, box));
}
