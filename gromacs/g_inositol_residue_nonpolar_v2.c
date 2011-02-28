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
	#include "physics.h"`
	#include "gstat.h"
	#include "tpxio.h"
}


//#define DEBUG_INOS
//#define DEBUG_COM
#define DEBUG_CONTACT
//#define EDGES
//#define DEBUG_KEY

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

int calculate_com ( t_topology *top, atom_id** index, rvec* x, int groupNum, int groupSize, vector<real*> &inositol_com, vector<int> &table ) {
    int g=0;

    int prevResnum =  top->atoms.atom[index[groupNum][0]].resnr;
    int currentResnum = prevResnum;
    real* com = new real[3];
    com[0]=0;com[1]=0;com[2]=0;
    //cerr<<prevResnum<<" "<<currentResnum<<endl;
    //cerr<<groupNum<<" "<<groupSize<<endl;
    int num_atoms_in_residue = 0;
    int atomIndex=0;
        for(int i=0; (i<groupSize); i++) {
            atomIndex = index[groupNum][i];
            currentResnum = top->atoms.atom[atomIndex].resnr;
            if(currentResnum != prevResnum) {

                for(int d=0;(d<DIM);d++) {
                    com[d] = com[d]/num_atoms_in_residue;
                }

#ifdef DEBUG_COM
                cerr<<"current = " << currentResnum<<" "<<num_atoms_in_residue<<endl;
                cerr<<com[0]*10<<" "<<com[1]*10<<" "<<com[2]*10<<endl;
#endif

                inositol_com.push_back(com);
                
                table.push_back(prevResnum);

                com = new real[3]; com[0]=0.0; com[1]=0.0; com[2]=0.0;
                for(int d=0;(d<DIM);d++) {
                    com[d] += x[atomIndex][d];
                }          

                //reset size of the residue
                num_atoms_in_residue = 1;
                //increment to next group
                g++;

            } else {
                for(int d=0;(d<DIM);d++) {
                    com[d] += x[atomIndex][d];
                }
                num_atoms_in_residue++;
            }
            prevResnum = currentResnum;
        }

        for(int d=0;(d<DIM);d++) {
            com[d] = com[d]/num_atoms_in_residue;
        }

#ifdef DEBUG_COM
//                cerr<<"current = " << currentResnum<<" "<<num_atoms_in_residue<<endl;
//                cerr<<com[0]*10<<" "<<com[1]*10<<" "<<com[2]*10<<endl;
#endif

        inositol_com.push_back(com);
        table.push_back(currentResnum);
        return g+1;
}

void delete_vector(vector<real*> &v) {
    for(int i=0; i<v.size(); i++) {
        delete v[i];
    }
}

// for each residue in an index group
// build an index of atomic counts for each residue
// if there is more than one peptide, then it would be accumulated across all peptides
void build_residue_atom_counts(t_topology* top, atom_id** index, int* isize, int groupNum, map<string,int> &residue_atom_count) {
    map <string, int>::iterator residue_atom_count_it;
    int total_num_atoms = isize[groupNum];

//cout<<total_num_atoms<<endl;

    int phe_index=0;
    for(int i=0; i< total_num_atoms; i++ ) {
        int atomResidueNum = top->atoms.atom[index[groupNum][i]].resnr;  //residue number corresponding 
        string residueName = *( top->atoms.resname[ atomResidueNum ] );
        ostringstream str_out;

        str_out<<residueName<<atomResidueNum;

        string residueNameId = str_out.str();

//cout<<residueNameId<<endl;


        residue_atom_count_it = residue_atom_count.find( residueNameId );
        if(residue_atom_count_it != residue_atom_count.end() ) {
            residue_atom_count[ residueNameId ]++;

//cout<<residueNameId<<" counted " <<endl;

        } else {
            residue_atom_count.insert(pair<string,int>(residueNameId,1));

//cout<<residueNameId<<" inserted " <<endl;

        }
    }
}

void print_map(map<string,int>& m, ofstream& f_residue_table) {
    map<string,int>::iterator iter;
    for(iter=m.begin(); iter!=m.end(); iter++){
        f_residue_table<<(*iter).first<<" "<<(*iter).second<<endl;
    }
}

//initialize a map with the same keys as a reference map
void initialize_maps(map<string, int>& uninit, map<string, int>&ref){
    map<string,int>::iterator ref_iter;
    for(ref_iter=ref.begin(); ref_iter!=ref.end(); ref_iter++) {
        uninit[ref_iter->first]=0;
    }
}

void initialize_vectors(vector<int> &uninit, map<string,int> &ref){
    uninit.resize(ref.size(),0);
}

bool is_in_contact(t_pbc* pbc, rvec p1, rvec p2, real cutoff, real &dist2){
    //function stub
    rvec dx;
    pbc_dx(pbc, p1, p2, dx);
    real dist = norm(dx);
    dist2=dist;
    if(dist < cutoff) {
        return true;
    }

#ifdef DEBUG_CONTACT
    cout<<dist<<" ";
#endif

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
    
    int     ngrps = 2;     /* the number of index groups */
    atom_id **index = new atom_id*[ngrps];
    int     *isize = new int[ngrps];    /* the size of each group */
    char    **grpname = new char*[ngrps]; /* the name of each group */

    atom_id **max;   /* the index for the atom numbers */

    rvec    *com;
    rvec    *com_phe;
    real    *mass;
    FILE    *fp=NULL;

    t_pbc   pbc;
    //static char *f_contact="residue_nonpolar_contact.dat";
    static char* perResidueContacts = "per_residue_contacts.dat";
    static char* residueTable = "table.dat";
    static char* ffInfo = "ff_vs_t.dat";
    static char* perInositolPheContacts = "per_inositol_phe_contacts.dat";
    static char* perInositolContacts = "per_inositol_contacts.dat";
    static char* perInositolPheComDists="per_inositol_phe_com_dists.dat";
    //static int NINS = 2; //total number of inositol groups
    //static int NSIDECHAINS = 7;  //the total number of sidechains to read in
    static int NPHE = 0;
    real CUTOFF=0.45;     /*the np cbeta-inositol cutoff in nm*/
    static int bComputeFirstCOM=0;
    static int bOutputDists=0;

    const char  *leg[4] = { "|d|","d\\sx\\N","d\\sy\\N","d\\sz\\N" };

    static t_pargs pa[] = {
        { "-dist",      FALSE, etREAL, {&CUTOFF},
        "the cutoff"},

        { "-per_residue_contacts", FALSE, etSTR, {&perResidueContacts},
        "description XYZ"},

        { "-per_residue_table", FALSE, etSTR, {&residueTable},
        "description XYZ"},

        { "-per_inositol_contacts", FALSE, etSTR, {&perInositolContacts},
         "description XYZ"},

        {"-per_inositol_phe_contacts", FALSE, etSTR, {&perInositolPheContacts},
        "description XYZ"}, 
        
        { "-first_com", FALSE, etINT, {&bComputeFirstCOM},
        "set to compute the COM of the first group"},

        { "-FF_info", FALSE, etSTR, {&ffInfo},
         "total FF contacts vs t"},

	{ "-com_dists", FALSE, etINT, {&bOutputDists},
	"print com distances between all inositol and phes"},

	{"-com_dist_xvg", FALSE, etSTR, {&perInositolPheComDists},
	"file name to print to"}

    };

    #define NPA asize(pa)
    
    char title[STRLEN];

    /* %%% print the cutoff to stderr %%%*/
    fprintf(stderr, "CUTOFF used is %f", CUTOFF);

    t_filenm fnm[] = {
        { efTPS, NULL, NULL, ffREAD },
        { efTRX, "-f", NULL, ffREAD },
        { efNDX, NULL, NULL, ffOPTRD }
    /* { efXVG, NULL, "dist", ffOPTWR },*/
    };

    #define NFILE asize(fnm)
    CopyRight(stderr,argv[0]);    
    parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_BE_NICE,
                        NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
    int ePBC=-1;
    top=read_top(ftp2fn(efTPS,NFILE,fnm), &ePBC);


    //read_tps_conf(ftp2fn(efTPS,NFILE,fnm), title, top, &xtop, NULL, box, TRUE);
    //sfree(xtop);

    get_index(&top->atoms,ftp2fn(efNDX,NFILE,fnm),ngrps,isize,index,grpname);
    int natoms = read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
    
    const int PROTEIN = 0;  //protein group index
    const int INOSITOL = 1; //inositol group index
    const int NUM_ATOMS_PROTEIN = isize [ PROTEIN ]; //number of atoms in the protein group
    const int NUM_ATOMS_INOSITOL = isize [ INOSITOL ]; //number of atoms in the inositol group

    map <string, int> residue_atom_count;
    map <string, int> per_residue_contacts_snapshot;
    vector <int> phe_table;
    vector <int> inositol_table;


    build_residue_atom_counts(top, index, isize, PROTEIN, residue_atom_count);
    initialize_maps(per_residue_contacts_snapshot, residue_atom_count);
    //print_map(residue_atom_count, residue_atom_count);

    string residueName;
    int residue_id;
    ofstream f_per_residue_contacts(perResidueContacts);
    ofstream f_residue_table(residueTable);
    ofstream f_per_inositol_contacts(perInositolContacts);
    ofstream f_inos_phe_contacts(perInositolPheContacts);
    ofstream f_phe_phe_contacts(ffInfo);
    ofstream f_phe_com_dists(perInositolPheComDists);

    print_map(residue_atom_count,f_residue_table);
    do {
        set_pbc(&pbc, -1, box);
	    /* GMX v 4.0.5 I think this is needed to make system whole*/
		rm_pbc(&top->idef,ePBC,natoms,box,x,x);

        vector<real*> inositol_com;
        vector<real*> phe_com;
        real dist;
        vector <int> per_inositol_contacts_snapshot;

#ifdef EDGES
    cout<<endl<<t<<" ";
#endif


        int num_groups_detected = calculate_com(top, index, x, INOSITOL, NUM_ATOMS_INOSITOL, inositol_com, inositol_table);

#ifdef DEBUG_INOS
	    cout<<"num_groups_detected="<<num_groups_detected<<endl;
#endif            

        if(per_inositol_contacts_snapshot.size() == 0){
            per_inositol_contacts_snapshot.resize(num_groups_detected, 0);
        }
		
		bool bInContact = false;	
        for( int i=0; i<inositol_com.size(); i++ ) {
			float min_dist = 10; 
            for( int atomIndex = 0; atomIndex < NUM_ATOMS_PROTEIN; atomIndex++ ) {
                //residue_id = top->atoms.atom[ index[PROTEIN][atomIndex] ].resnr;
                //residueName = *(top->atoms.resname[ residue_id ]);
				bInContact = is_in_contact(&pbc, x[ index[PROTEIN][atomIndex] ], inositol_com[i], CUTOFF, dist); 
				
				if (dist < min_dist) {
					min_dist = dist;
				}

                if(bInContact) { 
               		cout << "At t=" << t << " inositol " << i << " is in contact at dist=" << dist << endl; 

                    residue_id = top->atoms.atom[ index[PROTEIN][atomIndex] ].resnr;
                    residueName = *(top->atoms.resname[ residue_id ]);

#ifdef EDGES
    cout<<residue_id<<" ";
#endif

                    ostringstream key_ss;
                    key_ss<<residueName<<residue_id;
                    string key_str = key_ss.str();

                    per_residue_contacts_snapshot[key_str]++;                        
                    per_inositol_contacts_snapshot[i]++;

                }
            }

#ifdef DEBUG_CONTACT
			cout << "min_dist=" << min_dist << endl;
#endif
        }

        //output per residue contacts
        f_per_residue_contacts<<t<<" ";
        for(map<string,int>::iterator iter = per_residue_contacts_snapshot.begin(); iter != per_residue_contacts_snapshot.end(); iter++) {

#ifdef DEBUG_KEY
	    	cerr << iter->first<<" ";
        	cerr << iter->second<<" "<<endl;
#endif 
            f_per_residue_contacts<<iter->second<<" ";
        }

        f_per_residue_contacts<<endl;

        //output per inositol contacts        
        f_per_inositol_contacts << t << " ";
        for( int i=0; i<inositol_com.size(); i++) {
            f_per_inositol_contacts << per_inositol_contacts_snapshot[i] << " ";
        }
        f_per_inositol_contacts << endl;

        vector<int>per_inositol_COM_contacts_snapshot;
        per_inositol_COM_contacts_snapshot.resize(num_groups_detected,0);

        if(bComputeFirstCOM) {
            int total_FF_contact=0;
            int num_phe_detected = calculate_com(top, index, x, PROTEIN, NUM_ATOMS_PROTEIN, phe_com, phe_table);
            //output inositol phe contacts
            f_inos_phe_contacts<<t<<" ";
            f_phe_com_dists<<t<<" ";
            for(int i=0; i < inositol_com.size(); i++) {
                for(int j=0; j < phe_com.size(); j++) {
                    if( is_in_contact(&pbc, inositol_com[i], phe_com[j], CUTOFF, dist) ) {
                        per_inositol_COM_contacts_snapshot[i]++;

#ifdef DEBUG_COM
                        cout <<t<<" "<< inositol_table[i] << " ";
                        for(int d=0; d<3;d++){
                            cout<<inositol_com[i][d]<<" ";
                        }

                        cout << phe_table[j] << " ";
                        for(int d=0; d<3; d++){
                            cout<<phe_com[j][d]<<" ";
                        }
                        cout<<dist<<endl;
#endif
                    }
		    if(bOutputDists){
			f_phe_com_dists<<dist<<" ";
		    }
                }
                f_inos_phe_contacts<<per_inositol_COM_contacts_snapshot[i]<<" ";
            }

            if(bOutputDists){
                    f_phe_com_dists<<endl;
            }
            

            f_inos_phe_contacts<<endl;
            

            //output phe phe contacts
            f_phe_phe_contacts<<t<<" ";
            for(int j=0; j < phe_com.size(); j++) {
                for(int k=j+1; k < phe_com.size(); k++) {
                    if( is_in_contact( &pbc, phe_com[j], phe_com[k], CUTOFF, dist ) ) {
                        total_FF_contact++;
#ifdef DEBUG_COM
//                         cout << phe_table[j] << " ";
//                         for(int d=0; d<3;d++){
//                             cout<<phe_com[j][d]<<" ";
//                         }
// 
//                         cout << phe_table[k] << " ";
//                         for(int d=0; d<3; d++){
//                             cout<<phe_com[k][d]<<" ";
//                         }
//                         cout<<dist<<" ";
#endif 

                    }
                }
            }
            f_phe_phe_contacts<<total_FF_contact<<endl;
        }

        //do some clean up
        delete_vector(inositol_com);
        initialize_maps(per_residue_contacts_snapshot, residue_atom_count); //reinitialize
    } while(read_next_x(status, &t, natoms, x, box));

}

