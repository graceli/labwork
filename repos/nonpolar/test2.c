#define CPLUSPLUS

using namespace std;

//c++ libraries
#include <map>
#include <string>
#include <iostream>
#include <vector>

//gromacs c libraries
#include <typedefs.h>
#include "smalloc.h"
#include "macros.h"
#include "math.h"
#include "xvgr.h"
#include "copyrite.h"
#include "statutil.h"
#include "vec.h"
#include "index.h"
#include "pbc.h"
#include "futil.h"
#include "gstat.h"

//takes a topology structure
//and outputs atom index, atom name, and the residue the atom is in
void dump_top ( t_topology *top ) {
    int totalNumAtoms = top->atoms.nr;
    int totalNumResidues = top->atoms.nres;
    cout<<"Number of atoms = " << totalNumAtoms <<endl;
    cout<<"Number of residues = " << totalNumResidues <<endl;
    for(int i = 0; i < totalNumAtoms; i++){
        int atomResidueNum = top->atoms.atom[i].resnr;  //residue number corresponding 
        int atomIndex = i+1;
        cout<<atomIndex<<" "<<*(top->atoms.atomname[ atomIndex-1 ])<<" "<<*(top->atoms.resname[ atomResidueNum ])<<endl;
    }
}

void dump_index ( t_topology *top, atom_id ** index, int* isize, int ngrps ) {
    for(int g=0; g < ngrps; g++){
        for(int gIndex=0; gIndex < isize[g]; gIndex ++){    
            cout<<g<< " " << index[g][gIndex]<< " "<<*(top->atoms.atomname[ gIndex ])<<" "<<*(top->atoms.resname[ top->atoms.atom[gIndex].resnr ])<<endl;
        }
        cout<<endl;
    }
}

int calculate_com ( t_topology *top, atom_id** index, rvec* x, int groupNum, int groupSize, vector<real*> &inositol_com ) {
    int g=0;

    int prevResnum =  top->atoms.atom[index[groupNum][0]].resnr;
    int currentResnum = prevResnum;
    real* com = new real[3];
    com[0]=0;com[1]=0;com[2]=0;
    cerr<<prevResnum<<" "<<currentResnum<<endl;
    cerr<<groupNum<<" "<<groupSize<<endl;

        for(int i=0;(i<groupSize);i++) {
            int atomIndex = index[groupNum][i];
            currentResnum = top->atoms.atom[atomIndex].resnr;
            if(currentResnum != prevResnum) {
                cerr<<"current = " << currentResnum<<endl;
                cerr<<com[0]<<" "<<com[1]<<" "<<com[2]<<endl;
                inositol_com.push_back(com);
                //com[0]=0.0; com[1]=0.0; com[2]=0.0;
                com = new real[3];
                //increment to next group
                g++;
            }else{
                for(int d=0;(d<DIM);d++) {
                    com[d] += x[atomIndex][d]/groupSize;
                }
            }
            prevResnum = currentResnum;
        }

        inositol_com.push_back(com);
        return g+1;
}


int main(int argc,char *argv[]) {
        static char *desc[] = {
        "g_inos_np calculates the number of inositol to protein nonpolar contacts for each frame in the trajectory (num_contacts.dat)",
        "g_inos_np also outputs the minimum distance of inositol (ring) to protein for each inositol (min_dists.dat).",
        "This file is useful for debugging purposes."
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
    static char* perInositolContacts = "per_inositol_contacts.dat";

    static int NINS = 2; //total number of inositol groups
    static int NSIDECHAINS = 7;  //the total number of sidechains to read in
    static int NPHE = 0;
    static real CUTOFF=0.45;     /*the np cbeta-inositol cutoff in nm*/

    const char  *leg[4] = { "|d|","d\\sx\\N","d\\sy\\N","d\\sz\\N" };

    static t_pargs pa[] = {
        { "-dist",      FALSE, etREAL, {&CUTOFF},
        "the cutoff"},

        { "-per_residue_contacts", FALSE, etSTR, {&perResidueContacts},
        "Compute a total count of inositol-residue sidechain nonpolar contacts for each residue group inputted"},

        { "-per_inositol_contacts", FALSE, etSTR, {&perInositolContacts},
        "total number of NP contacts for each inositol"},

        { "-num_sidechains", FALSE, etINT, {&NSIDECHAINS},
        "total number of sidechains in the system (calculated as num_peptides*7)"
        },
    
        { "-num_inositols", FALSE, etINT, {&NINS},
        "total number of inositols in the system"
        },

        { "-phe", FALSE, etINT, {&NPHE},
        "number of phenylalanines in a system; if more than 0, will only output contacts to phenyalanine"
        }
    };

    #define NPA asize(pa)
    
    char title[STRLEN];

    /* %%% print the cutoff to stderr %%%*/
    fprintf(stderr, "CUTOFF used is %f", CUTOFF);

    t_filenm fnm[] = {
        { efTRX, "-f", NULL, ffREAD },
        { efTPS, NULL, NULL, ffREAD },
        { efNDX, NULL, NULL, ffOPTRD }
    /* { efXVG, NULL, "dist", ffOPTWR },*/
    };

    #define NFILE asize(fnm)
    CopyRight(stderr,argv[0]);    
    parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_BE_NICE,
                        NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
    top=read_top(ftp2fn(efTPS,NFILE,fnm));
    get_index(&top->atoms,ftp2fn(efNDX,NFILE,fnm),ngrps,isize,index,grpname);

    int natoms = read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
    // dump_index( top, index, isize, ngrps);
    //core algorithm
    map<string, int> per_residue_contacts;
    map<string, int> per_inositol_contacts;
    
    const int PROTEIN = 0;  //protein group index
    const int INOSITOL = 1; //inositol group index
    const int NUM_ATOMS_PROTEIN = isize[PROTEIN]; //number of atoms in the protein group
    const int NUM_ATOMS_INOSITOL = isize[INOSITOL]; //number of atoms in the inositol group


    do {
        set_pbc(&pbc, box);

        vector<real*> inositol_com;
        int num_groups_detected = calculate_com(top, index, x, INOSITOL, NUM_ATOMS_INOSITOL, inositol_com);
        cout<<"number of groups found = " << num_groups_detected<<endl;

        for( int i=0; i<inositol_com.size(); i++){
            cout<<inositol_com[i][0]<<" "<<inositol_com[i][1]<<" "<<inositol_com[i][2]<<endl;
        }
//        delete_vector(inositol_com);
    } while(read_next_x(status, &t, natoms, x, box));
}

