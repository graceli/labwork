#define CPLUSPLUS

#define DEBUG_PBC
#define DEBUG_RESID



//c++ libraries
#include <map>
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>

using namespace std;

//gromacs c libraries
extern "C" {
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
#include "tpxio.h"
}

/* 
This program is renamed from g_inositol_residue_nonpolar_v2.c -- which was written in September 2009 to 
look at peptide-peptide nonpolar contacts -- in response to Liu et al. paper on Trehalose and KLVFFAE
*/

/*
Be careful of number of residues use to determine the peptide id (that is, the chain number)
*/
bool is_in_contact(t_pbc* pbc, rvec p1, rvec p2, real cutoff, real &dist2){
    rvec dx;

#ifdef DEBUG_PBC
    //cout<<"pbc->ePBCDX "<<pbc->ePBCDX<<endl;
#endif

    pbc_dx(pbc, p1, p2, dx);
    real dist = norm(dx);
    dist2=dist;
    if(dist < cutoff) {
        return true;
    }

#ifdef DEBUG_CONTACT
    cout<<dist<<endl;
#endif

    return false;
}

void COMMENT(ofstream &out, const string &comment){
	out<<comment<<endl;
}

int main(int argc,char *argv[]) {
    static char *desc[] = {
        " Computes the total number of interpeptide nonpolar contacts "
    };
    
    t_topology *top=NULL;
    rvec *xtop; 
    real t,cut2,dist2;
    rvec *x=NULL,*v=NULL,dx;
    matrix box;
    int status;
    int g=0,d=0,i=0,j=0,res=0,teller=0;
    atom_id aid;
    
    int     ngrps = 1;     /* the number of index groups */
    atom_id **index = new atom_id*[ngrps];
    int     *isize = new int[ngrps];    /* the size of each group */
    char    **grpname = new char*[ngrps]; /* the name of each group */

    atom_id **max;   /* the index for the atom numbers */

    rvec    *com;
    rvec    *com_phe;
    real    *mass;
    FILE    *fp=NULL;

    t_pbc   pbc;
    int     ePBC=-1;  /* guess PBC from box*/
    cout<< "initialized ePBC ="<<ePBC<<endl;


    static char* contacts = "pp_nonpolar_vs_t.xvg";
    static char* prefix = "g_pp_nonpolar_";
    real CUTOFF=0.6;     /*the np cbeta-inositol cutoff in nm*/

    const char  *leg[4] = { "|d|","d\\sx\\N","d\\sy\\N","d\\sz\\N" };

    static t_pargs pa[] = {
	{ "-deffnm",      FALSE, etSTR, {&prefix},
	"name all files with this prefix"},
        { "-cutoff",      FALSE, etREAL, {&CUTOFF},
        "input: the cutoff"},
        { "-contacts", FALSE, etSTR, {&contacts},
        "output: file containing the time vs the total number of atomic nonpolar"}
    };

    #define NPA asize(pa)
    
    char title[STRLEN];

    /* %%% print the cutoff to stderr %%%*/
    fprintf(stderr, "CUTOFF used is %f", CUTOFF);

    t_filenm fnm[] = {
        { efTPX, NULL, NULL, ffREAD },
        { efTRX, "-f", NULL, ffREAD },
        { efNDX, NULL, NULL, ffOPTRD }
    /* { efXVG, NULL, "dist", ffOPTWR },*/
    };

    #define NFILE asize(fnm)
    CopyRight(stderr,argv[0]);    
    parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_BE_NICE,
                        NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
    top=read_top(ftp2fn(efTPX,NFILE,fnm), &ePBC);


#ifdef DEBUG_PBC	
    cout<<"ePBC = " << ePBC<<endl;
#endif

    get_index(&top->atoms,ftp2fn(efNDX,NFILE,fnm),ngrps,isize,index,grpname);
    int natoms = read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);

#ifdef DEBUG_PBC
    cout<<box[XX][XX]<<" "<<box[YY][YY]<<" "<<box[ZZ][ZZ]<<endl;
#endif

 
    const int PROTEIN = 0;  //protein group index
    const int NUM_ATOMS_PROTEIN = isize [ PROTEIN ]; //number of atoms in the protein group

    string contacts_full(prefix);  contacts_full.append(contacts);
    ofstream f_pep_nonpolar(contacts_full.c_str());
    
    COMMENT(f_pep_nonpolar, "#time inter-peptide-counts");

    do {
        set_pbc(&pbc, 0, box);
	/* make molecules whole*/
	rm_pbc(&top->idef, 0, natoms, box, x, x);

        real dist;
        int num_in_contact=0;

        f_pep_nonpolar<<t<<" ";

	/* loop over each atom in the index file group*/
        for( int i=0; i<NUM_ATOMS_PROTEIN; i++ ) {
            for( int j=0; j < NUM_ATOMS_PROTEIN; j++ ) {
		int atomi = index[PROTEIN][i];
		int atomj = index[PROTEIN][j];
                int residue_id_i = top->atoms.atom[ index[PROTEIN][i] ].resnr;
                int residue_id_j = top->atoms.atom[ index[PROTEIN][j] ].resnr;

                int peptide_i=ceil((double)residue_id_i/9);
                int peptide_j=ceil((double)residue_id_j/9);

#ifdef DEBUG_RESID
		cout<<"INTERACTION "<< atomi << " "<<atomj<<" "<<residue_id_i<<" "<<residue_id_j<<" "<<peptide_i << " " << peptide_j << endl; 
#endif
                string residueName_i = *(top->atoms.resname[ residue_id_i ]);
                string residueName_j = *(top->atoms.resname[ residue_id_j ]);

                if( i!=j && peptide_j >  peptide_i && is_in_contact(&pbc, x[ atomi ], x[ atomj ], CUTOFF, dist) ) {
#ifdef DEBUG_RESID
  		    	cout<<"IN CONTACT "<<atomi<< " "<<atomj<< " " << residue_id_i<<" "<<residue_id_j<<" "<<peptide_i << " " << peptide_j << " is in contact dist="<<dist<<endl;
#endif
               		num_in_contact++;
            	}
	    }
        }

        f_pep_nonpolar<<num_in_contact<<endl;

    } while(read_next_x(status, &t, natoms, x, box));

}

