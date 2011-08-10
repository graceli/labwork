#define CPLUSPLUS

using namespace std;
#include <string>

#ifdef HAVE_CONFIG_H
    #include <config.h>
#endif

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


/*
*    This program calculates the time evolution of the number of inositol-peptide nonpolar contacts
*
*    March 23 2009
*    This is the modified version to output the number of inositol NPcontacts for each residue individually 
*    as well as as a total sum over all frames in the xtc
*/

//#define DEBUG_CONTACT

int main(int argc,char *argv[])
{
    static char *desc[] = {
        "g_inos_np calculates the number of inositol to protein nonpolar contacts for each frame in the trajectory (num_contacts.dat)",
        "g_inos_np also outputs the minimum distance of inositol (ring) to protein for each inositol (min_dists.dat)."
        "This file is useful for debugging purposes."
    };
    
    t_topology *top=NULL;
    rvec *xtop; 
    real t,cut2,dist2;
    rvec *x=NULL,*v=NULL,dx;
    matrix box;
    int status;
    int natoms;
    
    int g=0,d=0,i=0,j=0,res=0,teller=0;
    atom_id aid;
    
    int     ngrps;     /* the number of index groups */
    atom_id **index,max;   /* the index for the atom numbers */

    int     *isize;    /* the size of each group */
    char    **grpname; /* the name of each group */
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
    
    ngrps = NINS + NSIDECHAINS; 

/*
    int     ngrps;     
    atom_id **index,max;   
    int     *isize;    
    char    **grpname; 
    rvec    *com;
    rvec    *com_phe;
    real    *mass;
    FILE    *fp=NULL;
*/

    com = new rvec[NINS];
    
    grpname = new char*[ngrps];
    index = new atom_id*[ngrps];
    isize = new int[ngrps];
    mass = new real[NINS];

    //allocate space for various things
//     snew(com, NINS);	 	// array to hold each of the com of inositol 
//     snew(com_phe, NPHE); 	// each of the com of PHE
//     snew(grpname, ngrps);
//     snew(index, ngrps);
//     snew(isize, ngrps);
//     snew(mass, ngrps);


    get_index(&top->atoms,ftp2fn(efNDX,NFILE,fnm),ngrps,isize,index,grpname);
    
    /* calculate mass for each inositol group */
    /* obsolete, calculate center of geometry instead as scyllo is symmetric, true for chiro though? */

//    snew(mass, NINS);int 
    for(g=0; (g<NINS); g++) {
        mass[g] = 0;
        for(i=0; (i<isize[g]); i++) {

            if (index[g][i] >= top->atoms.nr)
                gmx_fatal(FARGS,"Atom number %d, item %d of group %d, is larger than number of atoms in the topology (%d)\n", 
	    
	    index[g][i]+1, i+1, g+1, top->atoms.nr+1);
            mass[g] += top->atoms.atom[index[g][i]].m;
        }
    }
    

    FILE* fp_per_residue_contacts = fopen(perResidueContacts, "w");
    FILE* fp_per_inositol_contacts = fopen(perInositolContacts, "w");
    FILE* fp_phe_contacts = fopen("phetest.dat", "w");

    //calculate the minimum distances between center of masses of inositol and atoms of sidechain groups
    int s=0; 
    int** inositol_contact_count_matrix;
 
    inositol_contact_count_matrix = new int*[NINS];
//    snew(inositol_contact_count_matrix, NINS);

    for(g = 0; (g < NINS); g++) { // for each inositol group
        inositol_contact_count_matrix[g] = new int[NSIDECHAINS+NINS];
        //snew(inositol_contact_count_matrix[g], NSIDECHAINS+NINS);
        for(s = 0; (s < NSIDECHAINS+NINS); s++) { 
             inositol_contact_count_matrix[g][s]=0;
        }
    }

    for(s=NINS; s<NSIDECHAINS+NINS; s++){
	printf("size of sc # %d = %d\n", s, isize[s]);
    }

    real* mindist_phe_inos_com;
    //snew(mindist_phe_inos_com, NINS); //10 nm
    mindist_phe_inos_com = new real[NINS];
    for(i=0; i<NINS; i++){
        mindist_phe_inos_com[i]=10.0;
    }
    int* phe_number;
    //snew(phe_number, NINS);
    phe_number = new int[NINS];

    do {
        /* initialisation for correct distance calculations */
        set_pbc(&pbc,box);

        /* make molecules whole again */
//        rm_pbc(&top->idef,natoms,box,x,x);
    
        /* calculate center of mass for each inositol (the first two groups) */
        for(g=0;(g<NINS);g++) {

            for(d=0;(d<DIM);d++) {
    	       com[g][d] = 0;

    	       for(i=0;(i<isize[g]);i++) {
    	           com[g][d] += x[index[g][i]][d]*top->atoms.atom[index[g][i]].m;
    	       }

    	       com[g][d] /= mass[g];
            }
        }

        //the indices in index[][] is shifted down by 1 from the index.ndx
        real dist=0.0;
        real min_dists[2] = {3.0,3.0};
        real min_com_dists[3] = {3.0,3.0};

        int* per_inositol_contacts;
        //snew(per_inositol_contacts,NINS);
        per_inositol_contacts = new int[NINS];

        int* per_residue_contacts;
        //snew(per_residue_contacts, NSIDECHAINS);
        per_residue_contacts = new int[NSIDECHAINS];
    
        if(NPHE == 0) {
         //for each side chain groups count the number of nonpolar contacts by inositol to each 
          //sidechain AND average over all the chains in the system
          for(s=NINS; (s < NSIDECHAINS+NINS); s++) { // for each side chain group
    
                int sidechain_index = (s-2)%7+2;  //this formula wraps around indices to separate chains
    
                for(i = 0; (i < isize[s]); i++) { //for each atom in the sidechain group
                    for(g = 0; (g < NINS); g++) { // for each inositol group
                        pbc_dx(&pbc, com[g], x[ index[s][i] ], dx);
                        dist = norm(dx);
    
                        //if we found a distance that is less than CUTOFF
                        //count it as a contact between inositol g and sidechain s    
                        if( dist < CUTOFF ) {
    
#ifdef DEBUG_CONTACT
                            printf("sidechain_index = %i; s=%d; INS # %d\n", sidechain_index, s, g); 
                            printf("dist %f met the CUTOFF=%f\n", dist, CUTOFF);
#endif

                            inositol_contact_count_matrix[g][sidechain_index]++;
                            per_inositol_contacts[g]++;                 //sum across "columns" of above matrix
                            per_residue_contacts[sidechain_index]++;    //sum across "rows"
                        }
                    }
                }
            }
    
            //output per inositol contacts
            for(s=0; s<NINS; s++){
                fprintf(fp_per_inositol_contacts, "%d\n", per_inositol_contacts[0]);	    
            }
            fprintf(fp_per_inositol_contacts, "\n");
    
            //output per residue contacts
            for(s=NINS; (s<NSIDECHAINS+NINS); s++){
                fprintf(fp_per_residue_contacts, "%d ", per_residue_contacts[s]);
            }
            fprintf(fp_per_residue_contacts, "\n");

        } else {

            //if PHE flag is on, then compute only the center of geometry - cog contacts between benzenes and inositol
            //compute the center of geometry for each benzene ring 
            /* calculate center of mass for each inositol (the first two groups) */
            com_phe = new rvec[NPHE];
            for(g=NINS; (g < NINS + NPHE); g++) {
                for(d=0;(d<DIM);d++) {
                    com_phe[g][d] = 0;
    
                    for(i=0;(i<isize[g]);i++) {
                        com_phe[g][d] += x[index[g][i]][d]*top->atoms.atom[index[g][i]].m;
                    }

                    com_phe[g][d] /= mass[g];
                }
            }

            //for each inositol, compute the minimum distance between COM of inositol and benzene rings
            for(g=0; (g<NINS); g++) {
                for(s=NINS; (s<NINS+NPHE); s++) {
                    pbc_dx(&pbc, com[g], com_phe[s], dx);
                    double dist = norm(dx);
                    if(dist<mindist_phe_inos_com[g]) {  
                        mindist_phe_inos_com[g]=dist;
                        //phe_number[g]=s;
                    }
                }
            }

            for(g=0; (g<NINS); g++) {
                fprintf(fp_phe_contacts, "%f ", mindist_phe_inos_com[g]);
            }
            
            for(g=0; (g<NINS); g++) {
                if(mindist_phe_inos_com[g] < CUTOFF) {
                    fprintf(fp_phe_contacts, "%d ", 1);
                } else {
                    fprintf(fp_phe_contacts, "%d ", 0);
                }
            }

        }

        teller++;

    } while(read_next_x(status,&t,natoms,x,box));

    close_trj(status);
    //thanx(stderr);

    return 0;
}

