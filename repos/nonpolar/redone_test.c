#define CPLUSPLUS

using namespace std;

//c++ libraries
#include <map>
#include <string>

int main(int argc,char *argv[]) {
    //program description
    static char *desc[] = {
        "g_inos_np calculates the number of inositol to protein nonpolar contacts for each frame in the trajectory (num_contacts.dat)",
        "g_inos_np also outputs the minimum distance of inositol (ring) to protein for each inositol (min_dists.dat)."
        "This file is useful for debugging purposes."
    };    

    //various declarations for t_pargs
    t_pbc   pbc;
    static char* perResidueContacts = "per_residue_contacts.dat";
    static char* perInositolContacts = "per_inositol_contacts.dat";
    static int NSIDECHAINS = 7;  //the total number of sidechains to read in
    static int NINS = 2;         //total number of inositol groups
    static int NPHE = 0;
    static real CUTOFF=0.45;     /*the np cbeta-inositol cutoff in nm*/


    const char  *leg[4] = { "|d|","d\\sx\\N","d\\sy\\N","d\\sz\\N" };

    //optional args
    static t_pargs pa[] = {
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
        },

        { "-dist",      FALSE, etREAL, {&CUTOFF},
        "the cutoff"}
    };


    //set up for get_index
    int ngrps;
    
    get_index(&top->atoms,ftp2fn(efNDX,NFILE,fnm),ngrps,isize,index,grpname);

    int natoms = read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
    do {
        set_pbc(&pbc, box);
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


    } while (read_next_x(status,&t,natoms,x,box));

    //core algorithm
    map<string, int> per_residue_contacts;
    map<string, int> per_inositol_contacts;
    
    string resName = *(top->atoms.resname[resNum]);
    per_residue_contacts [ resname ] ++;