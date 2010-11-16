#define CPLUSPLUS

//#include "HB.h"
//#include "AminoAcid.h"
#include <fstream>
#include "peptide.h"

using namespace std;

static char *SRCID_template_c = "$Id: template.c,v 1.4 2001/07/23 15:28:29 lindahl Exp $";

extern "C" {
#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"
#include "index.h"
#include "pbc.h"
}

//#define DEBUG_protein_index
//#define DEBUG_INO
//#define DEBUG_GLU
//#define DEBUG_SHOWFRAME

int parseBackboneNH(atom_id** index, AminoAcid* aa, int g, int i){
        int N = i;
        int H = i+1;
	aa->setBackboneGroup("NH", index[g][H], index[g][N],"donor");
        return H+1;
}

int parseSidechain(atom_id** index, AminoAcid* aa, const string& resName, int g, int i){	
        if(resName == "LYSH"){
                for(int ep = i+1; ep <=i+3; ep++){
                    aa->setSidechainGroup("NH", index[g][ep], index[g][i], "donor");
                }
                return i+4;
        }else if(resName == "GLU"){
                //store sidechain specific to GLU
                for(int en = i+1; en <= i+2; en++){
                    aa->setSidechainGroup("CO", index[g][i], index[g][en], "acceptor");

                }
                return i+3;
        }else{
            //  do nothing
                return i;
        }
}

int parseBackboneCO(atom_id** index, AminoAcid* aa, int g, int i){
        int C = i;
        int O = i+1;
        //store CO
        aa->setBackboneGroup("CO", index[g][C], index[g][O], "acceptor"); 
        return O+1;
}

int parseAA(atom_id** index, AminoAcid* aa, const string& resName, int g,  int i){
        int firstSideAtom = parseBackboneNH(index, aa, g, i);

        int C = parseSidechain(index, aa, resName, g, firstSideAtom);

        int nextAtomPos = parseBackboneCO(index, aa, g, C);

        return nextAtomPos-1;
}

void debug_index(atom_id** index, int g, int gSize){
	cout<<"group="<<g<<endl;
	cout<<"group size="<<gSize<<endl;

	for(int i=0; i<gSize; i++){
		cout<<"index["<<g<<"]["<<i<<"] = "<<index[g][i]<<endl;
	}
}

void parse_protein_index(int totalNumChains, atom_id** index, int* isize, t_topology *top, vector<peptide*>&protein, const string& sequence){
	//parse the index for chain groups
	for(int g = 0; g < totalNumChains; g++){

#ifdef DEBUG_protein_index	
		debug_index(index, g, isize[g]);
		cout<<"DEBUG: group # = "<<g<<endl;
#endif

		peptide* aChain = new peptide(sequence);
		int groupSize = isize[g];

		for(int i=0; i < groupSize; i++){
                    //string atomName = *(top.atoms.atomname[atom_index]);
                    //double atomCharge = top.atoms.atom[atom_index].q;
                    int atomIndex = index[g][i];
                    int resNum = top->atoms.atom[atomIndex].resnr;	
                    string resName = *(top->atoms.resname[resNum]);

#ifdef DEBUG_protein_index
                    cout<<"DEBUG i="<<i<<" resName="<<resName<<endl;
                    cout<<"DEBUG i="<<i<<" atomIndex="<<atomIndex<<endl;
                    cout<<"DEBUG i="<<i<<" resNum="<<resNum<<endl;
#endif

                    AminoAcid* aa = new AminoAcid(resName,resNum+1);
                    i=parseAA(index, aa, resName, g, i);
                    aChain->attachAminoAcid(aa);

                    continue;
                }

		protein.push_back(aChain);
	}
}

void parse_inositol_index(int groupStart, int groupEnd, atom_id** index, int* isize,t_topology *top, vector<AminoAcid*> &inoMolecules){

        for(int g=groupStart; g<groupEnd; g++){
		int groupSize = isize[g];

		//represent an inositol with the AminoAcid type (note this is bad OO naming)
		AminoAcid* ino = new AminoAcid("INS");
		for(int i=0; i< groupSize; i++){
                    int atomIndex = index[g][i];
                    int resNum = top->atoms.atom[atomIndex].resnr;	
                    string resName = *(top->atoms.resname[resNum]);
                    int O=atomIndex;
                    int HA=atomIndex+1;
#ifdef DEBUG_INO
                    cout<<"DEBUG resname "<<resName<<endl;
                    cout<<"DEBUG O="<<O<<endl;
                    cout<<"DEBUG HA="<<HA<<endl;
#endif
                    ino->setSidechainGroup(resName, HA, O, "acceptor");
                    ino->setSidechainGroup(resName, HA, O, "donor");  
                    i+=1;
		}

		inoMolecules.push_back(ino);
	}
}

void peptide_inositol(t_pbc* pbc, rvec* x, vector<peptide*> &protein, vector<AminoAcid*> &inoMolecules, int totalNumChains, int totalNumInos, vector<int> &peptideBbToInosTotalHB, vector<int> &peptideScToInosTotalHB, vector<int> &inositolBack, vector<int> &inositolSideLys, vector<int> &inositolSideGlu,  vector<int> &backboneCounts, vector<int> &sidechainCounts) {
            
    for( int pi=0; pi < totalNumChains; pi++ ) {
        peptide* chaini = protein.at(pi);
        for( int ins=0; ins < totalNumInos; ins++ ) {
            AminoAcid* aInositol = inoMolecules.at(ins);
            int backbone = chaini->computeHBtoBackbone(pbc, x, aInositol, backboneCounts);
            chaini->computeHBtoSidechain(pbc, x, aInositol, sidechainCounts);
            int sidechainLys = chaini->computeHBtoSidechainAA(pbc, x, aInositol, 0);
            int sidechainGlu = chaini->computeHBtoSidechainAA(pbc, x, aInositol, 6);
            int sidechain = sidechainLys + sidechainGlu;
        
            inositolBack[ins] += backbone;
            inositolSideLys[ins] += sidechainLys;
            inositolSideGlu[ins] += sidechainGlu;
        
            peptideBbToInosTotalHB[pi] += backbone;
            peptideScToInosTotalHB[pi] += sidechain;        
        }
    }
}

int peptide_peptide(t_pbc* pbc, rvec* x, vector<peptide*> &protein, int totalNumChains, vector<int> &backboneCounts, vector<int> &sidechainCounts){
    int total=0;
    for( int pi=0; pi < totalNumChains; pi++ ) {
        peptide* chaini = protein.at(pi);
        for( int pj=pi+1; pj < totalNumChains; pj++ ) {
                peptide* chainj = protein.at(pj);

                total += chaini->computeHBtoBackboneChain(pbc,x, chainj, backboneCounts); //chaini backbone is hydrogen bonded to chainj
                total += chaini->computeHBtoSidechainChain(pbc,x, chainj, sidechainCounts); //chaini sidechain is hydrogen bonded to chainj                
        }
    }
    return total;
}

double peptide_peptide_intra(t_pbc* pbc, rvec* x, vector<peptide*> &protein, int totalNumChains, vector<int> &backboneCounts, vector<int> &sidechainCounts){
    int total=0;
    for( int pi=0; pi < totalNumChains; pi++ ) {
        peptide* chaini = protein.at(pi);
        for( int pj=0; pj < totalNumChains; pj++ ) {
            if(pi == pj){
                peptide* chainj = protein.at(pj);

                total += chaini->computeHBtoBackboneChain(pbc,x, chainj, backboneCounts); //chaini backbone is hydrogen bonded to chainj
                total += chaini->computeHBtoSidechainChain(pbc,x, chainj, sidechainCounts); //chaini sidechain is hydrogen bonded to chainj                
            }
        }
    }
    return total/2.0;
}


int main(int argc,char *argv[]) {

/*
    Example usage:
    g_parse_index -f traj_nosol.xtc -s traj_nosol.gro -n index_with_only_polar_parts.ndx [options ...]

    Notes:
    -f I've only ever passed in the nosol xtc because its faster to process;
       No waters are ever analyzed anyways 

    -s Gros worked for gromacs-3.3.1, not sure if they will with gromacs 4

    -n must contain a correct selection of backbone and lysine and glutamate sidechains.  If sidechains don't exist in selection, it will still work fine 
    for backbone.

    [options ..]  these are for setting files names for the output files. if
                  not specifically set, then the default hard coded file names
                  will be used
*/

	static char *desc[] = {
		"g_parse_index",
		"is a program to compute side chain polar and nonpolar contacts",
        "To use, need a index file selecting the polar groups of protein and OH groups of inositol",
        "Need a gro file or tpr file containing the system with no solvent (water)",
        "When running, select the protein first and each molecule of inositol second"
	};

	//set various variables to initialize the command line options
	static int totalNumChains = 1;
	static int totalNumInos = 2;
	static char* prefix="g_parse_index_";
	static char* countsName = "residue_counts.dat";
	static char* perPeptideBackboneName = "pep_bb.dat";
	static char* perPeptideSidechainName = "pep_side.dat";

        static char* perInositolBBName = "inos_bb.dat";
        static char* perInositolScLysName = "inos_lys.dat";
        static char* perInositolScGluName = "inos_glu.dat";
	static char* perInositolTotalName = "inos_total.dat";

     //   static char* perInositolBindMode = "per_inositol_bind_mode.dat";
        
        static char* perResidueBackboneName = "res_bb.dat";
        static char* perResidueSidechainName = "res_side.dat";
        static char* peptideSequence = "KLVFFAE";

        static char* pep2pepName = "p2p_vs_t.dat";

	t_pargs pa[] = {
	    {"-deffnm", FALSE, etSTR, {&prefix},
	    "defines a prefix for all the analysis files"
	    },

            { "-residue_counts", FALSE, etSTR, {&countsName},
            "Hydrogen bonding contact by backbone (first row) and residue (second row)"
            },
    
            { "-per_peptide_backbone", FALSE, etSTR, {&perPeptideBackboneName},
            "Plots total backbone hydrogen bonding contact for each peptide vs time"
            },
    
            { "-per_peptide_sidechain", FALSE, etSTR, {&perPeptideSidechainName},
                "Plots total sidechain hbonding contact for each peptide vs. time"
            },
    
            { "-per_inositol_bb", FALSE, etSTR, {&perInositolBBName},
            "total bb, sc-lys, sc-glu polar contacts made to sidechain for each inositol"
            },
    
            { "-per_inositol_sc_lys", FALSE, etSTR, {&perInositolScLysName},
            "total inositol-sidechain(lys) hb vs time"
            },

            { "-per_inositol_sc_glu", FALSE, etSTR, {&perInositolScGluName},
            "total inositol-sidechain(glu) hb vs time"
            },
   
	    { "-per_inositol_total", FALSE, etSTR, {&perInositolTotalName},
	    "total inositol-polar-protein hb vs time"
	    },
 
            { "-per_residue_backbone", FALSE, etSTR, {&perResidueBackboneName},
            "per residue binding of inositol to backbone"
            },
    
            { "-per_residue_sidechain", FALSE, etSTR, {&perResidueSidechainName},
            "per residue binding of inositol to sidechain"
            },
	 
            { "-p2p", FALSE, etSTR, {&pep2pepName},
                "total number of peptide to peptide hydrogen bonds vs. t"
            },
    
            { "-num_peptides", FALSE, etINT, {&totalNumChains},
                "Specify the total number of peptide chains in system"
            },
    
            { "-num_inositol", FALSE, etINT, {&totalNumInos},
            "Specify the total number of inositol"
            },

            { "-sequence", FALSE, etSTR, {&peptideSequence},
            "Residues of the protein analyzing"
            }
	};
  
	t_filenm fnm[] = {
            { efTPS,  NULL,  NULL, ffREAD },   /* this is for the topology */
            { efTRX, "-f", NULL, ffREAD },      /* and this for the trajectory */
            { efNDX, NULL, NULL, ffOPTRD }
	};
	#define NFILE asize(fnm)

	t_topology top;
	t_trxframe fr;
	rvec       *xtop;
	matrix     box;
	t_pbc 	   pbc;

	int        status;
	int        flags = TRX_READ_X;
	char       title[STRLEN];
	
	/* This is the routine responsible for adding default options,
   	 * calling the X/motif interface, etc. */
	parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,
				    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

	//read topology of system
	int ePBC=-1;
	read_tps_conf(ftp2fn(efTPS,NFILE,fnm), title, &top, &ePBC, &xtop, NULL, box, TRUE);
	sfree(xtop);

	int ngrps = totalNumChains + totalNumInos;
	int g=0, i=0;

	char** grpname = new char*[ngrps];
	atom_id** index = new atom_id*[ngrps];
	int* isize = new int[ngrps];

	get_index(&top.atoms, ftp2fn(efNDX,NFILE,fnm), ngrps, isize, index, grpname);

    	vector<peptide*> protein;
        vector<AminoAcid*> inoMolecules;
        string sequence(peptideSequence);
 
        parse_protein_index(totalNumChains, index, isize, &top, protein, sequence);
        parse_inositol_index(totalNumChains, totalNumChains+totalNumInos, index, isize, &top, inoMolecules);

#ifdef DEBUG
        cout<<"DEBUG number of peptides read: "<<protein.size()<<endl;
        cout<<"DEBUG inositol group starts at: "<<totalNumChains<<endl;
#endif

	real t;
	rvec *x=NULL,dx;
	int natoms = read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);

#ifdef DEBUG
	cout<<"DEBUG Number of atoms " << natoms<<endl;
    	cout<<"DEBUG Total number of inositol read " << totalNumInos<<endl;
#endif

        enum {K, L, V, F1, F2, A, E};
	
	string perInositolBBNameNew = string(prefix)+string(perInositolBBName);
	string perInositolScLysNameNew = string(prefix)+string(perInositolScLysName);
	string perInositolScGluNameNew = string(prefix)+string(perInositolScGluName);
	string perInositolTotal = string(prefix)+string(perInositolTotalName);

	string perPeptideBackboneNameNew = string(prefix) + string(perPeptideBackboneName);
	string perPeptideSidechainNameNew = string(prefix) + string(perPeptideSidechainName);
	string perResidueBackboneNameNew = string(prefix) + string(perResidueBackboneName);
	string perResidueSidechainNameNew = string(prefix) + string(perResidueSidechainName);
	string countsNameNew = string(prefix) + string(countsName);
	string pep2pepNameNew = string(prefix) + string(pep2pepName);
		
        ofstream PERINOSITOL_BB(perInositolBBNameNew.c_str());
        ofstream PERINOSITOL_SC_LYS(perInositolScLysNameNew.c_str());
        ofstream PERINOSITOL_SC_GLU(perInositolScGluNameNew.c_str());
	//added Jan 25 2010
	ofstream PERINOSITOL_TOT(perInositolTotal.c_str());

        ofstream PEPBB(perPeptideBackboneNameNew.c_str());
        ofstream PEPSC(perPeptideSidechainNameNew.c_str());
    
        ofstream PER_RESIDUE_BB(perResidueBackboneNameNew.c_str());
        ofstream PER_RESIDUE_SC(perResidueSidechainNameNew.c_str());

        ofstream PER_RESIDUE_T(countsNameNew.c_str());

        ofstream PEP2PEP(pep2pepNameNew.c_str());

        int peptideLength = protein.at(0)->chainLength();
        vector<int> backboneCountsTotal(peptideLength, 0);
        vector<int> sidechainCountsTotal(peptideLength, 0);

	//make this a COMMENT function
        PEP2PEP << "# pp_inter_total pp_intra_total"<<endl;
	do {
#ifdef DEBUG_SHOWFRAME
               //cout << "\n############  FRAME # " << t << " #############" << endl;
                cout<<"\n"<< t <<" ";
#endif
               set_pbc(&pbc, -1, box);
        //	rm_pbc(&top.idef, natoms, box, x, x);

               real my_box [] = {box[XX][XX],box[YY][YY], box[ZZ][ZZ]};
               
                //total number of backbone hb made for each inositol
               vector<int> inositolBack(totalNumInos,0);
               //total number of sidechain hbs made for each inositol
               vector<int> inositolSideLys(totalNumInos,0);
               vector<int> inositolSideGlu(totalNumInos,0);

                //total # of HB to inositol per peptide
               vector<int> peptideBbToInosTotalHB(totalNumChains,0);
               vector<int> peptideScToInosTotalHB(totalNumChains,0);

               vector<int> backboneCounts(peptideLength, 0);
               vector<int> backboneCountsDummy(peptideLength, 0);
               vector<int> sidechainCounts(peptideLength, 0);
               vector<int> sidechainCountsDummy(peptideLength, 0);

               peptide_inositol(&pbc, x, protein, inoMolecules, totalNumChains, totalNumInos, peptideBbToInosTotalHB, peptideScToInosTotalHB, inositolBack, inositolSideLys, inositolSideGlu, backboneCounts, sidechainCounts);
                
               int pp_total = peptide_peptide(&pbc, x, protein, totalNumChains, backboneCountsDummy, sidechainCountsDummy);
               int pp_total_intra = peptide_peptide_intra(&pbc, x, protein, totalNumChains, backboneCountsDummy, sidechainCountsDummy);
                PEP2PEP << t<<" "<<pp_total<<" "<<pp_total_intra<<endl;

                //output binding modes
                if(totalNumChains == 1) {
                    peptide* chaini = protein.at(0);
                    
                    chaini->getAminoAcid(K)->reset();
                    chaini->getAminoAcid(E)->reset();
                
                    for(int i=0; i<backboneCounts.size(); i++){
                        PER_RESIDUE_BB << backboneCounts[i] << " ";
                        PER_RESIDUE_SC << sidechainCounts[i] << " ";
                    }
                    PER_RESIDUE_BB<<endl;
                    PER_RESIDUE_SC<<endl;
                }


                PERINOSITOL_BB<< t <<" ";
                PERINOSITOL_SC_LYS<< t <<" ";
                PERINOSITOL_SC_GLU<< t <<" ";
                PERINOSITOL_TOT<< t <<" ";

                //sum per residue counts over all chains
               for(int i=0; i<backboneCounts.size(); i++){
                   backboneCountsTotal[i] += backboneCounts[i];
                   sidechainCountsTotal[i] += sidechainCounts[i];
               }
                
               for( int ins=0; ins < totalNumInos; ins++ ) {
                   PERINOSITOL_BB<<inositolBack[ins]<<" ";
               }
               PERINOSITOL_BB<<endl;
                
               
               for( int ins=0; ins < totalNumInos; ins++ ) {
                   PERINOSITOL_SC_LYS<<inositolSideLys[ins]<<" ";
               }
               PERINOSITOL_SC_LYS<<endl;

               for( int ins=0; ins < totalNumInos; ins++ ) {
                   PERINOSITOL_SC_GLU<<inositolSideGlu[ins]<<" ";
               }
               PERINOSITOL_SC_GLU<<endl;

	       for( int ins=0; ins < totalNumInos; ins++) {
		   PERINOSITOL_TOT << inositolBack[ins] + inositolSideLys[ins] + inositolSideGlu[ins] <<" ";
	       }
	       PERINOSITOL_TOT<<endl;

               for( int pi=0; pi < totalNumChains; pi++ ){
                   PEPBB<<peptideBbToInosTotalHB[pi]<< " ";
                   PEPSC<<peptideScToInosTotalHB[pi]<< " ";
               }
               PEPBB<<endl;
               PEPSC<<endl;
               
        } while (read_next_x(status,&t,natoms,x,box));

        for( int i=0; i<backboneCountsTotal.size(); i++ ) {
            PER_RESIDUE_T << backboneCountsTotal[i] << " ";
        }
        PER_RESIDUE_T << endl;

        for( int i=0; i<sidechainCountsTotal.size(); i++ ) {
            PER_RESIDUE_T << sidechainCountsTotal[i] << " ";
        }
        PER_RESIDUE_T << endl;
}

