#ifndef _PEPTIDE_H
#define _PEPTIDE_H

#include <fstream>
#include "AminoAcid.h"

using namespace std;

// #include "vec.h"
// #include "pbc.h"

class peptide {

	public:
		peptide();
		peptide(const string& seq);
		~peptide();
		void attachAminoAcid(AminoAcid* aa);
		AminoAcid* getAminoAcid(int i) const;
		int chainLength() const;
                int computeHBtoBackbone(t_pbc* pbc, rvec* x, AminoAcid* other, vector<int>&counts);
                int computeHBtoSidechain(t_pbc* pbc, rvec* x, AminoAcid* other, vector<int>&counts);
                int computeHBtoBackboneChain(t_pbc* pbc, rvec* x, peptide* otherChain, vector<int>&counts);
                int computeHBtoSidechainChain(t_pbc* pbc, rvec* x, peptide* otherChain, vector<int>&counts);

                int computeHBtoBackboneAA(t_pbc* pbc, rvec* x, AminoAcid* other, int aminoAcidThis);
                int computeHBtoSidechainAA(t_pbc* pbc, rvec* x, AminoAcid* other, int aminoAcidThis);
                string aminoAcidName(int aa);


	private:
		vector<AminoAcid*> chain;
		string sequence;
};

#endif
