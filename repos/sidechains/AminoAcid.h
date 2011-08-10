#ifndef _AMINO_ACID_H
#define _AMINO_ACID_H


#include <vector>
#include <string>

#include "hbgroup.h" //contains C libraries

class AminoAcid {

	public:
		AminoAcid();
		AminoAcid(const string& n, int rnum);
                AminoAcid(const string& n);

		~AminoAcid();
		string name() const;
                int resNum();
		int getBackboneDonorPos() const;
		int getBackboneDonorNeg() const;
		int getBackboneAcceptorPos() const;
		int getBackboneAcceptorNeg() const;

		int getSidechainDonorPos(int i) const;
		int getSidechainDonorNeg(int i) const;
		int getSidechainAcceptorPos(int i) const;
		int getSidechainAcceptorNeg(int i) const;

		void setBackboneGroup(const string& name, int pos, int neg, const string& type);
		void setSidechainGroup(const string& name, int pos, int neg, const string& type);
		
		int getNumDonors() const;  //only gets the sidechain donors 
		int getNumAcceptors() const; //only gets the sidechain acceptors

                bool hasBackbone();
                bool hasSidechain();
                int computeHBtoBackbone(t_pbc* pbc, rvec* x, AminoAcid* other);
                int computeHBtoSidechain(t_pbc* pbc, rvec* x, AminoAcid* other);
                void setNumGroupsBound(int nbound);
                void setNumGroupsBoundOther(int nbound);
                int getNumBound();
                int getNumBoundOther();
                void reset();

//                 int totalAcceptorGroupsBound();
//                 int totalDonorGroupsBound();
//                 void setSidechainAcceptorGroupBound(int i);
//                 void setSidechainDonorGroupBound(int i);
//                 void resetBound();

	private:
		string resname;  // three letter short name for an amino acid
                int resnum;      // the residue number read in from the gmx topology
		hbgroup* backboneDonor;  //backbone hydrogen bonding groups
		hbgroup* backboneAcceptor;  //backbone hydrogen bonding groups
		vector<hbgroup*> sidechainDonor; //hydrogen bonding groups belonging to sidechains, if any
		vector<hbgroup*> sidechainAcceptor; //hydrogen bonding groups belonging to sidechains, if any
		vector<bool> sidechainDonorBound;
		vector<bool> sidechainAcceptorBound;
                int totalNumHBGroupsBound;
                int totalNumHBGroupsBoundOther;
};
#endif




