#include "peptide.h"
using namespace std;

//#define DEBUG_PEPTIDE
//#define DEBUG_PEPTIDE_HB
#define EDGE

peptide::peptide()
:chain(0), sequence("sequence") {
}

peptide::peptide(const string& seq)
:chain(0), sequence(seq){ 

}

peptide::~peptide() {
	for(int i=0; i<chain.size(); i++){
		delete chain.at(i);
	}
}

void peptide::attachAminoAcid(AminoAcid* aa) {
	chain.push_back(aa);
}

//pre: a valid position on the sequence
//post: returns the AminoAcid at position i or NULL if chain is 0
AminoAcid* peptide::getAminoAcid(int i) const {
        if(i < 0 || chain.size() == 0)
            return 0;

	return chain.at(i);
}

int peptide::chainLength() const {
	return chain.size();
}

//pre: a periodic boundary condition, snapshot coordinates not NULL, counts size of length of peptide sequence
//post: returns the total number of hydrogen bonds between the backbone of the peptide represented by this class 
//and another AminoAcid object "other" (which can be inositol)
int peptide::computeHBtoBackbone(t_pbc* pbc, rvec* x, AminoAcid* other, vector<int> &counts) {
    int toChainBackTotal=0;
    for(int aa=0; aa<chain.size(); aa++){
        AminoAcid* aAmino = chain.at(aa);
        int count = aAmino->computeHBtoBackbone(pbc, x, other);

#ifdef EDGE
        if(count > 0 && other->name() == "INS"){
            cout<<aAmino->resNum()<<" ";
        }
#endif

        counts[aa] += count;
        toChainBackTotal += count;
    }
    return toChainBackTotal;
}

int peptide::computeHBtoBackboneChain(t_pbc* pbc, rvec* x, peptide* otherChain, vector<int>&counts){
    int total=0;
    for(int aai=0; aai<otherChain->chainLength(); aai++) {
        AminoAcid* otherAA = otherChain->getAminoAcid(aai);
        total += computeHBtoBackbone(pbc, x, otherAA, counts);
    }
    return total;
}

int peptide::computeHBtoBackboneAA(t_pbc* pbc, rvec* x, AminoAcid* other, int aminoAcidThis) {
    int toChainBackTotal=0;

    AminoAcid* aAmino = chain.at(aminoAcidThis);
    toChainBackTotal += aAmino->computeHBtoBackbone(pbc, x, other);

    return toChainBackTotal;
}

//pre: a periodic boundary condition, snapshot coordinates not NULL, counts size of length of peptide sequence
//post: returns the total number of hydrogen bonds between the sidechain of the peptide represented by this class 
//and another AminoAcid object "other" (which can be inositol)
int peptide::computeHBtoSidechain(t_pbc* pbc, rvec* x, AminoAcid* other, vector<int> &counts) {
    int toChainSideTotal=0;
    for(int aa=0; aa<chain.size(); aa++) {
        AminoAcid* aAmino = chain.at(aa);
        int count =  aAmino->computeHBtoSidechain(pbc,x,other);

#ifdef EDGE
        if(count > 0 && other->name() == "INS"){
            cout<<aAmino->resNum()<<" ";
        }
#endif

        counts[aa] += count;
        toChainSideTotal += count;
    }
    return toChainSideTotal;
}

int peptide::computeHBtoSidechainChain(t_pbc* pbc, rvec* x, peptide* otherChain, vector<int> &counts) {
    int total=0;
    for(int aai=0; aai<otherChain->chainLength(); aai++){
        AminoAcid* otherAA = otherChain->getAminoAcid(aai);
        total += computeHBtoSidechain(pbc, x, otherAA, counts);
    }

    return total;
}

int peptide::computeHBtoSidechainAA(t_pbc* pbc, rvec* x, AminoAcid* other, int aminoAcidThis) {
    int toChainSideTotal=0;

    AminoAcid* aAmino = chain.at(aminoAcidThis);
    toChainSideTotal += aAmino->computeHBtoSidechain(pbc,x,other);

    return toChainSideTotal;
}

string peptide::aminoAcidName(int aa) {

    return chain.at(aa)->name();

}


