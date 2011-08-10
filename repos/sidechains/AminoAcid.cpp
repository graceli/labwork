#include "AminoAcid.h"
using namespace std;

//#define DEBUG_HB_PP
//#define DEBUG_AA_BB
//#define DEBUG_AA_SC
//#define DEBUG_OTHER

AminoAcid::AminoAcid()
:resname("aa"), resnum(-1), sidechainDonor(0), sidechainAcceptor(0),totalNumHBGroupsBound(0), totalNumHBGroupsBoundOther(0) {
	backboneDonor=new hbgroup();
	backboneAcceptor=new hbgroup();

	sidechainDonorBound.resize(6,false);
	sidechainAcceptorBound.resize(6,false);
}

AminoAcid::AminoAcid(const string& n, int rnum)
:resname(n), resnum(rnum), sidechainDonor(0), sidechainAcceptor(0),totalNumHBGroupsBound(0),totalNumHBGroupsBoundOther(0) {
	backboneDonor=new hbgroup();
	backboneAcceptor=new hbgroup();
	sidechainDonorBound.resize(6,false);
	sidechainAcceptorBound.resize(6,false);
}

AminoAcid::AminoAcid(const string& n)
:resname(n), resnum(-1), sidechainDonor(0), sidechainAcceptor(0),totalNumHBGroupsBound(0),totalNumHBGroupsBoundOther(0) {
	backboneDonor=new hbgroup();
	backboneAcceptor=new hbgroup();
	sidechainDonorBound.resize(6,false);
	sidechainAcceptorBound.resize(6,false);
}


AminoAcid::~AminoAcid(){
	delete backboneDonor;
	delete backboneAcceptor;

	for(int d=0; d<sidechainDonor.size(); d++) {
		delete sidechainDonor.at(d);
	}

	for(int a=0; a<sidechainAcceptor.size(); a++) {
		delete sidechainAcceptor.at(a);
	}
}

int AminoAcid::computeHBtoBackbone(t_pbc* pbc, rvec* x, AminoAcid* other){

#ifdef DEBUG_OTHER
        cout<<"Computing hydrogen bonds to backbone of peptide"<<endl;
#endif

    int toAminoBackTotal = 0;

    //if "other" has backbone then compute number of HB of backbone (other) to backbone (this)
    if(other->hasBackbone()){

#ifdef DEBUG_OTHER
        cout<<"Other has a backbone"<<endl;
#endif

            int otherBackDonorPos = other->getBackboneDonorPos();
            int otherBackDonorNeg = other->getBackboneDonorNeg();
            int otherBackAcceptorPos = other->getBackboneAcceptorPos();
            int otherBackAcceptorNeg = other->getBackboneAcceptorNeg();

            
            int hbondAccept = backboneAcceptor->isHydrogenBondedTo(pbc, x, otherBackDonorPos, otherBackDonorNeg, "donor");
            
#ifdef DEBUG_HB_PP
        if(hbondAccept){
            cout<<backboneAcceptor->getName()<<" "<<backboneAcceptor->getElectronPos()<<" "<<backboneAcceptor->getElectronNeg()<<" "<<otherBackDonorPos<<" "<<otherBackDonorNeg<<" ";
        }
#endif
 
            int hbondDonor = backboneDonor->isHydrogenBondedTo(pbc,x, otherBackAcceptorPos, otherBackAcceptorNeg, "acceptor");

#ifdef DEBUG_HB_PP
        if(hbondDonor){
            cout<<backboneDonor->getName()<<" "<<backboneDonor->getElectronPos()<<" "<<backboneDonor->getElectronNeg()<<" "<<otherBackAcceptorPos<<" "<<otherBackAcceptorNeg<<endl;
        }
#endif
        

        toAminoBackTotal+=hbondAccept+hbondDonor;
    }
    
    if(other->hasSidechain()){

#ifdef DEBUG_OTHER
        cout<<other->name()<<" has a sidechain"<<endl;
#endif

        vector<int>hbondGroupsOther(other->getNumDonors(),0);
        
        //for each donor group in the sidechain count the hydrogen
        //bonds to the backbone of this
        for(int otherDonorIndex = 0; otherDonorIndex < other->getNumDonors(); otherDonorIndex++){
            int otherSideDonorPos=other->getSidechainDonorPos(otherDonorIndex);
            int otherSideDonorNeg=other->getSidechainDonorNeg(otherDonorIndex);

            int hbond = backboneAcceptor->isHydrogenBondedTo(pbc,x,otherSideDonorPos,otherSideDonorNeg, "donor");

#ifdef DEBUG_HB
    if(hbond){
        cout<<"### backbone ##"<<endl;
        cout<<backboneAcceptor->getName()<<" hbonded to "<<otherSideDonorPos<<" "<<otherSideDonorNeg<<endl;
    }
#endif

            toAminoBackTotal += hbond;
        }

        //for each acceptor group in the sidechain count the hydrogen bonds to the backbone of this
        for(int otherAcceptorIndex = 0; otherAcceptorIndex < other->getNumAcceptors(); otherAcceptorIndex++){
            int otherSideAcceptorPos=other->getSidechainAcceptorPos(otherAcceptorIndex);
            int otherSideAcceptorNeg=other->getSidechainAcceptorNeg(otherAcceptorIndex);
            
            int hbond = backboneDonor->isHydrogenBondedTo(pbc, x, otherSideAcceptorPos, otherSideAcceptorNeg, "acceptor");

#ifdef DEBUG_HB
    if(hbond){
        cout<<"### backbone ##"<<endl;
        cout<<backboneDonor->getName()<<" hbonded to "<<otherSideAcceptorPos<<" "<<otherSideAcceptorNeg<<endl;
    }
#endif

            toAminoBackTotal += hbond;
        }

    }

    return toAminoBackTotal;
}

//post: returns the total number of hydrogen made by this amino acid to aminoacid other
//      sets the total number (donor+acceptor) of sidechain hbgroups bound for this and for other
//      double counts the donors and acceptors for other because assumed to be inositol (hacky)
int AminoAcid::computeHBtoSidechain(t_pbc* pbc, rvec* x, AminoAcid* other){
#ifdef DEBUG_OTHER
        cout<<"Computing HB to sidechains of peptides"<<endl;
#endif


    int toAminoSideTotal=0;
    if(other->hasBackbone()){

#ifdef DEBUG_OTHER
        cout<<"Other has a backbone"<<endl;
#endif


        //other backbone donor
        int otherBackDonorPos = other->getBackboneDonorPos();
        int otherBackDonorNeg = other->getBackboneDonorNeg();
        for(int thisSideAcceptorIndex=0; thisSideAcceptorIndex < getNumAcceptors(); thisSideAcceptorIndex++){
            toAminoSideTotal += sidechainAcceptor.at(thisSideAcceptorIndex)->isHydrogenBondedTo(pbc, x, otherBackDonorPos, otherBackDonorNeg, "donor");
        }

        int otherBackAcceptorPos = other->getBackboneAcceptorPos();
        int otherBackAcceptorNeg = other->getBackboneAcceptorNeg();
        for(int thisSideDonorIndex=0; thisSideDonorIndex < getNumDonors(); thisSideDonorIndex++){
            toAminoSideTotal += sidechainDonor.at(thisSideDonorIndex)->isHydrogenBondedTo(pbc, x, otherBackAcceptorPos, otherBackAcceptorNeg, "acceptor");
        }
    }

    //compute the hydrogen bonds of the sidechain of this to the sidechain of "other"
    if(other->hasSidechain()) {

#ifdef DEBUG_OTHER
        cout<<"Other has a sidechain"<<endl;
#endif

        vector<int>hbondGroupsOther(other->getNumDonors(),0);
        vector<int>hbondGroupsThisDonors(getNumDonors(),0);
        vector<int>hbondGroupsThisAcceptors(getNumAcceptors(),0);

        for(int thisSideAcceptorIndex=0; thisSideAcceptorIndex < getNumAcceptors(); thisSideAcceptorIndex++) {
            for(int otherSideDonorIndex=0; otherSideDonorIndex<other->getNumDonors(); otherSideDonorIndex++){
                int otherSideDonorPos = other->getSidechainDonorPos(otherSideDonorIndex);
                int otherSideDonorNeg = other->getSidechainDonorNeg(otherSideDonorIndex);
                
                int hbonds = sidechainAcceptor.at(thisSideAcceptorIndex)->isHydrogenBondedTo(pbc, x, otherSideDonorPos, otherSideDonorNeg, "donor");
    
                //tabulate specific group has a hydrogen bond or not to this
                //note that otherDonor and otherAcceptor groups are overlayed because
                //each hbonding group on inositol is a donor and an acceptor
                //essentially we're counting n. hydrogen bonds made per group of inositol but we don't care 
                //about how many hbs, we just care about how many groups are bound
                if(hbonds){
#ifdef DEBUG_HB
    cout<<"### sidechain acceptor ##"<<endl;
    cout<<sidechainAcceptor.at(thisSideAcceptorIndex)->getName()<<" hbonded to "<<otherSideDonorPos<<" "<<otherSideDonorNeg<<endl;
#endif
                    hbondGroupsOther[otherSideDonorIndex]++;
                    hbondGroupsThisAcceptors[thisSideAcceptorIndex]++;
                }else{
#ifdef DEBUG_HB
    //cout<<sidechainAcceptor.at(thisSideAcceptorIndex)->getName()<<" not hbonded to "<<otherSideDonorPos<<" "<<otherSideDonorNeg<<endl;
#endif
                }
                toAminoSideTotal += hbonds;
            }
        }

        for(int thisSideDonorIndex=0; thisSideDonorIndex < getNumDonors(); thisSideDonorIndex++){
            for(int otherSideAcceptorIndex=0; otherSideAcceptorIndex < other->getNumAcceptors(); otherSideAcceptorIndex++){
                int otherSideAcceptorPos = other->getSidechainAcceptorPos(otherSideAcceptorIndex);
                int otherSideAcceptorNeg = other->getSidechainAcceptorNeg(otherSideAcceptorIndex);
         
                int hbonds = sidechainDonor.at(thisSideDonorIndex)->isHydrogenBondedTo(pbc, x, otherSideAcceptorPos, otherSideAcceptorNeg, "acceptor");
                if(hbonds){

#ifdef DEBUG_HB
    cout<<"### sidechain donor ##"<<endl;
    cout<<sidechainDonor.at(thisSideDonorIndex)->getName()<<" hbonded to "<<otherSideAcceptorPos<<" "<<otherSideAcceptorNeg<<endl;
#endif

                    hbondGroupsOther[otherSideAcceptorIndex]++;   
                    hbondGroupsThisDonors[thisSideDonorIndex]++;
                }else{

#ifdef DEBUG_HB
    //cout<<sidechainDonor.at(thisSideDonorIndex)->getName()<<" not hbonded to "<<otherSideAcceptorPos<<" "<<otherSideAcceptorNeg<<endl;
#endif

                }

                toAminoSideTotal += hbonds;
            }
        }

        //set number of inositol groups bound (this case is a bit of a hack for inositol aminoacids only)
        int totalGroupsOtherBound=0;
        for(int i=0; i<hbondGroupsOther.size(); i++){
            if(hbondGroupsOther[i]!=0){
                totalGroupsOtherBound++;
            }
        }
        //other->setNumGroupsBound(totalGroupsOtherBound);
        
        int totalGroupsThisBoundDonor=0;
        int totalGroupsThisBoundAcceptor=0;
        for(int i=0; i<hbondGroupsThisDonors.size(); i++){
                if(hbondGroupsThisDonors[i]){
                    totalGroupsThisBoundDonor++;
                }
        }
        for(int i=0; i<hbondGroupsThisAcceptors.size();i++){
                if(hbondGroupsThisAcceptors[i]){
                    totalGroupsThisBoundAcceptor++;
                }
        }
        setNumGroupsBound(totalGroupsThisBoundDonor+totalGroupsThisBoundAcceptor);
        setNumGroupsBoundOther(totalGroupsOtherBound);
    }

    return toAminoSideTotal;
}

//pre: this is not null
//post: return true if has a backbone (that is NH and CO) or false
//otherwise
bool AminoAcid::hasBackbone(){
    if(backboneDonor->nonempty() && backboneAcceptor->nonempty()){
        return true;
    }
    return false;
}

//pre: this is not null
//post: return true if has a sidechain or false otherwise
bool AminoAcid::hasSidechain(){
    if(sidechainDonor.size()!=0 || sidechainAcceptor.size()!=0){
        return true;
    }
    return false;
}

string AminoAcid::name() const {
    return resname;
}

int AminoAcid::resNum() {
    return resnum;
}

int AminoAcid::getBackboneDonorPos() const {
		return backboneDonor->getElectronPos();
}

int AminoAcid::getBackboneDonorNeg() const {
	return backboneDonor->getElectronNeg();
}

int AminoAcid::getBackboneAcceptorPos() const {
	return backboneAcceptor->getElectronPos();
}

int AminoAcid::getBackboneAcceptorNeg() const {
	return backboneAcceptor->getElectronNeg();
}

int AminoAcid::getSidechainDonorPos(int i) const {
	if(sidechainDonor.size()==0)
		return -1;
	return sidechainDonor.at(i)->getElectronPos();
}

int AminoAcid::getSidechainDonorNeg(int i) const {
	if(sidechainDonor.size()==0)
		return -1;
	return sidechainDonor.at(i)->getElectronNeg();
}

int AminoAcid::getSidechainAcceptorPos(int i) const {
	if(sidechainAcceptor.size()==0)
		return -1;
	return sidechainAcceptor.at(i)->getElectronPos();
}

int AminoAcid::getSidechainAcceptorNeg(int i) const {
	if(sidechainAcceptor.size()==0)
		return -1;
	return sidechainAcceptor.at(i)->getElectronNeg();
}

void AminoAcid::setBackboneGroup(const string& name, int pos, int neg, const string& type) {
	if(type == "donor") {

#ifdef DEBUG_AA_BB
        cout<<"attached backbone "<<name<<" "<< pos << " "<< neg <<" "<<type<<endl;
#endif

		backboneDonor=new hbgroup(name, pos, neg, type);
	}else{

#ifdef DEBUG_AA_BB
        cout<<"attached backbone "<<name<<" "<< pos << " "<< neg <<" "<<type<<endl;
#endif

		backboneAcceptor=new hbgroup(name, pos, neg, type);
	}
}

void AminoAcid::setSidechainGroup(const string& name, int pos, int neg, const string &type) {
	if(type == "donor"){
#ifdef DEBUG_AA_SC
        cout<<"attached sidechain "<<name<<" "<< pos << " "<< neg <<" "<<type<<endl;
#endif
		sidechainDonor.push_back(new hbgroup(name, pos, neg, type));
	}else{

#ifdef DEBUG_AA_SC
        cout<<"attached sidechain "<<name<<" "<< pos << " "<< neg <<" "<<type<<endl;
#endif

		sidechainAcceptor.push_back(new hbgroup(name, pos, neg, type));
	}
}

int AminoAcid::getNumDonors() const{
	return sidechainDonor.size();
}

int AminoAcid::getNumAcceptors() const {
	return sidechainAcceptor.size();
}

void AminoAcid::setNumGroupsBound(int nbound){
        totalNumHBGroupsBound=nbound;
}

void AminoAcid::setNumGroupsBoundOther(int nbound){
        totalNumHBGroupsBoundOther=nbound;
}

int AminoAcid::getNumBound(){
        return totalNumHBGroupsBound;
}

int AminoAcid::getNumBoundOther(){
        return totalNumHBGroupsBoundOther;
}

void AminoAcid::reset(){
        totalNumHBGroupsBound=0;
        totalNumHBGroupsBoundOther=0;
}
// int AminoAcid::totalAcceptorGroupsBound(){
// 
//         int totalNumBound=0;
// 
//         for(int i=0; i<sidechainAcceptorBound.size(); i++){
//             if(sidechainAcceptorBound.at(i) == true){
//                 totalNumBound++;
//             }
//         }
//         return totalNumBound;
// }
// 
// int AminoAcid::totalDonorGroupsBound(){
// 
//         int totalNumBound=0;
// 
//         for(int i=0; i<sidechainDonorBound.size(); i++){
//             if(sidechainDonorBound.at(i) == true){
//                 totalNumBound++;
//             }
//         }
//         return totalNumBound;
// }
