#include "hbgroup.h"

//#define DEBUG_HB

hbgroup::hbgroup()
:name("name"), indexPair(-1,-1), type("type") , hb(0.35, 0.25,120){
}

hbgroup::hbgroup(const string &n, int pos, int neg, const string &t)
:name(n), indexPair(pos,neg), type(t), hb(0.35, 0.25,120) {
}

hbgroup::~hbgroup() {
// #ifdef DEBUG
// 	cerr<<"DEBUG: hbgroup destroyed"<<endl;
// #endif
}

int hbgroup::isHydrogenBondedTo(t_pbc* pbc, rvec* x, int pos, int neg, const string &type){
        if(type == "donor") {
            if(hb.isHbonded(pbc, x[neg], x[pos], x[getElectronNeg()], x[getElectronPos()]) == true) {

#ifdef DEBUG_HB
        cout<<" hbonded to donor "<<endl;
#endif

                return 1;
            }
        }else if(type == "acceptor") {
            if(hb.isHbonded(pbc, x[getElectronNeg()], x[getElectronPos()], x[neg], x[pos]) == true) {

#ifdef DEBUG_HB
        cout<<" hbonded to acceptor "<<endl;
#endif

                return 1;
            }
        }

        return 0;

}

void hbgroup::setPair(int pos, int neg){
	indexPair=make_pair(pos,neg);
}

int hbgroup::getElectronPos() const {
	return indexPair.first;
}

int hbgroup::getElectronNeg() const {
	return indexPair.second;
}

string hbgroup::getName() const {
	return name;
}

string hbgroup::getType() const {
	return type;
}

bool hbgroup::nonempty() const {
    if(getElectronNeg() !=-1){
        return true;
    }
    return false;
}

