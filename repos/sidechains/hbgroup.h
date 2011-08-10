#ifndef _HBGROUP_H_
#define _HBGROUP_H_


#include <utility>
#include <string>
#include <iostream>

#include "HB.h"

using namespace std;

class hbgroup {

	public:
		hbgroup();
		~hbgroup();
		hbgroup(const string &name, int pos, int neg, const string &type); 
		void setPair(int pos, int neg); 
		int getElectronPos() const;
		int getElectronNeg() const;
		string getName() const;
		string getType() const;
                int isHydrogenBondedTo(t_pbc* pbc, rvec* x, int pos, int neg, const string &type);
                bool nonempty() const;


	private:
		string name;  //name of the group
		pair<int,int> indexPair; //the (electropos, electronneg) pair gromacs atomic index
		string type;   //"donor" or "acceptor"
                HB hb;
};

#endif 

