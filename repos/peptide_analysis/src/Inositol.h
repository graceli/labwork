#ifndef _INOSITOL_H
#define _INOSITOL_H

#include "PepGroup.h"
#include "HB.h"
#include <vector>
#include <iostream>
#include <cmath>

using namespace std;
//const double PI=3.141592;

class Inositol{
	public:
		Inositol(int id, double** coords);
		~Inositol();
		void getOCoords(int groupId, double O[3]);
		void getHCoords(int groupId, double H[3]);
		void setBoundGroup(int OHgroupId, PepGroup* boundGroup);
		PepGroup* getBoundGroup(int OHgroupId);
		void print();
		int totalNumHBs();		
		int numChainsBound(int numChains);
		int totalBoundGroups();

	private:
		int inosId;	  				//the id of the inositol
		double** OHcoords; 				//coordinates of the Oxygen of the OH groups of an inositol (12 by 3 matrix)
								//where O coordinates are at even indices, and H coords at odd indices
		vector<PepGroup*>boundGroups;   		//array of peptide groups bound to an inositol; 
								//the indices of the vector map to the OH group ID  of the inositol
};
#endif


