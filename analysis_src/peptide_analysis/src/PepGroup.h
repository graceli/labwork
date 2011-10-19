#ifndef _PEPGROUP_H_
#define _PEPGROUP_H_

#include <string>
#include <vector>

using namespace std; 	//should be avoided -- http://www.parashift.com/c++-faq-lite/coding-standards.html#faq-27.5

//class representing a bound peptide group to inositol
//note that the peptide group does not know which inositol, and which inositol OH group it is bound to
//this is a consequence of working from the inositol perspective
class PepGroup {
	public:
		PepGroup();
		//constructor initialized with (peptide id, residue id, group name (string))
		PepGroup(int pId, int rId, const string& group);
		void addPepGroup(int pId, int rId, const string& group);
		int getPepId(int bgroupNum);
		int getResId(int bgroupNum);
		string getGroupName(int bgroupNum);	//returns the group name of the first bound pep group
		int numGroups();

	private:
		vector<string> names;	//names of the group (CO or NH)
		vector<int> pepId;		//ids of the peptide chain the group sits on
		vector<int> resId; 		//ids of the residue the group sits on
};
#endif

