#include "PepGroup.h"

//default constructor
//initialize each vector to have size of zero
PepGroup::PepGroup()
:names(0), pepId(0), resId(0)
{
}

//constructor initializing with one peptide group
//(id of the peptide, id of the residue, name of the group)
PepGroup::PepGroup(int pId, int rId, const string& group)
:names(1), pepId(1), resId(1)
{
	names[0] = group;
	pepId[0] = pId;
	resId[0] = rId;
}
//adds a peptide group specified by (peptide Id, residue Id, name of group)
void PepGroup::addPepGroup(int pId, int rId, const string& group){
	names.push_back(group);
	pepId.push_back(pId);
	resId.push_back(rId);
}

//returns the first in the list of peptide backbone groups that the OH group
//of inositol is bound to
int PepGroup::getPepId(int bgroupNum){
	return pepId.at(bgroupNum);
}

//returns the first residue id belonging to the peptide the OH group is bound to
int PepGroup::getResId(int bgroupNum){
	return resId.at(bgroupNum);
}

//returns the first name of peptide group the OH group is bound to
string PepGroup::getGroupName(int bgroupNum){
	return names.at(bgroupNum);
}

//returns the number of bound peptide backbone groups
//the number of bound groups should be equal to names.size()==pepId.size()==resId.size()
int PepGroup::numGroups(){
	return names.size();
}

