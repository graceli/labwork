#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "Dipeptide.h"
#include "Inositol.h"
#include <cmath>

//#define DEBUG
using namespace std;

const int MAX=100;
const int NUM_LINES=61;
const int SIZE_RAMA = 60000;

//Note that combinations of peptide groups are declared as an enum in Inositol.h
enum pepOne {CO0, NH1, CO2, NH3};	//peptide groups declared as enum

//determines the stereochemistry of the OH (hydroxyl) groups on the inositol
char get_stereo(string isomer, int groupNum){
	if(isomer == "myo"){
		if(groupNum == 2){
			return 'A';
		}
		return 'E';
	}else if(isomer == "epi"){
		if(groupNum == 2 || groupNum ==4 ){
			return 'A';
		}
		return 'E';
	}else {
		return '-';
	}
}


bool readGroFile(ifstream& gro, double** dipepCoords, vector<Inos*>& inositols){
	string title, numAtoms;
	string line;

	//get a line and parse each line using the stringstream
	int lineNum=0;
	string resname, atomname;	
	int atomnum;
	double x,y,z;
	int inosIndex=0;
	double** inosCoords;
	inosCoords=new double*[12];
	for(int i=0; i<12; i++){
		inosCoords[i]=new double[3];
	}

	while(getline(gro,line)){
		//ignore comments the number of atoms
		if(lineNum  < 2){
#ifdef DEBUG
			cerr<<lineNum<<" ignored: "<<line<<endl;
#endif
			lineNum++;
			continue;
		}
		//read body: read in dipeptide coordinates
		istringstream ist(line);
		ist>>resname>>atomname>>atomnum>>x>>y>>z;

		//cerr<<resname<<" "<<atomname<<" "<<atomnum<<" "<<x<<" "<<y<<" "<<z<<endl;
		if(lineNum < 12) {
#ifdef DEBUG
			cerr<<lineNum<<":recording dipeptide coordinates"<<endl;
#endif
			dipepCoords[lineNum-2][0]=x; dipepCoords[lineNum-2][1]=y; dipepCoords[lineNum-2][2]=z;
		}else if(lineNum >= 12 && lineNum < 60){
			if(lineNum == 12){
#ifdef DEBUG
				cerr<<lineNum<<":read one dipeptide"<<endl;
#endif
			}
			//read in inositol coordinates
			inosCoords[inosIndex%12][0]=x; inosCoords[inosIndex%12][1]=y; inosCoords[inosIndex%12][2]=z;
			inosIndex++;
			if(inosIndex % 12 == 0){
#ifdef DEBUG
				cerr<<"read 1 inositol"<<endl;
#endif
				Inos* inos = new Inos(inosCoords);
				inositols.push_back(inos);
				
				inosCoords=new double*[12];
				for(int i=0; i<12; i++){
					inosCoords[i]=new double[3];
				}
			}
		}else if(lineNum == 60){
#ifdef DEBUG
			cerr<<lineNum<<" "<<line<<endl;
			cerr<<"ignored box dimensions.Read "<<lineNum<<". Finished reading one frame. "<<endl;
#endif
			return true;
		}
		lineNum++;
	}
	return false;
}

/* *.class file fields
0)NONE 1)CO0CO2_12 2)CO0NH3_12 3)NH1CO2_12 4)NH1NH3_12 5)CO0CO2_13 6)CO0NH3_13 7)NH1CO2_13 8)NH1NH3_13 9)CO0CO2_14 
10)CO0NH3_14 11)NH1CO2_14 12)NH1NH3_14 13) MONO 14) isDIID 15)Phi 16) Psi 
*/
void print_state(int time, int state[14], bool diid, double** rama, char stereo[]){
	for(int i=0; i<14; i++){
		cout<<state[i]<<"  ";
	}
	cout<<diid<<"   "<<rama[time][0]<<"  "<<rama[time][1]<<"  ";
	cout<<stereo[0]<<stereo[1]<<endl;
}
void readRamaFile(ifstream& rama, double** ramaHistogram){
	string line;
	int lineNum=0;
	double phi, psi;
	while(getline(rama,line)){
		string first = line.substr(0,1);
		if(first == "#" || first == "@"){
			continue;
		}
	//#ifdef DEBUG
		cerr<<line<<endl;
	//#endif 
		istringstream ist(line);
		ist>>phi>>psi;
		ramaHistogram[lineNum][0]=phi;
		ramaHistogram[lineNum][1]=psi;
		lineNum++;
	}
}

//function is called when numPepHbonds > 2 and classifies inositol as 
//bi-monodentate (2 to 1) or bi-monodentate (1 to 2) or one of the 12 bidentate states
//also sets the stereochemistry (that is, axial or equatorial) for each of the OH group participating in the hydrogen bonds
//in the bidentate states
bidstates classifyBidInos(string isomer, int dipepHbonds[], char stereo[]){
	int CONH_long=abs(dipepHbonds[NH3]-dipepHbonds[CO0]);
	int CONH_short=abs(dipepHbonds[CO2]-dipepHbonds[NH1]);
	int COCO=abs(dipepHbonds[CO2] - dipepHbonds[CO0]);
	int NHNH=abs(dipepHbonds[NH3]-dipepHbonds[NH1]);
	
	
	int index[4]={0,0,0,0};	 //stores dipepHbond[i]>0;
	int count=0;

	//test for bi-monodentate states
	for(int i=0; i<4; i++){
		if(dipepHbonds[i] > 0){
			count++;
			index[count-1]=dipepHbonds[i];
		}
	}

	if( (count==1) || ((count == 2) && (index[0] == index[1])) )  {
		return MONO_TWO;
	}
	
	//test for bidentate states
	if(dipepHbonds[NH3] && dipepHbonds[CO0]){
		stereo[0] = get_stereo(isomer, dipepHbonds[CO0]);
		stereo[1] = get_stereo(isomer, dipepHbonds[NH3]);

		if(CONH_long>3){
			CONH_long=6-CONH_long;
		}
		switch(CONH_long){
			case 1:
				return CO0NH3_12;
			case 2: 
				return CO0NH3_13;
			case 3:
				return CO0NH3_14;
			default:
				break;
		}
	}

	if(dipepHbonds[CO2] && dipepHbonds[NH1]){
		stereo[0] = get_stereo(isomer, dipepHbonds[NH1]);
		stereo[1] = get_stereo(isomer, dipepHbonds[CO2]);

		if(CONH_short>3){
			CONH_short=6-CONH_short;
		}
		switch(CONH_short){
			case 1:
				return NH1CO2_12;
			case 2: 
				return NH1CO2_13;
			case 3:
				return NH1CO2_14;
			default:
				break;
		}
	}

	if (dipepHbonds[CO2] && dipepHbonds[CO0]){
		stereo[0] = get_stereo(isomer, dipepHbonds[CO0]);
		stereo[1] = get_stereo(isomer, dipepHbonds[CO2]);
		
		if(COCO>3){
			COCO=6-COCO;
		}
		switch(COCO){
			case 1:
				return CO0CO2_12;
			case 2:
				return CO0CO2_13;
			case 3: 
				return CO0CO2_14;
			default:
				break;
		}
	}

	if(dipepHbonds[NH3] && dipepHbonds[NH1]){
		stereo[0] = get_stereo(isomer, dipepHbonds[NH1]);
		stereo[1] = get_stereo(isomer, dipepHbonds[NH3]);

		if(NHNH>3){
			NHNH=6-NHNH;
		}
		switch(NHNH){
			case 1:
				return NH1NH3_12;
			case 2:
				return NH1NH3_13;
			case 3:
				return NH1NH3_14;
			default:
				break;
		}
	}
	
	return MONO_TWO;
}

bool isDIID(vector<Inos*>&inositols){
	for(int i=0; i<inositols.size(); i++){
		//if i is hbonded to the dipeptide
		if(inositols[i]->isHbondedPep()){
			for(int j=0; j<inositols.size(); j++){
				//if i is hbonded to j and j is hbonded to the dipeptide, then we have a DIID snapshot
				if(i!=j && inositols[i]->isInosHbonded(j) && inositols[j]->isHbondedPep()){
					//make sure the inositols aren't bound to the same peptide groups (if they are bond monodentate)
					if(inositols[i]->getBindMode() != inositols[j]->getBindMode())
						return true;
				}
			}
		}
	}
	return false;
}


int main(int argc, char* argv[]){

	if(argc < 4){
		cerr<<"Usage: classify <grofile-name> <rama-file> <isomer> "<<endl;
		return 0;
	}
	
	vector<Inos*> inositols;
	string grofile(argv[1]);
	string rama(argv[2]);
	string isomer(argv[3]);
	double** dipep=new double*[10];
	for(int i=0; i<10; i++){
		dipep[i]=new double[3];
	}
	ifstream gro(grofile.c_str());
	ifstream rama_total(rama.c_str());

	double CO1coords[3], CO2coords[3];
	double OHcoords[3];
	//int inosHbonds[6]={0,0,0,0,0,0};	
	
	int time=0;

	double** ramaHistogram = new double*[SIZE_RAMA];
	for(int i=0; i<SIZE_RAMA; i++){
		ramaHistogram[i] = new double[2];
	}

	readRamaFile(rama_total,ramaHistogram);

	while(readGroFile(gro,dipep,inositols)){
		bool bid=false;
		bool mono=false;
		Dipeptide ala(dipep);
		ala.getCO(0,CO1coords);
		ala.getCO(1,CO2coords);

		#ifdef DEBUG
			cout<<"time/frame = "<<time<<endl; 
		#endif

		//snapshot level
		int snapshot[14]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};		//stores stats for each snapshot
		bool diid=false;
		char stereo[2]={'-','-'};

		for(int ins=0; ins<inositols.size(); ins++){
			int dipepHbonds[4]={0,0,0,0};	//stores the OH group of inositol ins that makes hbonds with the dipeptide
			int numPepHbonds=0;		//count of the number of hbonds the peptide makes with inositol ins

			bidstates inos_mode;
			//inositol level
			//determine if inositol ins is mono, bidentate, or none by determining which group(s) on the dipeptide it hbonds with
			//determine if inositol ins is hbonded to another inositol, if so, store the ID of the other inositol
			for(int g=0; g<6; g++){
				inositols[ins]->getOH(g,OHcoords);
				bool CO_0 = inositols[ins]->isHbonded(g,CO1coords);	//CO0
				bool CO_2 = inositols[ins]->isHbonded(g,CO2coords);	//CO2
				bool NH_1 = ala.isHbonded(0,OHcoords);			//NH1
				bool NH_3 = ala.isHbonded(1,OHcoords);			//NH3
				
				//label the dipeptide group with the # of the inositol-OH group that forms the hbond  
				//if more than one OH group forms an hbond with a dipeptide group, the # of the LAST hbonded inositol-OH 
				//is recorded
				if(CO_0 || CO_2 || NH_1 || NH_3){
					if(CO_0) {
						dipepHbonds[CO0]=g+1;
#ifdef DEBUG
						cout<<"CO-0 ";
#endif
					}
					if(CO_2) {
						dipepHbonds[CO2]=g+1;
#ifdef DEBUG
						cout<<"CO-2 ";
#endif
						//numPepHbonds++;
					}
					if(NH_1) {
						dipepHbonds[NH1]=g+1;
#ifdef DEBUG
						cout<<"NH-1 ";
#endif
					}
					if(NH_3) {
						dipepHbonds[NH3]=g+1;
#ifdef DEBUG
						cout<<"NH-3 ";
#endif
					}
					numPepHbonds++;
				}
#ifdef DEBUG
				else{
					cout<<"- ";
				}
#endif
				
				//determine if OH group g is hbonded to OH group of another inositol
				for(int k=0; k<inositols.size(); k++){
					if(k!=ins){	//assumed that inositols do not form hydrogen bonds with themselves
						//check for hbonds between inositols
						for(int j=0; j<6; j++){
							double OHcoordsk[3];
							inositols[k]->getOH(j,OHcoordsk);
							bool gDonate = inositols[ins]->isHbonded(g,OHcoordsk);	//is g of ins donating?
							bool gAccept = inositols[k]->isHbonded(j,OHcoords);	//is g of ins accepting (or is j of k donating?)
							if(gDonate || gAccept){
								inositols[ins]->setHbondedInos(k);	//store the inositol k that inositol ins is hbonded to
								break;					//stop looking for hbonds on the same inositol
							}
						}
					}
				}
			}
#ifdef DEBUG
			cout<<endl;
#endif
			//classify peptide binding mode of the inositol - is it mono, bidentate, or non-binding?

			switch(numPepHbonds){
				case 0: //non-bonding
					inos_mode=NONE;
					inositols[ins]->setBindMode(inos_mode);
					break;
				case 1: //monodentate
					if(dipepHbonds[CO0]) {
						inos_mode=MONO_CO0;
						//stereo[0] = get_stereo(isomer,dipepHbonds[CO0]);
					} else if (dipepHbonds[CO2]) {
						inos_mode=MONO_CO2;
						//stereo[0] = get_stereo(isomer,dipepHbonds[CO2]);
					} else if (dipepHbonds[NH1]) {
						inos_mode=MONO_NH1;
						//stereo[0] = get_stereo(isomer,dipepHbonds[NH1]);
					} else if (dipepHbonds[NH3]) {
						inos_mode=MONO_NH3;
						//stereo[0] = get_stereo(isomer,dipepHbonds[NH3]);
					}
					inositols[ins]->setBindMode(inos_mode);
					inositols[ins]->setHbondedPep();
					break;
				default:
					inos_mode=classifyBidInos(isomer,dipepHbonds,stereo);
					inositols[ins]->setBindMode(inos_mode);
					inositols[ins]->setHbondedPep();
					break;
			}
			if(inos_mode >= MONO_CO0)
				snapshot[13]++;	//increment count of mono states
			else {
				snapshot[inos_mode]++;
			}
		}
		//classify if snapshot consists of a Dipeptide -Inositol-Inositol-Dipeptide hydrogen bonding network
		//this condition is too loose..but will do for now
		if(snapshot[NONE] < 4 ){	
			diid = isDIID(inositols);
		}
#ifdef DEBUG
		cout<<endl;
		for(int ii=0; ii<4; ii++){
			cout<<inositols[ii]->getBindMode()<<"\t";
		}
		cout<<endl;
#endif
		print_state(time,snapshot,diid,ramaHistogram, stereo);
#ifdef DEBUG
		cout<<endl;
#endif
		inositols.resize(0);
		time++;
	}
}




