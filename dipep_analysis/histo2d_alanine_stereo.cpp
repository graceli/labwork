#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>

using namespace std;

#define DEBUG

const int MAX = 500000;
const double RT = 0.59616;

//returns the percentage
double percent(int num, int denom){
	return double(num)*100/denom;
}

void init_histogram(int numBinsX, int numBinsY, int** histo2d) {
	for(int i=0; i<numBinsX; i++){
		histo2d[i] = new int[numBinsY];
	}
}

void integrate_well_alanine(int wells[6][3], double phi, double psi, string stereo, int count){
	enum {R1, R2, R3, R4, R5, UNDEF};
	enum {EE,AE,AA};
	if(phi >= -180 && phi <=-125 && psi >= 100 && psi <= 180){
		//beta
		if(stereo == "EE"){
			wells[R1][EE]++;
		}else if(stereo == "EA" || stereo == "AE"){
			wells[R1][AE]++;
		}else{
			wells[R1][AA]++;
		}
	}else if(phi > -125 && phi <= -25 && psi >= 90 && psi <= 180){
		//beta
		if(stereo == "EE"){
			wells[R2][EE]++;
		}else if(stereo == "EA" || stereo == "AE"){
			wells[R2][AE]++;
		}else{
			wells[R2][AA]++;
		}
	}else if(phi >= -180 && phi <=-120 && psi >= -75 && psi <= 75){
		//alpha-R
		if(stereo == "EE"){
			wells[R3][EE]++;
		}else if(stereo == "EA" || stereo == "AE"){
			wells[R3][AE]++;
		}else{
			wells[R3][AA]++;
		}		
	}else if(phi >-120 && phi <= -25 && psi >= -75 && psi <= 50) {
		//alpha-R
		if(stereo == "EE"){
			wells[R4][EE]++;
		}else if(stereo == "EA" || stereo == "AE"){
			wells[R4][AE]++;
		}else{
			wells[R4][AA]++;
		}
	}else if(phi >= 50 && phi <= 100 && psi >=-40 && psi <= 50){
		//alpha-L
		if(stereo == "EE"){
			wells[R5][EE]++;
		}else if(stereo == "EA" || stereo == "AE"){
			wells[R5][AE]++;
		}else{
			wells[R5][AA]++;
		}
	}else{
		//points that do not go into either well
		if(stereo == "EE"){
			wells[UNDEF][EE]++;
		}else if(stereo == "EA" || stereo == "AE"){
			wells[UNDEF][AE]++;
		}else{
			wells[UNDEF][AA]++;
		}
	}
}

double binned_angle(double angle, double binWidth){
	return floor((angle+180)/binWidth)*binWidth - 180;
}


//this program reads in a 3 columned datafile.  the expected datafile is a datafile containing phi,psi angles and the stereochemistry of the inositol groups for
//some bidentate state, eg COCO-12
int main(int argc, char* argv[]){
	if(argc < 1){
		cerr<<"usage: histo2d_alanine_stereo <binWidth> < file > outfile"<<endl;
		return 0;
	}
	double binWidth = atof(argv[1]);

	string* stereo = new string[MAX];
	//allocate memory to store all the data points
	double **data = new double*[MAX];
	for(int d=0; d < MAX; d++) {
		data[d] = new double[2];
	}

	//table of PMF basins broken down in to AA, AE(EA), or EE populations (%)
	int wells[6][3];
	//initialize matrix
	for(int i=0; i<6; i++){
		for(int j=0; j<3; j++){
			wells[i][j]=0;
		}
	}

	double phi,psi;
	int totalNumData = 0;
	while(true){
		cin>>phi;
		if(cin.eof()){
			break;
		}
		cin>>psi;
		cin>>stereo[totalNumData];
		
		//readjust phi psi using the binniing approach so that I get consistent numbers with histo2d_alanine
		double phi_new = binned_angle(phi,binWidth);
		double psi_new = binned_angle(psi,binWidth);
		//update wells[][]	
		integrate_well_alanine(wells, phi_new,psi_new,stereo[totalNumData],1);
		//store the phi psi angles read in; don't think I really need this
		data[totalNumData][0]=phi;
		data[totalNumData][1]=psi;
		totalNumData++;
	}

	cout<<"totalNumData = "<<totalNumData<<endl;

	//get row totals , that is get the total number of points in each basin well
	int sum_rows[6]={0,0,0,0,0,0};
	for(int i=0; i<6; i++){
		for(int j=0; j<3; j++){
			sum_rows[i] += wells[i][j];
		}
	}		

	//print stereochemical analysis 
	int sum_wells=0;
	for(int i=0; i<6; i++){
		for(int j=0; j<3; j++){
			//compute the percentage of the total number of points in a certain basin (or well)
			// that is AA, AE(EA) or EE
			if(sum_rows[i]){
				cout<<setw(8)<<setfill(' ')<<fixed<<setprecision(2)<<wells[i][j]*100.0/sum_rows[i];
			}else{
				cout<<setw(8)<<setfill(' ')<<wells[i][j];
			}
		}
		//output the total number of points in a basin (row sum of the matrix wells[][])
		//as a percentage of the total num pts for the particular bidentate state (=sum over all entries of wells[][])
		cout<<setw(8)<<setfill(' ')<<fixed<<setprecision(2)<<sum_rows[i]*100.0/totalNumData;
		cout<<endl;
	}
	//print just the % EE,AE, AA - sum over the rows
	int sum_cols=0;
	for(int j=0; j<3; j++){
		sum_cols=0;
		for(int i=0; i<6; i++){
			sum_cols+=wells[i][j];		
		}
		cout<<setw(8)<<setfill(' ')<<fixed<<setprecision(2)<<double(sum_cols)*100.0/totalNumData;
	}
	cout<<endl;
	return 0;
}

