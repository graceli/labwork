#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;

#define DEBUG

const int MAX = 500000;
const double RT = 0.59616;
enum {R1, R2, R3, R4, R5, R6, UNDEF};

//returns the percentage
double percent(int num, int denom){
	return double(num)*100/denom;
}

void init_histogram(int numBinsX, int numBinsY, int** histo2d) {
	for(int i=0; i<numBinsX; i++){
		histo2d[i] = new int[numBinsY];
	}
}

//initial build histogram
//and store the phi, psi data read in
int build_histogram(double** data, int** histo2d, double binWidthX, double binWidthY){
	double phi,psi;
	int totalNumData = 0;
	while(true){
		cin>>phi;
		if(cin.eof()){
			break;
		}
		cin>>psi;
		//store the phi psi angles read in
		data[totalNumData][0]=phi;
		data[totalNumData][1]=psi;
		//calculate bin numbers
		int binIndX = (int)floor((phi)/binWidthX);
		int binIndY = (int)floor((psi)/binWidthY);
		//histogram (phi,psi)
		histo2d[binIndX][binIndY]++;
		totalNumData++;
	}
	//cerr<<"total points = "<<totalNumData<<endl;
	return totalNumData;
}

//prints the histogram
void histo(int** histo2d, int numBinsX, int numBinsY, double binWidthX, double binWidthY, int max){
	
	//int glycine_wells[6]={0,0,0,0,0,0};
	double total=0;
	for(int i=0; i<numBinsX; i++){
		for(int j=0; j<numBinsY; j++){
			double phi = i*binWidthX;
			double psi = j*binWidthY;
			cout<<phi<<" "<<psi<<" "<<(double) histo2d[i][j]/max<<endl;
		}
		cout<<endl;
	}
}

void integrate_well_alanine(int wells[], double phi, double psi, int count){
	if(phi >= -180 && phi <-125 && psi >= 100 && psi < 180){
		//beta
		wells[R1]+=count;
	}else if(phi >= -125 && phi < -25 && psi >= 90 && psi < 180){
		//beta
		wells[R2]+=count;
	}else if(phi >= -180 && phi < -120 && psi >= -75 && psi < 75){
		//alpha-R
		wells[R3]+=count;
	}else if(phi >= -120 && phi < -25 && psi >= -75 && psi < 50) {
		//alpha-R
		wells[R4]+=count;
	}else if(phi >= 50 && phi < 100 && psi >=-40 && psi < 50){
		//alpha-L
		wells[R5]+=count;
	}else{
		//points that do not go into either well
		wells[UNDEF]+=count;
	}
}

void integrate_well_glycine(int wells[], double phi, double psi, int count){
	if(phi >= -180 && phi < -50 && psi >= 75 && psi < 180){
		//beta
		wells[R1]+=count;
	}else if(phi >= -180 && phi < -50 && psi >= -75 && psi < 75){
		//beta
		wells[R2]+=count;
	}else if(phi >= -180 && phi < -50 && psi >= -180 && psi < -75){
		//alpha-R
		wells[R3]+=count;
	}else if(phi >= 50 && phi < 180 && psi >= 75 && psi < 180) {
		//alpha-R
		wells[R4]+=count;
	}else if(phi >= 50 && phi < 180 && psi >=-75 && psi < 75){
		//alpha-L
		wells[R5]+=count;
	}else if (phi >= 50 && phi < 180 && psi >=-180 && psi < -75){
		wells[R6]+=count;
	}else{
		//points that do not go into either well
		wells[UNDEF]+=count;
	}
}

// output phi-psi probability distribution as a PMF, with the x-axis as phi, and y as psi
// as well as integrate over the basins
// well 0 = R1, 1=R2, 2=R3, 3=R4...,5=UNDEFINED
void pmf(int** histo2d, int numBinsX, int numBinsY, double binWidthX, double binWidthY, int max, int min, const string& pmfName){
	int wells[7] = {0,0,0,0,0,0,0};
	ofstream pmfout(pmfName.c_str());
	for(int i=0; i<numBinsX; i++){
		double phi = i*binWidthX;
		for(int j=0; j<numBinsY; j++){
			double psi = j*binWidthY;
			pmfout <<phi<<" "<<psi<<" ";
			pmfout << double(histo2d[i][j])/max << endl;
		}
		pmfout<<endl;
	}
}

int find_max(int** histo2d, int numBinsX, int numBinsY, double binWidthX, double binWidthY){
	int max = histo2d[0][0];
	int first_nonzero = 0; 
	for(int i=0; i<numBinsX; i++){
		for(int j=0; j<numBinsY; j++){
			if(histo2d[i][j]>max){
				max = histo2d[i][j];
			}
		}
	}
	return max;
}

int first_nonzero_count(int** histo2d, int numBinsX, int numBinsY, double binWidth){
	//finds the first nonzero integer in the matrix
	//actually all I really need is some nonzero integer in the matrix
	for(int i=0; i<numBinsX; i++){
		for(int j=0; j<numBinsY; j++){	
			if(histo2d[i][j] > 0){
				return histo2d[i][j];
			}
		}
	}
}

//returns the nonzero minimum count in the bins
// for my purpose, I will be using something simple. First find the first nonzero integer in the matrix
// set that as the minimum, loop through entire matrix again to find the nonzero minimum
int find_min(int** histo2d, int numBinsX, int numBinsY, double binWidth){
	int min = first_nonzero_count(histo2d, numBinsX, numBinsY, binWidth);
	for(int i=0; i<numBinsX; i++) {
		for(int j=0; j<numBinsY; j++){
			if(histo2d[i][j] != 0 &&  histo2d[i][j] < min){
				min = histo2d[i][j];
			}
		}
	}
	return min;
}

//this program reads in a 2 columned (col x and col y) datafile and outputs a file H(x,y) where H(x,y) is the number of counts of x,y
int main(int argc, char* argv[]){
	if(argc < 3){
		cerr<<"usage: histo2d <bin-width-x> <bin-width-y> <pmf-name> < file"<<endl;
		cerr<<"eg. histo2d 0.5 0.5 analysis.hist < datafile"<<endl;
		return 0;
	}

	double binWidthX = atof(argv[1]);
	double binWidthY = atof(argv[2]);
	string pmfName(argv[3]);
	//calculate the number of bins needed based on the given bin widths
	int numBinsX = (int)ceil(361.0/binWidthX);
	int numBinsY = (int)ceil(361.0/binWidthY);

	//allocate memory to store all the data points
	double **data = new double*[MAX];
	for(int d=0; d < MAX; d++){
		data[d] = new double[2];
	}

	int** histo2d = new int*[numBinsX];

	//initialize histogram
	init_histogram(numBinsX,numBinsY, histo2d);

	//build phi psi histogram
	//total is the total number of snapshots read in
	int total = build_histogram(data,histo2d, binWidthX, binWidthY);
	int max = find_max(histo2d, numBinsX, numBinsY,binWidthX, binWidthY);
	//int min = find_min(histo2d, numBinsX, numBinsY,binWidthX, binWidthY);
	int min = -1;
	
	pmf(histo2d, numBinsX, numBinsY,binWidthX, binWidthY, max, min, pmfName);

	return 0;
}

