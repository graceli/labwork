#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

const int MAX = 500000;
const double RT = 0.59616;

//returns the percentage
double percent(int num, int denom){
	return double(num)*100/denom;
}

void init_histogram(int numBinsX, int numBinsY, int** histo2d, int** error, double** sum_block_var){
	for(int i=0; i<numBinsX; i++){
		histo2d[i] = new int[numBinsY];
		error[i] = new int[numBinsY];
		sum_block_var[i] = new double[numBinsY];
	}
}

//initial build histogram
//and store the phi, psi data read in
int build_histogram(double** data, int** histo2d, double binWidth){
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
		int binIndX = (int)floor((phi+180)/binWidth);
		int binIndY = (int)floor((psi+180)/binWidth);
		//histogram (phi,psi)
		histo2d[binIndX][binIndY]++;
		totalNumData++;
	}
	//cerr<<"total points = "<<totalNumData<<endl;
	return totalNumData;
}

void histo(int** histo2d, int numBinsX, int numBinsY, double binWidth, int max){
	
	int glycine_wells[6]={0,0,0,0,0,0};
	double total=0;
	for(int i=0; i<numBinsX; i++){
		for(int j=0; j<numBinsY; j++){
			double phi = i*binWidth-180;
			double psi = j*binWidth-180;
			cout<<phi<<" "<<psi<<" "<<(double) histo2d[i][j]/max<<endl;

			//integrate wells for glycine
			if(phi<=-50 && phi>=-180 && psi <= 180 && psi >=100){
				glycine_wells[0]+=histo2d[i][j];
			} else if (phi<=-50 && phi>=-180 && psi <= 50 && psi >= -50){
				glycine_wells[1]+=histo2d[i][j];
			} else if (phi<=-50 && phi>=-180 && psi <= -100 && psi >= -180){
				glycine_wells[2]+=histo2d[i][j];
			} else if(phi>=50 && phi <=180 && psi >= 100 && psi <= 180){
				glycine_wells[3]+=histo2d[i][j];
			} else if(phi>=50 && phi <=180 && psi >= -50 && psi <= 50){
				glycine_wells[4]+=histo2d[i][j];
			} else if(phi>=50 && phi <=180 && psi <= -100 && psi >= -180){
				glycine_wells[5]+=histo2d[i][j];
			}
			total +=histo2d[i][j];
		}
		cout<<endl;
	}
	cerr<<"total = "<<total<<endl;
	for(int i=0; i<3; i++){
		cerr<<glycine_wells[i]/total<<" "<<glycine_wells[i+3]/total<<" "<<(glycine_wells[i]-glycine_wells[i+3])/total<<endl;
	}
}

void pmf(int** histo2d, int numBinsX, int numBinsY, double binWidth, int max){
	for(int i=0; i<numBinsX; i++){
		for(int j=0; j<numBinsY; j++){
			double phi = i*binWidth-180;
			double psi = j*binWidth-180;
			cout<<phi<<" "<<psi<<" ";
			if(histo2d[i][j]!=0){
				cout<<-RT*log(double(histo2d[i][j])/max)<<endl;
			}else {
				cout<<0<<endl;
			}
		}
		cout<<endl;
	}
}
int find_max(int** histo2d, int numBinsX, int numBinsY, double binWidth){
	int max = histo2d[0][0];
	for(int i=0; i<numBinsX; i++){
		for(int j=0; j<numBinsY; j++){
			double phi = i*binWidth-180;
			double psi = j*binWidth-180;
			if(histo2d[i][j]>max){
				max = histo2d[i][j];
			}
		}
	}
	return max;
}
void integrate_basins(int** histo2d, int numBinsX, int numBinsY, double binWidth, int totalNumData){
	int beta=0, alphaR=0, alphaL=0, pass=0;
	for(int i=0; i<numBinsX; i++){
		for(int j=0; j<numBinsY; j++){
			double phi = i*binWidth-180;
			double psi = j*binWidth-180;
			//cerr<<phi<<" "<<psi<<" "<<histo2d[i][j]<<endl;
			//count beta 
			if(phi >=-180 && phi <=0){
				if(psi >=90 && psi <=180){
					beta += histo2d[i][j];
				}
			}

			//count pass
			if(phi >= -180 && phi<=0){
				if(psi>=30 && psi <90){
					pass += histo2d[i][j];
				}
			}

			//count alpha R
			if(phi >= -180 && phi <= 0){
				if(psi >=-120 && psi < 30){
					alphaR += histo2d[i][j];	
				}
			}

			//count alpha L
			if(phi > 0 && phi <= 120){
				if(psi >= -20 && psi <= 90){
					alphaL += histo2d[i][j];
				}
			}
		}
		cerr<<endl;
	}

	cout<<"beta="<<percent(beta,totalNumData)<<endl;
	cout<<"alphaR="<<percent(alphaR,totalNumData)<<endl;
	cout<<"alphaL="<<percent(alphaL,totalNumData)<<endl;
	cout<<"pass="<<percent(pass,totalNumData)<<endl;
}

//this program reads in a 2 columned (col x and col y) datafile and outputs a file H(x,y) where H(x,y) is the number of counts of x,y
int main(int argc, char* argv[]){
	if(argc < 2){
		cerr<<"usage: histo2d <bin-width> < file > outfile"<<endl;
		return 0;
	}

	double binWidth = atof(argv[1]);
	//calculate the number of bins needed based on the given bin widths
	int numBinsX = (int)ceil(361.0/binWidth);
	int numBinsY = (int)ceil(361.0/binWidth);


	//initialize histogram

	//allocate memory
	double **data = new double*[MAX];
	for(int d=0; d < MAX; d++){
		data[d] = new double[2];
	}
	int** histo2d = new int*[numBinsX];
	int** error = new int*[numBinsX];
	double** sum_block_var = new double*[numBinsX];
	init_histogram(numBinsX,numBinsY, histo2d, error, sum_block_var);

	//build phi psi histogram
	//total -- total number of snapshots read in
	int total=build_histogram(data,histo2d,binWidth);
	int max = find_max(histo2d, numBinsX, numBinsY,binWidth);

	pmf(histo2d, numBinsX, numBinsY,binWidth, max);
	//cout<<max<<endl;
	//histo(histo2d, numBinsX, numBinsY,binWidth, 1);
	
	//integrate_basins(histo2d, numBinsX, numBinsY, binWidth, total);
}

