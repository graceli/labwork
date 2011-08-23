#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

const int MAX = 500000;
//returns the percentage
double percent(int num, int denom){
	return double(num)*100/denom;
}

void init_histogram(int numBinsX, int numBinsY, int** histo2d, int** error, double** sum_block_var, double** run_var){
	for(int i=0; i<numBinsX; i++){
		histo2d[i] = new int[numBinsY];
		error[i] = new int[numBinsY];
		sum_block_var[i] = new double[numBinsY];
		run_var[i]=new double[numBinsY];
		for(int j=0; j<numBinsY; j++){
			histo2d[i][j]=0;
			error[i][j] = 0;
			sum_block_var[i][j] = 0;
			run_var[i][j] = 0;
		}
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
	cerr<<"total points = "<<totalNumData<<endl;
	return totalNumData;
}

double calc_run_var(double** data, int** histo2d, double** run_var, int numBinsX, int numBinsY, int total, double binWidth){
	int** temp = new int*[numBinsX];
	for(int j=0; j<numBinsX; j++){
		temp[j]=new int[numBinsY];
		for(int k=0; k<numBinsY; k++){
			temp[j][k] = 0;
		}
	}

	for(int i=0; i< total; i++){
		double phi = data[i][0];
		double psi = data[i][1];
		//accumulate data point
		int binIndX = (int)floor((phi+180)/binWidth);
		int binIndY = (int)floor((psi+180)/binWidth);
		temp[binIndX][binIndY] ++;
		//update run variances for each bin
		for(int r=0; r<numBinsX; r++){
			for(int c=0; c<numBinsY; c++){
				run_var[r][c]+=pow((temp[r][c] - double(histo2d[r][c])/total),2)/total;
			}
		}
	}

}

void integrate_basins(int** histo2d, int numBinsX, int numBinsY, double binWidth, int totalNumData){
	int beta=0, alphaR=0, alphaL=0, pass=0;
	for(int i=0; i<numBinsX; i++){
		for(int j=0; j<numBinsY; j++){
			double phi = i*binWidth-180;
			double psi = j*binWidth-180;
			cerr<<phi<<" "<<psi<<" "<<histo2d[i][j]<<endl;
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

void reset_histogram(int** histo, int numBinsX, int numBinsY){
	for(int x=0; x<numBinsX; x++){
		for(int y=0; y<numBinsY; y++){
			histo[x][y]=0;
		}
	}
}

//this program reads in a 2 columned (col x and col y) datafile and outputs a file H(x,y) where H(x,y) is the number of counts of x,y
int main(int argc, char* argv[]){
	if(argc < 2){
		cerr<<"usage: histo2d <bin-width>"<<endl;
		return 0;
	}

	double binWidth = atof(argv[1]);
	//calculate the number of bins needed based on the given bin widths
	int numBinsX = (int)ceil(361.0/binWidth);
	int numBinsY = (int)ceil(361.0/binWidth);

	//allocate memory
	double** data = new double*[MAX];
	for(int d=0; d < MAX; d++){
		data[d] = new double[2];
	}
	int** histo2d = new int*[numBinsX];
	int** error = new int*[numBinsX];	//use a matrix identical to histo2d for block accumulation
	double** sum_block_var = new double*[numBinsX];
	double** run_var = new double*[numBinsX];
	//initialize histogram
	init_histogram(numBinsX,numBinsY, histo2d, error, sum_block_var,run_var);

	//build phi psi histogram
	//total -- total number of snapshots read in
	int total=build_histogram(data,histo2d,binWidth);
	calc_run_var(data,histo2d,run_var, numBinsX, numBinsY, total, binWidth);

	for(int block_size=10000; block_size<100000; block_size+=10000){
		int num_blocks = total/50;	//number of whole blocks with size block_size
		//for a block of a fixed size block_size
		for(int i=0; i<total; i++){
			double phi = data[i][0];
			double psi = data[i][1];
			int binIndX = (int)floor((phi+180)/binWidth);
			int binIndY = (int)floor((psi+180)/binWidth);
			if(i % block_size==0 && i>0){
	
				//reached the beginning of a new block
				//update running block variance sums for each bin
				for(int r=0; r<numBinsX; r++){	
					for(int c=0; c<numBinsY; c++){
						//error[][] is the accumulated sums for the current block, and the 
						if(error[r][c]!=0) {
							sum_block_var[r][c] += pow(double(error[r][c])/block_size - histo2d[r][c]/total,2)/num_blocks;
						}
					}
				}
				//reset counts in each bin of interest
				reset_histogram(error,numBinsX,numBinsY);
				//don't accumulate data for blocks smaller than a whole block
				if(total-i < block_size){
					break;
				}
			}else{
				//otherwise, accumulate, but only for the interested regions
				//accumulate block sums for bins in beta region
				if(phi >=-180 && phi <=0){
					if(psi >=90 && psi <=180){
						error[binIndX][binIndY]++;
					}
				}
	
				//pass
				if(phi >= -180 && phi<=0){
					if(psi>=30 && psi <90){
						error[binIndX][binIndY]++;
					}
				}
	
				//alpha R
				if(phi >= -180 && phi <= 0){
					if(psi >=-120 && psi < 30){
						error[binIndX][binIndY]++;
					}
				}
	
				//alpha L
				if(phi > 0 && phi <= 120){
					if(psi >= -20 && psi <= 90){
						error[binIndX][binIndY]++;
					}
				}
			}
		}
		
		bool done=false;
		double s=0;
		//calculate s for each bin
		//as a test, will calculate and output value of s for one bin
		//cerr<<block_size<<endl;
		for(int r=0; r<numBinsX; r++){	
			for(int c=0; c<numBinsY; c++){
				//error[][] is the accumulated sums for the current block, and the 
				if(sum_block_var[r][c]!=0) {
					s=double(block_size)*sum_block_var[r][c]/run_var[r][c];
					cerr<<block_size<<"  "<<sum_block_var[r][c]<<"  "<<s<<endl;
					done=true;
					break;
				}
			}
			if(done) break;
		}
	}
	
	//integrate_basins(histo2d, numBinsX, numBinsY, binWidth, total);
}

