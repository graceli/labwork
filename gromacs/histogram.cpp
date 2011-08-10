#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;
const int MAX = 2500;
int main(int argc, char* argv[]){
	//computes a histogram given a column data file	
	//assumes that data file inputted has only a single column
	//assumes that there are only 100,000 lines of data

	if(argc < 3){
		cerr<<"usage: histo <data-file> <bin-width> <normalize>"<<endl;
		return 0;
	}	

	ifstream file(argv[1]);
	double binWidth = atof(argv[2]);	

	//cerr<< binWidth<<endl;
	
	string nonorm;
	if(argv[3]==0){
		nonorm = "-norm";
	}else{
		nonorm = argv[3];
	}

//	cerr<<nonorm<<endl;

	double* data = new double[MAX];
	double aValue=0;

	//initialize histogram to 0 -- ie each bin
	//is init. to 0
	for(int i=0; i<MAX; i++){
		data[i]=0;	
	}

	//read in data and bin data
	int count=0; //count the number of data points read
	while(file>>aValue){
		int binNum = int(floor(aValue/binWidth));
		data[binNum]+=1;
		count++;
	}
//	cerr<<"read in "<<count<<endl;
	file.close();
	
	//print histogram
	for(int i=0; i<MAX; i++){
		//normalize by default unless no normalization is specified
		if(nonorm == "-nonorm"){
			cout<<i*binWidth<<" "<<data[i]<<endl;
		}else{
			cout<<i*binWidth<<" "<<data[i]/count<<endl;
		}
	}
	
	return 0;
}
