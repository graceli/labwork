#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

//max number of bins 
const int MAX = 10000;

int main(int argc, char* argv[]){
	//computes a histogram given a column data file	
	//assumes that data file inputted has only a single column
	//assumes that there are only 100,000 lines of data
	if(argc < 3){
		cerr<<"usage: histo <data-file> <bin-width> <normalize> <nm|ang>"<<endl;
		return 0;
	}	

	ifstream file(argv[1]);

	double binWidth = atof(argv[2]);	

	//cerr<< binWidth<<endl;
	
	string nonorm;
	nonorm = argv[3];

    string scale = argv[4];

    //data holds the histogram
	double* data = new double[MAX];

	//initialize histogram to 0 -- ie each bin
	//is init. to 0
	for(int i=0; i<MAX; i++){
		data[i]=0;	
	}

	//read in data and bin data
	double aValue=0;
	double time=0.0;
	int count=0; //count the number of data points read
    int maxBinNum = 0;
	while(file>>time){
		file>>aValue;

        int binNum=0;
		//use sarah's bin assignment method
        if(scale == "nm"){
    		binNum = int(floor(aValue/binWidth+0.5));
		//int binNum = int(floor(aValue/binWidth));
    		data[binNum]+=1;
        }else{
            binNum = int(floor(aValue*10/binWidth+0.5));
            data[binNum]++;
        }
		count++;

        if(maxBinNum < binNum){
            maxBinNum=binNum;
        }
	}
	file.close();
	
	cerr<<count<<"number of data points read in."<<endl;
    cerr<<"bin width of "<<binWidth<<endl;

	//print histogram
	for(int i = 0; i <= maxBinNum+1; i++){
		//normalize by default unless no normalization is specified
		if(nonorm == "-nonorm"){
			cout<<i*binWidth<<" "<<data[i]<<endl;
		}else{
			cout<<i*binWidth<<" "<<data[i]/count<<endl;
		}
	}
	
	return 0;
}
