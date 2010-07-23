#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>

using namespace std;
const int MAX = 3000000;
int main(int argc, char* argv[]){
	/*
	 this program assumes the following:
	 -the first column is a time point unless num_cols == 1 
	 -data is for every picosecond 
	 -column data can be int or float
	 -output time point is the midpoint (2t+1)/2 of the time interval [t,t+1]
	*/

	if(argc < 4 ){
		cerr<<"moving <file> <num-cols> <period> <time>"<<endl;
		cerr<<"assumes that first column is time"<<endl;
		return 0;
	}
	
	ifstream infile(argv[1]);
	int num_cols = atoi(argv[2]);
	int period = atoi(argv[3]);
	int isFirstColTime = atoi(argv[4]);

	string line;
	double* tokens = new double[num_cols+1];	//array holding tokens belonging to one line
	double** data = new double*[MAX];		//array holding window sums
	for(int i=0; i<MAX; i++){
		data[i] = new double[num_cols];
	}

	int time=0;
	while(getline(infile, line)){
		if(line.size() == 0){
			break;
		}
		istringstream ist(line);
		for(int nToken = 0; nToken < num_cols; nToken++){
			ist>>data[time][nToken];
		}
		time++;
	}

	//compute and output moving average for each column of data
	double* window_sum = new double[num_cols];
	//initialize
	for(int w=0; w<num_cols; w++){
		window_sum[w]=0;
	}
	//cout<<"time = "<<time<<endl;

	for(int t=0; t<time-period+1; t++){
		cout<<double(t+t+period-1)/2<<" ";
		for(int c=1; c<num_cols; c++){
			for(int i=0; i<period; i++){
				if(isFirstColTime == 1){
					window_sum[c]+=data[t+i][c];
				}else{
					window_sum[c-1]+=data[t+i][c-1];
				}
			}
			
			cout<<double(window_sum[c])/period<<" ";
		}

		//reset column windo sums
		for(int w=0; w<num_cols; w++){
			window_sum[w]=0;
		}
		cout<<endl;
	}
	return 0;
}
