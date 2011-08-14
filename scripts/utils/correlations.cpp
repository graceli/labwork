#include <iostream>
#include <cmath>

using namespace std;

int main(){
	double sum_sq_x=0;
	double sum_sq_y=0;
	
	double sum_coproduct=0;
	double* dataX=new double[200000];
	double* dataY=new double[200000];
	//read in the X and Y variable data
	//from stdin
	double x,y;
	int total=0;
	while(cin>>x){
		cin>>y;
		dataX[total+1]=x;
		dataY[total+1]=y;	
		total++;
	}	
//	cout<<"total="<<total<<endl;
	double mean_x = dataX[1];
	double mean_y = dataY[1];

	for(int i=2; i<=total; i++){
		double sweep = (i-1.0)/i;
		double delta_x = dataX[i] - mean_x;
		double delta_y = dataY[i] - mean_y;
	//	cout<<"i="<<i<<" "<<delta_x<<" "<<delta_y<<endl;
		sum_sq_x += pow(delta_x,2)*sweep;
		sum_sq_y += pow(delta_y,2)*sweep;
		sum_coproduct += delta_x*delta_y*sweep;
		mean_x += delta_x/i;
		mean_y += delta_y/i;
	}

	double pop_sd_x = sqrt( sum_sq_x / total );
	double pop_sd_y = sqrt( sum_sq_y / total );

	//cout<<pop_sd_x<<" "<<pop_sd_y<<endl;
	double covXY = sum_coproduct/total;
	double correlation = covXY/(pop_sd_x*pop_sd_y);

	cout<<correlation<<endl;
	return 0;
}
