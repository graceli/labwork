#include <iostream>
#include <cmath>

using namespace std;

int main(){
	double v1[]={1,1};
	double v2[]={1,0};
	double v3[]={-1,0};
	double dot=0.0;
	cout<<"dot(v1,v2)=";
	for(int i=0; i<2; i++){
		cout<<v1[i]<<"."<<v3[i]<<"+";
		dot+=v1[i]*v3[i];
	}	
	cout<<"="<<dot<<endl;
	double norm_v1=sqrt(2);
	double norm_v3=1;
	cout<<"theta(v1,v3)="<<acos(dot/(norm_v1*norm_v3))*180/3.14159<<endl;
	return 0;
}
