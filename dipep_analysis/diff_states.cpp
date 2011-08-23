#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char* argv[]){

	ifstream file1(argv[1]);
	ifstream file2(argv[2]);
	string line1,line2;
	int lineIndex=0;
	while(getline(file1, line1)){
		getline(file2,line2);
		if(line1 != line2){
			cout<<"time "<<lineIndex<<"\t"<<line1<<"\t"<<line2<<endl;
		}
		lineIndex++;
	}
	cout<<"total lines read = "<<lineIndex<<endl;
	return 0;
}