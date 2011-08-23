#include <iostream>

using namespace std;

int main(){
	string a="../blah.ext";
	cout<<a<<endl;
 	int pos=a.rfind('.', a.size()-1);
	cout<<"found . at"<<pos<<endl;
	int slash=a.rfind('/',a.size()-1);
	cout<<"found / at"<<slash<<endl;
	string sub = a.substr(slash+1,pos-slash-1);
	cout<<"the substr of length "<< pos-slash<<" is "<<sub<<endl;
	return 0;
}
