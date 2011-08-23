#include <iostream>
using namespace std;
void foo(int a[]){
	a[0]=888;
}
int main(){
	int b[]={0};
	cout<<"pre-foo "<<b[0]<<endl;
	foo(b);
	cout<<"post-foo "<<b[0]<<endl;
	return 0;
}
