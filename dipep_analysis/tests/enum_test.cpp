#include <iostream>

using namespace std;
enum group {NONE,
		CO0CO2_12,CO0NH3_12,NH1CO2_12,NH1NH3_12,
		CO0CO2_13,CO0NH3_13,NH1CO2_13,NH1NH3_13,
		CO0CO2_14,CO0NH3_14,NH1CO2_14,NH1NH3_14,
		MONO_CO0,MONO_NH1,MONO_CO2,MONO_NH3,
		MONO_TWO};

group switch_test(group var){
	switch(var){
		case NONE:
			cout<<"NONE"<<endl;
			return NONE;
		case CO0CO2_12:
			cout<<"CO0CO2_12"<<endl;
			return CO0CO2_12;
	}
	return CO0NH3_12;
}

int main(){


group var = MONO_TWO;

group ret_val = switch_test(var);

cout<<"ret_val =" <<(group)ret_val<<endl;

}
