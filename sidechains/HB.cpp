#include "HB.h"

using namespace std;

//#define DEBUG_CASTING
//#define DEBUG_HB

HB::HB(real dDH, real dHdH, real alpha) 
:dDonorHeavy(dDH), dHdonorHeavy(dHdH), angle(alpha)
{
}

bool HB::isHbonded(t_pbc* pbc, rvec donor, rvec donorH, rvec acceptor, rvec acceptorH){
//H will be used in the general sense for the atom attached to the donor or acceptor heavy atom
//involved in the hydrogen bonding

#ifdef DEBUG_CASTING
        cout<<"In isHbonded:"<<endl;
        cout<<"DEBUG: "<<donor[0]<<" "<<donor[1]<<" "<<donor[2]<<endl;
        cout<<"DEBUG: "<<donorH[0]<<" "<<donorH[1]<<" "<<donorH[2]<<endl;
#endif

	real a, b, c, rAD, rHaHd, rAHd, rHaD;
        rvec dx;

	pbc_dx(pbc, donor, donorH, dx);
        a = norm(dx);

	pbc_dx(pbc, donorH, acceptor, dx);
        b = norm(dx);
	
        pbc_dx(pbc, donor, acceptor, dx);
        c = norm(dx);
	
	//calculate DSSP HB interaction energy
        pbc_dx(pbc, acceptor, donor, dx);
        rAD = 10*norm(dx);
        
        pbc_dx(pbc, acceptorH, donorH, dx);
        rHaHd=10*norm(dx);

        pbc_dx(pbc, acceptor, donorH, dx);
        rAHd=10*norm(dx);

        pbc_dx(pbc, acceptorH, donor, dx);	
        rHaD = 10*norm(dx);

	double dssp_energy = E_HB*((1.0/rAD) + (1.0/rHaHd) - (1.0/rAHd) - (1.0/rHaD));
	real alpha = cutoffAngle(a,b,c);

#ifdef DEBUG_HB
        cout<<"D-A = "<<c<<" DH-A = "<<b<<" D-H = "<<a<<" alpha = "<<alpha<<" dssp_energy="<<dssp_energy<<" ";
#endif

	if(c < dDonorHeavy && alpha > angle && b < dHdonorHeavy && dssp_energy < -0.5){

#ifdef DEBUG_HB
        cout<<"yes"<<endl;
#endif

		return true;
	} else {

#ifdef DEBUG_HB
        cout<<"no"<<endl;
#endif

        }

	return false;
}

real HB::cutoffAngle(real a, real b, real c){
	real alpha = acos((a*a+b*b-c*c)/(2*a*b));
	return alpha*180.0/PI;
}





