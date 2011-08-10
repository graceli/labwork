#include "HB.h"
using namespace std;

HB::HB(double dDH, double dHdH, double alpha) 
:dDonorHeavy(dDH), dHdonorHeavy(dHdH), angle(alpha)
{
}

bool HB::isHbonded(double donor[3], double donorH[3], double acceptor[3], double acceptorH[3],  double box[3]){
//H will be used in the general sense for the atom attached to the donor or acceptor heavy atom
//involved in the hydrogen bonding

	double a,b,c;
	a = dist(donor, donorH, box);
	b = dist(donorH, acceptor, box);
	c = dist(donor, acceptor, box);
	
	//calculate DSSP HB interaction energy
	double rAD = 10*dist(acceptor,donor,box);
	double rHaHd = 10*dist(acceptorH,donorH,box);
	double rAHd = 10*dist(acceptor, donorH,box);
	double rHaD = 10*dist(acceptorH,donor,box);	
	double dssp_energy = E_HB*((1.0/rAD) + (1.0/rHaHd) - (1.0/rAHd) - (1.0/rHaD));
	
	if(c < dDonorHeavy && cutoffAngle(a,b,c) > angle && b < dHdonorHeavy && dssp_energy < -0.5){
		#ifdef DEBUG_HB
			cout<<"dDonorHeavy = "<<c<<" dHdonorHeavy = "<<b<<" dDonorDonorH = "<<a<<endl;
		#endif
		return true;
	}

	return false;
}

inline double HB::dist(double a1[3], double a2[3], double box[3]) {
	double dx_img = image_dist(a1[0] - a2[0], box[0]);
	double dy_img = image_dist(a1[1] - a2[1], box[1]);
	double dz_img = image_dist(a1[2] - a2[2], box[2]); 

	return sqrt(dx_img*dx_img + dy_img*dy_img + dz_img*dz_img);
}

//calculates the norm of a vector
inline double HB::norm(double v[3]){
	return sqrt(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));
}

//calculates the dot product of two vectors
inline double HB::dot(double v1[3], double v2[3]){
	return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}

inline double HB::image_dist(double unimaged_dist, double box_side){
	double abs_dist = abs(unimaged_dist);
	int num_boxes = (int)floor(abs_dist/box_side + 0.5);

	return abs_dist-box_side*num_boxes;
}

inline double HB::cutoffAngle(double a, double b, double c){
	double alpha = acos((a*a+b*b-c*c)/(2*a*b));
	return alpha*180.0/PI;
}
