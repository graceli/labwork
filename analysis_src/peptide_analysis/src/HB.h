#include <cmath>
#include <iostream>

using namespace std;
#ifndef _HB_
#define _HB_
class HB {
	public:
		//constructor: sets the hydrogen bond criteria, where angle is in degrees, and distances are in nm
		HB(double dDonorHeavy, double dHdonorHeavy, double angle);
		bool isHbonded(double donor[3], double donorH[3], double acceptor[3], double acceptorH[3], double box[3]);

	private:
		double dist(double a1[3], double a2[3], double box[3]);
		double norm(double v[3]);
		double dot(double v1[3], double v2[3]);
		double image_dist(double unimaged_dist, double box_side);
		double cutoffAngle(double a, double b, double c);
		double dDonorHeavy;
		double dHdonorHeavy;
		double angle;
};

const double PI=3.141592;
const double E_HB=27.888;

#endif

