#ifndef _HB_
#define _HB_

#define CPLUSPLUS

#include <cmath>
#include <iostream>

using namespace std;

#include "typedefs.h"
#include "pbc.h"
#include "vec.h"

class HB {
	public:
		//constructor: sets the hydrogen bond criteria, where angle is in degrees, and distances are in nm
		HB(real dDonorHeavy, real dHdonorHeavy, real angle);
                bool isHbonded(t_pbc* pbc, rvec donor, rvec donorH, rvec acceptor, rvec acceptorH);

	private:
                real cutoffAngle(real a, real b, real c);
		real dDonorHeavy;
		real dHdonorHeavy;
		real angle;
};

const double PI=3.141592;
const double E_HB=27.888;

#endif

