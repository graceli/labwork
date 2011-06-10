#include "Vector.h"
#include <math.h>
double Distance(double *u, double *v)
{
    return sqrt(Distsq(u,v));
}
double DistPBC(double *u, double *v, double *boxk)
{
     if (boxk[0]>0)
     {
     int k = 0;
     double tmpdist = 0;
     for (k=0;k<3;k++)
     {
       double r = fabs(u[k]-v[k]);
       double nBoxes = floor(r/boxk[k] + 0.5) * boxk[k];
       double tmp1 = r - nBoxes;
       tmpdist += tmp1*tmp1;
     }
     return sqrt(tmpdist);
    }
    else
    {
       return Distance(u,v); /*Distance(u,v); */
    }
}
void Norm(double *x, double *xnorm)
{
  /* RETURNS INPUT VECTOR X NORMALIZED TO UNIT LENGTH.
     XNORM IS THE ORIGINAL LENGTH OF X.                         */
  double TEMP, TEMP1, TEMP2;

  TEMP = x[0];
  TEMP1 = x[1];
  TEMP2 = x[2];
  *xnorm = TEMP * TEMP + TEMP1 * TEMP1 + TEMP2 * TEMP2;
  if (*xnorm <= 0.0)
    return;
  *xnorm = sqrt(*xnorm);
  x[0] /= *xnorm;
  x[1] /= *xnorm;
  x[2] /= *xnorm;
}

double VLength(double *u)
{
    return sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
}
