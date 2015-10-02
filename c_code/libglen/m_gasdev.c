#include <math.h>

// Uses function ran1() to generate uniform random numbers on [0,1), then transforms
// to normal distribution with mean 0, var 1
// Use idum = -1 to initialize, then don't change idum between calls.
// modified to use seed value every time, not just every other.

float m_gasdev(long *idum)
{
	float ran1(long *idum);
	float fac,rsq,v1,v2;

	do {
	  v1=2.0*ran1(idum)-1.0;
	  v2=2.0*ran1(idum)-1.0;
	  rsq=v1*v1+v2*v2;
	} while (rsq >= 1.0 || rsq == 0.0);
	fac=sqrt(-2.0*log(rsq)/rsq);
	return v2*fac;
}
