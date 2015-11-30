#include "nr.h"
#include "parameters.h"


namespace NRquad2d
{
	DP yy1,yy2;
	DP xsav;
	Vec_DP coeffs(numunknowns+1);
	DP tol;
	DP (*nrfunc)(const DP, const DP,Vec_I_DP &coeffs);

	// The integrand f(x,y) evaluated at fixed x
	DP f2(const DP y)
	{
		return nrfunc(xsav,y,coeffs);
	}

	DP f1(const DP x)
	{
		xsav=x;
		return NR::adaptlob(f2,yy1,yy2,tol);
	}
}

DP NR::quad2d(DP func(const DP, const DP, Vec_I_DP &coeffs), const DP x1, const DP x2, const DP yy1, const DP yy2, Vec_I_DP &coeffs, const DP tol)
{
	// Returns the integral of a user-supplied function func over a three-dimensional
	// region specified by the limits x1,x2, and by the user-supplied functions 
	// yy1, yy2, z1 and z2. Integration is performed by calling adaptlob recusively.
	NRquad2d::nrfunc=func;
	NRquad2d::coeffs=coeffs;
	NRquad2d::tol=tol;
	NRquad2d::yy1=yy1;
	NRquad2d::yy2=yy2;
	return adaptlob(NRquad2d::f1,x1,x2,tol);
}

