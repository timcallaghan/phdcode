#include "nr.h"
#include "parameters.h"


namespace NRquad2dmod
{
	DP yy1,yy2;
	DP xsav;
	Vec_DP coeffs(numunknowns+1);
	DP tol;
	int mval,nval;
	DP (*nrfunc)(const DP, const DP, Vec_I_DP &coeffs, const int, const int);

	// The integrand f(x,y) evaluated at fixed x
	DP f2(const DP y)
	{
		return nrfunc(xsav,y,coeffs,mval,nval);
	}

	DP f1(const DP x)
	{
		xsav=x;
		return NR::adaptlob(f2,yy1,yy2,tol);
	}
}

DP NR::quad2dmod(DP func(const DP, const DP, Vec_I_DP &coeffs, const int, const int), const DP x1, const DP x2, const DP yy1, const DP yy2, Vec_I_DP &coeffs, const DP tol, const int mval, const int nval)
{
	// Returns the integral of a user-supplied function func over a three-dimensional
	// region specified by the limits x1,x2, and by the user-supplied functions 
	// yy1, yy2, z1 and z2. Integration is performed by calling adaptlob recusively.
	NRquad2dmod::nrfunc=func;
	NRquad2dmod::coeffs=coeffs;
	NRquad2dmod::tol=tol;
	NRquad2dmod::yy1=yy1;
	NRquad2dmod::yy2=yy2;
	NRquad2dmod::mval=mval;
	NRquad2dmod::nval=nval;
	return adaptlob(NRquad2dmod::f1,x1,x2,tol);
}

