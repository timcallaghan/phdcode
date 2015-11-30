#include "nr.h"
#include "fieldvars.h"
#include "parameters.h"

using namespace std;

// These functions evaluate the field variables
// at the grid point (kapeta, phi). They use the
// cached basis functions for faster evaluations


DP ulameval(Vec_I_DP &x, const int i, const int j, const DP phi, Mat_I_DP &Ceta, Mat_I_DP &Cphi1)
{
	DP value = 0.0;
	// Initialise to the zonal flow
    value = w*cos(phi);
    // Add in the rest of the series
    for (int m=1;m<=M;m++)
	{
	    for (int n=1;n<=N;n++)
		{
		    value = value + x[(m-1)*N+n-1]*Ceta[m-1][i-1]*Cphi1[n-1][j-1];
		}
	}
	return value;
}

DP uphieval(Vec_I_DP &x, const int i, const int j,Mat_I_DP &Seta, Mat_I_DP &Sphi2)
{
    // Initialise to the zonal flow
    DP value = 0.0;
    // Add in the rest of the series
    for (int m=1;m<=M;m++)
	{
	    for (int n=1;n<=N;n++)
		{
            value = value + x[M*N+(m-1)*N+n-1]*Seta[m-1][i-1]*Sphi2[n-1][j-1];
		}
	}
	return value;
}

DP heval(Vec_I_DP &x, const int i, const int j,Mat_I_DP &Ceta, Mat_I_DP &Cphi2)
{
    // Initialise to zero.
    DP value = 0.0;
    // Add in eta-independent series components
    // First add in constant term
    value = value + x[2*M*N];
    // Now add in the rest of the terms
	int n;
    for (n=1;n<=N;n++)
	{
        value = value + x[2*M*N+n]*Cphi2[n-1][j-1];
	}
    // Add in the case when n=1
	int m;
    for (m=1;m<=M-1;m++)
	{
        value = value - x[2*M*N+m*N+1]*Ceta[m-1][i-1]*(Cphi2[0][j-1]+1.0);
	}
    // Add in the rest of the series
    for (m=1;m<=M-1;m++)
	{
	    for (n=2;n<=N;n++)
		{
		    value = value + x[2*M*N+m*N+n]*Ceta[m-1][i-1]*pow(-1.0,(double) n)*(Cphi2[n-1][j-1]+Cphi2[n-2][j-1]);
		}
	}
	return value;
}

// Functions for evaluating the eta
// derivaties of the field variables at a given 
// point in the flow grid.

DP ulamevaleta(Vec_I_DP &x, const int i, const int j, Mat_I_DP &Seta, Mat_I_DP &Cphi1)
{
    // Initialise the sum
    DP value = 0.0;
    // Add in the rest of the series
    for (int m=1;m<=M;m++)
	{
	    for (int n=1;n<=N;n++)
		{
		    value = value - kappa*m*x[(m-1)*N+n-1]*Seta[m-1][i-1]*Cphi1[n-1][j-1];
		}
	}
	return value;
}

DP uphievaleta(Vec_I_DP &x, const int i, const int j,Mat_I_DP &Ceta, Mat_I_DP &Sphi2)
{
    // Initialise the sum
    DP value = 0.0;
    // Add in the rest of the series
    for (int m=1;m<=M;m++)
	{
	    for (int n=1;n<=N;n++)
		{
            value = value + kappa*m*x[M*N+(m-1)*N+n-1]*Ceta[m-1][i-1]*Sphi2[n-1][j-1];
		}
	}
	return value;
}

DP hevaleta(Vec_I_DP &x, const int i, const int j,Mat_I_DP &Seta, Mat_I_DP &Cphi2)
{
    // Initialise the sum
    DP value = 0.0;
    // Add in the case when n=1
	int m,n;
    for (m=1;m<=M-1;m++)
	{
        value = value + kappa*m*x[2*M*N+m*N+1]*Seta[m-1][i-1]*(Cphi2[0][j-1]+1.0);
	}
    // Add in the rest of the series
    for (m=1;m<=M-1;m++)
	{
	    for (n=2;n<=N;n++)
		{
		    value = value - kappa*m*x[2*M*N+m*N+n]*Seta[m-1][i-1]*pow(-1.0,(double) n)*(Cphi2[n-1][j-1]+Cphi2[n-2][j-1]);
		}
	}
	return value;
}

// Functions for evaluating the phi
// derivaties of the field variables at a given 
// point in the flow grid.


DP ulamevalphi(Vec_I_DP &x, const int i, const int j, const DP phi, Mat_I_DP &Ceta, Mat_I_DP &Sphi1)
{
    // Initialise to the zonal flow
    DP value = -w*sin(phi);
    // Add in the rest of the series
    for (int m=1;m<=M;m++)
	{
	    for (int n=1;n<=N;n++)
		{
		    value = value - (2*n-1)*x[(m-1)*N+n-1]*Ceta[m-1][i-1]*Sphi1[n-1][j-1];
		}
	}
	return value;
}

DP uphievalphi(Vec_I_DP &x, const int i, const int j,Mat_I_DP &Seta, Mat_I_DP &Cphi2)
{
    // Initialise to the zonal flow
    DP value = 0.0;
    // Add in the rest of the series
    for (int m=1;m<=M;m++)
	{
	    for (int n=1;n<=N;n++)
		{
            value = value + 2*n*x[M*N+(m-1)*N+n-1]*Seta[m-1][i-1]*Cphi2[n-1][j-1];
		}
	}
	return value;
}

DP hevalphi(Vec_I_DP &x, const int i, const int j,Mat_I_DP &Ceta, Mat_I_DP &Sphi2)
{
    // Initialise to the zonal flow
    DP value = 0.0;
    // Add in eta-independent series components
    // No constant terms so add in the rest of the phi terms
	int m,n;
    for (n=1;n<=N;n++)
	{
        value = value - 2*n*x[2*M*N+n]*Sphi2[n-1][j-1];
	}
    // Add in the case when n=1
    for (m=1;m<=M-1;m++)
	{
        value = value + x[2*M*N+m*N+1]*Ceta[m-1][i-1]*(2.0*Sphi2[0][j-1]+0.0);
	}
    // Add in the rest of the series
    for (m=1;m<=M-1;m++)
	{
	    for (n=2;n<=N;n++)
		{
		    value = value - x[2*M*N+m*N+n]*Ceta[m-1][i-1]*pow(-1.0,(double) n)*(2.0*n*Sphi2[n-1][j-1]+2.0*(n-1)*Sphi2[n-2][j-1]);
		}
	}
	return value;
}

// Function used in evaluating the volume integral of the free surface
DP surfint(const DP eta, const DP phi, Vec_I_DP &coeffs)
{
	// Compute the nonlinear series at phi and eta
	double h=0.0;
	// First compute the eta-independent terms
	int m,n;
	for (n=0;n<=N;n++)
	{
		h += coeffs[2*M*N+n]*cos(2.0*n*phi);
	}
	// Now compute the rest of the series
	for (m=1;m<=M-1;m++)
	{
		for (n=1;n<=N;n++)
		{
			h += coeffs[2*M*N+m*N+n]*cos(kappa*m*eta)*pow(-1.0,(double) n)*(cos(2.0*n*phi)+cos(2.0*(n-1)*phi));
		}
	}


	// Now that we've worked out the value of h we can make the
	// rest of the integrand
	DP value;
	value = (pow(h,3.0)/(3.0*arad*arad)+ h + pow(h,2.0)/arad)*cos(phi);

	return value;

}

// Functions used to evaluate the jacobian volume elements
DP I10(const int j)
{
	DP value;
	value = 2.0*pi*arad*arad*pow(-1.0,(double) j)/(Vzon*(4.0*j*j-1.0));
	return value;
}

DP I20(const int j, Vec_I_DP &coeffs)
{
	DP value = 0.0;
	DP term;
	for (int n=0;n<=N;n++)
	{
		term = pow(-1.0,(double) j-n)*(1.0-4.0*j*j-4.0*n*n)/(16.0*pow((double) j, 4.0)+pow((1.0-4.0*n*n),2.0)-8.0*j*j*(1.0+4.0*n*n));
		value += coeffs[2*M*N+n]*term;
	}
	value = -4.0*arad*pi*value/Vzon;	
	return value;
}


// Function used in evaluating the volume integral of the free surface
DP I3(const DP eta, const DP phi, Vec_I_DP &coeffs, const int mval, const int nval)
{
	// Compute the nonlinear series at phi and eta
	double h=0.0;
	// First compute the eta-independent terms
	int m,n;
	for (n=0;n<=N;n++)
	{
		h += coeffs[2*M*N+n]*cos(2.0*n*phi);
	}
	// Now compute the rest of the series
	for (m=1;m<=M-1;m++)
	{
		for (n=1;n<=N;n++)
		{
			h += coeffs[2*M*N+m*N+n]*cos(kappa*m*eta)*pow(-1.0,(double) n)*(cos(2.0*n*phi)+cos(2.0*(n-1)*phi));
		}
	}


	// Now that we've worked out the value of h we can make the
	// rest of the integrand
	DP value, dhdHval;
	if (mval == 0)
	{
		dhdHval = cos(phi*2.0*nval);
	}
	else
	{
		dhdHval = cos(kappa*mval*eta)*pow(-1.0,(double) nval)*(cos(2.0*nval*phi)+cos(2.0*(nval-1)*phi));
	}
	value = pow(h,2.0)*cos(phi)*dhdHval;
	return value;
}

DP I2i(const int i, const int j, Vec_I_DP &coeffs)
{
	DP value =0.0;
	DP numer,denom;
	for (int n=1;n<=N;n++)
	{
		numer = 96.0*(-1.0+2.0*j)*(-1.0+2.0*n)*(4.0*(-1.0+j)*j+(-3.0+2.0*n)*(1.0+2.0*n));
		denom = (-3.0+2.0*j-2.0*n)*(-1.0+2.0*j-2.0*n)*(1.0+2.0*j-2.0*n)*(3.0+2.0*j-2.0*n)*(-5.0+2.0*j+2.0*n)*(-3.0+2.0*j+2.0*n)*(-1.0+2.0*j+2.0*n)*(1.0+2.0*j+2.0*n);
		value += coeffs[2*M*N+i*N+n]*numer/denom;
	}
	value = -2.0*arad*pi*value/Vzon;
	return value;
}
