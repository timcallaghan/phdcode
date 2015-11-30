#include "nr.h"
#include "parameters.h"

using namespace std;

void NR::loadLinCoeffs(Vec_IO_DP & coeffs)
{
	// Routine for loading the linearised coefficients
	// from file into a vector for use in the nonlinear
	// code

	int i;

	// First check to see if we have valid truncation values
	if (Nlin > N)
	{
		nrerror("N must be >= Nlin...make N larger");
	}
	// Now that we have passed the truncation test we load the
	// linearised coeffs
	Vec_DP lincoeffs(3*Nlin+1);
	ifstream inLinCoeffs("lincoeffs.txt", ios::in);
	if (!inLinCoeffs)
	{
		cerr << "File lincoeffs.txt could not be opened" << endl ;
		exit(1);
	}
	for (i=0;i<3*Nlin+1;i++)
	{
		inLinCoeffs >> lincoeffs[i];
	}
	inLinCoeffs.close();

	// Assign all the linearized starting guess coeffs
	for (i=1;i<=Nlin;i++)
	{
		coeffs[i-1]=lincoeffs[i-1];
		coeffs[M*N+i-1]=lincoeffs[Nlin+i-1];
		coeffs[2*M*N+N+i]=lincoeffs[2*Nlin+i-1];
	}
	// Assign the zonal flow structure to h
	coeffs[2*M*N] = h0 + w*Fr*Fr*(1/Ro+w)/4.0;
	coeffs[2*M*N+1] = w*Fr*Fr*(1/Ro+w)/4.0;
	// Define the wavespeed c
	coeffs[numunknowns]=lincoeffs[3*Nlin];
}

