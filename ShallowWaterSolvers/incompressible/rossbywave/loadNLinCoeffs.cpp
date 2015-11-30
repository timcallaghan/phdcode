#include "nr.h"
#include "parameters.h"

using namespace std;

void NR::loadNLinCoeffs(Vec_IO_DP & coeffs)
{
	// Routine for loading the old nonlinear coefficients
	// from file into a vector for use in the nonlinear
	// code

	int i,m,n;

	// First check to see if we have valid truncation values
	if ((Nold > N) || (Mold > M))
	{
		nrerror("N and M must be >= Nold and Mold respectively...make N and/or M larger");
	}
	// Now that we have passed the truncation test we load the
	// old nonlinear coeffs
	Vec_DP nlincoeffs(3*Mold*Nold+2);
	ifstream inNLinCoeffs("current.txt", ios::in);
	if (!inNLinCoeffs)
	{
		cerr << "File current.txt could not be opened" << endl ;
		exit(1);
	}
	for (i=0;i<3*Mold*Nold+2;i++)
	{
		inNLinCoeffs >> nlincoeffs[i];
	}
	inNLinCoeffs.close();

	// Assign all the nonlinear starting guess coeffs
	// We need to take care of h separately so we do it later on...
	for (m=1;m<=Mold;m++)
	{
		for (n=1;n<=Nold;n++)
		{
			// This populates the ulam and uphi coeff components respectively
			coeffs[(m-1)*N+n-1] = nlincoeffs[(m-1)*Nold+n-1];
			coeffs[M*N+(m-1)*N+n-1] = nlincoeffs[Mold*Nold+(m-1)*Nold+n-1];
		}
	}
	// Now take care of the h initialisation
	for (m=1;m<=Mold-1;m++)
	{
		for (n=1;n<=Nold;n++)
		{
			coeffs[2*M*N+m*N+n] = nlincoeffs[2*Mold*Nold+m*Nold+n];
		}
	}
	// Since the h series contains 1 more coeff than the rest we need to 
	// initialise all those ete-independent coeffs here
	for (n=0;n<=Nold;n++)
	{
		coeffs[2*M*N+n] = nlincoeffs[2*Mold*Nold+n];
	}
	// Define the wavespeed c
	coeffs[numunknowns] = nlincoeffs[3*Mold*Nold+1];
}

