#include "nr.h"
#include "parameters.h"

void NR::cacheCphi1(Mat_O_DP &Cphi1, Vec_I_DP &phi)
{
	// Computes the values of the cosine phi1 basis
	// functions at each collocation point.

	for (int n=1;n<=N+1;n++)
	{
		for (int i=1;i<=N;i++)
		{
			Cphi1[n-1][i-1]=cos((2.0*n-1.0)*phi[i-1]);
		}
	}
}
