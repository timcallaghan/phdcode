#include "nr.h"
#include "parameters.h"

void NR::cacheSphi2(Mat_O_DP &Sphi2, Vec_I_DP &phi)
{
	// Computes the values of the cosine phi1 basis
	// functions at each collocation point.

	for (int n=1;n<=N+1;n++)
	{
		for (int i=1;i<=N;i++)
		{
			Sphi2[n-1][i-1]=sin(2.0*n*phi[i-1]);
		}
	}
}
