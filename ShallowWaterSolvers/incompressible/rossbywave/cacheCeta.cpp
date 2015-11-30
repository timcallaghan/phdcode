#include "nr.h"
#include "parameters.h"

void NR::cacheCeta(Mat_O_DP &Ceta)
{
	// Computes the values of the cosine eta basis
	// functions at each collocation point.

	DP deleta = pi/(M*kappa);
	DP eps = deleta/2.0;
	for (int m=1;m<=M;m++)
	{
		for (int i=1;i<=M;i++)
		{
			Ceta[m-1][i-1]=cos((deleta*(i-1)+eps)*m*kappa);
		}
	}
}
