#include "nr.h"

void NR::qrsolv(Mat_I_DP &a, Vec_I_DP &c, Vec_I_DP &d, Vec_IO_DP &b)
{
	// Solves the set of n linear equations A.x = b. a[0..n-1][0..n-1],
	// c[0..n-1] and d[0..n-1] are input as the output of the routine qrdcmp
	// and are not modified. b[0..n-1] is input as the right hand side
	// vector, and is overwritten with the solution vector on output.

	int i,j;
	DP sum,tau;

	int n=a.nrows();
	for (j=0;j<n-1;j++)
	{
		// Form Q'.b
		for (sum=0.0,i=j;i<n;i++)
		{
			sum += a[i][j]*b[i];
		}
		tau = sum/c[j];
		for (i=j;i<n;i++)
		{
			b[i] -= tau*a[i][j];
		}
	}
	// Solve R.x=Q'.b
	rsolv(a,d,b);
}

